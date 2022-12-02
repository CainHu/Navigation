#include "gen_utils.h"

namespace generator {
    using namespace std;
    using namespace Eigen;

    Matrix3d askew(const Vector3d &vec) {
        Matrix3d mat;

        mat << 0., -vec[2], vec[1],
            vec[2], 0., -vec[0],
            -vec[1], vec[0], 0.;

        return mat;       
    }

    Matrix3d euler2rot(const Vector3d &euler) {
        double cx = cosf(euler[0]), sx = sinf(euler[0]);
        double cy = cosf(euler[1]), sy = sinf(euler[1]);
        double cz = cosf(euler[2]), sz = sinf(euler[2]);

        Matrix3d rot;
        rot << cy*cz,cz*sx*sy - cx*sz,cx*cz*sy + sx*sz,
                cy*sz,cx*cz + sx*sy*sz,-(cz*sx) + cx*sy*sz,
                -sy,cy*sx,cx*cy;
        return rot;
    }

    Vector3d rot2euler(const Matrix3d &rot) {
        Vector3d euler;

        double sy = sqrtf(rot(0, 0) * rot(0, 0) + rot(1, 0) * rot(1, 0));
        if (sy < 1e-6f) {
            euler[0] = atan2f(-rot(1, 2), rot(1, 1));
            euler[1] = atan2f(-rot(2, 0), sy);
            euler[2] = 0.;
        } else {
            euler[0] = atan2f(rot(2, 1), rot(2, 2));
            euler[1] = atan2f(-rot(2, 0), sy);
            euler[2] = atan2f(rot(1, 0), rot(0, 0));
        }

        return euler;
    }

    Quaterniond euler2quat(const Vector3d &euler) {
        Quaterniond quat;

        double cx = cosf(0.5 * euler[0]), sx = sinf(0.5 * euler[0]);
        double cy = cosf(0.5 * euler[1]), sy = sinf(0.5 * euler[1]);
        double cz = cosf(0.5 * euler[2]), sz = sinf(0.5 * euler[2]);

        quat.w() = cz * cy * cx + sz * sy * sx;
        quat.x() = cz * cy * sx - sz * sy * cx;
        quat.y() = sz * cy * sx + cz * sy * cx;
        quat.z() = sz * cy * cx - cz * sy * sx;

        return quat;
    }

    Vector3d quat2euler(const Quaterniond &quat) {
        Vector3d euler;

        euler[0] = atan2f(2. * (quat.w() * quat.x() + quat.y() * quat.z()), 1. - 2. * (quat.x() * quat.x() + quat.y() * quat.y()));
        euler[1] = asinf(2. * (quat.w() * quat.y() - quat.x() * quat.z()));
        euler[2] = atan2f(2. * (quat.w() * quat.z() + quat.x() * quat.y()), 1. - 2. * (quat.y() * quat.y() + quat.z() * quat.z()));

        return euler;
    }

    Matrix3d quat2rot(const Quaterniond &quat) {
        Matrix3d rot;

        const double q00 = quat.w() * quat.w(), q01 = quat.w() * quat.x(), q02 = quat.w() * quat.y(), q03 = quat.w() * quat.z();
        const double q11 = quat.x() * quat.x(), q12 = quat.x() * quat.y(), q13 = quat.x() * quat.z();
        const double q22 = quat.y() * quat.y(), q23 = quat.y() * quat.z();
        const double q33 = quat.z() * quat.z();

        rot << q00 + q11 - q22 - q33, 2. * (q12 - q03), 2. * (q02 + q13),
            2. * (q12 + q03), q00 - q11 + q22 - q33, 2. * (q23 - q01),
            2. * (q13 - q02), 2. * (q01 + q23), q00 - q11 - q22 + q33;

        return rot;       
    }

    Quaterniond rot2quat(const Matrix3d &rot) {
        Quaterniond quat;

        if (rot(0, 0) >= rot(1, 1) + rot(2, 2)) {
            quat.x() = 0.5 * sqrtf(1. + rot(0, 0) - rot(1, 1) - rot(2, 2));
            quat.w() = (rot(2, 1) - rot(1, 2)) / (4. * quat.x());
            quat.y() = (rot(0, 1) + rot(1, 0)) / (4. * quat.x());
            quat.z() = (rot(0, 2) + rot(2, 0)) / (4. * quat.x());
        } else if (rot(1, 1) >= rot(0, 0) + rot(2, 2)) {
            quat.y() = 0.5 * sqrtf(1. - rot(0, 0) + rot(1, 1) - rot(2, 2));
            quat.w() = (rot(0, 2) - rot(2, 0)) / (4. * quat.y());
            quat.x() = (rot(0, 1) + rot(1, 0)) / (4. * quat.y());
            quat.z() = (rot(1, 2) + rot(2, 1)) / (4. * quat.y());
        } else if (rot(2, 2) >= rot(0, 0) + rot(1, 1)) {
            quat.z() = 0.5 * sqrtf(1. - rot(0, 0) - rot(1, 1) + rot(2, 2));
            quat.w() = (rot(1, 0) - rot(0, 1)) / (4. * quat.z());
            quat.x() = (rot(0, 2) + rot(2, 0)) / (4. * quat.z());
            quat.y() = (rot(1, 2) + rot(2, 1)) / (4. * quat.z());
        } else {
            quat.w() = 0.5 * sqrtf(1. + rot(0, 0) + rot(1, 1) + rot(2, 2));
            quat.x() = (rot(2, 1) - rot(1, 2)) / (4. * quat.w());
            quat.y() = (rot(0, 2) - rot(2, 0)) / (4. * quat.w());
            quat.z() = (rot(1, 0) - rot(0, 1)) / (4. * quat.w());
        }

        return quat;
    }

    Matrix3d rv2rot(const Vector3d &rv) {
        double a, b;
        const double rv2 = rv.squaredNorm();

        if (rv2 < 1e-8f) {
            a = 1 - rv2 * (1. / 6. - rv2 / 120.);
            b = 0.5 - rv2 * (1. / 24. - rv2 / 720.);
        } else {
            double norm = sqrtf(rv2);
            a = sinf(norm) / norm;
            b = (1. - cosf(norm)) / rv2;
        }

        Matrix3d rv_hat = askew(rv);

        return Matrix3d::Identity() + a * rv_hat + b * rv_hat * rv_hat;
    }

    Quaterniond rv2quat(const Vector3d &rv) {
        double q0, q1, q2, q3, s;
        const double rv2 = rv.squaredNorm();

        if (rv2 < 1e-8f) {
            q0 = 1. - rv2 * (1. / 8. - rv2 / 384.);
            s = 1. / 2. - rv2 * (1. / 48. - rv2 / 3840.);
        } else {
            double norm = sqrtf(rv2);
            q0 = cosf(0.5 * norm);
            s = sinf(0.5 * norm) / norm;
        }

        q1 = s * rv[0];
        q2 = s * rv[1];
        q3 = s * rv[2];

        return Quaterniond(q0, q1, q2, q3);
    }

    Vector3d quat2rv(Quaterniond quat) {
        if (quat.w() < 0) {
            quat.w() = -quat.w();
            quat.x() = -quat.x();
            quat.y() = -quat.y();
            quat.z() = -quat.z();
        }

        double b;
        double nm_half = acosf(quat.w());
        if (nm_half > 1e-12f) {
            b = 2. * nm_half / sinf(nm_half);
        } else {
            b = 2.;
        }

        return Vector3d(b * quat.x(), b * quat.y(), b * quat.z());
    }
}
