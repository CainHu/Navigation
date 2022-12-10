//
// Created by Cain on 2022/12/10.
//

#include "utils.h"

namespace eskf {
    void rotation_from_axis_angle(Matrix3f &r, const array<float, 3> &a) {
        /* Rodrigues's Formula:
        * a = n * θ
        * R = cosθ*I + (1 - cosθ)*n*n' + sinθ*n^
        * */
        const float a_norm_square = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
        if (a_norm_square < 1e-18f) {
            r << 1.f, -a[2], a[1],
                    a[2], 1.f, -a[0],
                    -a[1], a[0], 1.f;
        } else {
            const float a_norm = sqrtf(a_norm_square);
            const array<float, 3> a_unit = {a[0] / a_norm, a[1] / a_norm, a[2] / a_norm};
            const float theta = a_norm;
            const float cos_theta = cosf(theta), sin_theta = sinf(theta);
            const float tmp = 1.f - cos_theta;

            const float xx = a_unit[0] * a_unit[0] * tmp;
            const float xy = a_unit[0] * a_unit[1] * tmp;
            const float xz = a_unit[0] * a_unit[2] * tmp;
            const float yy = a_unit[1] * a_unit[1] * tmp;
            const float yz = a_unit[1] * a_unit[2] * tmp;
            const float zz = a_unit[2] * a_unit[2] * tmp;

            const float sx = sin_theta * a_unit[0];
            const float sy = sin_theta * a_unit[1];
            const float sz = sin_theta * a_unit[2];

            r(0, 0) = cos_theta + xx;
            r(0, 1) = xy - sz;
            r(0, 2) = xz + sy;
            r(1, 0) = xy + sz;
            r(1, 1) = cos_theta + yy;
            r(1, 2) = yz - sx;
            r(2, 0) = xz - sy;
            r(2, 1) = yz + sx;
            r(2, 2) = cos_theta + zz;
        }
    }

    void quaternion_from_axis_angle(Quaternionf &q, const array<float, 3> &a) {
        /*
        a = n * θ
        q = [cos(θ/2), n*sin(θ/2)]
        */

        const float a_norm_square = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
        if (a_norm_square < 1e-18f) {
            q.w() = sqrtf(1.f - 0.25f * a_norm_square);
            q.x() = 0.5f * a[0];
            q.y() = 0.5f * a[1];
            q.z() = 0.5f * a[2];
        } else {
            const float a_norm = sqrtf(a_norm_square);
            const array<float, 3> a_unit = {a[0] / a_norm, a[1] / a_norm, a[2] / a_norm};
            const float half_theta = 0.5f * a_norm;
            const float c = cosf(half_theta), s = sinf(half_theta);

            q.w() = c;
            q.x() = s * a_unit[0];
            q.y() = s * a_unit[1];
            q.z() = s * a_unit[2];
        }
    }

    void rotation_from_axis_angle(Matrix3f &r, const Vector3f &a) {
        /* Rodrigues's Formula:
        * a = n * θ
        * R = cosθ*I + (1 - cosθ)*n*n' + sinθ*n^
        * */
        const float a_norm_square = a.squaredNorm();
        if (a_norm_square < 1e-18f) {
            r << 1.f, -a[2], a[1],
                    a[2], 1.f, -a[0],
                    -a[1], a[0], 1.f;
        } else {
            const float a_norm = sqrtf(a_norm_square);
            const Vector3f a_unit = a / a_norm;
            const float theta = a_norm;
            const float cos_theta = cosf(theta), sin_theta = sinf(theta);
            const float tmp = 1.f - cos_theta;

            const float xx = a_unit[0] * a_unit[0] * tmp;
            const float xy = a_unit[0] * a_unit[1] * tmp;
            const float xz = a_unit[0] * a_unit[2] * tmp;
            const float yy = a_unit[1] * a_unit[1] * tmp;
            const float yz = a_unit[1] * a_unit[2] * tmp;
            const float zz = a_unit[2] * a_unit[2] * tmp;

            const float sx = sin_theta * a_unit[0];
            const float sy = sin_theta * a_unit[1];
            const float sz = sin_theta * a_unit[2];

            r(0, 0) = cos_theta + xx;
            r(0, 1) = xy - sz;
            r(0, 2) = xz + sy;
            r(1, 0) = xy + sz;
            r(1, 1) = cos_theta + yy;
            r(1, 2) = yz - sx;
            r(2, 0) = xz - sy;
            r(2, 1) = yz + sx;
            r(2, 2) = cos_theta + zz;
        }
    }

    void quaternion_from_axis_angle(Quaternionf &q, const Vector3f &a) {
        /*
        a = n * θ
        q = [cos(θ/2), n*sin(θ/2)]
        */

        const float a_norm_square = a.squaredNorm();
        if (a_norm_square < 1e-18f) {
            q.w() = sqrtf(1.f - 0.25f * a_norm_square);
            q.x() = 0.5f * a[0];
            q.y() = 0.5f * a[1];
            q.z() = 0.5f * a[2];
        } else {
            const float a_norm = sqrtf(a_norm_square);
            const array<float, 3> a_unit = {a[0] / a_norm, a[1] / a_norm, a[2] / a_norm};
            const float half_theta = 0.5f * a_norm;
            const float c = cosf(half_theta), s = sinf(half_theta);

            q.w() = c;
            q.x() = s * a_unit[0];
            q.y() = s * a_unit[1];
            q.z() = s * a_unit[2];
        }
    }

    void normalize_angle(float &angle) {
        while (angle > M_PI) {
            angle -= 2.f * M_PI;
        }
        while (angle < -M_PI) {
            angle += 2.f * M_PI;
        }
    }

    float normalize_angle(float angle) {
        while (angle > M_PI) {
            angle -= 2.f * M_PI;
        }
        while (angle < -M_PI) {
            angle += 2.f * M_PI;
        }
        return angle;
    }
}