#include "generator.h"
#include "utils.h"
#include "gvar.h"

using namespace std;
using namespace Eigen;
using namespace generator;

void Generator::traj_gen(const Vector3d &euler0, const Vector3d &vn0, const Vector3d &pos0, const vector<Matrix<double, 5, 1>> &wat, const double ts) {
    double t_sum = 0.;
    for (auto & var : wat) {
        t_sum += var[4];
    }

    int len = int(t_sum / ts);
    _euler.resize(len);
    _vn.resize(len);
    _pos.resize(len);

    int index = 0;
    _euler[index] = euler0;
    _vn[index] = vn0;
    _pos[index] = pos0;
    ++index;

    Vector3d vb = euler2rot(euler0) * vn0;
    double vbx = vb[0];
    double roll = euler0[0], pitch = euler0[1], yaw = euler0[2];

    Fir roll_filter(roll), pitch_filter(pitch), yaw_filter(yaw), vbx_filter(vbx);

    Earth eth;
    Vector3d vn01;
    for (auto & var : wat) {
        for (double t = ts; t < var[4] + 0.1f * ts; t += ts) {
            roll += var[0] * ts;
            pitch += var[1] * ts;
            yaw += var[2] * ts;
            vbx += var[3] * ts;

            _euler[index][0] = roll_filter(roll);
            _euler[index][1] = pitch_filter(pitch);
            _euler[index][2] = yaw_filter(yaw);
            _vn[index] = euler2rot(_euler[index]) * Vector3d(vbx_filter(vbx), 0., 0.);
            vn01 = 0.5 * (_vn[index - 1] + _vn[index]);
            earth(_pos[index - 1], vn01, eth);
            _pos[index] = _pos[index - 1] + Vector3d(vn01[0] / eth.RMh, vn01[1] / eth.clRNh, -vn01[2]) * ts;
            ++index;
        }
    }

    _euler.resize(index);
    _vn.resize(index);
    _pos.resize(index);
}

void Generator::ev2imu(const vector<Vector3d> &euler, const vector<Vector3d> &vn, const vector<Vector3d> &pos, const double ts) {
    Vector3d wm0 = Vector3d::Zero(), vm0 = Vector3d::Zero();
    Matrix3d I = Matrix3d::Identity();

    _wm.resize(euler.size() - 1);
    _vm.resize(_wm.size());

    Earth eth;
    Quaterniond qbb;
    Vector3d phim, wm1, vm1, dvnsf;
    Matrix3d Cnb0;
    for (unsigned int i = 1; i < euler.size(); ++i) {
        earth(0.5 * (pos[i - 1] + pos[i]), 0.5 * (vn[i - 1] + vn[i]), eth);
        qbb = euler2quat(euler[i - 1]).inverse() * rv2quat(eth.wnin * ts) * euler2quat(euler[i]);
        phim = quat2rv(qbb);
        wm1 = (I + askew(1./12. * wm0)).fullPivLu().solve(phim);
        dvnsf = vn[i] - vn[i - 1] - eth.gcc * ts;
        Cnb0 = euler2rot(euler[i - 1]);
        vm1 = (I + 0.5 * askew(1./6. * wm0 + wm1)).fullPivLu().solve(Cnb0.transpose() * (I + askew(0.5 * eth.wnin * ts))) * dvnsf - 1./12. * vm0.cross(wm1);
        _wm[i - 1] = wm1;
        _vm[i - 1] = vm1;
        wm0 = wm1;
        vm0 = vm1;
    }
}

void Generator::earth(const Vector3d &pos, const Vector3d &vn, Earth &eth) {
    eth.sl = sin(pos[0]);
    eth.cl = cos(pos[0]);
    eth.tl = eth.sl / eth.cl;
    eth.sl2 = eth.sl * eth.sl;

    double sl4 = eth.sl2 * eth.sl2;
    double sq = 1. - e2 * eth.sl2;
    double sq2 = sqrt(sq);

    eth.RMh = Re * (1. - e2) / sq / sq2 + pos[2];
    eth.RNh = Re / sq2 + pos[2];
    eth.clRNh = eth.cl * eth.RNh;
    eth.wnie = wie * Vector3d(eth.cl, 0., -eth.sl);
    eth.vn = vn;
    eth.wnen = Vector3d(vn[1] / eth.RNh, -vn[0] / eth.RMh, -vn[1] / eth.RNh * eth.tl);
    eth.wnin = eth.wnie + eth.wnen;
    eth.wnien = eth.wnie + eth.wnin;

    // grs80重力模型
    double gLh = g0 * (1. + 5.27094e-3 * eth.sl2 + 2.32718e-5 * sl4) - 3.086e-6 * pos[2];
    eth.gn = Vector3d(0., 0., gLh);
    eth.gcc = eth.gn - eth.wnien.cross(vn);

}

void Generator::ev2imu(const double ts) {
    ev2imu(_euler, _vn, _pos, ts);
}

Quaterniond Generator::quat_add_phi(const Quaterniond &quat, const Vector3d &phi) {
    return rv2quat(-phi) * quat;
}

Quaterniond Generator::quat_del_phi(const Quaterniond &quat, const Vector3d &phi) {
    return rv2quat(phi) * quat;
}

Vector3d Generator::qq2phi(const Quaterniond &qpb, const Quaterniond &qnb) {
    Quaterniond qe = qnb * qpb.inverse();
    return quat2rv(qe);
}

template<unsigned int N>
void Generator::cnscl(const array<Vector3d, N> &wm, const array<Vector3d, N> &vm, Vector3d &phim, Vector3d &dvbm) {
    static const array<array<double, 5>, 5> cs {
        {{2. / 3., 0., 0., 0., 0.},
         {9. / 20., 27. / 20., 0., 0., 0.},
         {54. / 105., 92. / 105., 214. / 105., 0., 0.},
         {250. / 504., 525. / 504., 650. / 504., 1375. / 504., 0.},
         {2315. / 4620., 4558. / 4620., 7296. / 4620., 7834. / 4620., 15797. / 4620.}}
    };

    Vector3d wmm, vmm, dphim, scullm;
    dphim.setZero();
    scullm.setZero();

    for (unsigned int i = 0; i < N; ++i) {
        wmm += wm[i];
        vmm += vm[i];
    }

    if (N > 1) {
        Vector3d csw, csv;
        csw.setZero();
        csv.setZero();

        for (unsigned int i = 0; i < N - 1; ++i) {
            csw += cs[N - 1][i] * wm[i];
            csv += cs[N - 1][i] * vm[i];
        }

        dphim = csw.cross(wm[N - 1]);
        scullm = csw.cross(vm[N - 1]) + csv.cross(wm[N - 1]);
    }

    phim = (wmm + dphim);
    dvbm = (vmm + 0.5 * wmm.cross(vmm) + scullm);
}

void Generator::cnscl(const vector<Vector3d> &wm, const vector<Vector3d> &vm, Vector3d &phim, Vector3d &dvbm) {
    static const array<array<double, 5>, 5> cs {
        {{2. / 3., 0., 0., 0., 0.},
         {9. / 20., 27. / 20., 0., 0., 0.},
         {54. / 105., 92. / 105., 214. / 105., 0., 0.},
         {250. / 504., 525. / 504., 650. / 504., 1375. / 504., 0.},
         {2315. / 4620., 4558. / 4620., 7296. / 4620., 7834. / 4620., 15797. / 4620.}}
    };

    unsigned int N = wm.size();

    Vector3d wmm, vmm, dphim, scullm;
    dphim.setZero();
    scullm.setZero();

    for (unsigned int i = 0; i < N; ++i) {
        wmm += wm[i];
        vmm += vm[i];
    }

    if (N > 1) {
        Vector3d csw, csv;
        csw.setZero();
        csv.setZero();

        for (unsigned int i = 0; i < N - 1; ++i) {
            csw += cs[N - 1][i] * wm[i];
            csv += cs[N - 1][i] * vm[i];
        }

        dphim = csw.cross(wm[N - 1]);
        scullm = csw.cross(vm[N - 1]) + csv.cross(wm[N - 1]);
    }

    phim = (wmm + dphim);
    dvbm = (vmm + 0.5 * wmm.cross(vmm) + scullm);
}

void Generator::imu_add_err(const vector<Vector3d> &wm, const vector<Vector3d> &vm, const Vector3d &eb, const Vector3d &web, const Vector3d &db, const Vector3d &wdb, const double ts,
                 vector<Vector3d> &wm_err, vector<Vector3d> &vm_err) {
    unsigned int m = wm.size();
    double sts = sqrtf(ts);

    for (unsigned int i = 0; i < m; ++i) {
        wm_err[i] = wm[i] + ts * eb + sts * web * Vector3d(_dist(_random_engine), _dist(_random_engine), _dist(_random_engine));
        vm_err[i] = vm[i] + ts * db + sts * wdb * Vector3d(_dist(_random_engine), _dist(_random_engine), _dist(_random_engine));
    }                        
}