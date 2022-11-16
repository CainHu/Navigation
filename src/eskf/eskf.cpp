//
// Created by Cain on 2022/11/11.
//

#include "eskf.h"
#include <cfloat>
#include <iostream>

using namespace std;
using namespace eskf;

void ESKF::initialize() {
    _rot.setIdentity();

    reset_state();
    reset_error_state();
    reset_covariance_matrix(0, 16, _q_cov);
    regular_covariance_to_symmetric();
    reset_accmulator_cov();
}

void ESKF::rotation_from_axis_angle(Matrix3f &r, const array<float, 3> &a) const {
    /* Rodrigues's Formula:
    * a = n * θ
    * R = cosθ*I + (1 - cosθ)*n*n' + sinθ*n^
    * */
    const float a_norm_square = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    if (a_norm_square < _dt2 * 1e-12f) {
        r.setIdentity();
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

void ESKF::quaternion_from_axis_angle(Quaternionf &q, const array<float, 3> &a) const {
    /*
    a = n * θ
    q = [cos(θ/2), n*sin(θ/2)]
    */

    const float a_norm_square = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    if (a_norm_square < _dt2 * 1e-12f) {
        q.setIdentity();
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

void ESKF::predict_state(const Vector3f &w, const Vector3f &a) {
    // state: [p, v, bg, ba, g]
    // error_state : [δp, δv, δθ, δbg, δba, δg]

    /* Nominal Model:
    * p = p + v * Δt + 0.5 * (R * (a - ba) + g) * Δt^2
    * v = v + (R * (a - ba) + g) * Δt
    * R = R * Exp((w - bg) * Δt)
    * bg = bg
    * ba = ba
    * g = g
    * */

    // a_body_corr = a - ba
    const Vector3f a_body_corr = {a[0] - _ba[0], a[1] - _ba[1], a[2] - _ba[2]};

    // a_world = R * (a - ba) + [0;0;g]
    const Vector3f a_world(_rot(0, 0) * a_body_corr[0] + _rot(0, 1) * a_body_corr[1] +_rot(0, 2) * a_body_corr[2],
                           _rot(1, 0) * a_body_corr[0] + _rot(1, 1) * a_body_corr[1] +_rot(1, 2) * a_body_corr[2],
                           _rot(2, 0) * a_body_corr[0] + _rot(2, 1) * a_body_corr[1] +_rot(2, 2) * a_body_corr[2] + _g);

    // v
    const Vector3f v_last = _v;

    // v' = v + a * Δt
    _v += a_world * _dt;

    // p' = p + 0.5 * (v' + v) * Δt
    _p += 0.5f * (_v + v_last) * _dt;

    // bg' = bg
    // ba' = ba
    // g = g

    // axis_angle = (w - bg) * Δt
    const array<float, 3> axis_angle = {(w[0] - _bg[0]) * _dt, (w[1] - _bg[1]) * _dt, (w[2] - _bg[2]) * _dt};

    // // ΔR = Exp((w - bg) * Δt)
    // Matrix3f delta_rot {};
    // rotation_from_axis_angle(delta_rot, axis_angle);

    // // R' = R * ΔR
    // const Matrix3f rot = _rot;
    // _rot = rot * delta_rot;

    // Δq = Exp((w - bg) * Δt)
    Quaternionf delta_q {};
    quaternion_from_axis_angle(delta_q, axis_angle);

    // q' = q * Δq
    const Quaternionf q = _q;
    _q = q * delta_q;
    _q.normalize();

    _rot = _q;
}

void ESKF::predict_covariance(const Vector3f &w, const Vector3f &a) {
    // R * dt               (乘法: 3, 加法: 0)
    const Matrix3f rot = _rot * _dt;

    // x = R * (a - ba)^ * dt       (乘法: 9, 加法: 12)
    const array<float, 3> a_corr = {a[0] - _ba[0], a[1] - _ba[1], a[2] - _ba[2]};
    const array<array<float, 3>, 3> x = {{{rot(0, 1)*a_corr[2] - rot(0, 2)*a_corr[1], rot(0, 2)*a_corr[0] - rot(0, 0)*a_corr[2], rot(0, 0)*a_corr[1] - rot(0, 1)*a_corr[0]},
                                          {rot(1, 1)*a_corr[2] - rot(1, 2)*a_corr[1], rot(1, 2)*a_corr[0] - rot(1, 0)*a_corr[2], rot(1, 0)*a_corr[1] - rot(1, 1)*a_corr[0]},
                                          {rot(2, 1)*a_corr[2] - rot(2, 2)*a_corr[1], rot(2, 2)*a_corr[0] - rot(2, 0)*a_corr[2], rot(2, 0)*a_corr[1] - rot(2, 1)*a_corr[0]}}};

    // y = Exp(-(w - bg) * Δt)      (乘法: 3, 加法: 3)
    const array<float, 3> axis_angle = {(_bg[0] - w[0]) * _dt, (_bg[1] - w[1]) * _dt, (_bg[2] - w[2]) * _dt};
    Matrix3f y {};
    rotation_from_axis_angle(y, axis_angle);

    // P11_ = P11 + (P21 + P12) * dt + P22 * dt^2   (乘法: 12, 加法: 18)
    _cov[0][0] += _dt*(_cov[0][3] + _cov[0][3] + _cov[3][3]*_dt);
    _cov[0][1] += _dt*(_cov[0][4] + _cov[1][3] + _cov[3][4]*_dt);
    _cov[0][2] += _dt*(_cov[0][5] + _cov[2][3] + _cov[3][5]*_dt);
    _cov[1][1] += _dt*(_cov[1][4] + _cov[1][4] + _cov[4][4]*_dt);
    _cov[1][2] += _dt*(_cov[1][5] + _cov[2][4] + _cov[4][5]*_dt);
    _cov[2][2] += _dt*(_cov[2][5] + _cov[2][5] + _cov[5][5]*_dt);

    // P16_ = P16 + P26 * dt    (乘法: 3, 加法: 3)
    _cov[0][15] += _cov[3][15] * _dt;
    _cov[1][15] += _cov[4][15] * _dt;
    _cov[2][15] += _cov[5][15] * _dt;

    // P15_ = P15 + P25 * dt    (乘法: 9, 加法: 9)
    _cov[0][12] += _cov[3][12] * _dt;
    _cov[0][13] += _cov[3][13] * _dt;
    _cov[0][14] += _cov[3][14] * _dt;
    _cov[1][12] += _cov[4][12] * _dt;
    _cov[1][13] += _cov[4][13] * _dt;
    _cov[1][14] += _cov[4][14] * _dt;
    _cov[2][12] += _cov[5][12] * _dt;
    _cov[2][13] += _cov[5][13] * _dt;
    _cov[2][14] += _cov[5][14] * _dt;

    // P14_ = P14 + P24 * dt    (乘法: 9, 加法: 9)
    _cov[0][9] += _cov[3][9] * _dt;
    _cov[0][10] += _cov[3][10] * _dt;
    _cov[0][11] += _cov[3][11] * _dt;
    _cov[1][9] += _cov[4][9] * _dt;
    _cov[1][10] += _cov[4][10] * _dt;
    _cov[1][11] += _cov[4][11] * _dt;
    _cov[2][9] += _cov[5][9] * _dt;
    _cov[2][10] += _cov[5][10] * _dt;
    _cov[2][11] += _cov[5][11] * _dt;

    // P13 + P23 * dt           (乘法: 9, 加法: 9)
    const array<array<float, 3>, 3> c1 = {{{_cov[0][6] + _cov[3][6]*_dt,_cov[0][7] + _cov[3][7]*_dt,_cov[0][8] + _cov[3][8]*_dt},
                                           {_cov[1][6] + _cov[4][6]*_dt,_cov[1][7] + _cov[4][7]*_dt,_cov[1][8] + _cov[4][8]*_dt},
                                           {_cov[2][6] + _cov[5][6]*_dt,_cov[2][7] + _cov[5][7]*_dt,_cov[2][8] + _cov[5][8]*_dt}}};

    // P13_ = (P13 + P23 * dt) * Y' - P14_ * dt     (乘法: 36, 加法: 27)
    _cov[0][6] = c1[0][0]*y(0, 0) + c1[0][1]*y(0, 1) + c1[0][2]*y(0, 2) - _cov[0][9]*_dt;
    _cov[0][7] = c1[0][0]*y(1, 0) + c1[0][1]*y(1, 1) + c1[0][2]*y(1, 2) - _cov[0][10]*_dt;
    _cov[0][8] = c1[0][0]*y(2, 0) + c1[0][1]*y(2, 1) + c1[0][2]*y(2, 2) - _cov[0][11]*_dt;
    _cov[1][6] = c1[1][0]*y(0, 0) + c1[1][1]*y(0, 1) + c1[1][2]*y(0, 2) - _cov[1][9]*_dt;
    _cov[1][7] = c1[1][0]*y(1, 0) + c1[1][1]*y(1, 1) + c1[1][2]*y(1, 2) - _cov[1][10]*_dt;
    _cov[1][8] = c1[1][0]*y(2, 0) + c1[1][1]*y(2, 1) + c1[1][2]*y(2, 2) - _cov[1][11]*_dt;
    _cov[2][6] = c1[2][0]*y(0, 0) + c1[2][1]*y(0, 1) + c1[2][2]*y(0, 2) - _cov[2][9]*_dt;
    _cov[2][7] = c1[2][0]*y(1, 0) + c1[2][1]*y(1, 1) + c1[2][2]*y(1, 2) - _cov[2][10]*_dt;
    _cov[2][8] = c1[2][0]*y(2, 0) + c1[2][1]*y(2, 1) + c1[2][2]*y(2, 2) - _cov[2][11]*_dt;

    // P12_ = (P12 + P22 * dt) - (P13 + P23 * dt) * X' * dt - P15_ * R' * dt + P16_ * ez' * dt      
    //      = P12 + (P22 - ((P13 + P23 * dt) * X' + P15_ * R') + P16_ * ez') * dt
    //      = P12 + (P22 - (c1 * X' + P15_ * R') + P16_ * ez') * dt         (乘法: 63, 加法: 66)
    _cov[0][3] += _cov[3][3] * _dt  - (c1[0][0] * x[0][0] + c1[0][1] * x[0][1] + c1[0][2] * x[0][2] + _cov[0][12] * rot(0, 0) + _cov[0][13] * rot(0, 1) + _cov[0][14] * rot(0, 2));
    _cov[0][4] += _cov[3][4] * _dt  - (c1[0][0] * x[1][0] + c1[0][1] * x[1][1] + c1[0][2] * x[1][2] + _cov[0][12] * rot(1, 0) + _cov[0][13] * rot(1, 1) + _cov[0][14] * rot(1, 2));
    _cov[0][5] += (_cov[3][5] + _cov[0][15]) * _dt - (c1[0][0] * x[2][0] + c1[0][1] * x[2][1] + c1[0][2] * x[2][2] + _cov[0][12] * rot(2, 0) + _cov[0][13] * rot(2, 1) + _cov[0][14] * rot(2, 2));
    _cov[1][3] += _cov[4][3] * _dt - (c1[1][0] * x[0][0] + c1[1][1] * x[0][1] + c1[1][2] * x[0][2] + _cov[1][12] * rot(0, 0) + _cov[1][13] * rot(0, 1) + _cov[1][14] * rot(0, 2));
    _cov[1][4] += _cov[4][4] * _dt - (c1[1][0] * x[1][0] + c1[1][1] * x[1][1] + c1[1][2] * x[1][2] + _cov[1][12] * rot(1, 0) + _cov[1][13] * rot(1, 1) + _cov[1][14] * rot(1, 2));
    _cov[1][5] += (_cov[4][5] + _cov[1][15]) * _dt - (c1[1][0] * x[2][0] + c1[1][1] * x[2][1] + c1[1][2] * x[2][2] + _cov[1][12] * rot(2, 0) + _cov[1][13] * rot(2, 1) + _cov[1][14] * rot(2, 2));
    _cov[2][3] += _cov[5][3] * _dt - (c1[2][0] * x[0][0] + c1[2][1] * x[0][1] + c1[2][2] * x[0][2] + _cov[2][12] * rot(0, 0) + _cov[2][13] * rot(0, 1) + _cov[2][14] * rot(0, 2));
    _cov[2][4] += _cov[5][4] * _dt - (c1[2][0] * x[1][0] + c1[2][1] * x[1][1] + c1[2][2] * x[1][2] + _cov[2][12] * rot(1, 0) + _cov[2][13] * rot(1, 1) + _cov[2][14] * rot(1, 2));
    _cov[2][5] += (_cov[5][5] + _cov[2][15]) * _dt - (c1[2][0] * x[2][0] + c1[2][1] * x[2][1] + c1[2][2] * x[2][2] + _cov[2][12] * rot(2, 0) + _cov[2][13] * rot(2, 1) + _cov[2][14] * rot(2, 2));

    // X * P32              (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> XP32 = {{{_cov[3][6]*x[0][0] + _cov[3][7]*x[0][1] + _cov[3][8]*x[0][2],
                                              _cov[4][6]*x[0][0] + _cov[4][7]*x[0][1] + _cov[4][8]*x[0][2],
                                              _cov[5][6]*x[0][0] + _cov[5][7]*x[0][1] + _cov[5][8]*x[0][2]},
                                             {_cov[3][6]*x[1][0] + _cov[3][7]*x[1][1] + _cov[3][8]*x[1][2],
                                              _cov[4][6]*x[1][0] + _cov[4][7]*x[1][1] + _cov[4][8]*x[1][2],
                                              _cov[5][6]*x[1][0] + _cov[5][7]*x[1][1] + _cov[5][8]*x[1][2]},
                                             {_cov[3][6]*x[2][0] + _cov[3][7]*x[2][1] + _cov[3][8]*x[2][2],
                                              _cov[4][6]*x[2][0] + _cov[4][7]*x[2][1] + _cov[4][8]*x[2][2],
                                              _cov[5][6]*x[2][0] + _cov[5][7]*x[2][1] + _cov[5][8]*x[2][2]}}};

    // R * P52              (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> RP52 = {{{_cov[3][12]*rot(0, 0) + _cov[3][13]*rot(0, 1) + _cov[3][14]*rot(0, 2),
                                              _cov[4][12]*rot(0, 0) + _cov[4][13]*rot(0, 1) + _cov[4][14]*rot(0, 2),
                                              _cov[5][12]*rot(0, 0) + _cov[5][13]*rot(0, 1) + _cov[5][14]*rot(0, 2)},
                                             {_cov[3][12]*rot(1, 0) + _cov[3][13]*rot(1, 1) + _cov[3][14]*rot(1, 2),
                                              _cov[4][12]*rot(1, 0) + _cov[4][13]*rot(1, 1) + _cov[4][14]*rot(1, 2),
                                              _cov[5][12]*rot(1, 0) + _cov[5][13]*rot(1, 1) + _cov[5][14]*rot(1, 2)},
                                             {_cov[3][12]*rot(2, 0) + _cov[3][13]*rot(2, 1) + _cov[3][14]*rot(2, 2),
                                              _cov[4][12]*rot(2, 0) + _cov[4][13]*rot(2, 1) + _cov[4][14]*rot(2, 2),
                                              _cov[5][12]*rot(2, 0) + _cov[5][13]*rot(2, 1) + _cov[5][14]*rot(2, 2)}}};

    // X * P33              (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> XP33 = {{{_cov[6][6]*x[0][0] + _cov[6][7]*x[0][1] + _cov[6][8]*x[0][2],
                                              _cov[6][7]*x[0][0] + _cov[7][7]*x[0][1] + _cov[7][8]*x[0][2],
                                              _cov[6][8]*x[0][0] + _cov[7][8]*x[0][1] + _cov[8][8]*x[0][2]},
                                             {_cov[6][6]*x[1][0] + _cov[6][7]*x[1][1] + _cov[6][8]*x[1][2],
                                              _cov[6][7]*x[1][0] + _cov[7][7]*x[1][1] + _cov[7][8]*x[1][2],
                                              _cov[6][8]*x[1][0] + _cov[7][8]*x[1][1] + _cov[8][8]*x[1][2]},
                                             {_cov[6][6]*x[2][0] + _cov[6][7]*x[2][1] + _cov[6][8]*x[2][2],
                                              _cov[6][7]*x[2][0] + _cov[7][7]*x[2][1] + _cov[7][8]*x[2][2],
                                              _cov[6][8]*x[2][0] + _cov[7][8]*x[2][1] + _cov[8][8]*x[2][2]}}};

    // R * P53              (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> RP53 = {{{_cov[6][12]*rot(0, 0) + _cov[6][13]*rot(0, 1) + _cov[6][14]*rot(0, 2),
                                              _cov[7][12]*rot(0, 0) + _cov[7][13]*rot(0, 1) + _cov[7][14]*rot(0, 2),
                                              _cov[8][12]*rot(0, 0) + _cov[8][13]*rot(0, 1) + _cov[8][14]*rot(0, 2)},
                                             {_cov[6][12]*rot(1, 0) + _cov[6][13]*rot(1, 1) + _cov[6][14]*rot(1, 2),
                                              _cov[7][12]*rot(1, 0) + _cov[7][13]*rot(1, 1) + _cov[7][14]*rot(1, 2),
                                              _cov[8][12]*rot(1, 0) + _cov[8][13]*rot(1, 1) + _cov[8][14]*rot(1, 2)},
                                             {_cov[6][12]*rot(2, 0) + _cov[6][13]*rot(2, 1) + _cov[6][14]*rot(2, 2),
                                              _cov[7][12]*rot(2, 0) + _cov[7][13]*rot(2, 1) + _cov[7][14]*rot(2, 2),
                                              _cov[8][12]*rot(2, 0) + _cov[8][13]*rot(2, 1) + _cov[8][14]*rot(2, 2)}}};

    // X * P34              (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> XP34 = {{{_cov[6][9]*x[0][0] + _cov[7][9]*x[0][1] + _cov[8][9]*x[0][2],
                                              _cov[6][10]*x[0][0] + _cov[7][10]*x[0][1] + _cov[8][10]*x[0][2],
                                              _cov[6][11]*x[0][0] + _cov[7][11]*x[0][1] + _cov[8][11]*x[0][2]},
                                             {_cov[6][9]*x[1][0] + _cov[7][9]*x[1][1] + _cov[8][9]*x[1][2],
                                              _cov[6][10]*x[1][0] + _cov[7][10]*x[1][1] + _cov[8][10]*x[1][2],
                                              _cov[6][11]*x[1][0] + _cov[7][11]*x[1][1] + _cov[8][11]*x[1][2]},
                                             {_cov[6][9]*x[2][0] + _cov[7][9]*x[2][1] + _cov[8][9]*x[2][2],
                                              _cov[6][10]*x[2][0] + _cov[7][10]*x[2][1] + _cov[8][10]*x[2][2],
                                              _cov[6][11]*x[2][0] + _cov[7][11]*x[2][1] + _cov[8][11]*x[2][2]}}};

    // R * P54              (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> RP54 = {{{_cov[9][12]*rot(0, 0) + _cov[9][13]*rot(0, 1) + _cov[9][14]*rot(0, 2),
                                              _cov[10][12]*rot(0, 0) + _cov[10][13]*rot(0, 1) + _cov[10][14]*rot(0, 2),
                                              _cov[11][12]*rot(0, 0) + _cov[11][13]*rot(0, 1) + _cov[11][14]*rot(0, 2)},
                                             {_cov[9][12]*rot(1, 0) + _cov[9][13]*rot(1, 1) + _cov[9][14]*rot(1, 2),
                                              _cov[10][12]*rot(1, 0) + _cov[10][13]*rot(1, 1) + _cov[10][14]*rot(1, 2),
                                              _cov[11][12]*rot(1, 0) + _cov[11][13]*rot(1, 1) + _cov[11][14]*rot(1, 2)},
                                             {_cov[9][12]*rot(2, 0) + _cov[9][13]*rot(2, 1) + _cov[9][14]*rot(2, 2),
                                              _cov[10][12]*rot(2, 0) + _cov[10][13]*rot(2, 1) + _cov[10][14]*rot(2, 2),
                                              _cov[11][12]*rot(2, 0) + _cov[11][13]*rot(2, 1) + _cov[11][14]*rot(2, 2)}}};

    // X * P35              (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> XP35 = {{{_cov[6][12]*x[0][0] + _cov[7][12]*x[0][1] + _cov[8][12]*x[0][2],
                                              _cov[6][13]*x[0][0] + _cov[7][13]*x[0][1] + _cov[8][13]*x[0][2],
                                              _cov[6][14]*x[0][0] + _cov[7][14]*x[0][1] + _cov[8][14]*x[0][2]},
                                             {_cov[6][12]*x[1][0] + _cov[7][12]*x[1][1] + _cov[8][12]*x[1][2],
                                              _cov[6][13]*x[1][0] + _cov[7][13]*x[1][1] + _cov[8][13]*x[1][2],
                                              _cov[6][14]*x[1][0] + _cov[7][14]*x[1][1] + _cov[8][14]*x[1][2]},
                                             {_cov[6][12]*x[2][0] + _cov[7][12]*x[2][1] + _cov[8][12]*x[2][2],
                                              _cov[6][13]*x[2][0] + _cov[7][13]*x[2][1] + _cov[8][13]*x[2][2],
                                              _cov[6][14]*x[2][0] + _cov[7][14]*x[2][1] + _cov[8][14]*x[2][2]}}};

    // R * P55              (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> RP55 = {{{_cov[12][12]*rot(0, 0) + _cov[12][13]*rot(0, 1) + _cov[12][14]*rot(0, 2),
                                              _cov[12][13]*rot(0, 0) + _cov[13][13]*rot(0, 1) + _cov[13][14]*rot(0, 2),
                                              _cov[12][14]*rot(0, 0) + _cov[13][14]*rot(0, 1) + _cov[14][14]*rot(0, 2)},
                                             {_cov[12][12]*rot(1, 0) + _cov[12][13]*rot(1, 1) + _cov[12][14]*rot(1, 2),
                                              _cov[12][13]*rot(1, 0) + _cov[13][13]*rot(1, 1) + _cov[13][14]*rot(1, 2),
                                              _cov[12][14]*rot(1, 0) + _cov[13][14]*rot(1, 1) + _cov[14][14]*rot(1, 2)},
                                             {_cov[12][12]*rot(2, 0) + _cov[12][13]*rot(2, 1) + _cov[12][14]*rot(2, 2),
                                              _cov[12][13]*rot(2, 0) + _cov[13][13]*rot(2, 1) + _cov[13][14]*rot(2, 2),
                                              _cov[12][14]*rot(2, 0) + _cov[13][14]*rot(2, 1) + _cov[14][14]*rot(2, 2)}}};

    // X * P36              (乘法: 9, 加法: 6)
    const array<float, 3> XP36 = {_cov[6][15]*x[0][0] + _cov[7][15]*x[0][1] + _cov[8][15]*x[0][2],
                                  _cov[6][15]*x[1][0] + _cov[7][15]*x[1][1] + _cov[8][15]*x[1][2],
                                  _cov[6][15]*x[2][0] + _cov[7][15]*x[2][1] + _cov[8][15]*x[2][2]};

    // R * P56              (乘法: 9, 加法: 6)
    const array<float, 3> RP56 = {_cov[12][15]*rot(0, 0) + _cov[13][15]*rot(0, 1) + _cov[14][15]*rot(0, 2),
                                  _cov[12][15]*rot(1, 0) + _cov[13][15]*rot(1, 1) + _cov[14][15]*rot(1, 2),
                                  _cov[12][15]*rot(2, 0) + _cov[13][15]*rot(2, 1) + _cov[14][15]*rot(2, 2)};

    // R * P53 * X'         (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> RP53XT = {{{RP53[0][0]*x[0][0] + RP53[0][1]*x[0][1] + RP53[0][2]*x[0][2],
                                                RP53[0][0]*x[1][0] + RP53[0][1]*x[1][1] + RP53[0][2]*x[1][2],
                                                RP53[0][0]*x[2][0] + RP53[0][1]*x[2][1] + RP53[0][2]*x[2][2]},
                                               {RP53[1][0]*x[0][0] + RP53[1][1]*x[0][1] + RP53[1][2]*x[0][2],
                                                RP53[1][0]*x[1][0] + RP53[1][1]*x[1][1] + RP53[1][2]*x[1][2],
                                                RP53[1][0]*x[2][0] + RP53[1][1]*x[2][1] + RP53[1][2]*x[2][2]},
                                               {RP53[2][0]*x[0][0] + RP53[2][1]*x[0][1] + RP53[2][2]*x[0][2],
                                                RP53[2][0]*x[1][0] + RP53[2][1]*x[1][1] + RP53[2][2]*x[1][2],
                                                RP53[2][0]*x[2][0] + RP53[2][1]*x[2][1] + RP53[2][2]*x[2][2]}}};

    // X * P33 * X'         (乘法: 18, 加法: 12)
    const float XP33XT00 = XP33[0][0]*x[0][0] + XP33[0][1]*x[0][1] + XP33[0][2]*x[0][2];
    const float XP33XT01 = XP33[0][0]*x[1][0] + XP33[0][1]*x[1][1] + XP33[0][2]*x[1][2];
    const float XP33XT02 = XP33[0][0]*x[2][0] + XP33[0][1]*x[2][1] + XP33[0][2]*x[2][2];
    const float XP33XT11 = XP33[1][0]*x[1][0] + XP33[1][1]*x[1][1] + XP33[1][2]*x[1][2];
    const float XP33XT12 = XP33[1][0]*x[2][0] + XP33[1][1]*x[2][1] + XP33[1][2]*x[2][2];
    const float XP33XT22 = XP33[2][0]*x[2][0] + XP33[2][1]*x[2][1] + XP33[2][2]*x[2][2];

    // R * P55 * R'         (乘法: 18, 加法: 12)
    const float RP55RT00 = RP55[0][0]*rot(0, 0) + RP55[0][1]*rot(0, 1) + RP55[0][2]*rot(0, 2);
    const float RP55RT01 = RP55[0][0]*rot(1, 0) + RP55[0][1]*rot(1, 1) + RP55[0][2]*rot(1, 2);
    const float RP55RT02 = RP55[0][0]*rot(2, 0) + RP55[0][1]*rot(2, 1) + RP55[0][2]*rot(2, 2);
    const float RP55RT11 = RP55[1][0]*rot(1, 0) + RP55[1][1]*rot(1, 1) + RP55[1][2]*rot(1, 2);
    const float RP55RT12 = RP55[1][0]*rot(2, 0) + RP55[1][1]*rot(2, 1) + RP55[1][2]*rot(2, 2);
    const float RP55RT22 = RP55[2][0]*rot(2, 0) + RP55[2][1]*rot(2, 1) + RP55[2][2]*rot(2, 2);

    // P22_ = P22 - (X * P32 + P23 * X') * dt - (R * P52 + P25 * R') * dt + (ez * P62 + P26 * ez') * dt
    //        + X * P33 * X' * dt^2 + (R * P53 * X' + X * P35 * R') * dt^2 + R * P55 * R' * dt^2
    //        - (ez * P63 * X' + X * P36 * ez') * dt^2 - (ez * P65 * R' + R * P56 * ez') * dt^2 
    //        + ez * P66 * ez' * dt^2          (乘法: 6, 加法: 55)
    _cov[3][3] += (XP33XT00 + RP53XT[0][0] + RP53XT[0][0] + RP55RT00) - (XP32[0][0] + XP32[0][0] + RP52[0][0] + RP52[0][0]);
    _cov[3][4] += (XP33XT01 + RP53XT[0][1] + RP53XT[1][0] + RP55RT01) - (XP32[0][1] + XP32[1][0] + RP52[0][1] + RP52[1][0]);
    _cov[3][5] += (XP33XT02 + RP53XT[0][2] + RP53XT[2][0] + RP55RT02 - (XP36[0] + RP56[0]) * _dt) + (_cov[15][3] * _dt - XP32[0][2] - XP32[2][0] - RP52[0][2] - RP52[2][0]);
    _cov[4][4] += (XP33XT11 + RP53XT[1][1] + RP53XT[1][1] + RP55RT11) - (XP32[1][1] + XP32[1][1] + RP52[1][1] + RP52[1][1]);
    _cov[4][5] += (XP33XT12 + RP53XT[1][2] + RP53XT[2][1] + RP55RT12 - (XP36[1] + RP56[1]) * _dt) + (_cov[15][4] * _dt - XP32[1][2] - XP32[2][1] - RP52[1][2] - RP52[2][1]);
    _cov[5][5] += (XP33XT22 + RP53XT[2][2] + RP53XT[2][2] + RP55RT22 + _cov[15][15] * _dt2 - (XP36[2] + XP36[2] + RP56[2] + RP56[2]) * _dt) + ((_cov[15][5] + _cov[15][5]) * _dt - XP32[2][2] - XP32[2][2] - RP52[2][2] - RP52[2][2]);

    // P26_ = P26 - XP36 * dt - R * P56 * dt + ez * P66 * dt        (乘法: 1, 加法: 7)
    _cov[3][15] -= (XP36[0] + RP56[0]);
    _cov[4][15] -= (XP36[1] + RP56[1]);
    _cov[5][15] += (_cov[15][15] * _dt - (XP36[2] + RP56[2]));

    // P25_ = P25 - XP35 * dt - R * P55 * dt + ez * P65 * dt        (乘法: 3, 加法: 21)
    _cov[3][12] -= (XP35[0][0] + RP55[0][0]);
    _cov[3][13] -= (XP35[0][1] + RP55[0][1]);
    _cov[3][14] -= (XP35[0][2] + RP55[0][2]);
    _cov[4][12] -= (XP35[1][0] + RP55[1][0]);
    _cov[4][13] -= (XP35[1][1] + RP55[1][1]);
    _cov[4][14] -= (XP35[1][2] + RP55[1][2]);
    _cov[5][12] += (_cov[12][15] * _dt - (XP35[2][0] + RP55[2][0]));
    _cov[5][13] += (_cov[13][15] * _dt - (XP35[2][1] + RP55[2][1]));
    _cov[5][14] += (_cov[14][15] * _dt - (XP35[2][2] + RP55[2][2]));

    // P24_ = P24 - XP34 * dt - R * P54 * dt + ez * P64 * dt        (乘法: 3, 加法: 21)
    _cov[3][9] -= (XP34[0][0] + RP54[0][0]);
    _cov[3][10] -= (XP34[0][1] + RP54[0][1]);
    _cov[3][11] -= (XP34[0][2] + RP54[0][2]);
    _cov[4][9] -= (XP34[1][0] + RP54[1][0]);
    _cov[4][10] -= (XP34[1][1] + RP54[1][1]);
    _cov[4][11] -= (XP34[1][2] + RP54[1][2]);
    _cov[5][9] += (_cov[9][15] * _dt - (XP34[2][0] + RP54[2][0]));
    _cov[5][10] += (_cov[10][15] * _dt - (XP34[2][1] + RP54[2][1]));
    _cov[5][11] += (_cov[11][15] * _dt - (XP34[2][2] + RP54[2][2]));

    // (P23 - X * P33 * dt - R * P53 * dt + ez * P63 * dt)          (乘法: 3, 加法: 21)
    const array<array<float, 3>, 3> c2 = {{{_cov[3][6] - (XP33[0][0] + RP53[0][0]),
                                            _cov[3][7] - (XP33[0][1] + RP53[0][1]),
                                            _cov[3][8] - (XP33[0][2] + RP53[0][2])},
                                           {_cov[4][6] - (XP33[1][0] + RP53[1][0]),
                                            _cov[4][7] - (XP33[1][1] + RP53[1][1]),
                                            _cov[4][8] - (XP33[1][2] + RP53[1][2])},
                                           {_cov[5][6] + (_cov[6][15] * _dt - XP33[2][0] - RP53[2][0]),
                                            _cov[5][7] + (_cov[7][15] * _dt - XP33[2][1] - RP53[2][1]),
                                            _cov[5][8] + (_cov[8][15] * _dt - XP33[2][2] - RP53[2][2])}}};

    // P23_ = (P23 - X * P33 * dt - R * P53 * dt + ez * P63 * dt) * Y' - P24_ * dt          (乘法: 36, 加法: 27)
    _cov[3][6] = c2[0][0]*y(0, 0) + c2[0][1]*y(0, 1) + c2[0][2]*y(0, 2) - _cov[3][9] * _dt;
    _cov[3][7] = c2[0][0]*y(1, 0) + c2[0][1]*y(1, 1) + c2[0][2]*y(1, 2) - _cov[3][10] * _dt;
    _cov[3][8] = c2[0][0]*y(2, 0) + c2[0][1]*y(2, 1) + c2[0][2]*y(2, 2) - _cov[3][11] * _dt;
    _cov[4][6] = c2[1][0]*y(0, 0) + c2[1][1]*y(0, 1) + c2[1][2]*y(0, 2) - _cov[4][9] * _dt;
    _cov[4][7] = c2[1][0]*y(1, 0) + c2[1][1]*y(1, 1) + c2[1][2]*y(1, 2) - _cov[4][10] * _dt;
    _cov[4][8] = c2[1][0]*y(2, 0) + c2[1][1]*y(2, 1) + c2[1][2]*y(2, 2) - _cov[4][11] * _dt;
    _cov[5][6] = c2[2][0]*y(0, 0) + c2[2][1]*y(0, 1) + c2[2][2]*y(0, 2) - _cov[5][9] * _dt;
    _cov[5][7] = c2[2][0]*y(1, 0) + c2[2][1]*y(1, 1) + c2[2][2]*y(1, 2) - _cov[5][10] * _dt;
    _cov[5][8] = c2[2][0]*y(2, 0) + c2[2][1]*y(2, 1) + c2[2][2]*y(2, 2) - _cov[5][11] * _dt;

    // Y * P36              (乘法: 9, 加法: 6)
    const array<float, 3> YP36 = {y(0, 0) * _cov[6][15] + y(0, 1) * _cov[7][15] + y(0, 2) * _cov[8][15],
                                  y(1, 0) * _cov[6][15] + y(1, 1) * _cov[7][15] + y(1, 2) * _cov[8][15],
                                  y(2, 0) * _cov[6][15] + y(2, 1) * _cov[7][15] + y(2, 2) * _cov[8][15]};

    // Y * P35              (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> YP35 = {{{_cov[6][12]*y(0, 0) + _cov[7][12]*y(0, 1) + _cov[8][12]*y(0, 2),
                                              _cov[6][13]*y(0, 0) + _cov[7][13]*y(0, 1) + _cov[8][13]*y(0, 2),
                                              _cov[6][14]*y(0, 0) + _cov[7][14]*y(0, 1) + _cov[8][14]*y(0, 2)},
                                             {_cov[6][12]*y(1, 0) + _cov[7][12]*y(1, 1) + _cov[8][12]*y(1, 2),
                                              _cov[6][13]*y(1, 0) + _cov[7][13]*y(1, 1) + _cov[8][13]*y(1, 2),
                                              _cov[6][14]*y(1, 0) + _cov[7][14]*y(1, 1) + _cov[8][14]*y(1, 2)},
                                             {_cov[6][12]*y(2, 0) + _cov[7][12]*y(2, 1) + _cov[8][12]*y(2, 2),
                                              _cov[6][13]*y(2, 0) + _cov[7][13]*y(2, 1) + _cov[8][13]*y(2, 2),
                                              _cov[6][14]*y(2, 0) + _cov[7][14]*y(2, 1) + _cov[8][14]*y(2, 2)}}};

    // Y * P34              (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> YP34 = {{{_cov[6][9]*y(0, 0) + _cov[7][9]*y(0, 1) + _cov[8][9]*y(0, 2),
                                              _cov[6][10]*y(0, 0) + _cov[7][10]*y(0, 1) + _cov[8][10]*y(0, 2),
                                              _cov[6][11]*y(0, 0) + _cov[7][11]*y(0, 1) + _cov[8][11]*y(0, 2)},
                                             {_cov[6][9]*y(1, 0) + _cov[7][9]*y(1, 1) + _cov[8][9]*y(1, 2),
                                              _cov[6][10]*y(1, 0) + _cov[7][10]*y(1, 1) + _cov[8][10]*y(1, 2),
                                              _cov[6][11]*y(1, 0) + _cov[7][11]*y(1, 1) + _cov[8][11]*y(1, 2)},
                                             {_cov[6][9]*y(2, 0) + _cov[7][9]*y(2, 1) + _cov[8][9]*y(2, 2),
                                              _cov[6][10]*y(2, 0) + _cov[7][10]*y(2, 1) + _cov[8][10]*y(2, 2),
                                              _cov[6][11]*y(2, 0) + _cov[7][11]*y(2, 1) + _cov[8][11]*y(2, 2)}}};

    // Y * P33              (乘法: 27, 加法: 18)
    const array<array<float, 3>, 3> YP33 = {{{_cov[6][6]*y(0, 0) + _cov[6][7]*y(0, 1) + _cov[6][8]*y(0, 2),
                                              _cov[6][7]*y(0, 0) + _cov[7][7]*y(0, 1) + _cov[7][8]*y(0, 2),
                                              _cov[6][8]*y(0, 0) + _cov[7][8]*y(0, 1) + _cov[8][8]*y(0, 2)},
                                             {_cov[6][6]*y(1, 0) + _cov[6][7]*y(1, 1) + _cov[6][8]*y(1, 2),
                                              _cov[6][7]*y(1, 0) + _cov[7][7]*y(1, 1) + _cov[7][8]*y(1, 2),
                                              _cov[6][8]*y(1, 0) + _cov[7][8]*y(1, 1) + _cov[8][8]*y(1, 2)},
                                             {_cov[6][6]*y(2, 0) + _cov[6][7]*y(2, 1) + _cov[6][8]*y(2, 2),
                                              _cov[6][7]*y(2, 0) + _cov[7][7]*y(2, 1) + _cov[7][8]*y(2, 2),
                                              _cov[6][8]*y(2, 0) + _cov[7][8]*y(2, 1) + _cov[8][8]*y(2, 2)}}};

    // Y * P33 * Y'         (乘法: 18, 加法: 12)
    const float YP33YT00 = YP33[0][0]*y(0, 0) + YP33[0][1]*y(0, 1) + YP33[0][2]*y(0, 2);
    const float YP33YT01 = YP33[0][0]*y(1, 0) + YP33[0][1]*y(1, 1) + YP33[0][2]*y(1, 2);
    const float YP33YT02 = YP33[0][0]*y(2, 0) + YP33[0][1]*y(2, 1) + YP33[0][2]*y(2, 2);
    const float YP33YT11 = YP33[1][0]*y(1, 0) + YP33[1][1]*y(1, 1) + YP33[1][2]*y(1, 2);
    const float YP33YT12 = YP33[1][0]*y(2, 0) + YP33[1][1]*y(2, 1) + YP33[1][2]*y(2, 2);
    const float YP33YT22 = YP33[2][0]*y(2, 0) + YP33[2][1]*y(2, 1) + YP33[2][2]*y(2, 2);

    // P33_ = Y * P33 * Y' - (P43 * Y' + Y * P34) * dt + P44 * dt^2         (乘法: 12, 加法: 18)
    _cov[6][6] = YP33YT00 - (YP34[0][0] + YP34[0][0]) * _dt + _cov[9][9] * _dt2;
    _cov[6][7] = YP33YT01 - (YP34[0][1] + YP34[1][0]) * _dt + _cov[9][10] * _dt2;
    _cov[6][8] = YP33YT02 - (YP34[0][2] + YP34[2][0]) * _dt + _cov[9][11] * _dt2;
    _cov[7][7] = YP33YT11 - (YP34[1][1] + YP34[1][1]) * _dt + _cov[10][10] * _dt2;
    _cov[7][8] = YP33YT12 - (YP34[1][2] + YP34[2][1]) * _dt + _cov[10][11] * _dt2;
    _cov[8][8] = YP33YT22 - (YP34[2][2] + YP34[2][2]) * _dt + _cov[11][11] * _dt2;

    // P36_ = Y * P36 - P46 * dt        (乘法: 3, 加法: 3)
    _cov[6][15] = YP36[0] - _cov[9][15] * _dt;
    _cov[7][15] = YP36[1] - _cov[10][15] * _dt;
    _cov[8][15] = YP36[2] - _cov[11][15] * _dt;

    // P35_ = Y * P35 - P45 * dt        (乘法: 9, 加法: 9)
    _cov[6][12] = YP35[0][0] - _cov[9][12] * _dt;
    _cov[6][13] = YP35[0][1] - _cov[9][13] * _dt;
    _cov[6][14] = YP35[0][2] - _cov[9][14] * _dt;
    _cov[7][12] = YP35[1][0] - _cov[10][12] * _dt;
    _cov[7][13] = YP35[1][1] - _cov[10][13] * _dt;
    _cov[7][14] = YP35[1][2] - _cov[10][14] * _dt;
    _cov[8][12] = YP35[2][0] - _cov[11][12] * _dt;
    _cov[8][13] = YP35[2][1] - _cov[11][13] * _dt;
    _cov[8][14] = YP35[2][2] - _cov[11][14] * _dt;

    // P34_ = Y * P34 - P44 * dt        (乘法: 9, 加法: 9)
    _cov[6][9] = YP34[0][0] - _cov[9][9] * _dt;
    _cov[6][10] = YP34[0][1] - _cov[9][10] * _dt;
    _cov[6][11] = YP34[0][2] - _cov[9][11] * _dt;
    _cov[7][9] = YP34[1][0] - _cov[9][10] * _dt;
    _cov[7][10] = YP34[1][1] - _cov[10][10] * _dt;
    _cov[7][11] = YP34[1][2] - _cov[10][11] * _dt;
    _cov[8][9] = YP34[2][0] - _cov[9][11] * _dt;
    _cov[8][10] = YP34[2][1] - _cov[10][11] * _dt;
    _cov[8][11] = YP34[2][2] - _cov[11][11] * _dt;

    // F *P *F' + Q
    _cov[0][0] = kahan_summation(_cov[0][0], _q_cov[0] * _dt2, _accumulator_cov[0]);
    _cov[1][1] = kahan_summation(_cov[1][1], _q_cov[1] * _dt2, _accumulator_cov[1]);
    _cov[2][2] = kahan_summation(_cov[2][2], _q_cov[2] * _dt2, _accumulator_cov[2]);
    if (_q_cov[3] == _q_cov[4]) {
        if (_q_cov[4] == _q_cov[5]) {
            _cov[3][3] = kahan_summation(_cov[3][3], _q_cov[3] * _dt2, _accumulator_cov[3]);
            _cov[4][4] = kahan_summation(_cov[4][4], _q_cov[4] * _dt2, _accumulator_cov[4]);
            _cov[5][5] = kahan_summation(_cov[5][5], _q_cov[5] * _dt2, _accumulator_cov[5]);
        } else {
            float sx = _q_cov[3] * _dt2;
            float sz = _q_cov[5] * _dt2;
            float r00_sx = _rot(0, 0) * sx;
            float r01_sx = _rot(0, 1) * sx;
            float r02_sz = _rot(0, 2) * sz;
            float r00_sx_r10_plus_r01_sx_r11_plus_r02_sz_r12 = r00_sx * _rot(1, 0) + r01_sx * _rot(1, 1) + r02_sz * _rot(1, 2);
            float r00_sx_r20_plus_r01_sx_r21_plus_r02_sz_r22 = r00_sx * _rot(2, 0) + r01_sx * _rot(2, 1) + r02_sz * _rot(2, 2);
            float r10_r20_sx_plus_r11_r21_sx_plus_r12_r22_sz = (_rot(1, 0) * _rot(2, 0) + _rot(1, 1) * _rot(2, 1)) * sx + _rot(1, 2) * _rot(2, 2) * sz;

            _cov[3][3] += kahan_summation(_cov[3][3], (_rot(0, 0) * _rot(0, 0) + _rot(0, 1) * _rot(0, 1)) * sx + _rot(0, 2) * _rot(0, 2) * sz, _accumulator_cov[3]);
            _cov[3][4] += r00_sx_r10_plus_r01_sx_r11_plus_r02_sz_r12;
            _cov[3][5] += r00_sx_r20_plus_r01_sx_r21_plus_r02_sz_r22;
            _cov[4][4] += kahan_summation(_cov[4][4], (_rot(1, 0) * _rot(1, 0) + _rot(1, 1) * _rot(1, 1)) * sx + _rot(1, 2) * _rot(1, 2) * sz, _accumulator_cov[4]);
            _cov[4][5] += r10_r20_sx_plus_r11_r21_sx_plus_r12_r22_sz;
            _cov[5][5] += kahan_summation(_cov[5][5], (_rot(2, 0) * _rot(2, 0) + _rot(2, 1) * _rot(2, 1)) * sx + _rot(2, 2) * _rot(2, 2) * sz, _accumulator_cov[5]);
        }
    } else {
        float sx = _q_cov[3] * _dt2;
        float sy = _q_cov[4] * _dt2;
        float sz = _q_cov[5] * _dt2;
        float r00_sx = _rot(0, 0) * sx;
        float r01_sy = _rot(0, 1) * sy;
        float r02_sz = _rot(0, 2) * sz;
        float r00_sx_r10_plus_r01_sy_r11_plus_r02_sz_r12 = r00_sx * _rot(1, 0) + r01_sy * _rot(1, 1) + r02_sz * _rot(1, 2);
        float r00_sx_r20_plus_r01_sy_r21_plus_r02_sz_r22 = r00_sx * _rot(2, 0) + r01_sy * _rot(2, 1) + r02_sz * _rot(2, 2);
        float r10_r20_sx_plus_r11_r21_sy_plus_r12_r22_sz = _rot(1, 0) * _rot(2, 0) * sx + _rot(1, 1) * _rot(2, 1) * sy + _rot(1, 2) * _rot(2, 2) * sz;
        
        _cov[3][3] += kahan_summation(_cov[3][3], _rot(0, 0) * _rot(0, 0) * sx + _rot(0, 1) * _rot(0, 1) * sy + _rot(0, 2) * _rot(0, 2) * sz, _accumulator_cov[3]);
        _cov[3][4] += r00_sx_r10_plus_r01_sy_r11_plus_r02_sz_r12;
        _cov[3][5] += r00_sx_r20_plus_r01_sy_r21_plus_r02_sz_r22;
        _cov[4][4] += kahan_summation(_cov[4][4], _rot(1, 0) * _rot(1, 0) * sx + _rot(1, 1) * _rot(1, 1) * sy + _rot(1, 2) * _rot(1, 2) * sz, _accumulator_cov[4]);
        _cov[4][5] += r10_r20_sx_plus_r11_r21_sy_plus_r12_r22_sz;
        _cov[5][5] += kahan_summation(_cov[5][5], _rot(2, 0) * _rot(2, 0) * sx + _rot(2, 1) * _rot(2, 1) * sy + _rot(2, 2) * _rot(2, 2) * sz, _accumulator_cov[5]);
    }
    _cov[6][6] = kahan_summation(_cov[6][6], _q_cov[6] * _dt2, _accumulator_cov[6]);
    _cov[7][7] = kahan_summation(_cov[7][7], _q_cov[7] * _dt2, _accumulator_cov[7]);
    _cov[8][8] = kahan_summation(_cov[8][8], _q_cov[8] * _dt2, _accumulator_cov[8]);
    _cov[9][9] = kahan_summation(_cov[9][9], _q_cov[9] * _dt2, _accumulator_cov[9]);
    _cov[10][10] = kahan_summation(_cov[10][10], _q_cov[10] * _dt2, _accumulator_cov[10]);
    _cov[11][11] = kahan_summation(_cov[11][11], _q_cov[11] * _dt2, _accumulator_cov[11]);
    _cov[12][12] = kahan_summation(_cov[12][12], _q_cov[12] * _dt2, _accumulator_cov[12]);
    _cov[13][13] = kahan_summation(_cov[13][13], _q_cov[13] * _dt2, _accumulator_cov[13]);
    _cov[14][14] = kahan_summation(_cov[14][14], _q_cov[14] * _dt2, _accumulator_cov[14]);
    _cov[15][15] = kahan_summation(_cov[15][15], _q_cov[15] * _dt2, _accumulator_cov[15]);

    regular_covariance_to_symmetric();
}

unsigned char ESKF::fuse_position(const Vector3f &pos, const Vector3f &w, const Vector3f &a, 
                                  const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) {
    /*
    s = [p, v, R, bg, ba, g]
    δs = [δp, δv, δθ, δbg, δba, δg]
    s + δs = [p+δp, v+δv, R*Exp(δθ), bg+δbg, ba+δba, g+δg]

    pos = p + R * dis

    δpos / δp = I
    δpos / δθ = -R * dis^

    H = [I, O, -R*dis^, O, O, O]
    */

    unsigned char info = 0;

    // array<float, 3> minus_rot_d_hat;   // -R * dis^
    // float cov_1_dim;    // _cov[1][dim]  or  _cov[dim][1]
    // array<float, 16> HP;    // H * P  or  P * H'
    // float HPHT_plus_R;    // H * P * H' + R
    // array<float, 16> K {};  // K = P * H' * (H * P * H' + R)^-1
    // float obs_error;    // e = y - h = pos - (p + R * dis)

    // Sequential Kalman Filter
    for (unsigned int dim = 0; dim < 3; ++dim) {
        if (!isfinite(noise_std[dim])) {
            continue;
        }

        /*
        K = P * H * (H * P * H' + R)^-1
        x = x + K * (y - h)
        P = (I - K * H) * P

        where, H = [I, O, -R*dis^, O, O, O]
               h = p + R * dis
        */

        // -R * dis^
        const array<float, 3> minus_rot_d_hat = {_rot(dim, 2)*dis[1] - _rot(dim, 1)*dis[2], _rot(dim, 0)*dis[2] - _rot(dim, 2)*dis[0], _rot(dim, 1)*dis[0] - _rot(dim, 0)*dis[1]};

        // H * P  or  P * H'
        const float cov_1_dim = (dim == 2) ? _cov[1][dim] : _cov[dim][1];
        const array<float, 16> HP = {_cov[0][dim] + _cov[0][6]*minus_rot_d_hat[0] + _cov[0][7]*minus_rot_d_hat[1] + _cov[0][8]*minus_rot_d_hat[2],
                                     cov_1_dim + _cov[1][6]*minus_rot_d_hat[0] + _cov[1][7]*minus_rot_d_hat[1] + _cov[1][8]*minus_rot_d_hat[2],
                                     _cov[dim][2] + _cov[2][6]*minus_rot_d_hat[0] + _cov[2][7]*minus_rot_d_hat[1] + _cov[2][8]*minus_rot_d_hat[2],
                                     _cov[dim][3] + _cov[3][6]*minus_rot_d_hat[0] + _cov[3][7]*minus_rot_d_hat[1] + _cov[3][8]*minus_rot_d_hat[2],
                                     _cov[dim][4] + _cov[4][6]*minus_rot_d_hat[0] + _cov[4][7]*minus_rot_d_hat[1] + _cov[4][8]*minus_rot_d_hat[2],
                                     _cov[dim][5] + _cov[5][6]*minus_rot_d_hat[0] + _cov[5][7]*minus_rot_d_hat[1] + _cov[5][8]*minus_rot_d_hat[2],
                                     _cov[dim][6] + _cov[6][6]*minus_rot_d_hat[0] + _cov[6][7]*minus_rot_d_hat[1] + _cov[6][8]*minus_rot_d_hat[2],
                                     _cov[dim][7] + _cov[6][7]*minus_rot_d_hat[0] + _cov[7][7]*minus_rot_d_hat[1] + _cov[7][8]*minus_rot_d_hat[2],
                                     _cov[dim][8] + _cov[6][8]*minus_rot_d_hat[0] + _cov[7][8]*minus_rot_d_hat[1] + _cov[8][8]*minus_rot_d_hat[2],
                                     _cov[dim][9] + _cov[6][9]*minus_rot_d_hat[0] + _cov[7][9]*minus_rot_d_hat[1] + _cov[8][9]*minus_rot_d_hat[2],
                                     _cov[dim][10] + _cov[6][10]*minus_rot_d_hat[0] + _cov[7][10]*minus_rot_d_hat[1] + _cov[8][10]*minus_rot_d_hat[2],
                                     _cov[dim][11] + _cov[6][11]*minus_rot_d_hat[0] + _cov[7][11]*minus_rot_d_hat[1] + _cov[8][11]*minus_rot_d_hat[2],
                                     _cov[dim][12] + _cov[6][12]*minus_rot_d_hat[0] + _cov[7][12]*minus_rot_d_hat[1] + _cov[8][12]*minus_rot_d_hat[2],
                                     _cov[dim][13] + _cov[6][13]*minus_rot_d_hat[0] + _cov[7][13]*minus_rot_d_hat[1] + _cov[8][13]*minus_rot_d_hat[2],
                                     _cov[dim][14] + _cov[6][14]*minus_rot_d_hat[0] + _cov[7][14]*minus_rot_d_hat[1] + _cov[8][14]*minus_rot_d_hat[2],
                                     _cov[dim][15] + _cov[6][15]*minus_rot_d_hat[0] + _cov[7][15]*minus_rot_d_hat[1] + _cov[8][15]*minus_rot_d_hat[2]};

        // H * P * H' + R
        const float HPHT_plus_R = HP[dim] + HP[6] * minus_rot_d_hat[0] + HP[7] * minus_rot_d_hat[1] + HP[8] * minus_rot_d_hat[2] + noise_std[dim] * noise_std[dim];

        // h = p + R * dis
        // e = y - h = pos - (p + R * dis)
        const float obs_error = pos[dim] - (_p[dim] + _rot(dim, 0) * dis[0] + _rot(dim, 1) * dis[1] + _rot(dim, 2) * dis[2]);

        /*
        K = P * H' * (H * P * H' + R)^-1
        P = P - K * H * P

        e = y - h = y - (v + R * (w - bg)^ * dis)
        x = x + K * (y - h)
        */
        switch (conservative_posteriori_estimate(HP, HPHT_plus_R, obs_error, gate[dim])) {
            case 0:
                break;
            case 1:
                info |= (1 << dim);
                break;
            case 2:
                info |= (8 << dim);  
                reset_covariance_matrix(dim, dim + 1, _q_cov);
                break;
        }
    }  
    regular_covariance_to_symmetric();

    return info;
}

unsigned char ESKF::fuse_velocity(const Vector3f &vel, const Vector3f &w, const Vector3f &a, 
                                  const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) {
    /*
    s = [p, v, R, bg, ba, g]
    δs = [δp, δv, δθ, δbg, δba, δg]
    s + δs = [p+δp, v+δv, R*Exp(δθ), bg+δbg, ba+δba, g+δg]

    vel = v + R * (w - bg)^ * dis

    δvel / δv = I
    δvel / δθ = -R * ((w - bg)^ * dis)^
    δvel / δbg = R * dis^

    H = [O, I, -R*((w-bg)^*dis)^, R*dis^, O, O]
    */ 

    unsigned char info = 0;

    // (w-bg)^ * dis
    const array<float, 3> w_corr = {w[0] - _bg[0], w[1] - _bg[1], w[2] - _bg[2]};
    const array<float, 3> w_corr_cross_d = {w_corr[1] * dis[2] - w_corr[2] * dis[1], 
                                            w_corr[2] * dis[0] - w_corr[0] * dis[2], 
                                            w_corr[0] * dis[1] - w_corr[1] * dis[0]};
    
    // Sequential Kalman Filter
    for (unsigned int dim = 0; dim < 3; ++dim) {
        if (!isfinite(noise_std[dim])) {
            continue;
        }

        /*
        K = P * H * (H * P * H' + R)^-1
        x = x + K * (y - h)
        P = (I - K * H) * P

        where, H = [O, I, R*(dis^*(w-bg))^, R*dis^, O, O]
               h = v + R * (w-bg)^ * dis
        */

        // R * dis^
        const array<float, 3> rot_d_hat = {_rot(dim, 1)*dis[2] - _rot(dim, 2)*dis[1], 
                                           _rot(dim, 2)*dis[0] - _rot(dim, 0)*dis[2], 
                                           _rot(dim, 0)*dis[1] - _rot(dim, 1)*dis[0]};

        // R * (dis^*(w-bg))^ = - R * ((w-bg)^*dis)^
        const array<float, 3> rot_d_cross_w_corr_hat = {_rot(dim, 2)*w_corr_cross_d[1] - _rot(dim, 1)*w_corr_cross_d[2], 
                                                        _rot(dim, 0)*w_corr_cross_d[2] - _rot(dim, 2)*w_corr_cross_d[0], 
                                                        _rot(dim, 1)*w_corr_cross_d[0] - _rot(dim, 0)*w_corr_cross_d[1]};

        // H * P  or  P * H'
        const unsigned int index = 3 + dim;
        const float cov_4_index = (dim == 2) ? _cov[4][index] : _cov[index][4];
        const array<float, 16> HP = {_cov[0][index] + _cov[0][6]*rot_d_cross_w_corr_hat[0] + _cov[0][7]*rot_d_cross_w_corr_hat[1] + _cov[0][8]*rot_d_cross_w_corr_hat[2] + _cov[0][9]*rot_d_hat[0] + _cov[0][10]*rot_d_hat[1] + _cov[0][11]*rot_d_hat[2],
                                     _cov[1][index] + _cov[1][6]*rot_d_cross_w_corr_hat[0] + _cov[1][7]*rot_d_cross_w_corr_hat[1] + _cov[1][8]*rot_d_cross_w_corr_hat[2] + _cov[1][9]*rot_d_hat[0] + _cov[1][10]*rot_d_hat[1] + _cov[1][11]*rot_d_hat[2],
                                     _cov[2][index] + _cov[2][6]*rot_d_cross_w_corr_hat[0] + _cov[2][7]*rot_d_cross_w_corr_hat[1] + _cov[2][8]*rot_d_cross_w_corr_hat[2] + _cov[2][9]*rot_d_hat[0] + _cov[2][10]*rot_d_hat[1] + _cov[2][11]*rot_d_hat[2],
                                     _cov[3][index] + _cov[3][6]*rot_d_cross_w_corr_hat[0] + _cov[3][7]*rot_d_cross_w_corr_hat[1] + _cov[3][8]*rot_d_cross_w_corr_hat[2] + _cov[3][9]*rot_d_hat[0] + _cov[3][10]*rot_d_hat[1] + _cov[3][11]*rot_d_hat[2],
                                     cov_4_index + _cov[4][6]*rot_d_cross_w_corr_hat[0] + _cov[4][7]*rot_d_cross_w_corr_hat[1] + _cov[4][8]*rot_d_cross_w_corr_hat[2] + _cov[4][9]*rot_d_hat[0] + _cov[4][10]*rot_d_hat[1] + _cov[4][11]*rot_d_hat[2],
                                     _cov[index][5] + _cov[5][6]*rot_d_cross_w_corr_hat[0] + _cov[5][7]*rot_d_cross_w_corr_hat[1] + _cov[5][8]*rot_d_cross_w_corr_hat[2] + _cov[5][9]*rot_d_hat[0] + _cov[5][10]*rot_d_hat[1] + _cov[5][11]*rot_d_hat[2],
                                     _cov[index][6] + _cov[6][6]*rot_d_cross_w_corr_hat[0] + _cov[6][7]*rot_d_cross_w_corr_hat[1] + _cov[6][8]*rot_d_cross_w_corr_hat[2] + _cov[6][9]*rot_d_hat[0] + _cov[6][10]*rot_d_hat[1] + _cov[6][11]*rot_d_hat[2],
                                     _cov[index][7] + _cov[6][7]*rot_d_cross_w_corr_hat[0] + _cov[7][7]*rot_d_cross_w_corr_hat[1] + _cov[7][8]*rot_d_cross_w_corr_hat[2] + _cov[7][9]*rot_d_hat[0] + _cov[7][10]*rot_d_hat[1] + _cov[7][11]*rot_d_hat[2],
                                     _cov[index][8] + _cov[6][8]*rot_d_cross_w_corr_hat[0] + _cov[7][8]*rot_d_cross_w_corr_hat[1] + _cov[8][8]*rot_d_cross_w_corr_hat[2] + _cov[8][9]*rot_d_hat[0] + _cov[8][10]*rot_d_hat[1] + _cov[8][11]*rot_d_hat[2],
                                     _cov[index][9] + _cov[6][9]*rot_d_cross_w_corr_hat[0] + _cov[7][9]*rot_d_cross_w_corr_hat[1] + _cov[8][9]*rot_d_cross_w_corr_hat[2] + _cov[9][9]*rot_d_hat[0] + _cov[9][10]*rot_d_hat[1] + _cov[9][11]*rot_d_hat[2],
                                     _cov[index][10] + _cov[6][10]*rot_d_cross_w_corr_hat[0] + _cov[7][10]*rot_d_cross_w_corr_hat[1] + _cov[8][10]*rot_d_cross_w_corr_hat[2] + _cov[9][10]*rot_d_hat[0] + _cov[10][10]*rot_d_hat[1] + _cov[10][11]*rot_d_hat[2],
                                     _cov[index][11] + _cov[6][11]*rot_d_cross_w_corr_hat[0] + _cov[7][11]*rot_d_cross_w_corr_hat[1] + _cov[8][11]*rot_d_cross_w_corr_hat[2] + _cov[9][11]*rot_d_hat[0] + _cov[10][11]*rot_d_hat[1] + _cov[11][11]*rot_d_hat[2],
                                     _cov[index][12] + _cov[6][12]*rot_d_cross_w_corr_hat[0] + _cov[7][12]*rot_d_cross_w_corr_hat[1] + _cov[8][12]*rot_d_cross_w_corr_hat[2] + _cov[9][12]*rot_d_hat[0] + _cov[10][12]*rot_d_hat[1] + _cov[11][12]*rot_d_hat[2],
                                     _cov[index][13] + _cov[6][13]*rot_d_cross_w_corr_hat[0] + _cov[7][13]*rot_d_cross_w_corr_hat[1] + _cov[8][13]*rot_d_cross_w_corr_hat[2] + _cov[9][13]*rot_d_hat[0] + _cov[10][13]*rot_d_hat[1] + _cov[11][13]*rot_d_hat[2],
                                     _cov[index][14] + _cov[6][14]*rot_d_cross_w_corr_hat[0] + _cov[7][14]*rot_d_cross_w_corr_hat[1] + _cov[8][14]*rot_d_cross_w_corr_hat[2] + _cov[9][14]*rot_d_hat[0] + _cov[10][14]*rot_d_hat[1] + _cov[11][14]*rot_d_hat[2],
                                     _cov[index][15] + _cov[6][15]*rot_d_cross_w_corr_hat[0] + _cov[7][15]*rot_d_cross_w_corr_hat[1] + _cov[8][15]*rot_d_cross_w_corr_hat[2] + _cov[9][15]*rot_d_hat[0] + _cov[10][15]*rot_d_hat[1] + _cov[11][15]*rot_d_hat[2]};                                                  
        // H * P * H' + R
        const float HPHT_plus_R = HP[index] + HP[6] * rot_d_cross_w_corr_hat[0] + HP[7] * rot_d_cross_w_corr_hat[1] + HP[8] * rot_d_cross_w_corr_hat[2] + HP[9] * rot_d_hat[0] + HP[10] * rot_d_hat[1] + HP[11] * rot_d_hat[2] + noise_std[dim] * noise_std[dim];

        // h = v + R * (w - bg)^ * dis
        // e = y - h = y - (v + R * (w - bg)^ * dis)
        const float obs_error = vel[dim] - (_v[dim] + _rot(dim, 0) * w_corr_cross_d[0] + _rot(dim, 1) * w_corr_cross_d[1] + _rot(dim, 2) * w_corr_cross_d[2]);

        /*
        K = P * H' * (H * P * H' + R)^-1
        P = P - K * H * P

        e = y - h = y - (v + R * (w - bg)^ * dis)
        x = x + K * (y - h)
        */
        switch (conservative_posteriori_estimate(HP, HPHT_plus_R, obs_error, gate[dim])) {
            case 0:
                break;
            case 1:
                info |= (1 << dim);
                break;
            case 2:
                info |= (8 << dim);    
                reset_covariance_matrix(index, index + 1, _q_cov);
                break;
        }
    }
    regular_covariance_to_symmetric();

    return info;
}

unsigned char ESKF::conservative_posteriori_estimate(const array<float, 16> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
    /*
    K = P * H' * (H * P * H' + R)^-1
    P = P - K * H * P

    e = y - h = y - (v + R * (w - bg)^ * dis)
    x = x + K * (y - h)
    */
    
    // Don't correct error state unless obs_error / √(H*P*H' + R) >= gate
    if (obs_error * obs_error  >= gate * gate * HPHT_plus_R) {
        return 1;
    }

    // Don't correct error state unless (I - KH) * P > 0 <=> P - PH'(HPH'+R)^-1HP > 0
    for (unsigned char i = 0; i < 16; ++i) {
        if (HP[i] * HP[i] >= HPHT_plus_R * _cov[i][i]) {
            return 2;    
        }
    }
    
    // K = P * H' * (H * P * H' + R)^-1
    array<float, 16> K {};
    for (unsigned char i = 0; i < 16; ++i) {
        K[i] = HP[i] / HPHT_plus_R;
    }

    // x = x + K * e
    for (unsigned char i = 0; i < 16; ++i) {
        _error_state[i] += K[i] * obs_error;
    }

    // P = P - K * H * P
    for (unsigned char i = 0; i < 16; ++i) {
        for (unsigned char j = i; j < 16; ++j) {
            _cov[i][j] -= K[i] * HP[j];
        }
    }
    return 0;
}

void ESKF::correct_state() {
    // for (float &es : _error_state) {
    //     cout << es << endl;
    // }
    // cout << "-------------------------" << endl;
    // state: [p, v, bg, ba, g], R
    // error_state : [δp, δv, δθ, δbg, δba, δg]

    /*
    p = p + δp
    v = v + δv
    bg = bg + δbg
    ba = ba + δba
    g = g + δg
    */
    for (unsigned char i = 0; i < 3; ++i) {
        _p[i] += _error_state[i];
        _v[i] += _error_state[3 + i];
        _bg[i] += _error_state[9 + i];
        _ba[i] += _error_state[12 + i];
    }
    _g += _error_state[15];
    
    // // R = R * Exp(δθ)
    // Matrix3f delta_rot;
    // array<float, 3> delta_theta = {_error_state[6], _error_state[7], _error_state[8]};
    // rotation_from_axis_angle(delta_rot, delta_theta);
    // Matrix3f rot = _rot;
    // _rot = rot * delta_rot;

    // q = q * Exp(δθ)
    Quaternionf delta_q;
    array<float, 3> delta_theta = {_error_state[6], _error_state[7], _error_state[8]};
    quaternion_from_axis_angle(delta_q, delta_theta);
    const Quaternionf q = _q;
    _q = q * delta_q;

    _rot = q;

    // [δp, δv, δθ, δbg, δba, δg] = 0
    for (float &es : _error_state) {
        es = 0.f;
    }
}

void ESKF::reset_state() {
    _p.setZero();
    _v.setZero();
    _bg.setZero();
    _ba.setZero();
    _g = _g_init;
    _rot.setIdentity();
    _q.setIdentity();
}

void ESKF::reset_error_state() { 
    for (float &es : _error_state) { 
        es = 0.f; 
    } 
}

void ESKF::reset_priori_covariance_matrix() {
    _q_cov.setZero();
}

void ESKF::reset_covariance_matrix(const unsigned int start_index, const unsigned int end_index, const Matrix<float, 16, 1> &diag_cov) {
    for (unsigned int i = start_index; i < end_index; ++i) {
        // Diaginal
        _cov[i][i] = _q_cov[i] * _dt2;

        // Upper triangular
        for (unsigned int j = start_index; j < i; ++j) {
            _cov[j][i] = 0.f;
        }

        // Columns
        for (unsigned int j = 0; j < start_index; ++j) {
            _cov[j][i] = 0.f;
        }

        // Rows
        for (unsigned int j = end_index; j < 16; ++j) {
            _cov[i][j] = 0.f;
        }
    }
}

void ESKF::reset_accmulator_cov() {
    for (float &a : _accumulator_cov) {
        a = 0.f;
    }
}

void ESKF::regular_covariance_to_symmetric(unsigned int start_index, unsigned int end_index) {
    for (unsigned int i = start_index + 1; i < end_index; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            _cov[i][j] = _cov[j][i];
        }
    }
}

float ESKF::kahan_summation(float sum_previous, float input, float &accumulator) {
    /*
    accumulator中记录了sum_previous + y中, y舍弃掉的部分
    */
    const float y = input - accumulator;
    const float t = sum_previous + y;
    accumulator = (t - sum_previous) - y;
    return t;
}