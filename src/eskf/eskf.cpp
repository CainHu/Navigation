//
// Created by Cain on 2022/11/11.
//

#include "eskf.h"

using namespace eskf;

void ESKF::rotation_from_axis_angle(array<array<float, 3>, 3> &r, const array<float, 3> &a) {
    /* Rodrigues's Formula:
        * a = n * θ
        * R = cosθ*I + (1 - cosθ)*n*n' + sinθ*n^
        * */
    float a_norm_square = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    if (a_norm_square < _dt2 * 1e-12f) {
        r[0][0] = r[1][1] = r[2][2] = 1.f;
        r[0][1] = r[0][2] = r[1][0] = r[1][2] = r[2][0] = r[2][1] = 0.f;
    } else {
        float a_norm = sqrtf(a_norm_square);
        array<float, 3> a_unit = {a[0] / a_norm, a[1] / a_norm, a[2] / a_norm};
        float theta = a_norm;
        float cos_theta = cosf(theta), sin_theta = sinf(theta);
        float tmp = 1.f - cos_theta;

        float xx = a_unit[0] * a_unit[0] * tmp;
        float xy = a_unit[0] * a_unit[1] * tmp;
        float xz = a_unit[0] * a_unit[2] * tmp;
        float yy = a_unit[1] * a_unit[1] * tmp;
        float yz = a_unit[1] * a_unit[2] * tmp;
        float zz = a_unit[2] * a_unit[2] * tmp;

        float sx = sin_theta * a_unit[0];
        float sy = sin_theta * a_unit[1];
        float sz = sin_theta * a_unit[2];

        r[0][0] = cos_theta + xx;
        r[0][1] = xy - sz;
        r[0][2] = xz + sy;
        r[1][0] = xy + sz;
        r[1][1] = cos_theta + yy;
        r[1][2] = yz - sx;
        r[2][0] = xz - sy;
        r[2][1] = yz + sx;
        r[2][2] = cos_theta + zz;
    }
}

void ESKF::predict_state(const array<float, 3> &w, const array<float, 3> &a) {
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
    array<float, 3> a_body_corr = {a[0] - _state[9], a[1] - _state[10], a[2] - _state[11]};

    // a_world = R * (a - ba) - [0;0;g]
    array<float, 3> a_world = {_rot[0][0] * a_body_corr[0] + _rot[0][1] * a_body_corr[1] +_rot[0][2] * a_body_corr[2],
                                _rot[1][0] * a_body_corr[0] + _rot[1][1] * a_body_corr[1] +_rot[1][2] * a_body_corr[2],
                                _rot[2][0] * a_body_corr[0] + _rot[2][1] * a_body_corr[1] +_rot[2][2] * a_body_corr[2] + _state[12]};

    // v
    array<float, 3> v_last = {_state[3], _state[4], _state[5]};

    // v' = v + a * Δt
    _state[3] += a_world[0] * _dt;
    _state[4] += a_world[1] * _dt;
    _state[5] += a_world[2] * _dt;

    // p' = p + 0.5 * (v' + v) * Δt
    _state[0] += 0.5f * (_state[3] + v_last[0]) * _dt;
    _state[1] += 0.5f * (_state[4] + v_last[1]) * _dt;
    _state[2] += 0.5f * (_state[5] + v_last[2]) * _dt;

    // bg' = bg
    // ba' = ba

    // axis_angle = (w - bg) * Δt
    array<float, 3> axis_angle = {(w[0] - _state[6]) * _dt, (w[1] - _state[7]) * _dt, (w[2] - _state[8]) * _dt};

    // ΔR = Exp((w - bg) * Δt)
    array<array<float, 3>, 3> delta_rot {};
    rotation_from_axis_angle(delta_rot, axis_angle);

    // R' = R * ΔR
    array<array<float, 3>, 3> rot = _rot;
    _rot[0][0] = rot[0][0] * delta_rot[0][0] + rot[0][1] * delta_rot[1][0] + rot[0][2] * delta_rot[2][0];
    _rot[0][1] = rot[0][0] * delta_rot[0][1] + rot[0][1] * delta_rot[1][1] + rot[0][2] * delta_rot[2][1];
    _rot[0][2] = rot[0][0] * delta_rot[0][2] + rot[0][1] * delta_rot[1][2] + rot[0][2] * delta_rot[2][2];
    _rot[1][0] = rot[1][0] * delta_rot[0][0] + rot[1][1] * delta_rot[1][0] + rot[1][2] * delta_rot[2][0];
    _rot[1][1] = rot[1][0] * delta_rot[0][1] + rot[1][1] * delta_rot[1][1] + rot[1][2] * delta_rot[2][1];
    _rot[1][2] = rot[1][0] * delta_rot[0][2] + rot[1][1] * delta_rot[1][2] + rot[1][2] * delta_rot[2][2];
    _rot[2][0] = rot[2][0] * delta_rot[0][0] + rot[2][1] * delta_rot[1][0] + rot[2][2] * delta_rot[2][0];
    _rot[2][1] = rot[2][0] * delta_rot[0][1] + rot[2][1] * delta_rot[1][1] + rot[2][2] * delta_rot[2][1];
    _rot[2][2] = rot[2][0] * delta_rot[0][2] + rot[2][1] * delta_rot[1][2] + rot[2][2] * delta_rot[2][2];

//    /* Error Model:
//     * δp = δp + δv * Δt
//     * δv = δv + (-R * (a - ba)^ * δθ - R * δba + δg) * Δt
//     * δθ = Exp(-(w - bg) * Δt) * δθ - δbg * Δt
//     * δbg = δbg
//     * δba = δba
//     * δg = δg
//     * */
//    // [p, v, bg, ba, g]
//    for (unsigned char i = 0; i < 3; ++i) {
//        _state[i] += _error_state[i];
//        _state[i + 3] += _error_state[i + 3];
//        _state[i + 6] += _error_state[i + 9];
//        _state[i + 9] += _error_state[i + 12];
//    }
//    _state[12] += _error_state[15];
//
//    // Exp(δθ)
//    // cosθ*I + (1 - cosθ)*n*n' + sinθ*n^
//    float w_corr[3] = {_error_state[6] / _dt, _error_state[7] / _dt, _error_state[8] / _dt};
//    float y[3][3] = {0};
//    float w_corr_norm_square = w_corr[0] * w_corr[0] + w_corr[1] * w_corr[1] + w_corr[2] * w_corr[2];
//    if (w_corr_norm_square < 1e-12f) {
//        y[0][0] = y[1][1] = y[2][2] = 1.f;
//    } else {
//        float w_corr_norm = sqrtf(w_corr_norm_square);
//        float w_corr_unit[3] = {w_corr[0] / w_corr_norm, w_corr[1] / w_corr_norm, w_corr[2] / w_corr_norm};
//        float theta = w_corr_norm * _dt;
//        float cos_theta = cosf(theta), sin_theta = sinf(theta);
//        float tmp = 1.f - cos_theta;
//
//        float xx = w_corr_unit[0] * w_corr_unit[0] * tmp;
//        float xy = w_corr_unit[0] * w_corr_unit[1] * tmp;
//        float xz = w_corr_unit[0] * w_corr_unit[2] * tmp;
//        float yy = w_corr_unit[1] * w_corr_unit[1] * tmp;
//        float yz = w_corr_unit[1] * w_corr_unit[2] * tmp;
//        float zz = w_corr_unit[2] * w_corr_unit[2] * tmp;
//
//        float sx = sin_theta * w_corr_unit[0];
//        float sy = sin_theta * w_corr_unit[1];
//        float sz = sin_theta * w_corr_unit[2];
//
//        y[0][0] = cos_theta + xx;
//        y[0][1] = xy - sz;
//        y[0][2] = xz + sy;
//        y[1][0] = xy + sz;
//        y[1][1] = cos_theta + yy;
//        y[1][2] = yz - sx;
//        y[2][0] = xz - sy;
//        y[2][1] = yz + sx;
//        y[2][2] = cos_theta + zz;
//    }
//
//    // R * Exp(δθ)
//    vector<vector<float>> rot = _rot;
//    _rot[0][0] = rot[0][0] * y[0][0] + rot[0][1] * y[1][0] + rot[0][2] * y[2][0];
//    _rot[0][1] = rot[0][0] * y[0][1] + rot[0][1] * y[1][1] + rot[0][2] * y[2][1];
//    _rot[0][2] = rot[0][0] * y[0][2] + rot[0][1] * y[1][2] + rot[0][2] * y[2][2];
//    _rot[1][0] = rot[1][0] * y[0][0] + rot[1][1] * y[1][0] + rot[1][2] * y[2][0];
//    _rot[1][1] = rot[1][0] * y[0][1] + rot[1][1] * y[1][1] + rot[1][2] * y[2][1];
//    _rot[1][2] = rot[1][0] * y[0][2] + rot[1][1] * y[1][2] + rot[1][2] * y[2][2];
//    _rot[2][0] = rot[2][0] * y[0][0] + rot[2][1] * y[1][0] + rot[2][2] * y[2][0];
//    _rot[2][1] = rot[2][0] * y[0][1] + rot[2][1] * y[1][1] + rot[2][2] * y[2][1];
//    _rot[2][2] = rot[2][0] * y[0][2] + rot[2][1] * y[1][2] + rot[2][2] * y[2][2];
}

void ESKF::predict_covariance(const array<float, 3> &w, const array<float, 3> &a) {
    const array<float, 13> &state = get_state();
    const array<array<float, 3>, 3> &r = get_rotation_matrix();

    // x = R * (a - ba)^
    array<float, 3> a_corr = {a[0] - state[9], a[1] - state[10], a[2] - state[11]};
    array<array<float, 3>, 3> x = {{{-r[0][2]*a_corr[1] + r[0][1]*a_corr[2], r[0][2]*a_corr[0] - r[0][0]*a_corr[2], -r[0][1]*a_corr[0] + r[0][0]*a_corr[1]},
                                    {-r[1][2]*a_corr[1] + r[1][1]*a_corr[2], r[1][2]*a_corr[0] - r[1][0]*a_corr[2], -r[1][1]*a_corr[0] + r[1][0]*a_corr[1]},
                                    {-r[2][2]*a_corr[1] + r[2][1]*a_corr[2], r[2][2]*a_corr[0] - r[2][0]*a_corr[2], -r[2][1]*a_corr[0] + r[2][0]*a_corr[1]}}};

    // y = Exp(-(w - bg) * Δt)
    array<float, 3> w_corr = {w[0] - state[6], w[1] - state[7], w[2] - state[8]};
    array<float, 3> axis_angle = {-w_corr[0] * _dt, -w_corr[1] * _dt, -w_corr[2] * _dt};
    array<array<float, 3>, 3> y {};
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
    array<array<float, 3>, 3> c1 = {{{_cov[0][6] + _cov[3][6]*_dt,_cov[0][7] + _cov[3][7]*_dt,_cov[0][8] + _cov[3][8]*_dt},
                                        {_cov[1][6] + _cov[4][6]*_dt,_cov[1][7] + _cov[4][7]*_dt,_cov[1][8] + _cov[4][8]*_dt},
                                        {_cov[2][6] + _cov[5][6]*_dt,_cov[2][7] + _cov[5][7]*_dt,_cov[2][8] + _cov[5][8]*_dt}}};

    // P13_ = (P13 + P23 * dt) * Y' - P14_ * dt     (乘法: 36, 加法: 27)
    _cov[0][6] = c1[0][0]*y[0][0] + c1[0][1]*y[0][1] + c1[0][2]*y[0][2] - _cov[0][9]*_dt;
    _cov[0][7] = c1[0][0]*y[1][0] + c1[0][1]*y[1][1] + c1[0][2]*y[1][2] - _cov[0][10]*_dt;
    _cov[0][8] = c1[0][0]*y[2][0] + c1[0][1]*y[2][1] + c1[0][2]*y[2][2] - _cov[0][11]*_dt;
    _cov[1][6] = c1[1][0]*y[0][0] + c1[1][1]*y[0][1] + c1[1][2]*y[0][2] - _cov[1][9]*_dt;
    _cov[1][7] = c1[1][0]*y[1][0] + c1[1][1]*y[1][1] + c1[1][2]*y[1][2] - _cov[1][10]*_dt;
    _cov[1][8] = c1[1][0]*y[2][0] + c1[1][1]*y[2][1] + c1[1][2]*y[2][2] - _cov[1][11]*_dt;
    _cov[2][6] = c1[2][0]*y[0][0] + c1[2][1]*y[0][1] + c1[2][2]*y[0][2] - _cov[2][9]*_dt;
    _cov[2][7] = c1[2][0]*y[1][0] + c1[2][1]*y[1][1] + c1[2][2]*y[1][2] - _cov[2][10]*_dt;
    _cov[2][8] = c1[2][0]*y[2][0] + c1[2][1]*y[2][1] + c1[2][2]*y[2][2] - _cov[2][11]*_dt;

    // P12_ = (P12 + P22 * dt) - (P13 + P23 * dt) * X' * dt - P15_ * R' * dt + P16_ * ez' * dt      (乘法: 63, 加法: 63)
    //      = P12 + (P22 - ((P13 + P23 * dt) * X' + P15_ * R') + P16_ * ez') * dt
    //      = P12 + (P22 - (c1 * X' + P15_ * R') + P16_ * ez') * dt
    _cov[0][3] += (_cov[3][3] - (c1[0][0] * x[0][0] + c1[0][1] * x[0][1] + c1[0][2] * x[0][2] +
                                    _cov[0][12] * r[0][0] + _cov[0][13] * r[0][1] + _cov[0][14] * r[0][2])) * _dt;
    _cov[0][4] += (_cov[3][4] - (c1[0][0] * x[1][0] + c1[0][1] * x[1][1] + c1[0][2] * x[1][2] +
                                    _cov[0][12] * r[1][0] + _cov[0][13] * r[1][1] + _cov[0][14] * r[1][2])) * _dt;
    _cov[0][5] += (_cov[3][5] - (c1[0][0] * x[2][0] + c1[0][1] * x[2][1] + c1[0][2] * x[2][2] +
                                    _cov[0][12] * r[2][0] + _cov[0][13] * r[2][1] + _cov[0][14] * r[2][2]) + _cov[0][15]) * _dt;
    _cov[1][3] += (_cov[4][3] - (c1[1][0] * x[0][0] + c1[1][1] * x[0][1] + c1[1][2] * x[0][2] +
                                    _cov[1][12] * r[0][0] + _cov[1][13] * r[0][1] + _cov[1][14] * r[0][2])) * _dt;
    _cov[1][4] += (_cov[4][4] - (c1[1][0] * x[1][0] + c1[1][1] * x[1][1] + c1[1][2] * x[1][2] +
                                    _cov[1][12] * r[1][0] + _cov[1][13] * r[1][1] + _cov[1][14] * r[1][2])) * _dt;
    _cov[1][5] += (_cov[4][5] - (c1[1][0] * x[2][0] + c1[1][1] * x[2][1] + c1[1][2] * x[2][2] +
                                    _cov[1][12] * r[2][0] + _cov[1][13] * r[2][1] + _cov[1][14] * r[2][2]) + _cov[1][15]) * _dt;
    _cov[2][3] += (_cov[5][3] - (c1[2][0] * x[0][0] + c1[2][1] * x[0][1] + c1[2][2] * x[0][2] +
                                    _cov[2][12] * r[0][0] + _cov[2][13] * r[0][1] + _cov[2][14] * r[0][2])) * _dt;
    _cov[2][4] += (_cov[5][4] - (c1[2][0] * x[1][0] + c1[2][1] * x[1][1] + c1[2][2] * x[1][2] +
                                    _cov[2][12] * r[1][0] + _cov[2][13] * r[1][1] + _cov[2][14] * r[1][2])) * _dt;
    _cov[2][5] += (_cov[5][5] - (c1[2][0] * x[2][0] + c1[2][1] * x[2][1] + c1[2][2] * x[2][2] +
                                    _cov[2][12] * r[2][0] + _cov[2][13] * r[2][1] + _cov[2][14] * r[2][2]) + _cov[2][15]) * _dt;

    // X * P32              (乘法: 27, 加法: 18)
    array<array<float, 3>, 3> XP32 = {{{_cov[3][6]*x[0][0] + _cov[3][7]*x[0][1] + _cov[3][8]*x[0][2],
                                        _cov[4][6]*x[0][0] + _cov[4][7]*x[0][1] + _cov[4][8]*x[0][2],
                                                _cov[5][6]*x[0][0] + _cov[5][7]*x[0][1] + _cov[5][8]*x[0][2]},
                                                {_cov[3][6]*x[1][0] + _cov[3][7]*x[1][1] + _cov[3][8]*x[1][2],
                                                        _cov[4][6]*x[1][0] + _cov[4][7]*x[1][1] + _cov[4][8]*x[1][2],
                                                        _cov[5][6]*x[1][0] + _cov[5][7]*x[1][1] + _cov[5][8]*x[1][2]},
                                                {_cov[3][6]*x[2][0] + _cov[3][7]*x[2][1] + _cov[3][8]*x[2][2],
                                                        _cov[4][6]*x[2][0] + _cov[4][7]*x[2][1] + _cov[4][8]*x[2][2],
                                                        _cov[5][6]*x[2][0] + _cov[5][7]*x[2][1] + _cov[5][8]*x[2][2]}}};

    // R * P52              (乘法: 27, 加法: 18)
    array<array<float, 3>, 3> RP52 = {{{_cov[3][12]*r[0][0] + _cov[3][13]*r[0][1] + _cov[3][14]*r[0][2],
                                _cov[4][12]*r[0][0] + _cov[4][13]*r[0][1] + _cov[4][14]*r[0][2],
                                _cov[5][12]*r[0][0] + _cov[5][13]*r[0][1] + _cov[5][14]*r[0][2]},
                        {_cov[3][12]*r[1][0] + _cov[3][13]*r[1][1] + _cov[3][14]*r[1][2],
                                _cov[4][12]*r[1][0] + _cov[4][13]*r[1][1] + _cov[4][14]*r[1][2],
                                _cov[5][12]*r[1][0] + _cov[5][13]*r[1][1] + _cov[5][14]*r[1][2]},
                        {_cov[3][12]*r[2][0] + _cov[3][13]*r[2][1] + _cov[3][14]*r[2][2],
                                _cov[4][12]*r[2][0] + _cov[4][13]*r[2][1] + _cov[4][14]*r[2][2],
                                _cov[5][12]*r[2][0] + _cov[5][13]*r[2][1] + _cov[5][14]*r[2][2]}}};

    // X * P33              (乘法: 27, 加法: 18)
    array<array<float, 3>, 3> XP33 = {{{_cov[6][6]*x[0][0] + _cov[6][7]*x[0][1] + _cov[6][8]*x[0][2],
                                _cov[6][7]*x[0][0] + _cov[7][7]*x[0][1] + _cov[7][8]*x[0][2],
                                _cov[6][8]*x[0][0] + _cov[7][8]*x[0][1] + _cov[8][8]*x[0][2]},
                        {_cov[6][6]*x[1][0] + _cov[6][7]*x[1][1] + _cov[6][8]*x[1][2],
                                _cov[6][7]*x[1][0] + _cov[7][7]*x[1][1] + _cov[7][8]*x[1][2],
                                _cov[6][8]*x[1][0] + _cov[7][8]*x[1][1] + _cov[8][8]*x[1][2]},
                        {_cov[6][6]*x[2][0] + _cov[6][7]*x[2][1] + _cov[6][8]*x[2][2],
                                _cov[6][7]*x[2][0] + _cov[7][7]*x[2][1] + _cov[7][8]*x[2][2],
                                _cov[6][8]*x[2][0] + _cov[7][8]*x[2][1] + _cov[8][8]*x[2][2]}}};

    // R * P53              (乘法: 27, 加法: 18)
    array<array<float, 3>, 3> RP53 = {{{_cov[6][12]*r[0][0] + _cov[6][13]*r[0][1] + _cov[6][14]*r[0][2],
                                _cov[7][12]*r[0][0] + _cov[7][13]*r[0][1] + _cov[7][14]*r[0][2],
                                _cov[8][12]*r[0][0] + _cov[8][13]*r[0][1] + _cov[8][14]*r[0][2]},
                        {_cov[6][12]*r[1][0] + _cov[6][13]*r[1][1] + _cov[6][14]*r[1][2],
                                _cov[7][12]*r[1][0] + _cov[7][13]*r[1][1] + _cov[7][14]*r[1][2],
                                _cov[8][12]*r[1][0] + _cov[8][13]*r[1][1] + _cov[8][14]*r[1][2]},
                        {_cov[6][12]*r[2][0] + _cov[6][13]*r[2][1] + _cov[6][14]*r[2][2],
                                _cov[7][12]*r[2][0] + _cov[7][13]*r[2][1] + _cov[7][14]*r[2][2],
                                _cov[8][12]*r[2][0] + _cov[8][13]*r[2][1] + _cov[8][14]*r[2][2]}}};

    // X * P34              (乘法: 27, 加法: 18)
    array<array<float, 3>, 3> XP34 = {{{_cov[6][9]*x[0][0] + _cov[7][9]*x[0][1] + _cov[8][9]*x[0][2],
                                _cov[6][10]*x[0][0] + _cov[7][10]*x[0][1] + _cov[8][10]*x[0][2],
                                _cov[6][11]*x[0][0] + _cov[7][11]*x[0][1] + _cov[8][11]*x[0][2]},
                        {_cov[6][9]*x[1][0] + _cov[7][9]*x[1][1] + _cov[8][9]*x[1][2],
                                _cov[6][10]*x[1][0] + _cov[7][10]*x[1][1] + _cov[8][10]*x[1][2],
                                _cov[6][11]*x[1][0] + _cov[7][11]*x[1][1] + _cov[8][11]*x[1][2]},
                        {_cov[6][9]*x[2][0] + _cov[7][9]*x[2][1] + _cov[8][9]*x[2][2],
                                _cov[6][10]*x[2][0] + _cov[7][10]*x[2][1] + _cov[8][10]*x[2][2],
                                _cov[6][11]*x[2][0] + _cov[7][11]*x[2][1] + _cov[8][11]*x[2][2]}}};

    // R * P54              (乘法: 27, 加法: 18)
    array<array<float, 3>, 3> RP54 = {{{_cov[9][12]*r[0][0] + _cov[9][13]*r[0][1] + _cov[9][14]*r[0][2],
                                _cov[10][12]*r[0][0] + _cov[10][13]*r[0][1] + _cov[10][14]*r[0][2],
                                _cov[11][12]*r[0][0] + _cov[11][13]*r[0][1] + _cov[11][14]*r[0][2]},
                        {_cov[9][12]*r[1][0] + _cov[9][13]*r[1][1] + _cov[9][14]*r[1][2],
                                _cov[10][12]*r[1][0] + _cov[10][13]*r[1][1] + _cov[10][14]*r[1][2],
                                _cov[11][12]*r[1][0] + _cov[11][13]*r[1][1] + _cov[11][14]*r[1][2]},
                        {_cov[9][12]*r[2][0] + _cov[9][13]*r[2][1] + _cov[9][14]*r[2][2],
                                _cov[10][12]*r[2][0] + _cov[10][13]*r[2][1] + _cov[10][14]*r[2][2],
                                _cov[11][12]*r[2][0] + _cov[11][13]*r[2][1] + _cov[11][14]*r[2][2]}}};

    // X * P35              (乘法: 27, 加法: 18)
    float XP35[3][3] = {{_cov[6][12]*x[0][0] + _cov[7][12]*x[0][1] + _cov[8][12]*x[0][2],
                                _cov[6][13]*x[0][0] + _cov[7][13]*x[0][1] + _cov[8][13]*x[0][2],
                                _cov[6][14]*x[0][0] + _cov[7][14]*x[0][1] + _cov[8][14]*x[0][2]},
                        {_cov[6][12]*x[1][0] + _cov[7][12]*x[1][1] + _cov[8][12]*x[1][2],
                                _cov[6][13]*x[1][0] + _cov[7][13]*x[1][1] + _cov[8][13]*x[1][2],
                                _cov[6][14]*x[1][0] + _cov[7][14]*x[1][1] + _cov[8][14]*x[1][2]},
                        {_cov[6][12]*x[2][0] + _cov[7][12]*x[2][1] + _cov[8][12]*x[2][2],
                                _cov[6][13]*x[2][0] + _cov[7][13]*x[2][1] + _cov[8][13]*x[2][2],
                                _cov[6][14]*x[2][0] + _cov[7][14]*x[2][1] + _cov[8][14]*x[2][2]}};

    // R * P55              (乘法: 27, 加法: 18)
    array<array<float, 3>, 3> RP55 = {{{_cov[12][12]*r[0][0] + _cov[12][13]*r[0][1] + _cov[12][14]*r[0][2],
                                _cov[12][13]*r[0][0] + _cov[13][13]*r[0][1] + _cov[13][14]*r[0][2],
                                _cov[12][14]*r[0][0] + _cov[13][14]*r[0][1] + _cov[14][14]*r[0][2]},
                        {_cov[12][12]*r[1][0] + _cov[12][13]*r[1][1] + _cov[12][14]*r[1][2],
                                _cov[12][13]*r[1][0] + _cov[13][13]*r[1][1] + _cov[13][14]*r[1][2],
                                _cov[12][14]*r[1][0] + _cov[13][14]*r[1][1] + _cov[14][14]*r[1][2]},
                        {_cov[12][12]*r[2][0] + _cov[12][13]*r[2][1] + _cov[12][14]*r[2][2],
                                _cov[12][13]*r[2][0] + _cov[13][13]*r[2][1] + _cov[13][14]*r[2][2],
                                _cov[12][14]*r[2][0] + _cov[13][14]*r[2][1] + _cov[14][14]*r[2][2]}}};

    // X * P36              (乘法: 9, 加法: 6)
    array<float, 3> XP36 = {_cov[6][15]*x[0][0] + _cov[7][15]*x[0][1] + _cov[8][15]*x[0][2],
                        _cov[6][15]*x[1][0] + _cov[7][15]*x[1][1] + _cov[8][15]*x[1][2],
                        _cov[6][15]*x[2][0] + _cov[7][15]*x[2][1] + _cov[8][15]*x[2][2]};

    // R * P56              (乘法: 9, 加法: 6)
    array<float, 3> RP56 = {_cov[12][15]*r[0][0] + _cov[13][15]*r[0][1] + _cov[14][15]*r[0][2],
                        _cov[12][15]*r[1][0] + _cov[13][15]*r[1][1] + _cov[14][15]*r[1][2],
                        _cov[12][15]*r[2][0] + _cov[13][15]*r[2][1] + _cov[14][15]*r[2][2]};

    // ez * P62
    array<float, 3> ezP62 = {_cov[3][15], _cov[4][15], _cov[5][15]};

    // R * P53 * X'         (乘法: 27, 加法: 18)
    array<array<float, 3>, 3> RP53XT = {{{RP53[0][0]*x[0][0] + RP53[0][1]*x[0][1] + RP53[0][2]*x[0][2],
                                    RP53[0][0]*x[1][0] + RP53[0][1]*x[1][1] + RP53[0][2]*x[1][2],
                                    RP53[0][0]*x[2][0] + RP53[0][1]*x[2][1] + RP53[0][2]*x[2][2]},
                            {RP53[1][0]*x[0][0] + RP53[1][1]*x[0][1] + RP53[1][2]*x[0][2],
                                    RP53[1][0]*x[1][0] + RP53[1][1]*x[1][1] + RP53[1][2]*x[1][2],
                                    RP53[1][0]*x[2][0] + RP53[1][1]*x[2][1] + RP53[1][2]*x[2][2]},
                            {RP53[2][0]*x[0][0] + RP53[2][1]*x[0][1] + RP53[2][2]*x[0][2],
                                    RP53[2][0]*x[1][0] + RP53[2][1]*x[1][1] + RP53[2][2]*x[1][2],
                                    RP53[2][0]*x[2][0] + RP53[2][1]*x[2][1] + RP53[2][2]*x[2][2]}}};

    // X * P33 * X'         (乘法: 18, 加法: 12)
    float XP33XT00 = XP33[0][0]*x[0][0] + XP33[0][1]*x[0][1] + XP33[0][2]*x[0][2];
    float XP33XT01 = XP33[0][0]*x[1][0] + XP33[0][1]*x[1][1] + XP33[0][2]*x[1][2];
    float XP33XT02 = XP33[0][0]*x[2][0] + XP33[0][1]*x[2][1] + XP33[0][2]*x[2][2];
    float XP33XT11 = XP33[1][0]*x[1][0] + XP33[1][1]*x[1][1] + XP33[1][2]*x[1][2];
    float XP33XT12 = XP33[1][0]*x[2][0] + XP33[1][1]*x[2][1] + XP33[1][2]*x[2][2];
    float XP33XT22 = XP33[2][0]*x[2][0] + XP33[2][1]*x[2][1] + XP33[2][2]*x[2][2];

    // R * P55 * R'         (乘法: 18, 加法: 12)
    float RP55RT00 = RP55[0][0]*r[0][0] + RP55[0][1]*r[0][1] + RP55[0][2]*r[0][2];
    float RP55RT01 = RP55[0][0]*r[1][0] + RP55[0][1]*r[1][1] + RP55[0][2]*r[1][2];
    float RP55RT02 = RP55[0][0]*r[2][0] + RP55[0][1]*r[2][1] + RP55[0][2]*r[2][2];
    float RP55RT11 = RP55[1][0]*r[1][0] + RP55[1][1]*r[1][1] + RP55[1][2]*r[1][2];
    float RP55RT12 = RP55[1][0]*r[2][0] + RP55[1][1]*r[2][1] + RP55[1][2]*r[2][2];
    float RP55RT22 = RP55[2][0]*r[2][0] + RP55[2][1]*r[2][1] + RP55[2][2]*r[2][2];

    // P22_ = P22 - (X * P32 + P23 * X') * dt - (R * P52 + P25 * R') * dt + (ez * P62 + P26 * ez') * dt
    //        + X * P33 * X' * dt^2 + (R * P53 * X' + X * P35 * R') * dt^2 - (ez * P63 * X' + X * P36 * ez') * dt^2
    //        - (ez * P65 * R' + R * P56 * ez') * dt^2 + ez * P66 * ez' * dt^2          (乘法: 12, 加法: 48)
    _cov[3][3] += -(XP32[0][0] + XP32[0][0] + RP52[0][0] + RP52[0][0]) * _dt
                    + (XP33XT00 + RP53XT[0][0] + RP53XT[0][0] + RP55RT00) * _dt2;
    _cov[3][4] += -(XP32[0][1] + XP32[1][0] + RP52[0][1] + RP52[1][0]) * _dt
                    + (XP33XT01 + RP53XT[0][1] + RP53XT[1][0] + RP55RT01) * _dt2;
    _cov[3][5] += (ezP62[0] - XP32[0][2] - XP32[2][0] - RP52[0][2] - RP52[2][0]) * _dt
                    + (XP33XT02 + RP53XT[0][2] + RP53XT[2][0] + RP55RT02 - XP36[0] - RP56[0]) * _dt2;
    _cov[4][4] += -(XP32[1][1] + XP32[1][1] + RP52[1][1] + RP52[1][1]) * _dt
                    + (XP33XT11 + RP53XT[1][1] + RP53XT[1][1] + RP55RT11) * _dt2;
    _cov[4][5] += (ezP62[1] - XP32[1][2] - XP32[2][1] - RP52[1][2] - RP52[2][1]) * _dt
                    + (XP33XT12 + RP53XT[1][2] + RP53XT[2][1] + RP55RT12 - XP36[1] - RP56[1]) * _dt2;
    _cov[5][5] += (ezP62[2] + ezP62[2] - XP32[2][2] - XP32[2][2] - RP52[2][2] - RP52[2][2]) * _dt
                    + (XP33XT22 + RP53XT[2][2] + RP53XT[2][2] + RP55RT22 + _cov[15][15] - XP36[2] - XP36[2] - RP56[2] - RP56[2]) * _dt2;

    // P26_ = P26 - XP36 * dt - R * P56 * dt + ez * P66 * dt        (乘法: 3, 加法: 6)
    _cov[3][15] += -(XP36[0] + RP56[0]) * _dt;
    _cov[4][15] += -(XP36[1] + RP56[1]) * _dt;
    _cov[5][15] += (_cov[15][15] - XP36[2] - RP56[2]) * _dt;

    // P25_ = P25 - XP35 * dt - R * P55 * dt + ez * P65 * dt        (乘法: 9, 加法: 18)
    _cov[3][12] += -(XP35[0][0] + RP55[0][0]) * _dt;
    _cov[3][13] += -(XP35[0][1] + RP55[0][1]) * _dt;
    _cov[3][14] += -(XP35[0][2] + RP55[0][2]) * _dt;
    _cov[4][12] += -(XP35[1][0] + RP55[1][0]) * _dt;
    _cov[4][13] += -(XP35[1][1] + RP55[1][1]) * _dt;
    _cov[4][14] += -(XP35[1][2] + RP55[1][2]) * _dt;
    _cov[5][12] += (_cov[12][15] - XP35[2][0] - RP55[2][0]) * _dt;
    _cov[5][13] += (_cov[13][15] - XP35[2][1] - RP55[2][1]) * _dt;
    _cov[5][14] += (_cov[14][15] - XP35[2][2] - RP55[2][2]) * _dt;

    // P24_ = P24 - XP34 * dt - R * P54 * dt + ez * P64 * dt        (乘法: 9, 加法: 18)
    _cov[3][9] += -(XP34[0][0] + RP54[0][0]) * _dt;
    _cov[3][10] += -(XP34[0][1] + RP54[0][1]) * _dt;
    _cov[3][11] += -(XP34[0][2] + RP54[0][2]) * _dt;
    _cov[4][9] += -(XP34[1][0] + RP54[1][0]) * _dt;
    _cov[4][10] += -(XP34[1][1] + RP54[1][1]) * _dt;
    _cov[4][11] += -(XP34[1][2] + RP54[1][2]) * _dt;
    _cov[5][9] += (_cov[9][15] - XP34[2][0] - RP54[2][0]) * _dt;
    _cov[5][10] += (_cov[10][15] - XP34[2][1] - RP54[2][1]) * _dt;
    _cov[5][11] += (_cov[11][15] - XP34[2][2] - RP54[2][2]) * _dt;

    // (P23 - X * P33 * dt - R * P53 * dt + ez * P63 * dt)          (乘法: 9, 加法: 18)
    float c2[3][3] = {{_cov[3][6] - (XP33[0][0] + RP53[0][0]) * _dt,
                                _cov[3][7] - (XP33[0][1] + RP53[0][1]) * _dt,
                                _cov[3][8] - (XP33[0][2] + RP53[0][2]) * _dt},
                        {_cov[4][6] - (XP33[1][0] + RP53[1][0]) * _dt,
                                _cov[4][7] - (XP33[1][1] + RP53[1][1]) * _dt,
                                _cov[4][8] - (XP33[1][2] + RP53[1][2]) * _dt},
                        {_cov[5][6] + (_cov[6][15] - XP33[2][0] - RP53[2][0]) * _dt,
                                _cov[5][7] + (_cov[7][15] - XP33[2][1] - RP53[2][1]) * _dt,
                                _cov[5][8] + (_cov[8][15] - XP33[2][2] - RP53[2][2]) * _dt}};

    // P23_ = (P23 - X * P33 * dt - R * P53 * dt + ez * P63 * dt) * Y' - P24_ * dt          (乘法: 36, 加法: 27)
    _cov[3][6] = c2[0][0]*y[0][0] + c2[0][1]*y[0][1] + c2[0][2]*y[0][2] - _cov[3][9] * _dt;
    _cov[3][7] = c2[0][0]*y[1][0] + c2[0][1]*y[1][1] + c2[0][2]*y[1][2] - _cov[3][10] * _dt;
    _cov[3][8] = c2[0][0]*y[2][0] + c2[0][1]*y[2][1] + c2[0][2]*y[2][2] - _cov[3][11] * _dt;
    _cov[4][6] = c2[1][0]*y[0][0] + c2[1][1]*y[0][1] + c2[1][2]*y[0][2] - _cov[4][9] * _dt;
    _cov[4][7] = c2[1][0]*y[1][0] + c2[1][1]*y[1][1] + c2[1][2]*y[1][2] - _cov[4][10] * _dt;
    _cov[4][8] = c2[1][0]*y[2][0] + c2[1][1]*y[2][1] + c2[1][2]*y[2][2] - _cov[4][11] * _dt;
    _cov[5][6] = c2[2][0]*y[0][0] + c2[2][1]*y[0][1] + c2[2][2]*y[0][2] - _cov[5][9] * _dt;
    _cov[5][7] = c2[2][0]*y[1][0] + c2[2][1]*y[1][1] + c2[2][2]*y[1][2] - _cov[5][10] * _dt;
    _cov[5][8] = c2[2][0]*y[2][0] + c2[2][1]*y[2][1] + c2[2][2]*y[2][2] - _cov[5][11] * _dt;

    // Y * P36              (乘法: 9, 加法: 6)
    array<float, 3> YP36 = {y[0][0] * _cov[6][15] + y[0][1] * _cov[7][15] + y[0][2] * _cov[8][15],
                        y[1][0] * _cov[6][15] + y[1][1] * _cov[7][15] + y[1][2] * _cov[8][15],
                        y[2][0] * _cov[6][15] + y[2][1] * _cov[7][15] + y[2][2] * _cov[8][15]};

    // Y * P35              (乘法: 27, 加法: 18)
    array<array<float, 3>, 3> YP35 = {{{_cov[6][12]*y[0][0] + _cov[7][12]*y[0][1] + _cov[8][12]*y[0][2],
                                _cov[6][13]*y[0][0] + _cov[7][13]*y[0][1] + _cov[8][13]*y[0][2],
                                _cov[6][14]*y[0][0] + _cov[7][14]*y[0][1] + _cov[8][14]*y[0][2]},
                        {_cov[6][12]*y[1][0] + _cov[7][12]*y[1][1] + _cov[8][12]*y[1][2],
                                _cov[6][13]*y[1][0] + _cov[7][13]*y[1][1] + _cov[8][13]*y[1][2],
                                _cov[6][14]*y[1][0] + _cov[7][14]*y[1][1] + _cov[8][14]*y[1][2]},
                        {_cov[6][12]*y[2][0] + _cov[7][12]*y[2][1] + _cov[8][12]*y[2][2],
                                _cov[6][13]*y[2][0] + _cov[7][13]*y[2][1] + _cov[8][13]*y[2][2],
                                _cov[6][14]*y[2][0] + _cov[7][14]*y[2][1] + _cov[8][14]*y[2][2]}}};

    // Y * P34              (乘法: 27, 加法: 18)
    array<array<float, 3>, 3> YP34 = {{{_cov[6][9]*y[0][0] + _cov[7][9]*y[0][1] + _cov[8][9]*y[0][2],
                                _cov[6][10]*y[0][0] + _cov[7][10]*y[0][1] + _cov[8][10]*y[0][2],
                                _cov[6][11]*y[0][0] + _cov[7][11]*y[0][1] + _cov[8][11]*y[0][2]},
                        {_cov[6][9]*y[1][0] + _cov[7][9]*y[1][1] + _cov[8][9]*y[1][2],
                                _cov[6][10]*y[1][0] + _cov[7][10]*y[1][1] + _cov[8][10]*y[1][2],
                                _cov[6][11]*y[1][0] + _cov[7][11]*y[1][1] + _cov[8][11]*y[1][2]},
                        {_cov[6][9]*y[2][0] + _cov[7][9]*y[2][1] + _cov[8][9]*y[2][2],
                                _cov[6][10]*y[2][0] + _cov[7][10]*y[2][1] + _cov[8][10]*y[2][2],
                                _cov[6][11]*y[2][0] + _cov[7][11]*y[2][1] + _cov[8][11]*y[2][2]}}};

    // Y * P33              (乘法: 27, 加法: 18)
    array<array<float, 3>, 3> YP33 = {{{_cov[6][6]*y[0][0] + _cov[6][7]*y[0][1] + _cov[6][8]*y[0][2],
                                _cov[6][7]*y[0][0] + _cov[7][7]*y[0][1] + _cov[7][8]*y[0][2],
                                _cov[6][8]*y[0][0] + _cov[7][8]*y[0][1] + _cov[8][8]*y[0][2]},
                        {_cov[6][6]*y[1][0] + _cov[6][7]*y[1][1] + _cov[6][8]*y[1][2],
                                _cov[6][7]*y[1][0] + _cov[7][7]*y[1][1] + _cov[7][8]*y[1][2],
                                _cov[6][8]*y[1][0] + _cov[7][8]*y[1][1] + _cov[8][8]*y[1][2]},
                        {_cov[6][6]*y[2][0] + _cov[6][7]*y[2][1] + _cov[6][8]*y[2][2],
                                _cov[6][7]*y[2][0] + _cov[7][7]*y[2][1] + _cov[7][8]*y[2][2],
                                _cov[6][8]*y[2][0] + _cov[7][8]*y[2][1] + _cov[8][8]*y[2][2]}}};

    // Y * P33 * Y'         (乘法: 18, 加法: 12)
    float YP33YT00 = YP33[0][0]*y[0][0] + YP33[0][1]*y[0][1] + YP33[0][2]*y[0][2];
    float YP33YT01 = YP33[0][0]*y[1][0] + YP33[0][1]*y[1][1] + YP33[0][2]*y[1][2];
    float YP33YT02 = YP33[0][0]*y[2][0] + YP33[0][1]*y[2][1] + YP33[0][2]*y[2][2];
    float YP33YT11 = YP33[1][0]*y[1][0] + YP33[1][1]*y[1][1] + YP33[1][2]*y[1][2];
    float YP33YT12 = YP33[1][0]*y[2][0] + YP33[1][1]*y[2][1] + YP33[1][2]*y[2][2];
    float YP33YT22 = YP33[2][0]*y[2][0] + YP33[2][1]*y[2][1] + YP33[2][2]*y[2][2];

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

    for (unsigned char i = 1; i < 16; ++i) {
        for (unsigned char j = 0; j < i; ++j) {
            _cov[i][j] = _cov[j][i];
        }
    }
}

void ESKF::fuse_position(const array<float, 3> &pos_left, const array<float, 3> &pos_right) {
    /*
    s = [p, v, R, bg, ba, g]
    δs = [δp, δv, δθ, δbg, δba, δg]
    s + δs = [p+δp, v+δv, R*Exp(δθ), bg+δbg, ba+δba, g+δg]

    pl = p + R * dl
    pr = p + R * dr
    vl = v + R * (w - bg)^ * dl
    vr = v + R * (w - bg)^ * dr

    δpl / δp = I
    δpl / δθ = -R * dl^
    δpr / δp = I
    δpl / δθ = -R * dr^

    δvl / δv = I
    δvl / δθ = R * (dl^ * (w - bg))^
    δvl / δbg = R * dl^
    δvr / δv = I
    δvr / δθ = R * (dr^ * (w - bg))^
    δvr / δbg = R * dr^

    H = [I, O, -R*dl^,          O,     O, O
         O, O, -R*dr^,          O,     O, O
         O, I, R*(dl^*(w-bg))^, R*dl^, O, O
         O, I, R*(dr^*(w-bg))^, R*dr^, O, 0]
    */

    // array<array<float, 3>, 3> rot_dl_hat = {{{_rot[0][1] * _dl[2] - _rot[2][0] * _dl[1]},
    //                                          {_rot[0][1]},
    //                                          {}}};
}
void fuse_velocity(const array<float, 3> &vel_left, const array<float, 3> &vel_right);

void fuse_position_left_only(const array<float, 3> pos);
void fuse_velocity_left_only(const array<float, 3> vel);
void fuse_position_right_only(const array<float, 3> pos);
void fuse_velocity_right_only(const array<float, 3> vel);

void fuse_position_z(float pos_z);
void fuse_velocity_xy(const array<float, 2> vel_xy);