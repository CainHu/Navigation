//#pragma GCC optimize(2)

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <Eigen/Dense>

#define P(n, m) _cov[n][m]

using namespace std;

class ESKF {
public:
    ESKF(float dt, float g=9.8f)
            : _error_state(16, 0.f), _cov(16, vector<float>(16, 0.f)),
              _state(13, 0.f), _rot(3, vector<float>(3, 0.f)), _dt(dt), _dt2(dt * dt) {
        _rot[0][0] = _rot[1][1] = _rot[2][2] = 1.f;
        for (unsigned char i = 0; i < 16; ++i) {
            _cov[i][i] = 1.f;
        }
        _state[12] = g;
    }

    void predict_covariance(const vector<float> &w, const vector<float> &a);
    void predict_state(const vector<float> &w, const vector<float> &a);
    void true_covariance(const vector<float> &w, const vector<float> &a);

    void rotation_from_axis_angle(vector<vector<float>> &r, const vector<float> &a);

    const vector<float> &get_state() { return _state; };
    const vector<float> &get_error_state() { return _error_state; };
    const vector<vector<float>> &get_covariance_matrix() { return _cov; };
    const vector<vector<float>> &get_rotation_matrix() { return _rot; };
    const float get_dt() { return _dt; }

    void set_dt(float dt) { _dt = dt; _dt2 = dt * dt; };
private:
    float _dt;                      // Sample time of IMU
    float _dt2;                     // Square of sample time

    vector<float> _state;           // [p, v, bg, ba, g]
    vector<vector<float>> _rot;     // Rotation matrix from body frame to world frame

    vector<float> _error_state;     // [δp, δv, δθ, δbg, δba, δg]
    vector<vector<float>> _cov;     // Covariance matrix of error state
};

void ESKF::rotation_from_axis_angle(vector<vector<float>> &r, const vector<float> &a) {
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
        float a_unit[3] = {a[0] / a_norm, a[1] / a_norm, a[2] / a_norm};
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

void ESKF::predict_state(const vector<float> &w, const vector<float> &a) {
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
    vector<float> a_body_corr(3);
    a_body_corr[0] = a[0] - _state[9];
    a_body_corr[1] = a[1] - _state[10];
    a_body_corr[2] = a[2] - _state[11];

    // a_world = R * (a - ba)
    vector<float> a_world(3);
    a_world[0] = _rot[0][0] * a_body_corr[0] + _rot[0][1] * a_body_corr[1] +_rot[0][2] * a_body_corr[2];
    a_world[1] = _rot[1][0] * a_body_corr[0] + _rot[1][1] * a_body_corr[1] +_rot[1][2] * a_body_corr[2];
    a_world[2] = _rot[2][0] * a_body_corr[0] + _rot[2][1] * a_body_corr[1] +_rot[2][2] * a_body_corr[2] + _state[12];

    // v
    vector<float> v_last(3);
    v_last[0] = _state[3];
    v_last[1] = _state[4];
    v_last[2] = _state[5];

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
    vector<float> axis_angle(3);
    axis_angle[0] = (w[0] - _state[6]) * _dt;
    axis_angle[1] = (w[1] - _state[7]) * _dt;
    axis_angle[2] = (w[2] - _state[8]) * _dt;

    // ΔR = Exp((w - bg) * Δt)
    vector<vector<float>> delta_rot(3, vector<float>(3));
    rotation_from_axis_angle(delta_rot, axis_angle);

    // R' = R * ΔR
    vector<vector<float>> rot = _rot;
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

void ESKF::predict_covariance(const vector<float> &w, const vector<float> &a) {
    const vector<float> &state = get_state();
    const vector<vector<float>> &r = get_rotation_matrix();

    // x = R * (a - ba)^
    float a_corr[3] = {a[0] - state[9], a[1] - state[10], a[2] - state[11]};
    float x[3][3] = {{-r[0][2]*a_corr[1] + r[0][1]*a_corr[2], r[0][2]*a_corr[0] - r[0][0]*a_corr[2], -r[0][1]*a_corr[0] + r[0][0]*a_corr[1]},
                     {-r[1][2]*a_corr[1] + r[1][1]*a_corr[2], r[1][2]*a_corr[0] - r[1][0]*a_corr[2], -r[1][1]*a_corr[0] + r[1][0]*a_corr[1]},
                     {-r[2][2]*a_corr[1] + r[2][1]*a_corr[2], r[2][2]*a_corr[0] - r[2][0]*a_corr[2], -r[2][1]*a_corr[0] + r[2][0]*a_corr[1]}};

    // y = Exp(-(w - bg) * Δt)
    float w_corr[3] = {w[0] - state[6], w[1] - state[7], w[2] - state[8]};
    vector<float> axis_angle = {-w_corr[0] * _dt, -w_corr[1] * _dt, -w_corr[2] * _dt};
    vector<vector<float>> y(3, vector<float>(3));
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
    float c1[3][3] = {{_cov[0][6] + _cov[3][6]*_dt,_cov[0][7] + _cov[3][7]*_dt,_cov[0][8] + _cov[3][8]*_dt},
                      {_cov[1][6] + _cov[4][6]*_dt,_cov[1][7] + _cov[4][7]*_dt,_cov[1][8] + _cov[4][8]*_dt},
                      {_cov[2][6] + _cov[5][6]*_dt,_cov[2][7] + _cov[5][7]*_dt,_cov[2][8] + _cov[5][8]*_dt}};

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
    //      = P12 + (P22 - ((P13 + P23 * dt) * X' + P15_ * R') + P16_ * ez) * dt
    //      = P12 + (P22 - (c1 * X' + P15_ * R') + P16_ * ez) * dt
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
    float XP32[3][3] = {{_cov[3][6]*x[0][0] + _cov[3][7]*x[0][1] + _cov[3][8]*x[0][2],
                                _cov[4][6]*x[0][0] + _cov[4][7]*x[0][1] + _cov[4][8]*x[0][2],
                                _cov[5][6]*x[0][0] + _cov[5][7]*x[0][1] + _cov[5][8]*x[0][2]},
                        {_cov[3][6]*x[1][0] + _cov[3][7]*x[1][1] + _cov[3][8]*x[1][2],
                                _cov[4][6]*x[1][0] + _cov[4][7]*x[1][1] + _cov[4][8]*x[1][2],
                                _cov[5][6]*x[1][0] + _cov[5][7]*x[1][1] + _cov[5][8]*x[1][2]},
                        {_cov[3][6]*x[2][0] + _cov[3][7]*x[2][1] + _cov[3][8]*x[2][2],
                                _cov[4][6]*x[2][0] + _cov[4][7]*x[2][1] + _cov[4][8]*x[2][2],
                                _cov[5][6]*x[2][0] + _cov[5][7]*x[2][1] + _cov[5][8]*x[2][2]}};

    // R * P52              (乘法: 27, 加法: 18)
    float RP52[3][3] = {{_cov[3][12]*r[0][0] + _cov[3][13]*r[0][1] + _cov[3][14]*r[0][2],
                                _cov[4][12]*r[0][0] + _cov[4][13]*r[0][1] + _cov[4][14]*r[0][2],
                                _cov[5][12]*r[0][0] + _cov[5][13]*r[0][1] + _cov[5][14]*r[0][2]},
                        {_cov[3][12]*r[1][0] + _cov[3][13]*r[1][1] + _cov[3][14]*r[1][2],
                                _cov[4][12]*r[1][0] + _cov[4][13]*r[1][1] + _cov[4][14]*r[1][2],
                                _cov[5][12]*r[1][0] + _cov[5][13]*r[1][1] + _cov[5][14]*r[1][2]},
                        {_cov[3][12]*r[2][0] + _cov[3][13]*r[2][1] + _cov[3][14]*r[2][2],
                                _cov[4][12]*r[2][0] + _cov[4][13]*r[2][1] + _cov[4][14]*r[2][2],
                                _cov[5][12]*r[2][0] + _cov[5][13]*r[2][1] + _cov[5][14]*r[2][2]}};

    // X * P33              (乘法: 27, 加法: 18)
    float XP33[3][3] = {{_cov[6][6]*x[0][0] + _cov[6][7]*x[0][1] + _cov[6][8]*x[0][2],
                                _cov[6][7]*x[0][0] + _cov[7][7]*x[0][1] + _cov[7][8]*x[0][2],
                                _cov[6][8]*x[0][0] + _cov[7][8]*x[0][1] + _cov[8][8]*x[0][2]},
                        {_cov[6][6]*x[1][0] + _cov[6][7]*x[1][1] + _cov[6][8]*x[1][2],
                                _cov[6][7]*x[1][0] + _cov[7][7]*x[1][1] + _cov[7][8]*x[1][2],
                                _cov[6][8]*x[1][0] + _cov[7][8]*x[1][1] + _cov[8][8]*x[1][2]},
                        {_cov[6][6]*x[2][0] + _cov[6][7]*x[2][1] + _cov[6][8]*x[2][2],
                                _cov[6][7]*x[2][0] + _cov[7][7]*x[2][1] + _cov[7][8]*x[2][2],
                                _cov[6][8]*x[2][0] + _cov[7][8]*x[2][1] + _cov[8][8]*x[2][2]}};

    // R * P53              (乘法: 27, 加法: 18)
    float RP53[3][3] = {{_cov[6][12]*r[0][0] + _cov[6][13]*r[0][1] + _cov[6][14]*r[0][2],
                                _cov[7][12]*r[0][0] + _cov[7][13]*r[0][1] + _cov[7][14]*r[0][2],
                                _cov[8][12]*r[0][0] + _cov[8][13]*r[0][1] + _cov[8][14]*r[0][2]},
                        {_cov[6][12]*r[1][0] + _cov[6][13]*r[1][1] + _cov[6][14]*r[1][2],
                                _cov[7][12]*r[1][0] + _cov[7][13]*r[1][1] + _cov[7][14]*r[1][2],
                                _cov[8][12]*r[1][0] + _cov[8][13]*r[1][1] + _cov[8][14]*r[1][2]},
                        {_cov[6][12]*r[2][0] + _cov[6][13]*r[2][1] + _cov[6][14]*r[2][2],
                                _cov[7][12]*r[2][0] + _cov[7][13]*r[2][1] + _cov[7][14]*r[2][2],
                                _cov[8][12]*r[2][0] + _cov[8][13]*r[2][1] + _cov[8][14]*r[2][2]}};

    // X * P34              (乘法: 27, 加法: 18)
    float XP34[3][3] = {{_cov[6][9]*x[0][0] + _cov[7][9]*x[0][1] + _cov[8][9]*x[0][2],
                                _cov[6][10]*x[0][0] + _cov[7][10]*x[0][1] + _cov[8][10]*x[0][2],
                                _cov[6][11]*x[0][0] + _cov[7][11]*x[0][1] + _cov[8][11]*x[0][2]},
                        {_cov[6][9]*x[1][0] + _cov[7][9]*x[1][1] + _cov[8][9]*x[1][2],
                                _cov[6][10]*x[1][0] + _cov[7][10]*x[1][1] + _cov[8][10]*x[1][2],
                                _cov[6][11]*x[1][0] + _cov[7][11]*x[1][1] + _cov[8][11]*x[1][2]},
                        {_cov[6][9]*x[2][0] + _cov[7][9]*x[2][1] + _cov[8][9]*x[2][2],
                                _cov[6][10]*x[2][0] + _cov[7][10]*x[2][1] + _cov[8][10]*x[2][2],
                                _cov[6][11]*x[2][0] + _cov[7][11]*x[2][1] + _cov[8][11]*x[2][2]}};

    // R * P54              (乘法: 27, 加法: 18)
    float RP54[3][3] = {{_cov[9][12]*r[0][0] + _cov[9][13]*r[0][1] + _cov[9][14]*r[0][2],
                                _cov[10][12]*r[0][0] + _cov[10][13]*r[0][1] + _cov[10][14]*r[0][2],
                                _cov[11][12]*r[0][0] + _cov[11][13]*r[0][1] + _cov[11][14]*r[0][2]},
                        {_cov[9][12]*r[1][0] + _cov[9][13]*r[1][1] + _cov[9][14]*r[1][2],
                                _cov[10][12]*r[1][0] + _cov[10][13]*r[1][1] + _cov[10][14]*r[1][2],
                                _cov[11][12]*r[1][0] + _cov[11][13]*r[1][1] + _cov[11][14]*r[1][2]},
                        {_cov[9][12]*r[2][0] + _cov[9][13]*r[2][1] + _cov[9][14]*r[2][2],
                                _cov[10][12]*r[2][0] + _cov[10][13]*r[2][1] + _cov[10][14]*r[2][2],
                                _cov[11][12]*r[2][0] + _cov[11][13]*r[2][1] + _cov[11][14]*r[2][2]}};

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
    float RP55[3][3] = {{_cov[12][12]*r[0][0] + _cov[12][13]*r[0][1] + _cov[12][14]*r[0][2],
                                _cov[12][13]*r[0][0] + _cov[13][13]*r[0][1] + _cov[13][14]*r[0][2],
                                _cov[12][14]*r[0][0] + _cov[13][14]*r[0][1] + _cov[14][14]*r[0][2]},
                        {_cov[12][12]*r[1][0] + _cov[12][13]*r[1][1] + _cov[12][14]*r[1][2],
                                _cov[12][13]*r[1][0] + _cov[13][13]*r[1][1] + _cov[13][14]*r[1][2],
                                _cov[12][14]*r[1][0] + _cov[13][14]*r[1][1] + _cov[14][14]*r[1][2]},
                        {_cov[12][12]*r[2][0] + _cov[12][13]*r[2][1] + _cov[12][14]*r[2][2],
                                _cov[12][13]*r[2][0] + _cov[13][13]*r[2][1] + _cov[13][14]*r[2][2],
                                _cov[12][14]*r[2][0] + _cov[13][14]*r[2][1] + _cov[14][14]*r[2][2]}};

    // X * P36              (乘法: 9, 加法: 6)
    float XP36[3] = {_cov[6][15]*x[0][0] + _cov[7][15]*x[0][1] + _cov[8][15]*x[0][2],
                     _cov[6][15]*x[1][0] + _cov[7][15]*x[1][1] + _cov[8][15]*x[1][2],
                     _cov[6][15]*x[2][0] + _cov[7][15]*x[2][1] + _cov[8][15]*x[2][2]};

    // R * P56              (乘法: 9, 加法: 6)
    float RP56[3] = {_cov[12][15]*r[0][0] + _cov[13][15]*r[0][1] + _cov[14][15]*r[0][2],
                     _cov[12][15]*r[1][0] + _cov[13][15]*r[1][1] + _cov[14][15]*r[1][2],
                     _cov[12][15]*r[2][0] + _cov[13][15]*r[2][1] + _cov[14][15]*r[2][2]};

    // ez * P62
    float ezP62[3] = {_cov[3][15], _cov[4][15], _cov[5][15]};

    // R * P53 * X'         (乘法: 27, 加法: 18)
    float RP53XT[3][3] = {{RP53[0][0]*x[0][0] + RP53[0][1]*x[0][1] + RP53[0][2]*x[0][2],
                                  RP53[0][0]*x[1][0] + RP53[0][1]*x[1][1] + RP53[0][2]*x[1][2],
                                  RP53[0][0]*x[2][0] + RP53[0][1]*x[2][1] + RP53[0][2]*x[2][2]},
                          {RP53[1][0]*x[0][0] + RP53[1][1]*x[0][1] + RP53[1][2]*x[0][2],
                                  RP53[1][0]*x[1][0] + RP53[1][1]*x[1][1] + RP53[1][2]*x[1][2],
                                  RP53[1][0]*x[2][0] + RP53[1][1]*x[2][1] + RP53[1][2]*x[2][2]},
                          {RP53[2][0]*x[0][0] + RP53[2][1]*x[0][1] + RP53[2][2]*x[0][2],
                                  RP53[2][0]*x[1][0] + RP53[2][1]*x[1][1] + RP53[2][2]*x[1][2],
                                  RP53[2][0]*x[2][0] + RP53[2][1]*x[2][1] + RP53[2][2]*x[2][2]}};

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

//    // P22_ = P22 - (X * P32 + P23 * X') * dt - (R * P52 + P25 * R') * dt + (ez * P62 + P26 * ez') * dt
//    //        + X * P33 * X' * dt^2 + (R * P53 * X' + X * P35 * R') * dt^2 - (ez * P63 * X' + X * P36 * ez') * dt^2
//    //        - (ez * P65 * R' + R * P56 * ez') * dt^2 + ez * P66 * ez' * dt^2          (乘法: 12, 加法: 48)
//    _cov[3][3] += -(XP32[0][0] + RP52[0][0] + c2[0][0] * x[0][0] + c2[0][1] * x[0][1] + c2[0][2] * x[0][2]
//                        + _cov[3][12] * r[0][0] + _cov[3][13] * r[0][1] + _cov[3][14] * r[0][2]) * _dt;
//    _cov[3][4] += -(XP32[0][1] + RP52[0][1] + c2[0][0] * x[1][0] + c2[0][1] * x[1][1] + c2[0][2] * x[1][2]
//                     + _cov[3][12] * r[1][0] + _cov[3][13] * r[1][1] + _cov[3][14] * r[1][2]) * _dt;
//    _cov[3][5] += (_cov[3][15] - (XP32[0][2] + RP52[0][2] + c2[0][0] * x[2][0] + c2[0][1] * x[2][1] + c2[0][2] * x[2][2]
//                                  + _cov[3][12] * r[2][0] + _cov[3][13] * r[2][1] + _cov[3][14] * r[2][2])) * _dt;
//    _cov[4][4] += -(XP32[1][1] + RP52[1][1] + c2[1][0] * x[1][0] + c2[1][1] * x[1][1] + c2[1][2] * x[1][2]
//                    + _cov[4][12] * r[1][0] + _cov[4][13] * r[1][1] + _cov[4][14] * r[1][2]) * _dt;
//    _cov[4][5] += (_cov[4][15] - (XP32[1][2] + RP52[1][2] + c2[1][0] * x[2][0] + c2[1][1] * x[2][1] + c2[1][2] * x[2][2]
//                                  + _cov[4][12] * r[2][0] + _cov[4][13] * r[2][1] + _cov[4][14] * r[2][2])) * _dt;
//    _cov[5][5] += (_cov[5][15] + _cov[5][15] - (XP32[2][2] + RP52[2][2] + c2[2][0] * x[2][0] + c2[2][1] * x[2][1] + c2[2][2] * x[2][2]
//                                  + _cov[5][12] * r[2][0] + _cov[5][13] * r[2][1] + _cov[5][14] * r[2][2])) * _dt;

    // Y * P36              (乘法: 9, 加法: 6)
    float YP36[3] = {y[0][0] * _cov[6][15] + y[0][1] * _cov[7][15] + y[0][2] * _cov[8][15],
                     y[1][0] * _cov[6][15] + y[1][1] * _cov[7][15] + y[1][2] * _cov[8][15],
                     y[2][0] * _cov[6][15] + y[2][1] * _cov[7][15] + y[2][2] * _cov[8][15]};

    // Y * P35              (乘法: 27, 加法: 18)
    float YP35[3][3] = {{_cov[6][12]*y[0][0] + _cov[7][12]*y[0][1] + _cov[8][12]*y[0][2],
                                _cov[6][13]*y[0][0] + _cov[7][13]*y[0][1] + _cov[8][13]*y[0][2],
                                _cov[6][14]*y[0][0] + _cov[7][14]*y[0][1] + _cov[8][14]*y[0][2]},
                        {_cov[6][12]*y[1][0] + _cov[7][12]*y[1][1] + _cov[8][12]*y[1][2],
                                _cov[6][13]*y[1][0] + _cov[7][13]*y[1][1] + _cov[8][13]*y[1][2],
                                _cov[6][14]*y[1][0] + _cov[7][14]*y[1][1] + _cov[8][14]*y[1][2]},
                        {_cov[6][12]*y[2][0] + _cov[7][12]*y[2][1] + _cov[8][12]*y[2][2],
                                _cov[6][13]*y[2][0] + _cov[7][13]*y[2][1] + _cov[8][13]*y[2][2],
                                _cov[6][14]*y[2][0] + _cov[7][14]*y[2][1] + _cov[8][14]*y[2][2]}};

    // Y * P34              (乘法: 27, 加法: 18)
    float YP34[3][3] = {{_cov[6][9]*y[0][0] + _cov[7][9]*y[0][1] + _cov[8][9]*y[0][2],
                                _cov[6][10]*y[0][0] + _cov[7][10]*y[0][1] + _cov[8][10]*y[0][2],
                                _cov[6][11]*y[0][0] + _cov[7][11]*y[0][1] + _cov[8][11]*y[0][2]},
                        {_cov[6][9]*y[1][0] + _cov[7][9]*y[1][1] + _cov[8][9]*y[1][2],
                                _cov[6][10]*y[1][0] + _cov[7][10]*y[1][1] + _cov[8][10]*y[1][2],
                                _cov[6][11]*y[1][0] + _cov[7][11]*y[1][1] + _cov[8][11]*y[1][2]},
                        {_cov[6][9]*y[2][0] + _cov[7][9]*y[2][1] + _cov[8][9]*y[2][2],
                                _cov[6][10]*y[2][0] + _cov[7][10]*y[2][1] + _cov[8][10]*y[2][2],
                                _cov[6][11]*y[2][0] + _cov[7][11]*y[2][1] + _cov[8][11]*y[2][2]}};

    // Y * P33              (乘法: 27, 加法: 18)
    float YP33[3][3] = {{_cov[6][6]*y[0][0] + _cov[6][7]*y[0][1] + _cov[6][8]*y[0][2],
                         _cov[6][7]*y[0][0] + _cov[7][7]*y[0][1] + _cov[7][8]*y[0][2],
                         _cov[6][8]*y[0][0] + _cov[7][8]*y[0][1] + _cov[8][8]*y[0][2]},
                        {_cov[6][6]*y[1][0] + _cov[6][7]*y[1][1] + _cov[6][8]*y[1][2],
                         _cov[6][7]*y[1][0] + _cov[7][7]*y[1][1] + _cov[7][8]*y[1][2],
                         _cov[6][8]*y[1][0] + _cov[7][8]*y[1][1] + _cov[8][8]*y[1][2]},
                        {_cov[6][6]*y[2][0] + _cov[6][7]*y[2][1] + _cov[6][8]*y[2][2],
                         _cov[6][7]*y[2][0] + _cov[7][7]*y[2][1] + _cov[7][8]*y[2][2],
                         _cov[6][8]*y[2][0] + _cov[7][8]*y[2][1] + _cov[8][8]*y[2][2]}};

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

//void ESKF::true_covariance(const vector<float> &w, const vector<float> &a) {
////    auto P = [&](int n, int m) -> float& {
////        return _cov[n][m];
////    };
//
//    const float dt = _dt;
//    const vector<float> &state = get_error_state();
//    const vector<vector<float>> &r = get_rotation_matrix();
//
//    const float R00 = r[0][0], R01 = r[0][1], R02 = r[0][2];
//    const float R10 = r[1][0], R11 = r[1][1], R12 = r[1][2];
//    const float R20 = r[2][0], R21 = r[2][1], R22 = r[2][2];
//
//    const float ax = a[0], ay = a[1], az = a[2];
//    const float bax =  state[9], bay = state[10], baz = state[11];
//
//    float w_corr[3] = {w[0] - state[6], w[1] - state[7], w[2] - state[8]};
//    vector<float> axis_angle = {-w_corr[0] * _dt, -w_corr[1] * _dt, -w_corr[2] * _dt};
//    vector<vector<float>> y(3, vector<float>(3));
//    rotation_from_axis_angle(y, axis_angle);
//
//    const float T00 = y[0][0], T01 = y[0][1], T02 = y[0][2];
//    const float T10 = y[1][0], T11 = y[1][1], T12 = y[1][2];
//    const float T20 = y[2][0], T21 = y[2][1], T22 = y[2][2];
//
//    const float PS0 = P(0,3) + P(3,3)*dt;
//    const float PS1 = P(3,4)*dt;
//    const float PS2 = P(0,4) + PS1;
//    const float PS3 = P(3,5)*dt;
//    const float PS4 = P(0,5) + PS3;
//    const float PS5 = P(3,12)*dt;
//    const float PS6 = P(0,12) + PS5;
//    const float PS7 = PS6*dt;
//    const float PS8 = P(3,13)*dt;
//    const float PS9 = P(0,13) + PS8;
//    const float PS10 = PS9*dt;
//    const float PS11 = P(3,14)*dt;
//    const float PS12 = P(0,14) + PS11;
//    const float PS13 = PS12*dt;
//    const float PS14 = az - baz;
//    const float PS15 = -ay + bay;
//    const float PS16 = PS14*R01 + PS15*R02;
//    const float PS17 = P(3,6)*dt;
//    const float PS18 = P(0,6) + PS17;
//    const float PS19 = PS18*dt;
//    const float PS20 = -az + baz;
//    const float PS21 = ax - bax;
//    const float PS22 = PS20*R00 + PS21*R02;
//    const float PS23 = P(3,7)*dt;
//    const float PS24 = P(0,7) + PS23;
//    const float PS25 = PS24*dt;
//    const float PS26 = ay - bay;
//    const float PS27 = -ax + bax;
//    const float PS28 = PS26*R00 + PS27*R01;
//    const float PS29 = P(3,8)*dt;
//    const float PS30 = P(0,8) + PS29;
//    const float PS31 = PS30*dt;
//    const float PS32 = PS14*R11 + PS15*R12;
//    const float PS33 = PS20*R10 + PS21*R12;
//    const float PS34 = PS26*R10 + PS27*R11;
//    const float PS35 = P(0,15) + P(3,15)*dt;
//    const float PS36 = PS14*R21 + PS15*R22;
//    const float PS37 = PS20*R20 + PS21*R22;
//    const float PS38 = PS26*R20 + PS27*R21;
//    const float PS39 = P(0,9) + P(3,9)*dt;
//    const float PS40 = P(0,10) + P(3,10)*dt;
//    const float PS41 = P(0,11) + P(3,11)*dt;
//    const float PS42 = P(1,4) + P(4,4)*dt;
//    const float PS43 = P(4,5)*dt;
//    const float PS44 = P(1,5) + PS43;
//    const float PS45 = P(4,12)*dt;
//    const float PS46 = P(1,12) + PS45;
//    const float PS47 = PS46*dt;
//    const float PS48 = P(4,13)*dt;
//    const float PS49 = P(1,13) + PS48;
//    const float PS50 = PS49*dt;
//    const float PS51 = P(4,14)*dt;
//    const float PS52 = P(1,14) + PS51;
//    const float PS53 = PS52*dt;
//    const float PS54 = P(4,6)*dt;
//    const float PS55 = P(1,6) + PS54;
//    const float PS56 = PS55*dt;
//    const float PS57 = P(4,7)*dt;
//    const float PS58 = P(1,7) + PS57;
//    const float PS59 = PS58*dt;
//    const float PS60 = P(4,8)*dt;
//    const float PS61 = P(1,8) + PS60;
//    const float PS62 = PS61*dt;
//    const float PS63 = P(1,15) + P(4,15)*dt;
//    const float PS64 = P(1,9) + P(4,9)*dt;
//    const float PS65 = P(1,10) + P(4,10)*dt;
//    const float PS66 = P(1,11) + P(4,11)*dt;
//    const float PS67 = P(2,5) + P(5,5)*dt;
//    const float PS68 = P(5,12)*dt;
//    const float PS69 = P(2,12) + PS68;
//    const float PS70 = PS69*dt;
//    const float PS71 = P(5,13)*dt;
//    const float PS72 = P(2,13) + PS71;
//    const float PS73 = PS72*dt;
//    const float PS74 = P(5,14)*dt;
//    const float PS75 = P(2,14) + PS74;
//    const float PS76 = PS75*dt;
//    const float PS77 = P(5,6)*dt;
//    const float PS78 = P(2,6) + PS77;
//    const float PS79 = PS78*dt;
//    const float PS80 = P(5,7)*dt;
//    const float PS81 = P(2,7) + PS80;
//    const float PS82 = PS81*dt;
//    const float PS83 = P(5,8)*dt;
//    const float PS84 = P(2,8) + PS83;
//    const float PS85 = PS84*dt;
//    const float PS86 = P(5,15)*dt;
//    const float PS87 = P(2,15) + PS86;
//    const float PS88 = P(2,9) + P(5,9)*dt;
//    const float PS89 = P(2,10) + P(5,10)*dt;
//    const float PS90 = P(2,11) + P(5,11)*dt;
//    const float PS91 = R00*dt;
//    const float PS92 = P(12,12)*PS91;
//    const float PS93 = R01*dt;
//    const float PS94 = P(12,13)*PS93;
//    const float PS95 = R02*dt;
//    const float PS96 = P(12,14)*PS95;
//    const float PS97 = PS16*dt;
//    const float PS98 = P(6,12)*PS97;
//    const float PS99 = PS22*dt;
//    const float PS100 = P(7,12)*PS99;
//    const float PS101 = PS28*dt;
//    const float PS102 = P(8,12)*PS101;
//    const float PS103 = -P(3,12) + PS100 + PS102 + PS92 + PS94 + PS96 + PS98;
//    const float PS104 = P(12,13)*PS91;
//    const float PS105 = P(13,13)*PS93;
//    const float PS106 = P(13,14)*PS95;
//    const float PS107 = P(6,13)*PS97;
//    const float PS108 = P(7,13)*PS99;
//    const float PS109 = P(8,13)*PS101;
//    const float PS110 = -P(3,13) + PS104 + PS105 + PS106 + PS107 + PS108 + PS109;
//    const float PS111 = P(12,14)*PS91;
//    const float PS112 = P(13,14)*PS93;
//    const float PS113 = P(14,14)*PS95;
//    const float PS114 = P(6,14)*PS97;
//    const float PS115 = P(7,14)*PS99;
//    const float PS116 = P(8,14)*PS101;
//    const float PS117 = -P(3,14) + PS111 + PS112 + PS113 + PS114 + PS115 + PS116;
//    const float PS118 = -P(3,8) + P(6,8)*PS97 + P(7,8)*PS99 + P(8,12)*PS91 + P(8,13)*PS93 + P(8,14)*PS95 + P(8,8)*PS101;
//    const float PS119 = -P(3,7) + P(6,7)*PS97 + P(7,12)*PS91 + P(7,13)*PS93 + P(7,14)*PS95 + P(7,7)*PS99 + P(7,8)*PS101;
//    const float PS120 = -P(3,6) + P(6,12)*PS91 + P(6,13)*PS93 + P(6,14)*PS95 + P(6,6)*PS97 + P(6,7)*PS99 + P(6,8)*PS101;
//    const float PS121 = PS103*dt;
//    const float PS122 = PS110*dt;
//    const float PS123 = PS117*dt;
//    const float PS124 = PS118*dt;
//    const float PS125 = PS119*dt;
//    const float PS126 = PS120*dt;
//    const float PS127 = P(12,15)*dt;
//    const float PS128 = PS127*R00;
//    const float PS129 = P(13,15)*dt;
//    const float PS130 = PS129*R01;
//    const float PS131 = P(14,15)*dt;
//    const float PS132 = PS131*R02;
//    const float PS133 = P(6,15)*dt;
//    const float PS134 = PS133*PS16;
//    const float PS135 = P(7,15)*dt;
//    const float PS136 = PS135*PS22;
//    const float PS137 = P(8,15)*dt;
//    const float PS138 = PS137*PS28;
//    const float PS139 = P(9,12)*dt;
//    const float PS140 = PS139*R00;
//    const float PS141 = P(9,13)*dt;
//    const float PS142 = PS141*R01;
//    const float PS143 = P(9,14)*dt;
//    const float PS144 = PS143*R02;
//    const float PS145 = P(6,9)*dt;
//    const float PS146 = PS145*PS16;
//    const float PS147 = P(7,9)*dt;
//    const float PS148 = PS147*PS22;
//    const float PS149 = P(8,9)*dt;
//    const float PS150 = PS149*PS28;
//    const float PS151 = P(10,12)*dt;
//    const float PS152 = PS151*R00;
//    const float PS153 = P(10,13)*dt;
//    const float PS154 = PS153*R01;
//    const float PS155 = P(10,14)*dt;
//    const float PS156 = PS155*R02;
//    const float PS157 = P(6,10)*dt;
//    const float PS158 = PS157*PS16;
//    const float PS159 = P(7,10)*dt;
//    const float PS160 = PS159*PS22;
//    const float PS161 = P(8,10)*dt;
//    const float PS162 = PS161*PS28;
//    const float PS163 = P(11,12)*dt;
//    const float PS164 = PS163*R00;
//    const float PS165 = P(11,13)*dt;
//    const float PS166 = PS165*R01;
//    const float PS167 = P(11,14)*dt;
//    const float PS168 = PS167*R02;
//    const float PS169 = P(6,11)*dt;
//    const float PS170 = PS16*PS169;
//    const float PS171 = P(7,11)*dt;
//    const float PS172 = PS171*PS22;
//    const float PS173 = P(8,11)*dt;
//    const float PS174 = PS173*PS28;
//    const float PS175 = R10*dt;
//    const float PS176 = P(12,12)*PS175;
//    const float PS177 = R11*dt;
//    const float PS178 = P(12,13)*PS177;
//    const float PS179 = R12*dt;
//    const float PS180 = P(12,14)*PS179;
//    const float PS181 = PS32*dt;
//    const float PS182 = P(6,12)*PS181;
//    const float PS183 = PS33*dt;
//    const float PS184 = P(7,12)*PS183;
//    const float PS185 = PS34*dt;
//    const float PS186 = P(8,12)*PS185;
//    const float PS187 = -P(4,12) + PS176 + PS178 + PS180 + PS182 + PS184 + PS186;
//    const float PS188 = P(12,13)*PS175;
//    const float PS189 = P(13,13)*PS177;
//    const float PS190 = P(13,14)*PS179;
//    const float PS191 = P(6,13)*PS181;
//    const float PS192 = P(7,13)*PS183;
//    const float PS193 = P(8,13)*PS185;
//    const float PS194 = -P(4,13) + PS188 + PS189 + PS190 + PS191 + PS192 + PS193;
//    const float PS195 = P(12,14)*PS175;
//    const float PS196 = P(13,14)*PS177;
//    const float PS197 = P(14,14)*PS179;
//    const float PS198 = P(6,14)*PS181;
//    const float PS199 = P(7,14)*PS183;
//    const float PS200 = P(8,14)*PS185;
//    const float PS201 = -P(4,14) + PS195 + PS196 + PS197 + PS198 + PS199 + PS200;
//    const float PS202 = -P(4,8) + P(6,8)*PS181 + P(7,8)*PS183 + P(8,12)*PS175 + P(8,13)*PS177 + P(8,14)*PS179 + P(8,8)*PS185;
//    const float PS203 = -P(4,7) + P(6,7)*PS181 + P(7,12)*PS175 + P(7,13)*PS177 + P(7,14)*PS179 + P(7,7)*PS183 + P(7,8)*PS185;
//    const float PS204 = -P(4,6) + P(6,12)*PS175 + P(6,13)*PS177 + P(6,14)*PS179 + P(6,6)*PS181 + P(6,7)*PS183 + P(6,8)*PS185;
//    const float PS205 = PS127*R10;
//    const float PS206 = PS129*R11;
//    const float PS207 = PS131*R12;
//    const float PS208 = PS133*PS32;
//    const float PS209 = PS135*PS33;
//    const float PS210 = PS137*PS34;
//    const float PS211 = R20*dt;
//    const float PS212 = R21*dt;
//    const float PS213 = R22*dt;
//    const float PS214 = PS38*dt;
//    const float PS215 = PS37*dt;
//    const float PS216 = PS36*dt;
//    const float PS217 = PS139*R10;
//    const float PS218 = PS141*R11;
//    const float PS219 = PS143*R12;
//    const float PS220 = PS145*PS32;
//    const float PS221 = PS147*PS33;
//    const float PS222 = PS149*PS34;
//    const float PS223 = PS151*R10;
//    const float PS224 = PS153*R11;
//    const float PS225 = PS155*R12;
//    const float PS226 = PS157*PS32;
//    const float PS227 = PS159*PS33;
//    const float PS228 = PS161*PS34;
//    const float PS229 = PS163*R10;
//    const float PS230 = PS165*R11;
//    const float PS231 = PS167*R12;
//    const float PS232 = PS169*PS32;
//    const float PS233 = PS171*PS33;
//    const float PS234 = PS173*PS34;
//    const float PS235 = P(15,15)*dt;
//    const float PS236 = PS127*R20;
//    const float PS237 = PS129*R21;
//    const float PS238 = PS131*R22;
//    const float PS239 = PS133*PS36;
//    const float PS240 = PS135*PS37;
//    const float PS241 = PS137*PS38;
//    const float PS242 = P(12,12)*PS211;
//    const float PS243 = P(12,13)*PS212;
//    const float PS244 = P(12,14)*PS213;
//    const float PS245 = P(6,12)*PS216;
//    const float PS246 = P(7,12)*PS215;
//    const float PS247 = P(8,12)*PS214;
//    const float PS248 = P(12,13)*PS211;
//    const float PS249 = P(13,13)*PS212;
//    const float PS250 = P(13,14)*PS213;
//    const float PS251 = P(6,13)*PS216;
//    const float PS252 = P(7,13)*PS215;
//    const float PS253 = P(8,13)*PS214;
//    const float PS254 = P(12,14)*PS211;
//    const float PS255 = P(13,14)*PS212;
//    const float PS256 = P(14,14)*PS213;
//    const float PS257 = P(6,14)*PS216;
//    const float PS258 = P(7,14)*PS215;
//    const float PS259 = P(8,14)*PS214;
//    const float PS260 = -P(5,8) + P(6,8)*PS216 + P(7,8)*PS215 + P(8,12)*PS211 + P(8,13)*PS212 + P(8,14)*PS213 + P(8,8)*PS214 - PS137;
//    const float PS261 = -P(5,7) + P(6,7)*PS216 + P(7,12)*PS211 + P(7,13)*PS212 + P(7,14)*PS213 + P(7,7)*PS215 + P(7,8)*PS214 - PS135;
//    const float PS262 = -P(5,6) + P(6,12)*PS211 + P(6,13)*PS212 + P(6,14)*PS213 + P(6,6)*PS216 + P(6,7)*PS215 + P(6,8)*PS214 - PS133;
//    const float PS263 = P(9,15)*dt;
//    const float PS264 = -PS263;
//    const float PS265 = PS139*R20;
//    const float PS266 = PS141*R21;
//    const float PS267 = PS143*R22;
//    const float PS268 = PS145*PS36;
//    const float PS269 = PS147*PS37;
//    const float PS270 = PS149*PS38;
//    const float PS271 = P(10,15)*dt;
//    const float PS272 = -PS271;
//    const float PS273 = PS151*R20;
//    const float PS274 = PS153*R21;
//    const float PS275 = PS155*R22;
//    const float PS276 = PS157*PS36;
//    const float PS277 = PS159*PS37;
//    const float PS278 = PS161*PS38;
//    const float PS279 = P(11,15)*dt;
//    const float PS280 = -PS279;
//    const float PS281 = PS163*R20;
//    const float PS282 = PS165*R21;
//    const float PS283 = PS167*R22;
//    const float PS284 = PS169*PS36;
//    const float PS285 = PS171*PS37;
//    const float PS286 = PS173*PS38;
//    const float PS287 = P(6,6)*T00 + P(6,7)*T01 + P(6,8)*T02 - PS145;
//    const float PS288 = P(6,7)*T00 + P(7,7)*T01 + P(7,8)*T02 - PS147;
//    const float PS289 = P(6,8)*T00 + P(7,8)*T01 + P(8,8)*T02 - PS149;
//    const float PS290 = P(6,9)*T00 + P(7,9)*T01 + P(8,9)*T02 - P(9,9)*dt;
//    const float PS291 = -P(9,10)*dt;
//    const float PS292 = P(6,10)*T00 + P(7,10)*T01 + P(8,10)*T02 + PS291;
//    const float PS293 = -P(9,11)*dt;
//    const float PS294 = P(6,11)*T00 + P(7,11)*T01 + P(8,11)*T02 + PS293;
//    const float PS295 = P(6,6)*T10 + P(6,7)*T11 + P(6,8)*T12 - PS157;
//    const float PS296 = P(6,7)*T10 + P(7,7)*T11 + P(7,8)*T12 - PS159;
//    const float PS297 = P(6,8)*T10 + P(7,8)*T11 + P(8,8)*T12 - PS161;
//    const float PS298 = -P(10,10)*dt + P(6,10)*T10 + P(7,10)*T11 + P(8,10)*T12;
//    const float PS299 = -P(10,11)*dt;
//    const float PS300 = P(6,11)*T10 + P(7,11)*T11 + P(8,11)*T12 + PS299;
//    const float PS301 = -P(11,11)*dt + P(6,11)*T20 + P(7,11)*T21 + P(8,11)*T22;
//
//    vector<vector<float>> nextP(16, vector<float>(16, 0.f));
//    nextP[0][0] = P(0,0) + P(0,3)*dt + PS0*dt;
//    nextP[0][1] = P(0,1) + P(1,3)*dt + PS2*dt;
//    nextP[1][1] = P(1,1) + P(1,4)*dt + PS42*dt;
//    nextP[0][2] = P(0,2) + P(2,3)*dt + PS4*dt;
//    nextP[1][2] = P(1,2) + P(2,4)*dt + PS44*dt;
//    nextP[2][2] = P(2,2) + P(2,5)*dt + PS67*dt;
//    nextP[0][3] = PS0 - PS10*R01 - PS13*R02 - PS16*PS19 - PS22*PS25 - PS28*PS31 - PS7*R00;
//    nextP[1][3] = P(1,3) + PS1 - PS16*PS56 - PS22*PS59 - PS28*PS62 - PS47*R00 - PS50*R01 - PS53*R02;
//    nextP[2][3] = P(2,3) - PS16*PS79 - PS22*PS82 - PS28*PS85 + PS3 - PS70*R00 - PS73*R01 - PS76*R02;
//    nextP[3][3] = P(3,3) + PS101*PS118 + PS103*PS91 - PS11*R02 + PS110*PS93 + PS117*PS95 + PS119*PS99 + PS120*PS97 - PS16*PS17 - PS22*PS23 - PS28*PS29 - PS5*R00 - PS8*R01;
//    nextP[0][4] = -PS10*R11 - PS13*R12 - PS19*PS32 + PS2 - PS25*PS33 - PS31*PS34 - PS7*R10;
//    nextP[1][4] = -PS32*PS56 - PS33*PS59 - PS34*PS62 + PS42 - PS47*R10 - PS50*R11 - PS53*R12;
//    nextP[2][4] = P(2,4) - PS32*PS79 - PS33*PS82 - PS34*PS85 + PS43 - PS70*R10 - PS73*R11 - PS76*R12;
//    nextP[3][4] = P(3,4) + PS121*R10 + PS122*R11 + PS123*R12 + PS124*PS34 + PS125*PS33 + PS126*PS32 - PS16*PS54 - PS22*PS57 - PS28*PS60 - PS45*R00 - PS48*R01 - PS51*R02;
//    nextP[4][4] = P(4,4) + PS175*PS187 + PS177*PS194 + PS179*PS201 + PS181*PS204 + PS183*PS203 + PS185*PS202 - PS32*PS54 - PS33*PS57 - PS34*PS60 - PS45*R10 - PS48*R11 - PS51*R12;
//    nextP[0][5] = -PS10*R21 - PS13*R22 - PS19*PS36 - PS25*PS37 - PS31*PS38 + PS35*dt + PS4 - PS7*R20;
//    nextP[1][5] = -PS36*PS56 - PS37*PS59 - PS38*PS62 + PS44 - PS47*R20 - PS50*R21 - PS53*R22 + PS63*dt;
//    nextP[2][5] = -PS36*PS79 - PS37*PS82 - PS38*PS85 + PS67 - PS70*R20 - PS73*R21 - PS76*R22 + PS87*dt;
//    nextP[3][5] = P(3,5) + PS121*R20 + PS122*R21 + PS123*R22 + PS124*PS38 + PS125*PS37 + PS126*PS36 - PS16*PS77 - PS22*PS80 - PS28*PS83 - PS68*R00 - PS71*R01 - PS74*R02 - dt*(-P(3,15) + PS128 + PS130 + PS132 + PS134 + PS136 + PS138);
//    nextP[4][5] = P(4,5) + PS187*PS211 + PS194*PS212 + PS201*PS213 + PS202*PS214 + PS203*PS215 + PS204*PS216 - PS32*PS77 - PS33*PS80 - PS34*PS83 - PS68*R10 - PS71*R11 - PS74*R12 - dt*(-P(4,15) + PS205 + PS206 + PS207 + PS208 + PS209 + PS210);
//    nextP[5][5] = P(5,5) + PS211*(-P(5,12) - PS127 + PS242 + PS243 + PS244 + PS245 + PS246 + PS247) + PS212*(-P(5,13) - PS129 + PS248 + PS249 + PS250 + PS251 + PS252 + PS253) + PS213*(-P(5,14) - PS131 + PS254 + PS255 + PS256 + PS257 + PS258 + PS259) + PS214*PS260 + PS215*PS261 + PS216*PS262 - PS36*PS77 - PS37*PS80 - PS38*PS83 - PS68*R20 - PS71*R21 - PS74*R22 + PS86 - dt*(-P(5,15) - PS235 + PS236 + PS237 + PS238 + PS239 + PS240 + PS241);
//    nextP[0][6] = PS18*T00 + PS24*T01 + PS30*T02 - PS39*dt;
//    nextP[1][6] = PS55*T00 + PS58*T01 + PS61*T02 - PS64*dt;
//    nextP[2][6] = PS78*T00 + PS81*T01 + PS84*T02 - PS88*dt;
//    nextP[3][6] = -PS118*T02 - PS119*T01 - PS120*T00 + dt*(-P(3,9) + PS140 + PS142 + PS144 + PS146 + PS148 + PS150);
//    nextP[4][6] = -PS202*T02 - PS203*T01 - PS204*T00 + dt*(-P(4,9) + PS217 + PS218 + PS219 + PS220 + PS221 + PS222);
//    nextP[5][6] = -PS260*T02 - PS261*T01 - PS262*T00 + dt*(-P(5,9) + PS264 + PS265 + PS266 + PS267 + PS268 + PS269 + PS270);
//    nextP[6][6] = PS287*T00 + PS288*T01 + PS289*T02 - PS290*dt;
//    nextP[0][7] = PS18*T10 + PS24*T11 + PS30*T12 - PS40*dt;
//    nextP[1][7] = PS55*T10 + PS58*T11 + PS61*T12 - PS65*dt;
//    nextP[2][7] = PS78*T10 + PS81*T11 + PS84*T12 - PS89*dt;
//    nextP[3][7] = -PS118*T12 - PS119*T11 - PS120*T10 + dt*(-P(3,10) + PS152 + PS154 + PS156 + PS158 + PS160 + PS162);
//    nextP[4][7] = -PS202*T12 - PS203*T11 - PS204*T10 + dt*(-P(4,10) + PS223 + PS224 + PS225 + PS226 + PS227 + PS228);
//    nextP[5][7] = -PS260*T12 - PS261*T11 - PS262*T10 + dt*(-P(5,10) + PS272 + PS273 + PS274 + PS275 + PS276 + PS277 + PS278);
//    nextP[6][7] = PS287*T10 + PS288*T11 + PS289*T12 - PS292*dt;
//    nextP[7][7] = PS295*T10 + PS296*T11 + PS297*T12 - PS298*dt;
//    nextP[0][8] = PS18*T20 + PS24*T21 + PS30*T22 - PS41*dt;
//    nextP[1][8] = PS55*T20 + PS58*T21 + PS61*T22 - PS66*dt;
//    nextP[2][8] = PS78*T20 + PS81*T21 + PS84*T22 - PS90*dt;
//    nextP[3][8] = -PS118*T22 - PS119*T21 - PS120*T20 + dt*(-P(3,11) + PS164 + PS166 + PS168 + PS170 + PS172 + PS174);
//    nextP[4][8] = -PS202*T22 - PS203*T21 - PS204*T20 + dt*(-P(4,11) + PS229 + PS230 + PS231 + PS232 + PS233 + PS234);
//    nextP[5][8] = -PS260*T22 - PS261*T21 - PS262*T20 + dt*(-P(5,11) + PS280 + PS281 + PS282 + PS283 + PS284 + PS285 + PS286);
//    nextP[6][8] = PS287*T20 + PS288*T21 + PS289*T22 - PS294*dt;
//    nextP[7][8] = PS295*T20 + PS296*T21 + PS297*T22 - PS300*dt;
//    nextP[8][8] = -PS301*dt + T20*(P(6,6)*T20 + P(6,7)*T21 + P(6,8)*T22 - PS169) + T21*(P(6,7)*T20 + P(7,7)*T21 + P(7,8)*T22 - PS171) + T22*(P(6,8)*T20 + P(7,8)*T21 + P(8,8)*T22 - PS173);
//    nextP[0][9] = PS39;
//    nextP[1][9] = PS64;
//    nextP[2][9] = PS88;
//    nextP[3][9] = P(3,9) - PS140 - PS142 - PS144 - PS146 - PS148 - PS150;
//    nextP[4][9] = P(4,9) - PS217 - PS218 - PS219 - PS220 - PS221 - PS222;
//    nextP[5][9] = P(5,9) + PS263 - PS265 - PS266 - PS267 - PS268 - PS269 - PS270;
//    nextP[6][9] = PS290;
//    nextP[7][9] = P(6,9)*T10 + P(7,9)*T11 + P(8,9)*T12 + PS291;
//    nextP[8][9] = P(6,9)*T20 + P(7,9)*T21 + P(8,9)*T22 + PS293;
//    nextP[9][9] = P(9,9);
//    nextP[0][10] = PS40;
//    nextP[1][10] = PS65;
//    nextP[2][10] = PS89;
//    nextP[3][10] = P(3,10) - PS152 - PS154 - PS156 - PS158 - PS160 - PS162;
//    nextP[4][10] = P(4,10) - PS223 - PS224 - PS225 - PS226 - PS227 - PS228;
//    nextP[5][10] = P(5,10) + PS271 - PS273 - PS274 - PS275 - PS276 - PS277 - PS278;
//    nextP[6][10] = PS292;
//    nextP[7][10] = PS298;
//    nextP[8][10] = P(6,10)*T20 + P(7,10)*T21 + P(8,10)*T22 + PS299;
//    nextP[9][10] = P(9,10);
//    nextP[10][10] = P(10,10);
//    nextP[0][11] = PS41;
//    nextP[1][11] = PS66;
//    nextP[2][11] = PS90;
//    nextP[3][11] = P(3,11) - PS164 - PS166 - PS168 - PS170 - PS172 - PS174;
//    nextP[4][11] = P(4,11) - PS229 - PS230 - PS231 - PS232 - PS233 - PS234;
//    nextP[5][11] = P(5,11) + PS279 - PS281 - PS282 - PS283 - PS284 - PS285 - PS286;
//    nextP[6][11] = PS294;
//    nextP[7][11] = PS300;
//    nextP[8][11] = PS301;
//    nextP[9][11] = P(9,11);
//    nextP[10][11] = P(10,11);
//    nextP[11][11] = P(11,11);
//    nextP[0][12] = PS6;
//    nextP[1][12] = PS46;
//    nextP[2][12] = PS69;
//    nextP[3][12] = P(3,12) - PS100 - PS102 - PS92 - PS94 - PS96 - PS98;
//    nextP[4][12] = P(4,12) - PS176 - PS178 - PS180 - PS182 - PS184 - PS186;
//    nextP[5][12] = P(5,12) + PS127 - PS242 - PS243 - PS244 - PS245 - PS246 - PS247;
//    nextP[6][12] = P(6,12)*T00 + P(7,12)*T01 + P(8,12)*T02 - PS139;
//    nextP[7][12] = P(6,12)*T10 + P(7,12)*T11 + P(8,12)*T12 - PS151;
//    nextP[8][12] = P(6,12)*T20 + P(7,12)*T21 + P(8,12)*T22 - PS163;
//    nextP[9][12] = P(9,12);
//    nextP[10][12] = P(10,12);
//    nextP[11][12] = P(11,12);
//    nextP[12][12] = P(12,12);
//    nextP[0][13] = PS9;
//    nextP[1][13] = PS49;
//    nextP[2][13] = PS72;
//    nextP[3][13] = P(3,13) - PS104 - PS105 - PS106 - PS107 - PS108 - PS109;
//    nextP[4][13] = P(4,13) - PS188 - PS189 - PS190 - PS191 - PS192 - PS193;
//    nextP[5][13] = P(5,13) + PS129 - PS248 - PS249 - PS250 - PS251 - PS252 - PS253;
//    nextP[6][13] = P(6,13)*T00 + P(7,13)*T01 + P(8,13)*T02 - PS141;
//    nextP[7][13] = P(6,13)*T10 + P(7,13)*T11 + P(8,13)*T12 - PS153;
//    nextP[8][13] = P(6,13)*T20 + P(7,13)*T21 + P(8,13)*T22 - PS165;
//    nextP[9][13] = P(9,13);
//    nextP[10][13] = P(10,13);
//    nextP[11][13] = P(11,13);
//    nextP[12][13] = P(12,13);
//    nextP[13][13] = P(13,13);
//    nextP[0][14] = PS12;
//    nextP[1][14] = PS52;
//    nextP[2][14] = PS75;
//    nextP[3][14] = P(3,14) - PS111 - PS112 - PS113 - PS114 - PS115 - PS116;
//    nextP[4][14] = P(4,14) - PS195 - PS196 - PS197 - PS198 - PS199 - PS200;
//    nextP[5][14] = P(5,14) + PS131 - PS254 - PS255 - PS256 - PS257 - PS258 - PS259;
//    nextP[6][14] = P(6,14)*T00 + P(7,14)*T01 + P(8,14)*T02 - PS143;
//    nextP[7][14] = P(6,14)*T10 + P(7,14)*T11 + P(8,14)*T12 - PS155;
//    nextP[8][14] = P(6,14)*T20 + P(7,14)*T21 + P(8,14)*T22 - PS167;
//    nextP[9][14] = P(9,14);
//    nextP[10][14] = P(10,14);
//    nextP[11][14] = P(11,14);
//    nextP[12][14] = P(12,14);
//    nextP[13][14] = P(13,14);
//    nextP[14][14] = P(14,14);
//    nextP[0][15] = PS35;
//    nextP[1][15] = PS63;
//    nextP[2][15] = PS87;
//    nextP[3][15] = P(3,15) - PS128 - PS130 - PS132 - PS134 - PS136 - PS138;
//    nextP[4][15] = P(4,15) - PS205 - PS206 - PS207 - PS208 - PS209 - PS210;
//    nextP[5][15] = P(5,15) + PS235 - PS236 - PS237 - PS238 - PS239 - PS240 - PS241;
//    nextP[6][15] = P(6,15)*T00 + P(7,15)*T01 + P(8,15)*T02 + PS264;
//    nextP[7][15] = P(6,15)*T10 + P(7,15)*T11 + P(8,15)*T12 + PS272;
//    nextP[8][15] = P(6,15)*T20 + P(7,15)*T21 + P(8,15)*T22 + PS280;
//    nextP[9][15] = P(9,15);
//    nextP[10][15] = P(10,15);
//    nextP[11][15] = P(11,15);
//    nextP[12][15] = P(12,15);
//    nextP[13][15] = P(13,15);
//    nextP[14][15] = P(14,15);
//    nextP[15][15] = P(15,15);
//
//    _cov = nextP;
//}

int main() {
    Eigen::MatrixXf matrix(3, 3);
    matrix << 1, 2, 3,
              4, 5, 6,
              7, 8, 9;

    cout << matrix << endl;

    ESKF eskf(0.001f);
    ESKF eskf1(0.001f);

    vector<float> w(3, 0.f), a(3, 0.f);
    w[0] = 1.f;
    w[1] = 2.f;
    w[2] = 3.f;
    a[0] = -10.f;
    a[1] = -7.f;
    a[2] = -4.f;

    eskf.predict_covariance(w, a);
//    eskf.predict_covariance(w, a);
//    eskf.predict_covariance(w, a);
//    eskf.predict_covariance(w, a);
//    eskf.predict_covariance(w, a);

//    eskf1.true_covariance(w, a);
//    eskf1.true_covariance(w, a);
//    eskf1.true_covariance(w, a);
//    eskf1.true_covariance(w, a);
//    eskf1.true_covariance(w, a);

    const vector<vector<float>> & s = eskf.get_covariance_matrix();
    const vector<vector<float>> & s1 = eskf1.get_covariance_matrix();
    const vector<float> state = eskf.get_state();
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            cout << s[i][j] << ", ";
        }
        cout << endl;
    }
//    for (int i = 0; i < 13; ++i) {
//        cout << state[i] << endl;
//    }

//    eskf.true_covariance(w, a);

//    time_t t1, t2;
//    t1 = clock();
//    for (int i = 0; i < 10000; ++i) {
//        eskf.predict_covariance(w, a);
////        eskf.true_covariance(w, a);
//    }
//    t2 = clock();
//    cout << t2 - t1 << endl;

    Eigen::VectorXf x = Eigen::VectorXf::Zero(16);
    cout << x << endl;
    for (int i = 0; i < x.size(); ++i) {
        x(i) = i;
    }
    Eigen::Vector3f y = x.segment(3, 3);
    x.segment(3, 3) *= 2;
    cout << x << endl;
    cout << y << endl;

    return 0;
}
