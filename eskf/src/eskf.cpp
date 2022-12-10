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
    reset_covariance_matrix(0, ESKF::DIM);
    regular_covariance_to_symmetric<ESKF::DIM>(0);
    reset_accmulator_cov();
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
    // h' = h
    // dec' = dec
    // bm' = bm
    // w' = w

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

unsigned char ESKF::conservative_posteriori_estimate(const array<float, ESKF::DIM> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
    /*
    K = P * H' * (H * P * H' + R)^-1
    P = P - K * H * P

    e = y - h = y - (v + R * (w - bg)^ * dis)
    x = x + K * (y - h)
    */
    
    // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
    if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
        return 1;
    }

    // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
    for (unsigned char i = 0; i < 16; ++i) {
        if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
            return 2;    
        }
    }

    // K = P * H' * (H * P * H' + R)^-1
    array<float, ESKF::DIM> K {};
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

    if (_control_status.flags.mag) {
        // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
        if (HP[16] * HP[16] > HPHT_plus_R * _cov[16][16]) {
            return 2;    
        }

        // K = P * H' * (H * P * H' + R)^-1
        K[16] = HP[16] / HPHT_plus_R;

        // x = x + K * e
        _error_state[16] += K[16] * obs_error;

        // P = P - K * H * P
        for (unsigned char j = 0; j <= 16; ++j) {
            _cov[j][16] -= K[j] * HP[16];
        }
    }

    if (_control_status.flags.dec) {
        // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
        for (unsigned char i = 17; i < 19; ++i) {
            if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
                return 2;    
            }
        }

        // K = P * H' * (H * P * H' + R)^-1
        for (unsigned char i = 17; i < 19; ++i) {
            K[i] = HP[i] / HPHT_plus_R;
        }

        // x = x + K * e
        for (unsigned char i = 17; i < 19; ++i) {
            _error_state[i] += K[i] * obs_error;
        }

        // P = P - K * H * P
        for (unsigned char i = 17; i < 19; ++i) {
            for (unsigned char j = 0; j <= i; ++j) {
                _cov[j][i] -= K[j] * HP[i];
            }
        }
    }

    if (_control_status.flags.mag_bias) {
        // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
        for (unsigned char i = 19; i < 22; ++i) {
            if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
                return 2;    
            }
        }

        // K = P * H' * (H * P * H' + R)^-1
        for (unsigned char i = 19; i < 22; ++i) {
            K[i] = HP[i] / HPHT_plus_R;
        }

        // x = x + K * e
        for (unsigned char i = 19; i < 22; ++i) {
            _error_state[i] += K[i] * obs_error;
        }

        // P = P - K * H * P
        for (unsigned char i = 19; i < 22; ++i) {
            for (unsigned char j = 0; j <= i; ++j) {
                _cov[j][i] -= K[j] * HP[i];
            }
        }
    }

    if (_control_status.flags.wind) {
        // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
        for (unsigned char i = 22; i < 24; ++i) {
            if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
                return 2;    
            }
        }

        // K = P * H' * (H * P * H' + R)^-1
        for (unsigned char i = 22; i < 24; ++i) {
            K[i] = HP[i] / HPHT_plus_R;
        }

        // x = x + K * e
        for (unsigned char i = 22; i < 24; ++i) {
            _error_state[i] += K[i] * obs_error;
        }

        // P = P - K * H * P
        for (unsigned char i = 22; i < 24; ++i) {
            for (unsigned char j = 0; j <= i; ++j) {
                _cov[j][i] -= K[j] * HP[i];
            }
        }
    }
    
    return 0;
}