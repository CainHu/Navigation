//
// Created by Cain on 2022/12/9.
//

#include "geskf.h"
#include <cfloat>
#include <iostream>

namespace geskf {
    using namespace std;

    void GESKF::correct_state() {
        // state: [p, v, bg, ba, g], R
        // error_state : [δp, δv, δθ, δbg, δba, δg]

        /*
        p = p + δp
        v = v + δv
        bg = bg + δbg
        ba = ba + δba
        g = g + δg
        m = m + δm
        bm = bm + δbm
        */
        for (unsigned char i = 0; i < 3; ++i) {
            _p[i] += _error_state[i];
            _v[i] += _error_state[3 + i];
            _bg[i] += _error_state[9 + i];
            _ba[i] += _error_state[12 + i];
            _bm[i] += _error_state[19 + i];
        }
        _g += _error_state[15];
        _h += _error_state[16];
        _dec[0] += _error_state[17];
        _dec[1] += _error_state[18];
        _w[0] += _error_state[22];
        _w[1] += _error_state[23];
        
        // q = Exp(δθ) * q
        Quaternionf delta_q;
        array<float, 3> delta_theta = {_error_state[6], _error_state[7], _error_state[8]};
        quaternion_from_axis_angle(delta_q, delta_theta);
        const Quaternionf q = _q;
        _q = delta_q * q;
        _rot = q;

        // [δp, δv, δθ, δbg, δba, δg] = 0
        for (float &es : _error_state) {
            es = 0.f;
        }
    }

    void GESKF::correct_output_states(const eskf::ImuSample &imu_sample) {
        eskf::PreIntegralSample pre_sample;

        pre_sample.dt = 0.5f * (imu_sample.delta_vel_dt + imu_sample.delta_ang_dt);

        const Vector3f axis_angle = imu_sample.delta_ang - _bg * imu_sample.delta_ang_dt;
        rotation_from_axis_angle(pre_sample.dr, axis_angle);

        pre_sample.dv = imu_sample.delta_vel - _ba * imu_sample.delta_vel_dt;
        pre_sample.dp = 0.5f * pre_sample.dv * pre_sample.dt;

        pre_buffer.push(pre_sample);
        for (unsigned char i = 0; i < pre_buffer.get_length() - 1; ++i) {
            int index = pre_buffer.get_oldest_index() + i;
            index = (index == pre_buffer.get_length()) ? (index - pre_buffer.get_length()) : index;
            if (index == pre_buffer.get_newest_index()) {
                break;
            }

            const Vector3f drdv = pre_buffer[index].dr * pre_sample.dv;
            pre_buffer[index].dp += (pre_buffer[index].dv + 0.5f * drdv) * pre_sample.dt;
            pre_buffer[index].dv += drdv;
            pre_buffer[index].dr *= pre_sample.dr;
            pre_buffer[index].dt += pre_sample.dt;
        }

        const eskf::PreIntegralSample &pre_oldest = pre_buffer.get_oldest();
        _output_state.q = _q * pre_oldest.dr;
        _output_state.q.normalize();
        _output_state.r = _output_state.q;

        _output_state.v = _v + _rot * pre_oldest.dv;
        _output_state.v[2] += _g * pre_oldest.dt;

        _output_state.p = _p + _rot * pre_oldest.dp + _v * pre_oldest.dt;
        _output_state.p[2] += 0.5f * _g * pre_oldest.dt * pre_oldest.dt;

        _output_state.time_us = imu_sample.time_us;
    }
}