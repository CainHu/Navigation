//
// Created by Cain on 2022/11/11.
//

#ifndef NAVIGATION_ESKF_ESKF_RUNNER_H
#define NAVIGATION_ESKF_ESKF_RUNNER_H

#include "eskf.h"

namespace eskf {
    template <unsigned char DELAYS>
    class ESKFRunner {
    public:
        ESKFRunner(ESKF &eskf) : _eskf(eskf) {

        };

        void update(const ImuSample &sample);
    
        const OutputSample &get_output_state();

        
    private:
        void update_output_state();

        ESKF &_eskf;
        // 用于滞后补偿
        ImuSample _newest_high_rate_imu_sample {};
        ImuSample _imu_sample_delayed {};
        OutputSample _output_state {};

        Queue<ImuSample, DELAYS> imu_buffer;
        Queue<GpsSample, DELAYS> gps_buffer;
        Queue<MagSample, DELAYS> mag_buffer;
        Queue<DecSample, DELAYS> dec_buffer;
        Queue<AirspeedSample, DELAYS> airspeed_buffer;
        Queue<PreIntegralSample, DELAYS> pre_buffer;
    }

    template<unsigned char DELAYS>
    void ESKFRunner::ESKFR<DELAYS> update_output_state(const ImuSample &sample) {
        PreIntegralSample pre_sample;

        pre_sample.dt = 0.5f * (sample.delta_vel_dt + sample.delta_ang_dt);

        const Vector3f axis_angle = sample.delta_ang - _bg * sample.delta_ang_dt;
        rotation_from_axis_angle(pre_sample.dr, axis_angle);

        pre_sample.dv = sample.delta_vel - _ba * sample.delta_vel_dt;
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
        _output_state.q = _eskf.get_quaternion() * pre_oldest.dr;
        _output_state.q.normalize();
        _output_state.r = _output_state.q;

        _output_state.v = _eskf.get_velocity() + _eskf.get_rotation_matrix() * pre_oldest.dv;
        _output_state.v[2] += _eskf.get_gravity() * pre_oldest.dt;

        _output_state.p = _eskf.get_position() + _eskf.get_rotation_matrix() * pre_oldest.dp + _eskf.get_velocity() * pre_oldest.dt;
        _output_state.p[2] += 0.5f * _eskf.get_gravity() * pre_oldest.dt * pre_oldest.dt;

        _output_state.time_us = sample.time_us;
    }
}

#endif //NAVIGATION_ESKF_ESKF_RUNNER_H