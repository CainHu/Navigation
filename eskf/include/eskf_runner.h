//
// Created by Cain on 2022/11/11.
//

#ifndef NAVIGATION_ESKF_ESKF_RUNNER_H
#define NAVIGATION_ESKF_ESKF_RUNNER_H

#include "eskf.h"
#include "common.h"
#include "utils.h"

namespace eskf {
    template <unsigned char DELAYS>
    class ESKFRunner {
    public:
        ESKFRunner(ESKF &eskf) : _eskf(eskf) {

        };

        void update();

        void set_imu_data(const Vector3f &w, const Vector3f &a, const unsigned long &time_us);
        void set_gps_data(const Vector3f &pl, const Vector3f &vl, const Vector3f &pr, const Vector3f &vr, const unsigned long &time_us);
        void set_magnet_data(const Vector3f &m, const unsigned long &time_us);
        void set_declination_data(const Vector2f &d, const unsigned long &time_us);
        void set_airspeed_data(const float &true_airspeed, const float &eas2tas, const unsigned long &time_us);
    
        const OutputSample &get_output_state() { return _output_state; };

        // 卡尔曼滤波器的参数
        parameters _params {};

    private:
        void predict_state_from_delay(const ImuSample &sample);

        // 卡尔曼滤波器
        ESKF &_eskf;

        // imu
        bool _imu_updated {false};  ///< imu是否更新
        ImuSample _imu_sample_last {};  ///< 上一次的imu采样

        // imu采样时间
        float _dt_imu {};

        // 用于滞后补偿
        OutputSample _output_state {};

        Queue<ImuSample, DELAYS> _imu_buffer;
        Queue<GpsSample, DELAYS> _gps_buffer;
        Queue<MagSample, DELAYS> _mag_buffer;
        Queue<DecSample, DELAYS> _dec_buffer;
        Queue<AirspeedSample, DELAYS> _airspeed_buffer;
        Queue<PreIntegralSample, DELAYS> _pre_buffer;

    };

    template<unsigned char DELAYS>
    void ESKFRunner<DELAYS>::update() {
        if (_imu_updated) {
            const ImuSample &imu_sample = _imu_buffer.get_oldest();
            GpsSample gps_sample;
            MagSample mag_sample;
            DecSample dec_sample;
            AirspeedSample airspeed_sample;

            // 计算角速度与加速度
            const Vector3f w = imu_sample.delta_ang / imu_sample.delta_ang_dt;
            const Vector3f a = imu_sample.delta_vel / imu_sample.delta_vel_dt;

            // 计算时间间隔
            // _dt_imu = 1e-6f * float(imu_sample.time_us - _imu_sample_delayed.time_us);
            _dt_imu = 0.5f * (imu_sample.delta_ang_dt + imu_sample.delta_ang_dt);

            // 对时间间隔进行低通滤波
//            _eskf.set_dt(0.99f * _eskf._dt + 0.01f * _dt_imu);

            // 在滞后时刻进行先验估计
            _eskf.predict_state(w, a);
            _eskf.predict_covariance(w, a);

            // 在滞后时刻进行后验估计
            if (_gps_buffer.pop_first_older_than(imu_sample.time_us, gps_sample)){
                _eskf.fuse_position(gps_sample.pos_l, w, a, _params.d_gps_left, _params.noise_std_pos_rtk, _params.gate_gps_pos);
                _eskf.fuse_velocity(gps_sample.vel_l, w, a, _params.d_gps_left, _params.noise_std_vel_gps, _params.gate_gps_vel);
                _eskf.fuse_position(gps_sample.pos_r, w, a, _params.d_gps_right, _params.noise_std_pos_rtk, _params.gate_gps_pos);
                _eskf.fuse_velocity(gps_sample.vel_r, w, a, _params.d_gps_right, _params.noise_std_vel_gps, _params.gate_gps_vel);
            }

            if (_mag_buffer.pop_first_older_than(imu_sample.time_us, mag_sample)) {
                _eskf.fuse_magnet(mag_sample.mag, w, a, _params.noise_std_mag, _params.gate_mag);
            }

            if (_dec_buffer.pop_first_older_than(imu_sample.time_us, dec_sample)) {
                _eskf.fuse_declination(dec_sample.dec, w, a, _params.noise_std_dec, _params.gate_dec);
            }

            //TODO: 速度融合还没完成
            // if (dec_buffer.pop_first_older_than(imu_sample.time_us, airspeed_sample)) {

            // }

            // 修正滞后时刻的状态与协方差
            _eskf.correct_state();
            // _eskf.correct_covariance();
            
            // 利用滞后时刻的后验估计值对当前时刻进行预测
            predict_state_from_delay(_imu_buffer.get_newest());

            _imu_updated = false;
        }
    }

    template<unsigned char DELAYS>
    void ESKFRunner<DELAYS>::predict_state_from_delay(const ImuSample &sample) {
//        PreIntegralSample pre_sample;
//
//        pre_sample.dt = 0.5f * (sample.delta_vel_dt + sample.delta_ang_dt);
//
//        const Vector3f axis_angle = sample.delta_ang - _eskf._bg * sample.delta_ang_dt;
//        rotation_from_axis_angle(pre_sample.dr, axis_angle);
//
//        pre_sample.dv = sample.delta_vel - _eskf._ba * sample.delta_vel_dt;
//        pre_sample.dp = 0.5f * pre_sample.dv * pre_sample.dt;
//
//        _pre_buffer.push(pre_sample);
//        for (unsigned char i = 0; i < _pre_buffer.get_length() - 1; ++i) {
//            int index = _pre_buffer.get_oldest_index() + i;
//            index = (index == _pre_buffer.get_length()) ? (index - _pre_buffer.get_length()) : index;
//            if (index == _pre_buffer.get_newest_index()) {
//                break;
//            }
//
//            const Vector3f drdv = _pre_buffer[index].dr * pre_sample.dv;
//            _pre_buffer[index].dp += (_pre_buffer[index].dv + 0.5f * drdv) * pre_sample.dt;
//            _pre_buffer[index].dv += drdv;
//            _pre_buffer[index].dr *= pre_sample.dr;
//            _pre_buffer[index].dt += pre_sample.dt;
//        }

        const eskf::PreIntegralSample &pre_oldest = _pre_buffer.get_oldest();
        _output_state.q = _eskf._q * pre_oldest.dr;
        _output_state.q.normalize();
        _output_state.r = _output_state.q;

        _output_state.v = _eskf._v + _eskf._rot * pre_oldest.dv;
        _output_state.v[2] += _eskf._g * pre_oldest.dt;

        _output_state.p = _eskf._p + _eskf._rot * pre_oldest.dp + _eskf._v * pre_oldest.dt;
        _output_state.p[2] += 0.5f * _eskf._g * pre_oldest.dt * pre_oldest.dt;

        _output_state.time_us = sample.time_us;
    }

    template<unsigned char DELAYS>
    void ESKFRunner<DELAYS>::set_imu_data(const Vector3f &w, const Vector3f &a, const unsigned long &time_us) {
        const float dt = 1e-6f * float(time_us - _imu_sample_last.time_us);

        if (dt <= 0.f) {
            _imu_updated = false;
            return;
        }

        Vector3f delta_ang = w * dt;
        Vector3f delta_vel = a * dt;

        // 划桨补偿
        _imu_sample_last.delta_vel = delta_vel + 0.5f * delta_ang.cross(delta_vel) + 1.f/12.f * (_imu_sample_last.delta_ang.cross(delta_vel) + _imu_sample_last.delta_vel.cross(delta_ang));
        _imu_sample_last.delta_vel_dt = dt;

        // 圆锥补偿
        _imu_sample_last.delta_ang = delta_ang + 1.f/12.f * _imu_sample_last.delta_ang.cross(delta_ang);
        _imu_sample_last.delta_ang_dt = dt;
        
        _imu_sample_last.delta_vel_clipping = {false, false, false};
        _imu_sample_last.time_us = time_us;

        _imu_buffer.push(_imu_sample_last);

        // 临时预积分
        PreIntegralSample pre_sample;

        pre_sample.dt = 0.5f * (_imu_sample_last.delta_vel_dt + _imu_sample_last.delta_ang_dt);

        const Vector3f axis_angle = _imu_sample_last.delta_ang - _eskf._bg * _imu_sample_last.delta_ang_dt;
        rotation_from_axis_angle(pre_sample.dr, axis_angle);

        pre_sample.dv = _imu_sample_last.delta_vel - _eskf._ba * _imu_sample_last.delta_vel_dt;
        pre_sample.dp = 0.5f * pre_sample.dv * pre_sample.dt;

        _pre_buffer.push(pre_sample);
        for (unsigned char i = 0; i < _pre_buffer.get_length() - 1; ++i) {
            int index = _pre_buffer.get_oldest_index() + i;
            index = (index < _pre_buffer.get_length()) ? index : (index - _pre_buffer.get_length());
            if (index == _pre_buffer.get_newest_index()) {
                break;
            }

            const Vector3f drdv = _pre_buffer[index].dr * pre_sample.dv;
            _pre_buffer[index].dp += (_pre_buffer[index].dv + 0.5f * drdv) * pre_sample.dt;
            _pre_buffer[index].dv += drdv;
            _pre_buffer[index].dr *= pre_sample.dr;
            _pre_buffer[index].dt += pre_sample.dt;
        }

        _imu_updated = true;
    }

    template<unsigned char DELAYS>
    void ESKFRunner<DELAYS>::set_gps_data(const Vector3f &pl, const Vector3f &vl, const Vector3f &pr, const Vector3f &vr, const unsigned long &time_us) {
        GpsSample sample;
        sample.pos_l = pl;
        sample.vel_l = vl;
        sample.pos_r = pr;
        sample.vel_r = vr;
        sample.time_us = time_us;

        _gps_buffer.push(sample);
    }

    template<unsigned char DELAYS>
    void ESKFRunner<DELAYS>::set_magnet_data(const Vector3f &m, const unsigned long &time_us) {
        MagSample sample;
        sample.mag = m;
        sample.time_us = time_us;

        _mag_buffer.push(sample);
    }

    template<unsigned char DELAYS>
    void ESKFRunner<DELAYS>::set_declination_data(const Vector2f &d, const unsigned long &time_us) {
        DecSample sample;
        sample.dec = d;
        sample.time_us = time_us;

        _dec_buffer.push(sample);
    }

    template<unsigned char DELAYS>
    void ESKFRunner<DELAYS>::set_airspeed_data(const float &true_airspeed, const float &eas2tas, const unsigned long &time_us) {
        AirspeedSample sample;
        sample.true_airspeed = true_airspeed;
        sample.eas2tas = eas2tas;
        sample.time_us = time_us;

        _airspeed_buffer.push(sample);
    }
}

#endif //NAVIGATION_ESKF_ESKF_RUNNER_H