//
// Created by Cain on 2022/11/11.
//

#ifndef NAVIGATION_ESKF_H
#define NAVIGATION_ESKF_H

#include <Eigen/Dense>
#include <cmath>
#include <array>

namespace eskf {
    using namespace std;
    using namespace Eigen;

    class ESKF {
    public:
        ESKF(float dt, float dl_x, float dl_y, float dl_z, 
             float dr_x, float dr_y, float dr_z, float g=9.8f) 
             : _dt(dt), _dt2(dt * dt), _dl{dl_x, dl_y, dl_z}, _dr{dr_x, dr_y, dr_z} {
            _rot[0][0] = _rot[1][1] = _rot[2][2] = 1.f;
            for (unsigned char i = 0; i < 16; ++i) {
                _cov[i][i] = 1.f;
            }
            _state[12] = g;
        }

        void predict_covariance(const array<float, 3> &w, const array<float, 3> &a);
        void predict_state(const array<float, 3> &w, const array<float, 3> &a);

        void fuse_position(const array<float, 3> &pos_left, const array<float, 3> &pos_right);
        void fuse_velocity(const array<float, 3> &vel_left, const array<float, 3> &vel_right);

        void fuse_position_left_only(const array<float, 3> pos);
        void fuse_velocity_left_only(const array<float, 3> vel);
        void fuse_position_right_only(const array<float, 3> pos);
        void fuse_velocity_right_only(const array<float, 3> vel);

        void fuse_position_z(float pos_z);
        void fuse_velocity_xy(const array<float, 2> vel_xy);

        void rotation_from_axis_angle(array<array<float, 3>, 3> &r, const array<float, 3> &a);

        const array<float, 13> &get_state() { return _state; };
        const array<float, 16> &get_error_state() { return _error_state; };
        const array<array<float, 16>, 16> &get_covariance_matrix() { return _cov; };
        const array<array<float, 3>, 3> &get_rotation_matrix() { return _rot; };
        const float get_dt() { return _dt; };

        void set_dt(float dt) { _dt = dt; _dt2 = dt * dt; };
        void set_g(float g) { _state[12] = g; };
        void set_left_antenna(float x, float y, float z) { _dl[0] = x, _dl[1] = y, _dl[2] = z; };
        void set_right_antenna(float x, float y, float z) { _dr[0] = x, _dr[1] = y, _dr[2] = z; };
        void set_accelerometer_standard_deviation(float std) { _q_cov[3] = _q_cov[4] = _q_cov[5] = std * std; };
        void set_gyroscope_standard_deviation(float std) { _q_cov[6] = _q_cov[7] = _q_cov[8] = std * std; };
        void set_accelerometer_drift_standard_deviation(float std) { _q_cov[9] = _q_cov[10] = _q_cov[11] = std * std; }
        void set_gyroscope_drift_standard_deviation(float std) { _q_cov[12] = _q_cov[13] = _q_cov[14] = std * std; }
        void set_gravity_standard_deviation(float std) { _q_cov[15] = std * std; };
        void set_processing_standard_deviation(float std) { for (unsigned char i = 0; i < _q_cov.size(); ++i) _q_cov[i] += std * std; };
        void set_position_standard_deviation(float std) { _r_cov[0] = _r_cov[1] = _r_cov[2] = _r_cov[6] = _r_cov[7] = _r_cov[8] = std * std; };
        void set_velocity_standard_deviation(float std) { _r_cov[3] = _r_cov[4] = _r_cov[5] = _r_cov[9] = _r_cov[10] = _r_cov[11] = std * std; };
        void set_barometer_standard_deviation(float std) { _r_cov[12] = std * std; };
        void set_optical_flow_standard_deviation(float std) { _r_cov[13] = _r_cov[14] = std * std; }

    private:
        float _dt;                      // Sample time of IMU
        float _dt2;                     // Square of sample time
        
        array<float, 3> _dl;     // Direction vector of GPS left antenna
        array<float, 3> _dr;     // Direction vector of GPS right antenna

        array<float, 16> _q_cov {};
        array<float, 15> _r_cov;

        array<float, 13> _state {};           // [p, v, bg, ba, g]
        array<array<float, 3>, 3> _rot {};     // Rotation matrix from body frame to world frame

        array<float, 16> _error_state {};     // [δp, δv, δθ, δbg, δba, δg]
        array<array<float, 16>, 16> _cov {};     // Covariance matrix of error state
    };
}

#endif //NAVIGATION_ESKF_H
