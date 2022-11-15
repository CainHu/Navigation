//
// Created by Cain on 2022/11/11.
//

#ifndef NAVIGATION_ESKF_ESKF_H
#define NAVIGATION_ESKF_ESKF_H

#include <Eigen/Dense>
#include <cmath>
#include <array>

namespace eskf {
    using namespace std;
    using namespace Eigen;

    class ESKF {
    public:
        ESKF(float dt, const float g=9.8f) 
             : _g(g), _dt(dt), _dt2(dt * dt) {
            _rot[0][0] = _rot[1][1] = _rot[2][2] = 1.f;
            for (unsigned char i = 0; i < 16; ++i) {
                _cov[i][i] = _q_cov[i] * _dt2;
            }
            _state[12] = g;
        }

        void initialize();

        // Priori
        void predict_covariance(const array<float, 3> &w, const array<float, 3> &a);
        void predict_state(const array<float, 3> &w, const array<float, 3> &a);

        // Posteriori
        unsigned char fuse_position(const array<float, 3> &pos, const array<float, 3> &w, const array<float, 3> &a, const array<float, 3> &dis, const array<float, 3> &noise_std, const array<float, 3> &gate);
        unsigned char fuse_velocity(const array<float, 3> &vel, const array<float, 3> &w, const array<float, 3> &a, const array<float, 3> &dis, const array<float, 3> &noise_std, const array<float, 3> &gate);
        unsigned char conservative_posteriori_estimate(const array<float, 16> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate);
        void correct_state();

        void regular_covariance_to_symmetric(unsigned int start_index=0, unsigned int end_index=16);
        void rotation_from_axis_angle(array<array<float, 3>, 3> &r, const array<float, 3> &a) const;

        // Resetters
        void reset_rot(const float roll=0.f, const float pitch=0.f, const float yaw=0.f);
        void reset_rot(const array<float, 3> &axis_angle);
        void reset_rot(const array<float, 4> &quaternion);
        void reset_rot(const array<array<float, 3>, 3> &rotation_maxtrix);
        void reset_state();
        void reset_error_state();
        void reset_covariance_matrix(const unsigned int start_index, const unsigned int end_index, const array<float, 16> &diag_cov);
        void reset_accmulator_cov();

        // Getters
        const array<float, 13> &get_state() const { return _state; };
        const array<float, 16> &get_error_state() const { return _error_state; };
        const array<array<float, 16>, 16> &get_covariance_matrix() const { return _cov; };
        const array<array<float, 3>, 3> &get_rotation_matrix() const { return _rot; };
        const float get_dt() const { return _dt; };

        // Setters
        void set_dt(const float dt) { _dt = dt; _dt2 = dt * dt; };
        void set_g(const float g) { _state[12] = g; };
        void set_accelerometer_standard_deviation(const array<float, 3> std) { _q_cov[3] = std[0] * std[0]; _q_cov[4] = std[1] * std[1]; _q_cov[5] = std[2] * std[2]; };
        void set_gyroscope_standard_deviation(const array<float, 3> std) { _q_cov[6] = std[0] * std[0]; _q_cov[7] = std[1] * std[1]; _q_cov[8] = std[2] * std[2]; };
        void set_drift_saccelerometer_tandard_deviation(const array<float, 3> std) { _q_cov[9] = std[0] * std[0]; _q_cov[10] = std[1] * std[1]; _q_cov[11] = std[2] * std[2]; }
        void set_drift_gyroscope_standard_deviation(const array<float, 3> std) { _q_cov[12] = std[0] * std[0]; _q_cov[13] = std[1] * std[1]; _q_cov[14] = std[2] * std[2]; }
        void set_gravity_standard_deviation(const float std) { _q_cov[15] = std * std; };
        void set_processing_standard_deviation(const float std) { for (float &cov : _q_cov) cov += std * std; };

    private:
        const float _g;                 // gravity

        float _dt;                      // Sample time of IMU
        float _dt2;                     // Square of sample time

        array<float, 16> _q_cov {};

        array<float, 13> _state {};             // [p, v, bg, ba, g]
        array<array<float, 3>, 3> _rot {};      // Rotation matrix from body frame to world frame

        array<float, 16> _error_state {};       // [δp, δv, δθ, δbg, δba, δg]
        array<array<float, 16>, 16> _cov {};    // Covariance matrix of error state

        array<float, 16> _accumulator_cov {};   // Accumulator of covariance matrix

        float kahan_summation(float sum_previous, float input, float &accumulator);
    };
}

#endif //NAVIGATION_ESKF_H
