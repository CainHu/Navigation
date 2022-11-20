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
             : _g_init(g), _g(g), _dt(dt), _dt2(dt * dt) {
            _rot.setIdentity();
            reset_priori_covariance_matrix();
            reset_state();
            reset_error_state();
            reset_accmulator_cov();
            for (unsigned int i = 0; i < 16; ++i) {
                _cov[i][i] = 1.f;
            }
        }

        void initialize();

        // Priori
        void predict_state(const Vector3f &w, const Vector3f &a);
        virtual void predict_covariance(const Vector3f &w, const Vector3f &a) = 0;

        // Posteriori
        virtual unsigned char fuse_position(const Vector3f &pos, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) = 0;
        virtual unsigned char fuse_velocity(const Vector3f &vel, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) = 0;
        virtual void correct_state() = 0;

        // Resetters
        void reset_state();
        void reset_error_state();
        void reset_priori_covariance_matrix();
        void reset_covariance_matrix(const unsigned int start_index, const unsigned int end_index, const Matrix<float, 16, 1> &diag_cov);
        void reset_accmulator_cov();

        // Getters
        const Vector3f &get_position() const { return _p; };
        const Vector3f &get_velocity() const { return _v; };
        const Vector3f &get_drift_gyro() const { return _bg; };
        const Vector3f &get_drift_acc() const { return _ba; };
        const Matrix3f &get_rotation_matrix() const { return _rot; };
        const Quaternionf &get_quaternion() const { return _q; };
        const float get_gravity() const { return _g; };
        const float get_dt() const { return _dt; };
        const array<float, 16> &get_error_state() const { return _error_state; };
        const array<array<float, 16>, 16> &get_covariance_matrix() const { return _cov; };

        // Setters
        void set_dt(const float dt) { _dt = dt; _dt2 = dt * dt; };
        void set_gravity(const float g) { _g_init = g; _g = g; };
        void set_position(const Vector3f &p) { _p = p; };
        void set_velocity(const Vector3f &v) { _v = v; };
        void set_drift_gyro(const Vector3f &bg) { _bg = bg; };
        void set_drift_acc(const Vector3f &ba) { _ba = ba; };
        void set_attitude(const Matrix3f &r) { _rot = r; _q = _rot; };
        void set_attitude(const Quaternionf &q) { _q = q; _rot = _q; };
        void set_attitude(const AngleAxisf &a) { _q = a; _rot = _q;};    
        void set_accelerometer_standard_deviation(const Vector3f std) { _q_cov[3] = std[0] * std[0]; _q_cov[4] = std[1] * std[1]; _q_cov[5] = std[2] * std[2]; };
        void set_gyroscope_standard_deviation(const Vector3f std) { _q_cov[6] = std[0] * std[0]; _q_cov[7] = std[1] * std[1]; _q_cov[8] = std[2] * std[2]; };
        void set_drift_saccelerometer_tandard_deviation(const Vector3f std) { _q_cov[12] = std[0] * std[0]; _q_cov[13] = std[1] * std[1]; _q_cov[14] = std[2] * std[2]; }
        void set_drift_gyroscope_standard_deviation(const Vector3f std) { _q_cov[9] = std[0] * std[0]; _q_cov[10] = std[1] * std[1]; _q_cov[11] = std[2] * std[2]; }
        void set_gravity_standard_deviation(const float std) { _q_cov[15] = std * std; };
        void set_processing_standard_deviation(const float std) { for (unsigned char i = 0; i < 16; ++i) { _q_cov[i] += std * std; } };
        void set_priori_covariance_matrix(Matrix<float, 16, 1> &q) { _q_cov = q; };

    protected:
        float _g_init;                 // Initial gravity constant

        float _dt;                      // Sample time of IMU
        float _dt2;                     // Square of sample time

        Matrix<float, 16, 1> _q_cov {};         // Priori covariance matrix
  
        Vector3f _p, _v, _bg, _ba;              // Translational state: p, v, bg, ba
        float _g;                               // Gravity constant
        Matrix3f _rot;                          // Rotation matrix from body frame to world frame
        Quaternionf _q;                         // Quaternion matrix from body frame to world frame

        array<float, 16> _error_state {};       // [δp, δv, δθ, δbg, δba, δg]
        array<array<float, 16>, 16> _cov {};    // Covariance matrix of error state

        array<float, 16> _accumulator_cov {};   // Accumulator of covariance matrix

        // Posteriori
        unsigned char conservative_posteriori_estimate(const array<float, 16> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate);

        void regular_covariance_to_symmetric(unsigned int start_index=0, unsigned int end_index=16);
        void rotation_from_axis_angle(Matrix3f &r, const array<float, 3> &a) const;
        void quaternion_from_axis_angle(Quaternionf &q, const array<float, 3> &a) const;
        float kahan_summation(float sum_previous, float input, float &accumulator);
    };
}

#endif //NAVIGATION_ESKF_ESKF_H
