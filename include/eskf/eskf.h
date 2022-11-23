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

    union filter_control_status {
        struct {
            bool acc_x_bias : 1;
            bool acc_y_bias : 1;
            bool acc_z_bias : 1;
            bool grav : 1;
            bool mag : 1;
            bool mag_bias : 1;
            bool wind : 1;
        } flags;
        unsigned char status;
    };

    union innovation_fault_status {
        struct {
            bool reject_pos : 1;
            bool reject_vel : 1;
            bool reject_mag : 1;
            bool reject_wind : 1;
        } flags;
        unsigned char status;
    };

    class ESKF {
    public:
        constexpr static unsigned char dim {24};

        ESKF(float dt, const float g=9.8f) 
             : _g_init(g), _g(g), _dt(dt), _dt2(dt * dt) {
            _rot.setIdentity();
            reset_priori_covariance_matrix();
            reset_state();
            reset_error_state();
            reset_accmulator_cov();
            for (unsigned int i = 0; i < dim; ++i) {
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
        virtual void correct_covariance() = 0;

        // Resetters
        void reset_state() {
            _p.setZero();
            _v.setZero();
            _bg.setZero();
            _ba.setZero();
            _g = _g_init;
            _m.setZero();
            _bm.setZero();
            _w.setZero();
            _rot.setIdentity();
            _q.setIdentity();
        }

        void reset_error_state() {
            for (float &es : _error_state) { 
                es = 0.f; 
            } 
        }

        void reset_priori_covariance_matrix() {
            _q_cov.setZero();
        }

        void reset_covariance_matrix(const unsigned char start_index, const unsigned char end_index, const Matrix<float, dim, 1> &diag_cov) {
            for (unsigned char i = start_index; i < end_index; ++i) {
                // Diaginal
                _cov[i][i] = diag_cov[i] * _dt2;

                // Upper triangular
                for (unsigned char j = start_index; j < i; ++j) {
                    _cov[j][i] = 0.f;
                }

                // Columns
                for (unsigned char j = 0; j < start_index; ++j) {
                    _cov[j][i] = 0.f;
                }

                // Rows
                for (unsigned char j = end_index; j < ESKF::dim; ++j) {
                    _cov[i][j] = 0.f;
                }
            }
        }

        void reset_covariance_matrix(const unsigned char start_index, const unsigned char end_index) {
            reset_covariance_matrix(start_index, end_index, _q_cov);
        }

        template <unsigned char N>
        void reset_covariance_matrix(const unsigned char start_index, const Matrix<float, N, 1> &diag_cov) {
            unsigned char end_index = start_index + N;
            for (unsigned char i = start_index; i < end_index; ++i) {
                // Diaginal
                _cov[i][i] = diag_cov[i - start_index] * _dt2;

                // Upper triangular
                for (unsigned char j = start_index; j < i; ++j) {
                    _cov[j][i] = 0.f;
                }

                // Columns
                for (unsigned char j = 0; j < start_index; ++j) {
                    _cov[j][i] = 0.f;
                }

                // Rows
                for (unsigned char j = end_index; j < ESKF::dim; ++j) {
                    _cov[i][j] = 0.f;
                }
            }
        }

        template <unsigned char N>
        void reset_covariance_matrix(const unsigned char start_index, const float diag_cov) {
            unsigned char end_index = start_index + N;
            for (unsigned char i = start_index; i < end_index; ++i) {
                // Diaginal
                _cov[i][i] = diag_cov * _dt2;

                // Upper triangular
                for (unsigned char j = start_index; j < i; ++j) {
                    _cov[j][i] = 0.f;
                }

                // Columns
                for (unsigned char j = 0; j < start_index; ++j) {
                    _cov[j][i] = 0.f;
                }

                // Rows
                for (unsigned char j = end_index; j < ESKF::dim; ++j) {
                    _cov[i][j] = 0.f;
                }
            }
        }

        template <unsigned char N>
        void reset_accmulator_cov(const unsigned char start_index) {
            unsigned char end_index = start_index + N;
            for (unsigned char i = start_index; i < end_index; ++i) {
                _accumulator_cov[i] = 0.f;
            }
        }

        void reset_accmulator_cov() {
            for (float &a : _accumulator_cov) {
                a = 0.f;
            }
        }

        // Switches
        void enable_estimation_acc_x_bias() { _control_status.flags.acc_x_bias = true; };
        void enable_estimation_acc_y_bias() { _control_status.flags.acc_y_bias = true; };
        void enable_estimation_acc_z_bias() { _control_status.flags.acc_z_bias = true; };
        void enable_estimation_acc_bias () { 
            enable_estimation_acc_x_bias(); 
            enable_estimation_acc_y_bias();
            enable_estimation_acc_z_bias(); 
        }
        void enable_estimation_gravity() { _control_status.flags.grav = true; };
        void enable_estimation_magnet() { _control_status.flags.mag = true; };
        void enable_estimation_magnet_bias() { 
            _control_status.flags.mag_bias = true; 
            enable_estimation_magnet();
        };
        void enable_estimation_wind() { _control_status.flags.wind = true; } ;

        void disable_estimation_acc_x_bias() { 
            _control_status.flags.acc_x_bias = false; 
            reset_covariance_matrix<1>(12, 0.f);
            reset_accmulator_cov<1>(12);
            regular_covariance_to_symmetric<1>(12);
        };
        void disable_estimation_acc_y_bias() { 
            _control_status.flags.acc_y_bias = false; 
            reset_covariance_matrix<1>(13, 0.f);
            reset_accmulator_cov<1>(13);
            regular_covariance_to_symmetric<1>(13);
        };
        void disable_estimation_acc_z_bias() { 
            _control_status.flags.acc_z_bias = false; 
            reset_covariance_matrix<1>(14, 0.f);
            reset_accmulator_cov<1>(14);
            regular_covariance_to_symmetric<1>(14);
        };
        void disable_estimation_acc_bias() {
            disable_estimation_acc_x_bias();
            disable_estimation_acc_y_bias();
            disable_estimation_acc_z_bias();
        }
        void disable_estimation_gravity() { 
            _control_status.flags.grav = false; 
            reset_covariance_matrix<1>(15, 0.f);
            reset_accmulator_cov<1>(15);
            regular_covariance_to_symmetric<1>(15);
        };
        void disable_estimation_magnet_bias() { 
            _control_status.flags.mag_bias = false; 
            reset_covariance_matrix<3>(19, 0.f);
            reset_accmulator_cov<3>(19);
            regular_covariance_to_symmetric<3>(19);
        };
        void disable_estimation_magnet() { 
            _control_status.flags.mag = false; 
            reset_covariance_matrix<3>(16, 0.f);
            reset_accmulator_cov<3>(16);
            regular_covariance_to_symmetric<3>(16);
            disable_estimation_magnet_bias();
        };
        void disable_estimation_wind() { 
            _control_status.flags.wind = false; 
            reset_covariance_matrix<2>(22, 0.f);
            reset_accmulator_cov<2>(22);
            regular_covariance_to_symmetric<2>(22);
        };


        // Getters
        const Vector3f &get_position() const { return _p; };
        const Vector3f &get_velocity() const { return _v; };
        const Vector3f &get_drift_gyro() const { return _bg; };
        const Vector3f &get_drift_acc() const { return _ba; };
        const Vector3f &get_magnet() const { return _m; };
        const Vector3f &get_drift_magnet() const { return _bm; };
        const Vector2f &get_wind() const { return _w; };
        const Matrix3f &get_rotation_matrix() const { return _rot; };
        const Quaternionf &get_quaternion() const { return _q; };
        const float get_gravity() const { return _g; };
        const float get_dt() const { return _dt; };
        const array<float, dim> &get_error_state() const { return _error_state; };
        const array<array<float, dim>, dim> &get_covariance_matrix() const { return _cov; };


        // Setters
        void set_dt(const float dt) { _dt = dt; _dt2 = dt * dt; };
        void set_gravity(const float g) { _g_init = g; _g = g; };
        void set_position(const Vector3f &p) { _p = p; };
        void set_velocity(const Vector3f &v) { _v = v; };
        void set_drift_gyro(const Vector3f &bg) { _bg = bg; };
        void set_drift_acc(const Vector3f &ba) { _ba = ba; };
        void set_magnet(const Vector3f &m) { _m = m; };
        void set_drift_magnet(const Vector3f &bm) { _bm = bm; };
        void set_wind(const Vector2f &w) { _w = w; };
        void set_attitude(const Matrix3f &r) { _rot = r; _q = _rot; };
        void set_attitude(const Quaternionf &q) { _q = q; _rot = _q; };
        void set_attitude(const AngleAxisf &a) { _q = a; _rot = _q;};    
        void set_accelerometer_standard_deviation(const Vector3f &std) { _q_cov[3] = std[0] * std[0]; _q_cov[4] = std[1] * std[1]; _q_cov[5] = std[2] * std[2]; };
        void set_gyroscope_standard_deviation(const Vector3f &std) { _q_cov[6] = std[0] * std[0]; _q_cov[7] = std[1] * std[1]; _q_cov[8] = std[2] * std[2]; };
        void set_drift_saccelerometer_tandard_deviation(const Vector3f &std) { _q_cov[12] = std[0] * std[0]; _q_cov[13] = std[1] * std[1]; _q_cov[14] = std[2] * std[2]; }
        void set_drift_gyroscope_standard_deviation(const Vector3f &std) { _q_cov[9] = std[0] * std[0]; _q_cov[10] = std[1] * std[1]; _q_cov[11] = std[2] * std[2]; }
        void set_gravity_standard_deviation(const float std) { _q_cov[15] = std * std; };
        void set_magnetometer_standard_deviation(const Vector3f &std) { _q_cov[16] = std[0] * std[0]; _q_cov[17] = std[1] * std[1]; _q_cov[18] = std[2] * std[2]; }
        void set_drift_magnetometer_standard_deviation(const Vector3f &std) { _q_cov[19] = std[0] * std[0]; _q_cov[20] = std[1] * std[1]; _q_cov[21] = std[2] * std[2]; }
        void set_wind_standard_deviation(const Vector2f &std) { _q_cov[22] = std[0] * std[0]; _q_cov[23] = std[1] * std[1]; }
        void set_processing_standard_deviation(const float std) { for (unsigned char i = 0; i < dim; ++i) { _q_cov[i] += std * std; } };
        void set_priori_covariance_matrix(Matrix<float, dim, 1> &q) { _q_cov = q; };

    protected:
        filter_control_status _control_status {0};
        innovation_fault_status _innovation_fault_status {0};

        float _g_init;                 // Initial gravity constant

        float _dt;                      // Sample time of IMU
        float _dt2;                     // Square of sample time

        Matrix<float, dim, 1> _q_cov {};         // Priori covariance matrix
  
        Vector3f _p, _v, _bg, _ba;              // Translational state: p, v, bg, ba
        float _g;                               // Gravity constant
        Vector3f _m, _bm;                       // Translational state: m, bm
        Vector2f _w;                            // Translational state: w
        Matrix3f _rot;                          // Rotation matrix from body frame to world frame
        Quaternionf _q;                         // Quaternion matrix from body frame to world frame

        array<float, dim> _error_state {};       // [δp, δv, δθ, δbg, δba, δg, δm, δmb, δw]
        array<array<float, dim>, dim> _cov {};    // Covariance matrix of error state

        array<float, dim> _accumulator_cov {};   // Accumulator of covariance matrix

        // Posteriori
        unsigned char conservative_posteriori_estimate(const array<float, dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate);

        // Utils
        void rotation_from_axis_angle(Matrix3f &r, const array<float, 3> &a) const;
        void quaternion_from_axis_angle(Quaternionf &q, const array<float, 3> &a) const;

        template <unsigned char N>
        void regular_covariance_to_symmetric(unsigned char start_index) {
            unsigned char end_index = start_index + N;
            for (unsigned char i = start_index; i < end_index; ++i) {
                for (unsigned char j = 0; j < i; ++j) {
                    _cov[i][j] = _cov[j][i];
                }
            }
            for (unsigned char i = end_index; i < ESKF::dim; ++i) {
                for (unsigned char j = start_index; j < end_index; ++j) {
                    _cov[i][j] = _cov[j][i];
                }
            }
        }

        float kahan_summation(float sum_previous, float input, float &accumulator) {
            /*
            accumulator中记录了sum_previous + y中, y舍弃掉的部分
            */
            const float y = input - accumulator;
            const float t = sum_previous + y;
            accumulator = (t - sum_previous) - y;
            return t;
        }

        template <unsigned char N>
        void add_processing_covariance(unsigned char start_index) {
            unsigned char end_index = start_index + N;
            for (unsigned char i = start_index; i < end_index; ++i) {
                _cov[i][i] = kahan_summation(_cov[i][i], _q_cov[i] * _dt2, _accumulator_cov[i]);
            }
        }
    };
}

#endif //NAVIGATION_ESKF_ESKF_H
