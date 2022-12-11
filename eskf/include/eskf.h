//
// Created by Cain on 2022/11/11.
//

#ifndef NAVIGATION_ESKF_ESKF_H
#define NAVIGATION_ESKF_ESKF_H

#include <Eigen/Dense>
#include <cmath>
#include <array>
#include "utils.h"
#include "common.h"
//#include "eskf_runner.h"

namespace eskf {
    using namespace std;
    using namespace Eigen;

    template <unsigned char DELAYS>
    class ESKFRunner;

    class ESKF {
    public:
        constexpr static unsigned char DIM {24};

        ESKF(float dt, const float g=9.8f, const float h=1.f, const float dec_y=0.f, const float dec_z=0.f) 
             : _g_init(g), _g(g), _h_init(h), _h(h), _dec_init(dec_y, dec_z), _dec(dec_y, dec_z), _dt(dt), _dt2(dt * dt) {   
            _rot.setIdentity();
            reset_priori_covariance_matrix();
            reset_state();
            reset_error_state();
            reset_accmulator_cov();
            for (unsigned int i = 0; i < DIM; ++i) {
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
        virtual unsigned char fuse_magnet(const Vector3f &mag, const Vector3f &w, const Vector3f &a, const Vector3f &noise_std, const Vector3f &gate) = 0;
        virtual unsigned char fuse_declination(const Vector2f &dec, const Vector3f &w, const Vector3f &a, const Vector2f &noise_std, const Vector2f &gate) = 0;

        virtual void correct_state() = 0;
        virtual void correct_covariance() = 0;

        // Resetters
        void reset_state() {
            _p.setZero();
            _v.setZero();
            _bg.setZero();
            _ba.setZero();
            _g = _g_init;
            _h = _h_init;
            _dec = _dec_init;
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

        void reset_covariance_matrix(const unsigned char start_index, const unsigned char end_index, const Matrix<float, DIM, 1> &diag_cov) {
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
                for (unsigned char j = end_index; j < ESKF::DIM; ++j) {
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
                for (unsigned char j = end_index; j < ESKF::DIM; ++j) {
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
                for (unsigned char j = end_index; j < ESKF::DIM; ++j) {
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
        void enable_estimation_declination() { _control_status.flags.dec = true; };
        void enable_estimation_magnet_bias() { _control_status.flags.mag_bias = true; };
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
        void disable_estimation_magnet() { 
            _control_status.flags.mag = false; 
            reset_covariance_matrix<1>(16, 0.f);
            reset_accmulator_cov<1>(16);
            regular_covariance_to_symmetric<1>(16);
        };
        void disable_estimation_declination() { 
            _control_status.flags.mag = false; 
            reset_covariance_matrix<2>(17, 0.f);
            reset_accmulator_cov<2>(17);
            regular_covariance_to_symmetric<2>(17);
        };
        void disable_estimation_magnet_bias() { 
            _control_status.flags.mag_bias = false; 
            reset_covariance_matrix<3>(19, 0.f);
            reset_accmulator_cov<3>(19);
            regular_covariance_to_symmetric<3>(19);
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
        const float get_magnet() const { return _h; };
        const Vector2f &get_declination() const { return _dec; };
        const Vector3f &get_drift_magnet() const { return _bm; };
        const Vector2f &get_wind() const { return _w; };
        const Matrix3f &get_rotation_matrix() const { return _rot; };
        const Quaternionf &get_quaternion() const { return _q; };
        const float get_gravity() const { return _g; };
        const float get_dt() const { return _dt; };
        const array<float, DIM> &get_error_state() const { return _error_state; };
        const array<array<float, DIM>, DIM> &get_covariance_matrix() const { return _cov; };


        // Setters
        void set_dt(const float dt) { _dt = dt; _dt2 = dt * dt; };
        void set_gravity(const float g) { _g_init = g; _g = g; };
        void set_position(const Vector3f &p) { _p = p; };
        void set_velocity(const Vector3f &v) { _v = v; };
        void set_bias_gyro(const Vector3f &bg) { _bg = bg; };
        void set_bias_acc(const Vector3f &ba) { _ba = ba; };
        void set_mag_norm(const float h) { _h_init = h; _h = h; };
        void set_mag_dec(const Vector2f &dec) { _dec_init = dec; _dec = dec; };
        void set_bias_magnet(const Vector3f &bm) { _bm = bm; };
        void set_wind(const Vector2f &w) { _w = w; };
        void set_attitude(const Matrix3f &r) { _rot = r; _q = _rot; };
        void set_attitude(const Quaternionf &q) { _q = q; _rot = _q; };
        void set_attitude(const AngleAxisf &a) { _q = a; _rot = _q;};    

        void set_std_proc_pos(const Vector3f &std) { _q_cov[0] = std[0] * std[0]; _q_cov[1] = std[1] * std[1]; _q_cov[2] = std[2] * std[2]; };
        void set_std_proc_vel(const Vector3f &std) { _q_cov[3] = std[0] * std[0]; _q_cov[4] = std[1] * std[1]; _q_cov[5] = std[2] * std[2]; };
        void set_std_proc_ang(const Vector3f &std) { _q_cov[6] = std[0] * std[0]; _q_cov[7] = std[1] * std[1]; _q_cov[8] = std[2] * std[2]; };
        void set_std_proc_bias_gyro(const Vector3f &std) { _q_cov[9] = std[0] * std[0]; _q_cov[10] = std[1] * std[1]; _q_cov[11] = std[2] * std[2]; };
        void set_std_proc_bias_acc(const Vector3f &std) { _q_cov[12] = std[0] * std[0]; _q_cov[13] = std[1] * std[1]; _q_cov[14] = std[2] * std[2]; };
        void set_std_proc_grav(const float &std) { _q_cov[15] = std * std; };
        void set_std_proc_mag_norm(const float &std) { _q_cov[16] = std * std; };
        void set_std_proc_mag_dec(const Vector2f &std) { _q_cov[17] = std[0] * std[0]; _q_cov[18] = std[1] * std[1]; };
        void set_std_proc_bias_mag(const Vector3f &std) { _q_cov[19] = std[0] * std[0]; _q_cov[20] = std[1] * std[1]; _q_cov[21] = std[2] * std[2]; };
        void set_std_proc_wind(const Vector2f &std) { _q_cov[22] = std[0] * std[0]; _q_cov[23] = std[1] * std[1]; };

        void set_std_acc(const Vector3f &std) { _q_cov[3] = std[0] * std[0]; _q_cov[4] = std[1] * std[1]; _q_cov[5] = std[2] * std[2]; };
        void set_std_gyro(const Vector3f &std) { _q_cov[6] = std[0] * std[0]; _q_cov[7] = std[1] * std[1]; _q_cov[8] = std[2] * std[2]; };

        void set_priori_covariance_matrix(Matrix<float, DIM, 1> &q) { _q_cov = q; };

    protected:
        filter_control_status _control_status {0};
        innovation_fault_status _innovation_fault_status {0};

        // State _state;
        // ErrorState _error_state;

        float _g_init;                  // Initial gravity constant
        float _h_init;                  // Initial magnetic field intensity
        Vector2f _dec_init;              // Initial Declination

        float _dt;                      // Sample time of IMU
        float _dt2;                     // Square of sample time

        Matrix<float, DIM, 1> _q_cov {};        // Priori covariance matrix
  
        Vector3f _p, _v, _bg, _ba;              // Translational state: p, v, bg, ba
        float _g;                               // Gravity constant
        float _h;                               // Magnetic field intensity
        Vector2f _dec;                          // declination              
        Vector3f _bm;                           // Translational state: bh
        Vector2f _w;                            // Translational state: w
        Matrix3f _rot;                          // Rotation matrix from body frame to world frame
        Quaternionf _q;                         // Quaternion matrix from body frame to world frame

        array<float, DIM> _error_state {};       // [δp, δv, δθ, δbg, δba, δg, δm, δmb, δw]
        array<array<float, DIM>, DIM> _cov {};    // Covariance matrix of error state

        array<float, DIM> _accumulator_cov {};   // Accumulator of covariance matrix

        bool _imu_updated {false};

        // Posteriori
        unsigned char conservative_posteriori_estimate(const array<float, DIM> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate);

        template <unsigned char N>
        void regular_covariance_to_symmetric(unsigned char start_index) {
            unsigned char end_index = start_index + N;
            for (unsigned char i = start_index; i < end_index; ++i) {
                for (unsigned char j = 0; j < i; ++j) {
                    _cov[i][j] = _cov[j][i];
                }
            }
            for (unsigned char i = end_index; i < ESKF::DIM; ++i) {
                for (unsigned char j = start_index; j < end_index; ++j) {
                    _cov[i][j] = _cov[j][i];
                }
            }
        }

        template <unsigned char N>
        void copy_covariance_cols_to_rows(unsigned char start_index) {
            unsigned char end_index = start_index + N;
            for (unsigned char i = start_index; i < end_index; ++i) {
                for (unsigned char j = 0; j < i; ++j) {
                    _cov[i][j] = _cov[j][i];
                }
            }
        }

        template <unsigned char N>
        void copy_covariance_rows_to_cols(unsigned char start_index) {
            unsigned char end_index = start_index + N;
            for (unsigned char i = start_index; i < end_index; ++i) {
                for (unsigned char j = i + 1; j < DIM; ++j) {
                    _cov[j][i] = _cov[i][j];
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

    private:
        template<unsigned char DELAYS>
        friend class ESKFRunner;
    };
}

#endif //NAVIGATION_ESKF_ESKF_H