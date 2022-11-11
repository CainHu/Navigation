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
        ESKF(float dt, float g=9.8f) : _dt(dt), _dt2(dt * dt) {
            _rot[0][0] = _rot[1][1] = _rot[2][2] = 1.f;
            for (unsigned char i = 0; i < 16; ++i) {
                _cov[i][i] = 1.f;
            }
            _state[12] = g;
        }

        void predict_covariance(const array<float, 3> &w, const array<float, 3> &a);
        void predict_state(const array<float, 3> &w, const array<float, 3> &a);

        void rotation_from_axis_angle(array<array<float, 3>, 3> &r, const array<float, 3> &a);

        const array<float, 13> &get_state() { return _state; };
        const array<float, 16> &get_error_state() { return _error_state; };
        const array<array<float, 16>, 16> &get_covariance_matrix() { return _cov; };
        const array<array<float, 3>, 3> &get_rotation_matrix() { return _rot; };
        const float get_dt() { return _dt; }

        void set_dt(float dt) { _dt = dt; _dt2 = dt * dt; };
    private:
        float _dt;                      // Sample time of IMU
        float _dt2;                     // Square of sample time

        array<float, 13> _state {};           // [p, v, bg, ba, g]
        array<array<float, 3>, 3> _rot {};     // Rotation matrix from body frame to world frame

        array<float, 16> _error_state {};     // [δp, δv, δθ, δbg, δba, δg]
        array<array<float, 16>, 16> _cov {};     // Covariance matrix of error state
    };
}

#endif //NAVIGATION_ESKF_H
