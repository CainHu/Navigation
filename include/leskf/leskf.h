//
// Created by Cain on 2022/11/19.
//

/*
 Local Error State Kalman Filter
*/

#ifndef NAVIGATION_LESKF_LESKF_H
#define NAVIGATION_LESKF_LESKF_H

#include <Eigen/Dense>
#include <cmath>
#include <array>
#include "eskf.h"

namespace leskf {
    using namespace std;
    using namespace Eigen;

    class LESKF : public eskf::ESKF {
    public:
        LESKF(float dt, const float g=9.8f) : eskf::ESKF(dt, g) {};

        // Priori
        void predict_covariance(const Vector3f &w, const Vector3f &a);

        // Posteriori
        unsigned char fuse_position(const Vector3f &pos, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate);
        unsigned char fuse_velocity(const Vector3f &vel, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate);
        unsigned char fuse_magnet(const Vector3f &vel, const Vector3f &w, const Vector3f &a, const Vector3f &noise_std, const Vector3f &gate);
        unsigned char fuse_declination(const Vector2f &dec, const Vector3f &w, const Vector3f &a, const Vector2f &noise_std, const Vector2f &gate);
        void correct_state();
        void correct_covariance();
    };
}

#endif //NAVIGATION_LESKF_LESKF_H
