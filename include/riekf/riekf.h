//
// Created by Cain on 2022/11/19.
//

#ifndef NAVIGATION_RIEKF_RIEKF_H
#define NAVIGATION_RIEKF_RIEKF_H

#include <Eigen/Dense>
#include <cmath>
#include <array>
#include "eskf.h"

namespace riekf {
    using namespace std;
    using namespace Eigen;

    class RIEKF : public eskf::ESKF {
    public:
        RIEKF(float dt, const float g=9.8f) : eskf::ESKF(dt, g) {};

        // Priori
        void predict_covariance(const Vector3f &w, const Vector3f &a);

        // Posteriori
        unsigned char fuse_position(const Vector3f &pos, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate);
        unsigned char fuse_velocity(const Vector3f &vel, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate);
        void correct_state();
        void correct_covariance();
    };
}

#endif //NAVIGATION_RIEKF_RIEKF_H
