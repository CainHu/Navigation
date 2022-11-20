//
// Created by Cain on 2022/11/18.
//

#ifndef NAVIGATION_LIEKF_LIEKF_H
#define NAVIGATION_LIEKF_LIEKF_H

#include <Eigen/Dense>
#include <cmath>
#include <array>
#include "eskf.h"

namespace liekf {
    using namespace std;
    using namespace Eigen;

    class LIEKF : public eskf::ESKF {
    public:
        LIEKF(float dt, const float g=9.8f) : eskf::ESKF(dt, g) {};

        // Priori
        void predict_covariance(const Vector3f &w, const Vector3f &a);

        // Posteriori
        unsigned char fuse_position(const Vector3f &pos, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate);
        unsigned char fuse_velocity(const Vector3f &vel, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate);
        void correct_state();
    };
}

#endif //NAVIGATION_LIEKF_LIEKF_H
