//
// Created by Cain on 2022/11/19.
//

#ifndef NAVIGATION_GESKF_GESKF_H
#define NAVIGATION_GESKF_GESKF_H

#include <Eigen/Dense>
#include <cmath>
#include <array>
#include "eskf.h"

namespace geskf {
    using namespace std;
    using namespace Eigen;

    class GESKF : public eskf::ESKF {
    public:
        GESKF(float dt, const float g=9.8f) : eskf::ESKF(dt, g) {};

        // Priori
        void predict_covariance(const Vector3f &w, const Vector3f &a);

        // Posteriori
        unsigned char fuse_position(const Vector3f &pos, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate);
        unsigned char fuse_velocity(const Vector3f &vel, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate);
        void correct_state();
    };
}

#endif //NAVIGATION_GESKF_GESKF_H