//
// Created by Cain on 2022/11/19.
//

#ifndef NAVIGATION_GLOBAL_ESKF_GESKF_H
#define NAVIGATION_GLOBAL_ESKF_GESKF_H

#include <Eigen/Dense>
#include <cmath>
#include <array>
#include "eskf.h"

namespace eskf {
    using namespace std;
    using namespace Eigen;

    class GESKF : public eskf::ESKF {
    public:
        explicit GESKF(float dt, const float g=9.8f) : eskf::ESKF(dt, g) {};

        // Priori
        void predict_covariance(const Vector3f &w, const Vector3f &a) override;

        // Posteriori
        unsigned char fuse_position(const Vector3f &pos, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) override;
        unsigned char fuse_velocity(const Vector3f &vel, const Vector3f &w, const Vector3f &a, const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) override;
        unsigned char fuse_magnet(const Vector3f &mag, const Vector3f &w, const Vector3f &a, const Vector3f &noise_std, const Vector3f &gate) override;
        unsigned char fuse_declination(const Vector2f &dec, const Vector3f &w, const Vector3f &a, const Vector2f &noise_std, const Vector2f &gate) override;
        void correct_state() override;
        void correct_covariance() override;
    };
}

#endif //NAVIGATION_GLOBAL_ESKF_GESKF_H
