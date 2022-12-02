#ifndef NAVIGATION_GENERATOR_UTILS_H
#define NAVIGATION_GENERATOR_UTILS_H

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <deque>
#include <array>

namespace generator {
    using namespace std;
    using namespace Eigen;
    
    Matrix3d askew(const Vector3d &vec);
    Matrix3d euler2rot(const Vector3d &euler);
    Vector3d rot2euler(const Matrix3d &rot);
    Quaterniond euler2quat(const Vector3d &euler);
    Vector3d quat2euler(const Quaterniond &quat);
    Matrix3d quat2rot(const Quaterniond &quat);
    Quaterniond rot2quat(const Matrix3d &rot);
    Matrix3d rv2rot(const Vector3d &rv);
    Quaterniond rv2quat(const Vector3d &rv);
    Vector3d quat2rv(Quaterniond quat);
}


#endif // NAVIGATION_GENERATOR_UTILS_H