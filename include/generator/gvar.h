//
// Created by Cain on 2022/12/1.
//

#ifndef NAVIGATION_GENERATOR_GVAR_H
#define NAVIGATION_GENERATOR_GVAR_H

#include <cmath>
#include <Eigen/Dense>

extern const double GM;
extern const double Re;
extern const double wie;

extern const double ff;
extern const double ee;
extern const double e2;
extern const double Rp;
extern const double ge;
extern const double gp;
extern const double g0;
extern const double ug;
extern const double arcdeg;
extern const double arcmin;
extern const double arcsec;
extern const double hur;
extern const double dph;
extern const double dpsh;
extern const double ugpsHz;

Eigen::Vector3d geo2earth(const Eigen::Vector3d &pos);
Eigen::Vector3d earth2nav(const Eigen::Vector3d &earth, const Eigen::Vector3d &geo);
Eigen::Vector3d diff_geo(const Eigen::Vector3d &pos2, const Eigen::Vector3d &pos1);

#endif // NAVIGATION_GENERATOR_GVAR_H