//
// Created by Cain on 2022/12/1.
//

#ifndef NAVIGATION_GENERATOR_GENERATOR_H
#define NAVIGATION_GENERATOR_GENERATOR_H

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <deque>
#include <array>
#include <algorithm>
#include <random>

namespace generator {
    using namespace std;
    using namespace Eigen;

    class Fir {
    public:
        Fir(const double elem=0.) : buff(DIM, elem) { };

        double operator() (const double elem) {
            buff.push_back(elem);

            double res = 0.;
            for (unsigned int i = 0; i < DIM; ++i) {
                res += PARAMS[i] * buff[i];
            }

            return res;
        }

    private:
        constexpr static unsigned int DIM {21};

        constexpr static array<double, DIM> PARAMS {
            7.252481e-03f, 9.322776e-03f, 1.530768e-02f, 2.464950e-02f, 3.645113e-02f,  4.956447e-02f, 6.270454e-02f, 
            7.457793e-02f, 8.401243e-02f, 9.007478e-02f, 9.216459e-02f, 9.007478e-02f, 8.401243e-02f, 7.457793e-02f, 
            6.270454e-02f, 4.956447e-02f, 3.645113e-02f, 2.464950e-02f, 1.530768e-02f, 9.322776e-03f, 7.252481e-03f
        };

        deque<double> buff;
    };

    typedef struct Earth {
        double sl, cl, tl, sl2;
        double RMh, RNh, clRNh;
        Vector3d wnie, vn, wnen, wnin, wnien;
        Vector3d gn, gcc;
    } Earth;

    class Generator {
    public:
        Generator() : _dist(0., 1.) {  };
        void traj_gen(const Vector3d &euler0, const Vector3d &vn0, const Vector3d &pos0, const vector<Matrix<double, 5, 1>> &wat, const double ts);
        void ev2imu(const vector<Vector3d> &euler, const vector<Vector3d> &vn, const vector<Vector3d> &pos, const double ts);
        void ev2imu(const double ts);
        

    private:
        default_random_engine _random_engine;
        normal_distribution<double> _dist;

        vector<Vector3d> _euler, _vn, _pos;
        vector<Vector3d> _wm, _vm;

        Quaterniond quat_add_phi(const Quaterniond &quat, const Vector3d &phi);
        Quaterniond quat_del_phi(const Quaterniond &quat, const Vector3d &phi);
        Vector3d qq2phi(const Quaterniond &qpb, const Quaterniond &qnb);

        void earth(const Vector3d &pos, const Vector3d &vn, Earth &eth);

        template<unsigned int N>
        void cnscl(const array<Vector3d, N> &wm, const array<Vector3d, N> &vm, Vector3d &phim, Vector3d &dvbm);
        void cnscl(const vector<Vector3d> &wm, const vector<Vector3d> &vm, Vector3d &phim, Vector3d &dvbm);
        void imu_add_err(const vector<Vector3d> &wm, const vector<Vector3d> &vm, const Vector3d &eb, const Vector3d &web, const Vector3d &db, const Vector3d &wdb, const double ts,
                         vector<Vector3d> &wm_err, vector<Vector3d> &vm_err);

    };
}

#endif // NAVIGATION_GENERATOR_GENERATOR_H