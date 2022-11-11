//#pragma GCC optimize(2)

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <ctime>
#include <Eigen/Dense>
#include "eskf.h"

using namespace std;

int main() {
    eskf::ESKF eskf(0.001f, 0.f, -0.1f, 0.f, 0.f, 0.1f, 0.f);

    array<float, 3> w = {1.f, 2.f, 3.f}, a = {-10.f, -7.f, -4.f};

    eskf.predict_covariance(w, a);

    const array<array<float, 16>, 16> & cov = eskf.get_covariance_matrix();
    const array<float, 13> & state = eskf.get_state();
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            cout << cov[i][j] << ", ";
        }
        cout << endl;
    }

//    time_t t1, t2;
//    t1 = clock();
//    for (int i = 0; i < 10000; ++i) {
//        eskf.predict_covariance(w, a);
//    }
//    t2 = clock();
//    cout << t2 - t1 << endl;

    return 0;
}
