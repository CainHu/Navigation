//#pragma GCC optimize(2)

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <ctime>
#include <Eigen/Dense>
#include <algorithm>
#include <random>
#include "eskf.h"
#include "param.h"

using namespace std;
using namespace Eigen;

int main() {
    // eskf::ESKF eskf_test(0.004f);
    // eskf_test.set_gyroscope_standard_deviation(eskf::noise_std_gyro);
    // eskf_test.set_accelerometer_standard_deviation(eskf::noise_std_acc);
    // eskf_test.set_drift_gyroscope_standard_deviation(eskf::noise_std_drift_gyro);
    // eskf_test.set_drift_saccelerometer_tandard_deviation(eskf::noise_std_drift_acc);
    // eskf_test.set_gravity_standard_deviation(eskf::noise_std_grav);
    // eskf_test.set_processing_standard_deviation(eskf::noise_std_proc);
    // eskf_test.initialize();

    // array<float, 3> w = {1.f, 2.f, 3.f};
    // array<float, 3> a = {-3.f, -2.f, -9.8f - 1.f};

    
    // for (int i = 0; i < 50; ++i) {
    //     eskf_test.predict_state(w, a);
    //     eskf_test.predict_covariance(w, a);
    // }
    // const array<array<float, 16>, 16> cov = eskf_test.get_covariance_matrix();
    // for (int i = 0; i < 16; ++i) {
    //     for (int j = 0; j < 16; ++j) {
    //         cout << cov[i][j] * 1000000.f << ", ";
    //     }
    //     cout << endl;
    // }









    static const unsigned int num_step = 100000;
    static const float g = 9.80665f;

    static const float dt = 0.0001f;
    float ts = 0.001f;

    Vector3f G(0.f, 0.f, g);

    Vector3f p(0.f, 0.f, 0.f), v(0.f, 0.f, 0.f), bg(0.f, 0.f, 0.f), ba(0.f, 0.f, 0.f);
    Quaternionf q(1.f, 0.f, 0.f, 0.f);
    Vector3f w(0.f, 0.f, 0.f), a(0.f, 0.f, -g);

    static array<Vector3f, num_step> p_true, v_true, bg_true, ba_true;
    static array<Quaternionf, num_step> q_true;
    static array<Vector3f, num_step> w_true, a_true;

    static array<Vector3f, num_step> pl_meas, vl_meas, pr_meas, vr_meas;
    static array<Vector3f, num_step> w_meas, a_meas;

    static array<Vector3f, num_step> p_hat, v_hat, bg_hat, ba_hat;
    static array<Quaternionf, num_step> q_hat;
    static array<float, num_step> g_hat;

    default_random_engine random_engine;
    normal_distribution<float> dist(0.f, 1.f);

    Vector3f dl(eskf::d_gps_left[0], eskf::d_gps_left[1], eskf::d_gps_left[2]);
    Vector3f dr(eskf::d_gps_right[0], eskf::d_gps_right[1], eskf::d_gps_right[2]);

     // Simulation
    unsigned int j = 0;
    w += 3.f * Vector3f(dist(random_engine), dist(random_engine), dist(random_engine)) * ts;
    a += 10.f * Vector3f(dist(random_engine), dist(random_engine), dist(random_engine)) * ts;
    for (unsigned int i = 0; i < num_step * (unsigned int)(ts / dt); ++i) {
        Vector3f v_last = v;
        Matrix3f R = q.toRotationMatrix();
        v += (R * a + G) * dt;
        p += 0.5f * (v + v_last) * dt;

        Vector3f d_ang = w * dt;
        AngleAxisf angle_axis(d_ang.norm(), d_ang.normalized());
        Quaternionf dq(angle_axis);
        q = q * dq;

        Vector3f nbg(eskf::noise_std_drift_gyro[0] * dist(random_engine),
                    eskf::noise_std_drift_gyro[1] * dist(random_engine),
                    eskf::noise_std_drift_gyro[2] * dist(random_engine));

        Vector3f nba(eskf::noise_std_drift_acc[0] * dist(random_engine),
                    eskf::noise_std_drift_acc[1] * dist(random_engine),
                    eskf::noise_std_drift_acc[2] * dist(random_engine));

        bg += nbg * dt;
        ba += nba * dt;

        if ((i + 1) % (unsigned int)(ts / dt) == 0) {
            // Real State
            p_true[j] = p;
            v_true[j] = v;
            bg_true[j] = bg;
            ba_true[j] = ba;

            q_true[j] = q;

            w_true[j] = w;
            a_true[j] = a;

            // Measurement
            Vector3f npl(dist(random_engine) * eskf::noise_std_rtk_pos[0],
                         dist(random_engine) * eskf::noise_std_rtk_pos[1],
                         dist(random_engine) * eskf::noise_std_rtk_pos[2]);
            Vector3f nvl(dist(random_engine) * eskf::noise_std_gps_vel[0],
                         dist(random_engine) * eskf::noise_std_gps_vel[1],
                         dist(random_engine) * eskf::noise_std_gps_vel[2]);
            Vector3f npr(dist(random_engine) * eskf::noise_std_rtk_pos[0],
                         dist(random_engine) * eskf::noise_std_rtk_pos[1],
                         dist(random_engine) * eskf::noise_std_rtk_pos[2]);
            Vector3f nvr(dist(random_engine) * eskf::noise_std_gps_vel[0],
                         dist(random_engine) * eskf::noise_std_gps_vel[1],
                         dist(random_engine) * eskf::noise_std_gps_vel[2]);
            Matrix3f R = q_true[j].toRotationMatrix();            
            pl_meas[j] = p_true[j] + R * dl + npl;
            vl_meas[j] = v_true[j] - R * (dl.cross(w_true[j])) + nvl;
            pr_meas[j] = p_true[j] + R * dr + npr;
            vr_meas[j] = v_true[j] - R * (dr.cross(w_true[j])) + nvr;

            Vector3f ng(dist(random_engine) * eskf::noise_std_gyro[0],
                        dist(random_engine) * eskf::noise_std_gyro[1],
                        dist(random_engine) * eskf::noise_std_gyro[2]);
            Vector3f na(dist(random_engine) * eskf::noise_std_acc[0],
                        dist(random_engine) * eskf::noise_std_acc[1],
                        dist(random_engine) * eskf::noise_std_acc[2]);
            w_meas[j] = w_true[j] + bg + ng;
            a_meas[j] = a_true[j] + ba + na;
        
            // Simulation
            w = 1.f * Vector3f(dist(random_engine), dist(random_engine), dist(random_engine));
            a = 3.f * Vector3f(dist(random_engine), dist(random_engine), dist(random_engine)) - 3.f * v_true[j];
            ++j;
        }
    }

    eskf::ESKF eskf_rtk(ts);
    eskf_rtk.set_gyroscope_standard_deviation(eskf::noise_std_gyro);
    eskf_rtk.set_accelerometer_standard_deviation(eskf::noise_std_acc);
    eskf_rtk.set_drift_gyroscope_standard_deviation(eskf::noise_std_drift_gyro);
    eskf_rtk.set_drift_saccelerometer_tandard_deviation(eskf::noise_std_drift_acc);
    eskf_rtk.set_gravity_standard_deviation(eskf::noise_std_grav);
    eskf_rtk.set_processing_standard_deviation(eskf::noise_std_proc);
    eskf_rtk.initialize();

    for (unsigned i = 0; i < num_step; ++i) {
        array<float, 3> eskf_w, eskf_a;
        eskf_w[0] = w_meas[i](0);
        eskf_w[1] = w_meas[i](1);
        eskf_w[2] = w_meas[i](2);
        eskf_a[0] = a_meas[i](0);
        eskf_a[1] = a_meas[i](1);
        eskf_a[2] = a_meas[i](2);

        array<float, 3> eskf_pl, eskf_vl, eskf_pr, eskf_vr;
        eskf_pl[0] = pl_meas[i](0);
        eskf_pl[1] = pl_meas[i](1);
        eskf_pl[2] = pl_meas[i](2);
        eskf_vl[0] = vl_meas[i](0);
        eskf_vl[1] = vl_meas[i](1);
        eskf_vl[2] = vl_meas[i](2);
        eskf_pr[0] = pr_meas[i](0);
        eskf_pr[1] = pr_meas[i](1);
        eskf_pr[2] = pr_meas[i](2);
        eskf_vr[0] = vr_meas[i](0);
        eskf_vr[1] = vr_meas[i](1);
        eskf_vr[2] = vr_meas[i](2);

        eskf_rtk.predict_state(eskf_w, eskf_a);
        eskf_rtk.predict_covariance(eskf_w, eskf_a);
        unsigned char info;
        info = eskf_rtk.fuse_position(eskf_pl, eskf_w, eskf_a, eskf::d_gps_left, eskf::noise_std_rtk_pos, eskf::gate_rtk_pos);
        if (info != 0) {
            cout << "pl: " << int(info) << endl;
        }
        info = eskf_rtk.fuse_velocity(eskf_vl, eskf_w, eskf_a, eskf::d_gps_left, eskf::noise_std_gps_vel, eskf::gate_gps_vel);
        if (info != 0) {
            cout << "vl: " << int(info) << endl;
        }
        info = eskf_rtk.fuse_position(eskf_pr, eskf_w, eskf_a, eskf::d_gps_right, eskf::noise_std_rtk_pos, eskf::gate_rtk_pos);
        if (info != 0) {
            cout << "pr: " << int(info) << endl;
        }
        info = eskf_rtk.fuse_velocity(eskf_vr, eskf_w, eskf_a, eskf::d_gps_right, eskf::noise_std_gps_vel, eskf::gate_gps_vel);
        if (info != 0) {
            cout << "vr: " << int(info) << endl;
        }
        eskf_rtk.correct_state();

        auto state = eskf_rtk.get_state();
        p_hat[i](0) = state[0];
        p_hat[i](1) = state[1];
        p_hat[i](2) = state[2];
        v_hat[i](0) = state[3];
        v_hat[i](1) = state[4];
        v_hat[i](2) = state[5];
        bg_hat[i](0) = state[6];
        bg_hat[i](1) = state[7];
        bg_hat[i](2) = state[8];
        ba_hat[i](0) = state[9];
        ba_hat[i](1) = state[10];
        ba_hat[i](2) = state[11];
        g_hat[i] = state[12];

        auto rot = eskf_rtk.get_rotation_matrix();
        Matrix3f r;
        r << rot[0][0], rot[0][1], rot[0][2],
             rot[1][0], rot[1][1], rot[1][2],
             rot[2][0], rot[2][1], rot[2][2];
        Quaternionf q(r);
        q_hat[i] = q;     

        cout << g_hat[i] << endl;

        // cout << v_true[i].transpose() << ", " << v_hat[i].transpose() << endl;
    }

    return 0;
}
