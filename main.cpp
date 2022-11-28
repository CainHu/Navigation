//#pragma GCC optimize(2)

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <ctime>
#include <Eigen/Dense>
#include <algorithm>
#include <random>
#include "leskf.h"
#include "geskf.h"
#include "liekf.h"
#include "riekf.h"
#include "param.h"

using namespace std;
using namespace Eigen;

int main() {
    // eskf::ESKF eskf_test(0.001f);
    // eskf_test.reset_priori_covariance_matrix();
    // eskf_test.set_gyroscope_standard_deviation(eskf::noise_std_gyro);
    // eskf_test.set_accelerometer_standard_deviation(eskf::noise_std_acc);
    // eskf_test.set_drift_gyroscope_standard_deviation(eskf::noise_std_drift_gyro);
    // eskf_test.set_drift_saccelerometer_tandard_deviation(eskf::noise_std_drift_acc);
    // eskf_test.set_gravity_standard_deviation(eskf::noise_std_grav);
    // eskf_test.set_processing_standard_deviation(eskf::noise_std_proc);
    // eskf_test.initialize();

    // Vector3f w = {1.f, 2.f, 3.f};
    // Vector3f a = {-3.f, -2.f, -9.8f - 1.f};

    
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
    Vector3f M(0.9f, sqrtf(1.f - 0.9f * 0.9f), 0.f);
    float declination = atan2f(M[1], M[0]);

    Vector3f p(0.f, 0.f, 0.f), v(0.f, 0.f, 0.f), bg(0.f, 0.f, 0.f), ba(0.f, 0.f, 0.f), bm(0.f, 0.f, 0.f);
    Quaternionf q(1.f, 0.f, 0.f, 0.f);
    Vector3f w(0.f, 0.f, 0.f), a(0.f, 0.f, -g);

    static array<Vector3f, num_step> p_true, v_true, bg_true, ba_true, mb_true, bm_true;
    static array<Quaternionf, num_step> q_true;
    static array<Vector3f, num_step> w_true, a_true;

    static array<Vector3f, num_step> pl_meas, vl_meas, pr_meas, vr_meas, mb_meas;
    static array<Vector3f, num_step> w_meas, a_meas;

    static array<Vector3f, num_step> p_hat, v_hat, bg_hat, ba_hat, mb_hat, bm_hat;
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
        q.normalize();

        Vector3f nbg(eskf::noise_std_drift_gyro[0] * dist(random_engine),
                     eskf::noise_std_drift_gyro[1] * dist(random_engine),
                     eskf::noise_std_drift_gyro[2] * dist(random_engine));

        Vector3f nba(eskf::noise_std_drift_acc[0] * dist(random_engine),
                     eskf::noise_std_drift_acc[1] * dist(random_engine),
                     eskf::noise_std_drift_acc[2] * dist(random_engine));

        Vector3f nbm(eskf::noise_std_drift_mag[0] * dist(random_engine),
                     eskf::noise_std_drift_mag[1] * dist(random_engine),
                     eskf::noise_std_drift_mag[2] * dist(random_engine));

        bg += nbg * dt;
        ba += nba * dt;
        bm += nbm * dt;

        if ((i + 1) % (unsigned int)(ts / dt) == 0) {
            // Real State
            p_true[j] = p;
            v_true[j] = v;
            bg_true[j] = bg;
            ba_true[j] = ba;
            mb_true[j] = q.toRotationMatrix().transpose() * M;
            bm_true[j] = bm;

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
            Vector3f nmb(dist(random_engine) * eskf::noise_std_mag[0],
                         dist(random_engine) * eskf::noise_std_mag[1],
                         dist(random_engine) * eskf::noise_std_mag[2]);
            Matrix3f R = q_true[j].toRotationMatrix();            
            pl_meas[j] = p_true[j] + R * dl + npl;
            vl_meas[j] = v_true[j] - R * (dl.cross(w_true[j])) + nvl;
            pr_meas[j] = p_true[j] + R * dr + npr;
            vr_meas[j] = v_true[j] - R * (dr.cross(w_true[j])) + nvr;
            mb_meas[j] = mb_true[j] + nmb;

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

    // leskf::LESKF eskf_rtk(ts);
    geskf::GESKF eskf_rtk(ts);
    // liekf::LIEKF eskf_rtk(ts);
    // riekf::RIEKF eskf_rtk(ts);
    // eskf_rtk.set_magnet(M);
    eskf_rtk.set_gyroscope_standard_deviation(eskf::noise_std_gyro);
    eskf_rtk.set_accelerometer_standard_deviation(eskf::noise_std_acc);
    eskf_rtk.set_drift_gyroscope_standard_deviation(eskf::noise_std_drift_gyro);
    eskf_rtk.set_drift_saccelerometer_tandard_deviation(eskf::noise_std_drift_acc);
    eskf_rtk.set_gravity_standard_deviation(eskf::noise_std_proc_grav);
    eskf_rtk.set_magnet_standard_deviation(eskf::noise_std_proc_mag);
    eskf_rtk.set_drift_magnetometer_standard_deviation(eskf::noise_std_drift_mag);
    eskf_rtk.set_processing_standard_deviation(eskf::noise_std_proc);
    eskf_rtk.enable_estimation_acc_bias();
    eskf_rtk.enable_estimation_gravity();
    eskf_rtk.enable_estimation_magnet();
    eskf_rtk.enable_estimation_magnet_bias();
    // eskf_rtk.disable_estimation_magnet();
    // eskf_rtk.disable_estimation_magnet_bias();
    eskf_rtk.disable_estimation_wind();
    eskf_rtk.initialize();

    clock_t t1 = clock();
    for (unsigned i = 0; i < num_step; ++i) {
        eskf_rtk.predict_state(w_meas[i], a_meas[i]);
        eskf_rtk.predict_covariance(w_meas[i], a_meas[i]);
        unsigned char info;
        info = eskf_rtk.fuse_position(pl_meas[i], w_meas[i], a_meas[i], eskf::d_gps_left, eskf::noise_std_rtk_pos, eskf::gate_rtk_pos);
        if (info != 0) {
            cout << "pl: " << int(info) << endl;
        }
        info = eskf_rtk.fuse_velocity(vl_meas[i], w_meas[i], a_meas[i], eskf::d_gps_left, eskf::noise_std_gps_vel, eskf::gate_gps_vel);
        if (info != 0) {
            cout << "vl: " << int(info) << endl;
        }
        info = eskf_rtk.fuse_position(pr_meas[i], w_meas[i], a_meas[i], eskf::d_gps_right, eskf::noise_std_rtk_pos, eskf::gate_rtk_pos);
        if (info != 0) {
            cout << "pr: " << int(info) << endl;
        }
        info = eskf_rtk.fuse_velocity(vr_meas[i], w_meas[i], a_meas[i], eskf::d_gps_right, eskf::noise_std_gps_vel, eskf::gate_gps_vel);
        if (info != 0) {
            cout << "vr: " << int(info) << endl;
        }
        info = eskf_rtk.fuse_magnet(mb_meas[i], w_meas[i], a_meas[i], eskf::noise_std_mag, eskf::gate_mag);
        if (info != 0) {
            cout << "m: " << int(info) << endl;
        }
        // info = eskf_rtk.fuse_declination(declination, mb_meas[i], w_meas[i], a_meas[i], eskf::noise_std_dec, eskf::gate_dec);
        // if (info != 0) {
        //     cout << "m: " << int(info) << endl;
        // }
        // info = eskf_rtk.fuse_magnet_1D(declination, mb_meas[i], w_meas[i], a_meas[i], eskf::noise_std_mag_1d, eskf::gate_mag_1d);
        // if (info != 0) {
        //     cout << "m: " << int(info) << endl;
        // }
        // eskf_rtk.correct_covariance();
        eskf_rtk.correct_state();

        q_hat[i] = eskf_rtk.get_quaternion();   
        p_hat[i] = eskf_rtk.get_position();
        v_hat[i] = eskf_rtk.get_velocity();
        bg_hat[i] = eskf_rtk.get_drift_gyro();
        ba_hat[i] = eskf_rtk.get_drift_acc();
        g_hat[i] = eskf_rtk.get_gravity();
        mb_hat[i] = q_hat[i].toRotationMatrix().transpose() * eskf_rtk.get_magnet(); 
        bm_hat[i] = eskf_rtk.get_drift_magnet();

        // cout << g_hat[i] << endl;

        // cout << eskf_rtk.get_wind() << endl;

        // cout << p_true[i].transpose() << ", " << p_hat[i].transpose() << endl;

        // cout << mb_true[i].transpose() << ", " << mb_hat[i].transpose() << endl;

        cout << M.transpose() << ", " << eskf_rtk.get_magnet().transpose() << endl;
    }
    // clock_t t2 = clock();
    // cout << "time: " << t2 - t1 << endl;

    return 0;
}
