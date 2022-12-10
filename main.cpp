//#pragma GCC optimize(2)

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <ctime>
#include <Eigen/Dense>
#include <algorithm>
#include <random>
// #include "eskf/include/common.h"
//#include "leskf.h"
//#include "liekf.h"
//#include "riekf.h"

// #include "param.h"

#include "eskf_runner.h"
#include "geskf.h"
#include "common.h"
#include "generator.h"
#include "gvar.h"
#include "gen_utils.h"

using namespace std;
using namespace Eigen;

generator::Generator env;
vector<unsigned long> time_us_seq;

int main() {
    // // eskf::ESKF eskf_test(0.001f);
    // // eskf_test.reset_priori_covariance_matrix();
    // // eskf_test.set_gyroscope_standard_deviation(eskf::noise_std_gyro);
    // // eskf_test.set_accelerometer_standard_deviation(eskf::noise_std_acc);
    // // eskf_test.set_drift_gyroscope_standard_deviation(eskf::noise_std_drift_gyro);
    // // eskf_test.set_drift_saccelerometer_tandard_deviation(eskf::noise_std_drift_acc);
    // // eskf_test.set_gravity_standard_deviation(eskf::noise_std_grav);
    // // eskf_test.set_processing_standard_deviation(eskf::noise_std_proc);
    // // eskf_test.initialize();

    // // Vector3f w = {1.f, 2.f, 3.f};
    // // Vector3f a = {-3.f, -2.f, -9.8f - 1.f};

    
    // // for (int i = 0; i < 50; ++i) {
    // //     eskf_test.predict_state(w, a);
    // //     eskf_test.predict_covariance(w, a);
    // // }
    // // const array<array<float, 16>, 16> cov = eskf_test.get_covariance_matrix();
    // // for (int i = 0; i < 16; ++i) {
    // //     for (int j = 0; j < 16; ++j) {
    // //         cout << cov[i][j] * 1000000.f << ", ";
    // //     }
    // //     cout << endl;
    // // }









    // static const unsigned int num_step = 100000;
    // static const float g = 9.80665f;

    // static const float dt = 0.0001f;
    // float ts = 0.001f;

    // Vector3f G(0.f, 0.f, g);
    // Vector3f M(0.9f, sqrtf(1.f - 0.9f * 0.9f), 0.f);
    // Vector2f declination(-asinf(M[2]), atan2f(M[1], M[0]));

    // Vector3f p(0.f, 0.f, 0.f), v(0.f, 0.f, 0.f), bg(0.f, 0.f, 0.f), ba(0.f, 0.f, 0.f), bm(0.f, 0.f, 0.f);
    // Quaternionf q(1.f, 0.f, 0.f, 0.f);
    // Vector3f w(0.f, 0.f, 0.f), a(0.f, 0.f, -g);

    // static array<Vector3f, num_step> p_true, v_true, bg_true, ba_true, mb_true, bm_true;
    // static array<Quaternionf, num_step> q_true;
    // static array<Vector3f, num_step> w_true, a_true;

    // static array<Vector3f, num_step> pl_meas, vl_meas, pr_meas, vr_meas, mb_meas;
    // static array<Vector3f, num_step> w_meas, a_meas;

    // static array<Vector3f, num_step> p_hat, v_hat, bg_hat, ba_hat, mb_hat, bm_hat;
    // static array<Quaternionf, num_step> q_hat;
    // static array<float, num_step> g_hat;

    // default_random_engine random_engine;
    // normal_distribution<float> dist(0.f, 1.f);

    // Vector3f dl(eskf::d_gps_left[0], eskf::d_gps_left[1], eskf::d_gps_left[2]);
    // Vector3f dr(eskf::d_gps_right[0], eskf::d_gps_right[1], eskf::d_gps_right[2]);

    //  // Simulation
    // unsigned int j = 0;
    // // w += 3.f * Vector3f(dist(random_engine), dist(random_engine), dist(random_engine)) * ts;
    // // a += 10.f * Vector3f(dist(random_engine), dist(random_engine), dist(random_engine)) * ts;
    // w = Vector3f::Zero();
    // a = Vector3f(1.f, 0.f, -g);
    // for (unsigned int i = 0; i < num_step * (unsigned int)(ts / dt); ++i) {
    //     Vector3f v_last = v;
    //     Matrix3f R = q.toRotationMatrix();
    //     v += (R * a + G) * dt;
    //     p += 0.5f * (v + v_last) * dt;

    //     Vector3f d_ang = w * dt;
    //     AngleAxisf angle_axis(d_ang.norm(), d_ang.normalized());
    //     Quaternionf dq(angle_axis);
    //     q = q * dq;
    //     q.normalize();

    //     Vector3f nbg(eskf::noise_std_drift_gyro[0] * dist(random_engine),
    //                  eskf::noise_std_drift_gyro[1] * dist(random_engine),
    //                  eskf::noise_std_drift_gyro[2] * dist(random_engine));

    //     Vector3f nba(eskf::noise_std_drift_acc[0] * dist(random_engine),
    //                  eskf::noise_std_drift_acc[1] * dist(random_engine),
    //                  eskf::noise_std_drift_acc[2] * dist(random_engine));

    //     Vector3f nbm(eskf::noise_std_drift_mag[0] * dist(random_engine),
    //                  eskf::noise_std_drift_mag[1] * dist(random_engine),
    //                  eskf::noise_std_drift_mag[2] * dist(random_engine));

    //     bg += nbg * dt;
    //     ba += nba * dt;
    //     bm += nbm * dt;

    //     if ((i + 1) % (unsigned int)(ts / dt) == 0) {
    //         // Real State
    //         p_true[j] = p;
    //         v_true[j] = v;
    //         bg_true[j] = bg;
    //         ba_true[j] = ba;
    //         mb_true[j] = q.toRotationMatrix().transpose() * M;
    //         bm_true[j] = bm;

    //         q_true[j] = q;

    //         w_true[j] = w;
    //         a_true[j] = a;

    //         // Measurement
    //         Vector3f npl(dist(random_engine) * eskf::noise_std_rtk_pos[0],
    //                      dist(random_engine) * eskf::noise_std_rtk_pos[1],
    //                      dist(random_engine) * eskf::noise_std_rtk_pos[2]);
    //         Vector3f nvl(dist(random_engine) * eskf::noise_std_gps_vel[0],
    //                      dist(random_engine) * eskf::noise_std_gps_vel[1],
    //                      dist(random_engine) * eskf::noise_std_gps_vel[2]);
    //         Vector3f npr(dist(random_engine) * eskf::noise_std_rtk_pos[0],
    //                      dist(random_engine) * eskf::noise_std_rtk_pos[1],
    //                      dist(random_engine) * eskf::noise_std_rtk_pos[2]);
    //         Vector3f nvr(dist(random_engine) * eskf::noise_std_gps_vel[0],
    //                      dist(random_engine) * eskf::noise_std_gps_vel[1],
    //                      dist(random_engine) * eskf::noise_std_gps_vel[2]);
    //         Vector3f nmb(dist(random_engine) * eskf::noise_std_mag[0],
    //                      dist(random_engine) * eskf::noise_std_mag[1],
    //                      dist(random_engine) * eskf::noise_std_mag[2]);
    //         Matrix3f R = q_true[j].toRotationMatrix();            
    //         pl_meas[j] = p_true[j] + R * dl + npl;
    //         vl_meas[j] = v_true[j] - R * (dl.cross(w_true[j])) + nvl;
    //         pr_meas[j] = p_true[j] + R * dr + npr;
    //         vr_meas[j] = v_true[j] - R * (dr.cross(w_true[j])) + nvr;
    //         mb_meas[j] = mb_true[j] + nmb;

    //         Vector3f ng(dist(random_engine) * eskf::noise_std_gyro[0],
    //                     dist(random_engine) * eskf::noise_std_gyro[1],
    //                     dist(random_engine) * eskf::noise_std_gyro[2]);
    //         Vector3f na(dist(random_engine) * eskf::noise_std_acc[0],
    //                     dist(random_engine) * eskf::noise_std_acc[1],
    //                     dist(random_engine) * eskf::noise_std_acc[2]);
    //         w_meas[j] = w_true[j] + bg + ng;
    //         a_meas[j] = a_true[j] + ba + na;
        
    //         // Simulation
    //         // w = 1.f * Vector3f(dist(random_engine), dist(random_engine), dist(random_engine));
    //         // a = 3.f * Vector3f(dist(random_engine), dist(random_engine), dist(random_engine)) - 3.f * v_true[j];
    //         ++j;
    //     }
    // }

    // // leskf::LESKF eskf_rtk(ts);
    // geskf::GESKF eskf_rtk(ts);
    // // liekf::LIEKF eskf_rtk(ts);
    // // riekf::RIEKF eskf_rtk(ts);
    // // eskf_rtk.set_magnet(M);
    // eskf_rtk.set_gyroscope_standard_deviation(eskf::noise_std_gyro);
    // eskf_rtk.set_accelerometer_standard_deviation(eskf::noise_std_acc);
    // eskf_rtk.set_drift_gyroscope_standard_deviation(eskf::noise_std_drift_gyro);
    // eskf_rtk.set_drift_saccelerometer_tandard_deviation(eskf::noise_std_drift_acc);
    // eskf_rtk.set_gravity_standard_deviation(eskf::noise_std_proc_grav);
    // eskf_rtk.set_magnet_standard_deviation(eskf::noise_std_proc_mag);
    // eskf_rtk.set_declination_standard_deviation(eskf::noise_std_dec);
    // eskf_rtk.set_drift_magnetometer_standard_deviation(eskf::noise_std_drift_mag);
    // eskf_rtk.set_processing_standard_deviation(eskf::noise_std_proc);
    // eskf_rtk.enable_estimation_acc_bias();
    // eskf_rtk.enable_estimation_gravity();
    // eskf_rtk.enable_estimation_magnet();
    // eskf_rtk.enable_estimation_declination();
    // eskf_rtk.enable_estimation_magnet_bias();
    // // eskf_rtk.disable_estimation_magnet();
    // // eskf_rtk.disable_estimation_magnet_bias();
    // eskf_rtk.disable_estimation_wind();
    // eskf_rtk.initialize();

    // clock_t t1 = clock();
    // for (unsigned i = 0; i < num_step; ++i) {
    //     eskf_rtk.predict_state(w_meas[i], a_meas[i]);
    //     eskf_rtk.predict_covariance(w_meas[i], a_meas[i]);
    //     unsigned char info;
    //     info = eskf_rtk.fuse_position(pl_meas[i], w_meas[i], a_meas[i], eskf::d_gps_left, eskf::noise_std_rtk_pos, eskf::gate_rtk_pos);
    //     if (info != 0) {
    //         cout << "pl: " << int(info) << endl;
    //     }
    //     info = eskf_rtk.fuse_velocity(vl_meas[i], w_meas[i], a_meas[i], eskf::d_gps_left, eskf::noise_std_gps_vel, eskf::gate_gps_vel);
    //     if (info != 0) {
    //         cout << "vl: " << int(info) << endl;
    //     }
    //     info = eskf_rtk.fuse_position(pr_meas[i], w_meas[i], a_meas[i], eskf::d_gps_right, eskf::noise_std_rtk_pos, eskf::gate_rtk_pos);
    //     if (info != 0) {
    //         cout << "pr: " << int(info) << endl;
    //     }
    //     info = eskf_rtk.fuse_velocity(vr_meas[i], w_meas[i], a_meas[i], eskf::d_gps_right, eskf::noise_std_gps_vel, eskf::gate_gps_vel);
    //     if (info != 0) {
    //         cout << "vr: " << int(info) << endl;
    //     }
    //     info = eskf_rtk.fuse_magnet(mb_meas[i], w_meas[i], a_meas[i], eskf::noise_std_mag, eskf::gate_mag);
    //     if (info != 0) {
    //         cout << "mag: " << int(info) << endl;
    //     }
    //     info = eskf_rtk.fuse_declination(declination, w_meas[i], a_meas[i], eskf::noise_std_dec, eskf::gate_dec);
    //     if (info != 0) {
    //         cout << "dec: " << int(info) << endl;
    //     }
    //     // eskf_rtk.correct_covariance();
    //     eskf_rtk.correct_state();

    //     q_hat[i] = eskf_rtk.get_quaternion();   
    //     p_hat[i] = eskf_rtk.get_position();
    //     v_hat[i] = eskf_rtk.get_velocity();
    //     bg_hat[i] = eskf_rtk.get_drift_gyro();
    //     ba_hat[i] = eskf_rtk.get_drift_acc();
    //     g_hat[i] = eskf_rtk.get_gravity();
    //     // mb_hat[i] = q_hat[i].toRotationMatrix().transpose() * eskf_rtk.get_magnet(); 
    //     bm_hat[i] = eskf_rtk.get_drift_magnet();

    //     float h = eskf_rtk.get_magnet();
    //     const Vector2f &dec = eskf_rtk.get_declination();

    //     // cout << g_hat[i] << endl;

    //     // cout << eskf_rtk.get_wind() << endl;

    //     // cout << p_true[i].transposce() << ", " << p_hat[i].transpose() << endl;

    //     // cout << mb_true[i].transpose() << ", " << mb_hat[i].transpose() << endl;

    //     // cout << M.transpose() << ", " << eskf_rtk.get_magnet().transpose() << endl;

    //     // cout << "H: (" << 1.f << ", " << h << "), dec: (" << declination[0] << ", " << dec[0] << "), (" << declination[1] << ", " << dec[1] << ")" << endl;

    //     // cout << q_true[i].y() << ", " << q_hat[i].y() << endl;
    //     cout << generator::quat2euler(q_hat[i].cast<double>()).transpose() << endl;
    // }
    // // clock_t t2 = clock();
    // // cout << "time: " << t2 - t1 << endl;



    // 飞行数据生成
    float ts = 0.001;
    double ts_sim = 0.001;

    // leskf::LESKF eskf_rtk(ts);
    eskf::GESKF eskf_rtk(ts);
    eskf::ESKFRunner<100> eskf_runner(eskf_rtk);
    // liekf::LIEKF eskf_rtk(ts);
    // riekf::RIEKF eskf_rtk(ts);
    // eskf_rtk.set_magnet(M);

    Vector3d euler0(0., 0., 0.), vn0(0., 0., 0.), pos0(23.1659394 * arcdeg, 113.4522718 * arcdeg, 20.);
    vector<array<double, 5>> wat(13);
    wat[0][0] = 0., wat[0][1] = 0., wat[0][2] = 0., wat[0][3] = 0., wat[0][4] = 10.;            // 静止
    wat[1][0] = 0., wat[1][1] = 0., wat[1][2] = 0., wat[1][3] = 1., wat[1][4] = 10.;            // 加速
    wat[2][0] = 0., wat[2][1] = 0., wat[2][2] = 0., wat[2][3] = 0., wat[2][4] = 10.;            // 匀速
    wat[3][0] = 0., wat[3][1] = 5., wat[3][2] = 0., wat[3][3] = 0., wat[3][4] = 4.;             // 抬头
    wat[4][0] = 0., wat[4][1] = 0., wat[4][2] = 0., wat[4][3] = 0., wat[4][4] = 10.;            // 匀速
    wat[5][0] = 0., wat[5][1] = -5, wat[5][2] = 0., wat[5][3] = 0., wat[5][4] = 4.;             // 低头
    wat[6][0] = 0., wat[6][1] = 0., wat[6][2] = 0., wat[6][3] = 0., wat[6][4] = 10.;            // 匀速
    wat[7][0] = 10., wat[7][1] = 0., wat[7][2] = 0., wat[7][3] = 0., wat[7][4] = 1.;            // 横滚
    wat[8][0] = 0., wat[8][1] = 0., wat[8][2] = 9., wat[8][3] = 0., wat[8][4] = 10.;            // 转弯
    wat[9][0] = -10., wat[9][1] = 0., wat[9][2] = 0., wat[9][3] = 0., wat[9][4] = 1.;           // 横滚
    wat[10][0] = 0., wat[10][1] = 0., wat[10][2] = 0., wat[10][3] = 0., wat[10][4] = 10.;       // 匀速
    wat[11][0] = 0., wat[11][1] = 0., wat[11][2] = 0., wat[11][3] = -1., wat[11][4] = 10.;      // 减速
    wat[12][0] = 0., wat[12][1] = 0., wat[12][2] = 0., wat[12][3] = 0., wat[12][4] = 10.;       // 静止

    for (unsigned int i = 0; i < wat.size(); ++i) {
        wat[i][0] *= arcdeg;
        wat[i][1] *= arcdeg;
        wat[i][2] *= arcdeg;
    }

    // pve数据
    env.traj_gen(euler0, vn0, pos0, wat, ts_sim);
    const vector<Vector3d> &pos = env.get_pos();
    const vector<Vector3d> &vn = env.get_vn();
    const vector<Vector3d> &euler = env.get_euler();

    // imu数据
    env.ev2imu(euler, vn, pos, ts_sim);
    const vector<Vector3d> &wm = env.get_wm();
    const vector<Vector3d> &vm = env.get_vm();

    // 量测imu数据
    const Vector3d bias_gyro(0.001, 0.001, 0.001);
    const Vector3d bias_acc(0.0001, 0.0001, 0.0001);
    const Vector3d std_gyro(0.01, 0.01, 0.01);
    const Vector3d std_acc(0.1, 0.1, 0.1);
    static vector<Vector3d> wm_meas;
    static vector<Vector3d> vm_meas;
    wm_meas.resize(wm.size());
    vm_meas.resize(vm.size());
    env.imu_add_err(wm, vm, bias_gyro, std_gyro, bias_acc, std_acc, ts_sim, wm_meas, vm_meas);

    // 时间序列
    time_us_seq.resize(pos.size());
    time_us_seq[0] = 0;
    unsigned long dt_us = (unsigned long)(ts_sim * 1e6f);
    for (unsigned int i = 1; i < pos.size(); ++i) {
        time_us_seq[i] = time_us_seq[i - 1] + dt_us;
    }

    // 地球坐标位置转导航坐标
    static vector<Vector3d> pn;
    pn.resize(pos.size());
    for (unsigned int i = 0; i < pos.size(); ++i) {
        pn[i] = diff_geo(pos[i], pos0);
        // if (i > 0) {
        //     pn[i] = pn[i - 1] + 0.5 * (vn[i - 1] + vn[i]) * ts_sim;
        // } else {
        //     pn[i].setZero();
        // }
        // // cout << pn[i].transpose() << ", " << diff_geo(pos[i], pos0).transpose() << endl;
    }

    // 位置和速度量测
    static vector<Vector3f> pl_meas, vl_meas, pr_meas, vr_meas, mb_meas;
    static vector<Vector3f> w_meas, a_meas;
    pl_meas.resize(pos.size()), vl_meas.resize(pos.size()), pr_meas.resize(pos.size()), vr_meas.resize(pos.size()), mb_meas.resize(pos.size());
    w_meas.resize(pos.size()), a_meas.resize(pos.size());

    static vector<Vector3f> p_hat, v_hat, bg_hat, ba_hat, mb_hat, bm_hat;
    static vector<Vector3f> euler_hat;
    static vector<float> g_hat;
    p_hat.resize(pos.size()), v_hat.resize(pos.size()), bg_hat.resize(pos.size()), ba_hat.resize(pos.size()), mb_hat.resize(pos.size()), bm_hat.resize(pos.size());
    euler_hat.resize(pos.size());
    g_hat.resize(pos.size());

    default_random_engine random_engine;
    normal_distribution<float> dist(0.f, 1.f);
    const Vector3f dl(eskf_runner._params.d_gps_left[0], eskf_runner._params.d_gps_left[1], eskf_runner._params.d_gps_left[2]);
    const Vector3f dr(eskf_runner._params.d_gps_right[0], eskf_runner._params.d_gps_right[1], eskf_runner._params.d_gps_right[2]);
    for (unsigned int i = 0; i < pos.size(); ++i) {
        Vector3f npl(dist(random_engine) * eskf_runner._params.noise_std_pos_rtk[0],
                        dist(random_engine) * eskf_runner._params.noise_std_pos_rtk[1],
                        dist(random_engine) * eskf_runner._params.noise_std_pos_rtk[2]);
        Vector3f nvl(dist(random_engine) * eskf_runner._params.noise_std_vel_gps[0],
                        dist(random_engine) * eskf_runner._params.noise_std_vel_gps[1],
                        dist(random_engine) * eskf_runner._params.noise_std_vel_gps[2]);
        Vector3f npr(dist(random_engine) * eskf_runner._params.noise_std_pos_rtk[0],
                        dist(random_engine) * eskf_runner._params.noise_std_pos_rtk[1],
                        dist(random_engine) * eskf_runner._params.noise_std_pos_rtk[2]);
        Vector3f nvr(dist(random_engine) * eskf_runner._params.noise_std_vel_gps[0],
                        dist(random_engine) * eskf_runner._params.noise_std_vel_gps[1],
                        dist(random_engine) * eskf_runner._params.noise_std_vel_gps[2]);

        Matrix3f R = generator::euler2rot(euler[i]).cast<float>();                

        pl_meas[i] = pn[i].cast<float>() + R * dl + npl;
        vl_meas[i] = vn[i].cast<float>() - R * (dl.cross(wm[i].cast<float>() / ts)) + nvl;
        pr_meas[i] = pn[i].cast<float>() + R * dr + npr;
        vr_meas[i] = vn[i].cast<float>() - R * (dr.cross(wm[i].cast<float>() / ts)) + nvr;
        w_meas[i] = wm_meas[i].cast<float>() / ts;
        a_meas[i] = vm_meas[i].cast<float>() / ts;
    }

    // 数据融合
//    eskf_rtk.set_std_proc_pos(eskf_runner._params.noise_std_proc);
    eskf_rtk.set_std_proc_vel(eskf_runner._params.noise_std_acc);
    eskf_rtk.set_std_proc_ang(eskf_runner._params.noise_std_gyro);
    eskf_rtk.set_std_proc_bias_gyro(eskf_runner._params.noise_std_proc_bias_gyro);
    eskf_rtk.set_std_proc_bias_acc(eskf_runner._params.noise_std_proc_bias_acc);
    eskf_rtk.set_std_proc_grav(eskf_runner._params.noise_std_proc_grav);
    eskf_rtk.set_std_proc_mag_norm(eskf_runner._params.noise_std_proc_mag_norm);
    eskf_rtk.set_std_proc_mag_dec(eskf_runner._params.noise_std_proc_mag_dec);
    eskf_rtk.set_std_proc_bias_mag(eskf_runner._params.noise_std_proc_bias_mag);
    eskf_rtk.set_std_proc_wind(eskf_runner._params.noise_std_proc_wind);

    eskf_rtk.enable_estimation_acc_bias();
    eskf_rtk.enable_estimation_gravity();
    // eskf_rtk.enable_estimation_magnet();
    // eskf_rtk.enable_estimation_declination();
    // eskf_rtk.enable_estimation_magnet_bias();
    eskf_rtk.disable_estimation_magnet();
    eskf_rtk.disable_estimation_declination();
    eskf_rtk.disable_estimation_magnet_bias();
    eskf_rtk.disable_estimation_wind();
    eskf_rtk.initialize();

    clock_t t1 = clock();
    for (unsigned i = 0; i < pos.size(); ++i) {
        Vector3f w, a;
        w = wm[i].cast<float>() / ts;
        a = vm[i].cast<float>() / ts;

        eskf_runner.set_imu_data(w, a, time_us_seq[i]);
        if (i > 99) {
            eskf_runner.set_gps_data(pl_meas[i - 100], vl_meas[i - 100], pr_meas[i - 100], vr_meas[i - 100], time_us_seq[i - 100]);
            eskf_runner.update();
//            cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        }
        const eskf::OutputSample &state = eskf_runner.get_output_state();

        euler_hat[i] = generator::quat2euler(state.q.cast<double>()).cast<float>();
        p_hat[i] = eskf_rtk.get_position();
        v_hat[i] = eskf_rtk.get_velocity();
        bg_hat[i] = eskf_rtk.get_drift_gyro();
        ba_hat[i] = eskf_rtk.get_drift_acc();
        g_hat[i] = eskf_rtk.get_gravity();
        // mb_hat[i] = q_hat[i].toRotationMatrix().transpose() * eskf_rtk.get_magnet();
        // bm_hat[i] = eskf_rtk.get_drift_magnet();

        // float h = eskf_rtk.get_magnet();
        // const Vector2f &dec = eskf_rtk.get_declination();

        // cout << g_hat[i] << endl;

        // cout << eskf_rtk.get_wind() << endl;

         cout << pn[i].transpose() << ", " << p_hat[i].transpose() << ", " << state.p.transpose() << endl;
        // cout << vn[i].transpose() << ", " << v_hat[i].transpose() << endl;
//        cout << euler[i].transpose() << ", " << euler_hat[i].transpose() << endl;

    }


//    clock_t t1 = clock();
//    Vector3f accum_w = Vector3f::Zero();
//    for (unsigned i = 0; i < pos.size(); ++i) {
//        Vector3f w, a;
//        if (i > 0) {
//            w = (wm[i] + 1./12. * wm[i - 1].cross(wm[i])).cast<float>() / ts;
//            a = (vm[i] + 0.5 * wm[i].cross(vm[i]) + 1./12. * (wm[i - 1].cross(vm[i]) + vm[i - 1].cross(wm[i]))).cast<float>() / ts;
//        } else {
//            w = wm[i].cast<float>() / ts;
//            a = vm[i].cast<float>() / ts;
//        }
//
//
//        eskf_rtk.predict_state(w, a);
//        eskf_rtk.predict_covariance(w, a);
//        unsigned char info;
//        info = eskf_rtk.fuse_position(pl_meas[i], w, a, eskf_runner._params.d_gps_left, eskf_runner._params.noise_std_pos_rtk, eskf_runner._params.gate_rtk_pos);
//        if (info != 0) {
//            cout << "pl: " << int(info) << endl;
//        }
//        info = eskf_rtk.fuse_velocity(vl_meas[i], w, a, eskf_runner._params.d_gps_left, eskf_runner._params.noise_std_vel_gps, eskf_runner._params.gate_gps_vel);
//        if (info != 0) {
//            cout << "vl: " << int(info) << endl;
//        }
//        info = eskf_rtk.fuse_position(pr_meas[i], w, a, eskf_runner._params.d_gps_right, eskf_runner._params.noise_std_pos_rtk, eskf_runner._params.gate_rtk_pos);
//        if (info != 0) {
//            cout << "pr: " << int(info) << endl;
//        }
//        info = eskf_rtk.fuse_velocity(vr_meas[i], w, a, eskf_runner._params.d_gps_right, eskf_runner._params.noise_std_vel_gps, eskf_runner._params.gate_gps_vel);
//        if (info != 0) {
//            cout << "vr: " << int(info) << endl;
//        }
//        // info = eskf_rtk.fuse_magnet(mb_meas[i], w_meas[i], a_meas[i], eskf::noise_std_mag, eskf::gate_mag);
//        // if (info != 0) {
//        //     cout << "mag: " << int(info) << endl;
//        // }
//        // info = eskf_rtk.fuse_declination(declination, w_meas[i], a_meas[i], eskf::noise_std_dec, eskf::gate_dec);
//        // if (info != 0) {
//        //     cout << "dec: " << int(info) << endl;
//        // }
//        // eskf_rtk.correct_covariance();
//        eskf_rtk.correct_state();
//
//        euler_hat[i] = generator::quat2euler(eskf_rtk.get_quaternion().cast<double>()).cast<float>();
//        p_hat[i] = eskf_rtk.get_position();
//        v_hat[i] = eskf_rtk.get_velocity();
//        bg_hat[i] = eskf_rtk.get_drift_gyro();
//        ba_hat[i] = eskf_rtk.get_drift_acc();
//        g_hat[i] = eskf_rtk.get_gravity();
//        // mb_hat[i] = q_hat[i].toRotationMatrix().transpose() * eskf_rtk.get_magnet();
//        // bm_hat[i] = eskf_rtk.get_drift_magnet();
//
//        // float h = eskf_rtk.get_magnet();
//        // const Vector2f &dec = eskf_rtk.get_declination();
//
//        // cout << g_hat[i] << endl;
//
//        // cout << eskf_rtk.get_wind() << endl;
//
//        // cout << pn[i].transpose() << ", " << p_hat[i].transpose() << endl;
//        // cout << vn[i].transpose() << ", " << v_hat[i].transpose() << endl;
//        cout << euler[i].transpose() << ", " << euler_hat[i].transpose() << endl;
//        // cout << bg_hat[i].transpose() << ", " << ba_hat[i].transpose() << ", " << g_hat[i] << endl;
//        // if (i > 0){
//        //     cout << (wm[i] + 1./12. * wm[i - 1].cross(wm[i])) / ts << endl;
//        // }
//        // accum_w += w * ts;
//        // cout << accum_w.transpose() << endl;
//
//        // cout << pn[i].transpose() << ", " << vn[i].transpose() << endl;
//
//        // cout << wm[i].transpose() / ts_sim << ", " << vm[i].transpose() / ts_sim << endl;
//
//        // cout << mb_true[i].transpose() << ", " << mb_hat[i].transpose() << endl;
//
//        // cout << M.transpose() << ", " << eskf_rtk.get_magnet().transpose() << endl;
//
//        // cout << "H: (" << 1.f << ", " << h << "), dec: (" << declination[0] << ", " << dec[0] << "), (" << declination[1] << ", " << dec[1] << ")" << endl;
//
//        // cout << q_true[i].z() << ", " << q_hat[i].z() << endl;
//    }
//    // clock_t t2 = clock();
//    // cout << "time: " << t2 - t1 << endl;


    return 0;
}
