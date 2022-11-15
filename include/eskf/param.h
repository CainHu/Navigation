#ifndef NAVIGATION_ESKF_PARAM_H
#define NAVIGATION_ESKF_PARAM_H

#include <array>

namespace eskf {
    using namespace std;

    // Standard deviation of noise
    array<float, 3> noise_std_rtk_pos {0.01f, 0.01f, 0.01f};
    array<float, 3> noise_std_gps_pos {0.1f, 0.1f, 0.1f};
    array<float, 3> noise_std_gps_vel {0.01f, 0.01f, 0.01f};
    array<float, 3> noise_std_vision_pos {0.01f, 0.01f, 0.01f};
    array<float, 3> noise_std_vision_vel {0.01f, 0.01f, 0.01f};
    array<float, 3> noise_std_gyro {0.01f, 0.01f, 0.01f};
    array<float, 3> noise_std_acc {0.1f, 0.1f, 0.1f};
    array<float, 3> noise_std_sat_gyro {1.f, 1.f, 1.f};
    array<float, 3> noise_std_sat_acc {10.f, 10.f, 10.f};
    array<float, 3> noise_std_drift_gyro {0.001f, 0.001f, 0.001f};
    array<float, 3> noise_std_drift_acc {0.0001f, 0.0001f, 0.0001f};
    float noise_std_baro {0.2f};
    float noise_std_grav {0.0001f};
    float noise_std_proc {0.00001f};

    // Direction vector
    array<float, 3> d_gps_left {-0.2f, -0.15f, -0.02f};
    array<float, 3> d_gps_right {-0.2f, 0.15f, -0.02f};
    array<float, 3> d_barometer {0.f, 0.f, 0.01f};
    array<float, 3> d_vision {0.2f, 0.f, 0.02f};

    // Gate of observation's error
    array<float, 3> gate_rtk_pos {1000.f, 1000.f, 1000.f};
    array<float, 3> gate_gps_pos {10000.f, 10000.f, 10000.f};
    array<float, 3> gate_gps_vel {1000.f, 1000.f, 1000.f};
    array<float, 3> gate_vision_pos {1000.f, 1000.f, 1000.f};
    array<float, 3> gate_vision_vel {1000.f, 1000.f, 1000.f};
    float gate_baro {100.f};
}

#endif // NAVIGATION_ESKF_PARAM_H