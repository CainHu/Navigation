#ifndef NAVIGATION_ESKF_COMMON_H
#define NAVIGATION_ESKF_COMMON_H

#include <Eigen/Dense>
#include <array>

namespace eskf {
    using namespace std;
    using namespace Eigen;

    struct parameters {
        // Standard deviation of noise
        Vector3f noise_std_rtk_pos {0.01f, 0.01f, 0.01f};
        Vector3f noise_std_gps_pos {0.1f, 0.1f, 0.1f};
        Vector3f noise_std_gps_vel {0.01f, 0.01f, 0.01f};
        Vector3f noise_std_vision_pos {0.01f, 0.01f, 0.01f};
        Vector3f noise_std_vision_vel {0.01f, 0.01f, 0.01f};
        Vector3f noise_std_gyro {0.01f, 0.01f, 0.01f};
        Vector3f noise_std_acc {0.1f, 0.1f, 0.1f};
        Vector3f noise_std_sat_gyro {1.f, 1.f, 1.f};
        Vector3f noise_std_sat_acc {10.f, 10.f, 10.f};
        Vector3f noise_std_drift_gyro {0.001f, 0.001f, 0.001f};
        Vector3f noise_std_drift_acc {0.0001f, 0.0001f, 0.0001f};
        Vector3f noise_std_mag {0.1f, 0.1f, 0.1f};
        float noise_std_proc_mag {0.01f};
        Vector2f noise_std_dec {0.01f, 0.01f};
        Vector3f noise_std_drift_mag {0.0001f, 0.0001f, 0.0001f};
        float noise_std_baro {0.2f};
        float noise_std_proc_grav {0.005f};
        float noise_std_proc {0.00001f};

        // Direction vector
        Vector3f d_gps_left {-0.2f, -0.15f, -0.02f};
        Vector3f d_gps_right {-0.2f, 0.15f, -0.02f};
        Vector3f d_barometer {0.f, 0.f, 0.01f};
        Vector3f d_vision {0.2f, 0.f, 0.02f};

        // Gate of observation's error
        Vector3f gate_rtk_pos {1000.f, 1000.f, 1000.f};
        Vector3f gate_gps_pos {10000.f, 10000.f, 10000.f};
        Vector3f gate_gps_vel {1000.f, 1000.f, 1000.f};
        Vector3f gate_vision_pos {1000.f, 1000.f, 1000.f};
        Vector3f gate_vision_vel {1000.f, 1000.f, 1000.f};
        Vector3f gate_mag {1000.f, 1000.f, 1000.f};
        Vector2f gate_dec {1000.f, 1000.f};
        float gate_baro {100.f};
    };

    union filter_control_status {
        struct {
            bool acc_x_bias : 1;
            bool acc_y_bias : 1;
            bool acc_z_bias : 1;
            bool grav : 1;
            bool mag : 1;
            bool dec : 1;
            bool mag_bias : 1;
            bool wind : 1;
        } flags;
        unsigned char status;
    };

    union innovation_fault_status {
        struct {
            bool reject_pos : 1;
            bool reject_vel : 1;
            bool reject_mag : 1;
            bool reject_dec : 1;
            bool reject_wind : 1;
        } flags;
        unsigned char status;
    };

    struct State {
        Vector3f pos;
        Vector3f vel;
        Quaternionf quat;
        Vector3f bias_gyro;
        Vector3f bias_acc;
        float grav;
        float mag;
        Vector2f dec;
        Vector3f bias_mag;
        Vector2f wind;

        Matrix3f rot;
    };

    struct ErrorState {
        Vector3f pos;
        Vector3f vel;
        Vector3f ang;
        Vector3f bias_gyro;
        Vector3f bias_acc;
        float g;
        float mag;
        Vector2f dec;
        Vector3f bias_mag;
        Vector2f wind;
    };

    struct Sample {
        unsigned long time_us {0};
    };

    struct OutputSample : Sample {
        Quaternionf q;
        Matrix3f r;
        Vector3f v;
        Vector3f p;
    };

    struct ImuSample : Sample {
        Vector3f delta_ang {};
        Vector3f delta_vel {};
        float delta_ang_dt {0.f};
        float delta_vel_dt {0.f};
        array<bool, 3> delta_vel_clipping;
    };

    struct GpsSample : Sample {
        Vector3f pos_l;
        Vector3f vel_l;    
        Vector3f pos_r;
        Vector3f vel_r;    
    };

    struct MagSample : Sample {
        Vector3f mag;    
    };

    struct DecSample : Sample {
        Vector2f dec;
    };

    struct AirspeedSample : Sample {
        float true_airspeed;
        float eas2tas;
    };

    struct PreIntegralSample : Sample {
        float dt;
        Matrix3f dr;
        Vector3f dv;
        Vector3f dp;
    };
    
}

#endif //NAVIGATION_ESKF_COMMON_H