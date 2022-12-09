#include <Eigen/Dense>
#include <array>

namespace eskf {
    using namespace Eigen;

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