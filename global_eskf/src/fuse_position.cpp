//
// Created by Cain on 2022/12/9.
//

#include "geskf.h"
#include <cfloat>
#include <iostream>

namespace eskf {
    using namespace std;

    unsigned char GESKF::fuse_position(const Vector3f &pos, const Vector3f &w, const Vector3f &a, 
                                    const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) {
        /*
        s = [p, v, R, bg, ba, g]
        δs = [δp, δv, δθ, δbg, δba, δg]
        s + δs = [p+δp, v+δv, Exp(δθ)*R, bg+δbg, ba+δba, g+δg]

        pos = p + R * dis

        δpos / δp = I
        δpos / δθ = -(R * dis)^

        H = [I, O, -(R*dis)^, O, O, O]
        */

        const Vector3f rot_d = _rot * dis;
        const array<array<float, 3>, 3> minus_rot_d_hat {
            {{0.f, rot_d[2], -rot_d[1]},
            {-rot_d[2], 0.f, rot_d[0]},
            {rot_d[1], -rot_d[0], 0.f}}
        };

        unsigned char info = 0;

        // Sequential Kalman Filter
        for (unsigned int dim = 0; dim < 3; ++dim) {
            if (!isfinite(noise_std[dim])) {
                continue;
            }

            /*
            K = P * H * (H * P * H' + R)^-1
            x = x + K * (y - h)
            P = (I - K * H) * P

            where, H = [I, O, -(R*dis)^, O, O, O]
                h = p + R * dis
            */

            // H * P  or  P * H'
            const float cov_1_dim = (dim == 2) ? _cov[1][dim] : _cov[dim][1];
            array<float, ESKF::DIM> HP = {_cov[0][dim] + _cov[0][6]*minus_rot_d_hat[dim][0] + _cov[0][7]*minus_rot_d_hat[dim][1] + _cov[0][8]*minus_rot_d_hat[dim][2],
                                                cov_1_dim + _cov[1][6]*minus_rot_d_hat[dim][0] + _cov[1][7]*minus_rot_d_hat[dim][1] + _cov[1][8]*minus_rot_d_hat[dim][2],
                                                _cov[dim][2] + _cov[2][6]*minus_rot_d_hat[dim][0] + _cov[2][7]*minus_rot_d_hat[dim][1] + _cov[2][8]*minus_rot_d_hat[dim][2],
                                                _cov[dim][3] + _cov[3][6]*minus_rot_d_hat[dim][0] + _cov[3][7]*minus_rot_d_hat[dim][1] + _cov[3][8]*minus_rot_d_hat[dim][2],
                                                _cov[dim][4] + _cov[4][6]*minus_rot_d_hat[dim][0] + _cov[4][7]*minus_rot_d_hat[dim][1] + _cov[4][8]*minus_rot_d_hat[dim][2],
                                                _cov[dim][5] + _cov[5][6]*minus_rot_d_hat[dim][0] + _cov[5][7]*minus_rot_d_hat[dim][1] + _cov[5][8]*minus_rot_d_hat[dim][2],
                                                _cov[dim][6] + _cov[6][6]*minus_rot_d_hat[dim][0] + _cov[6][7]*minus_rot_d_hat[dim][1] + _cov[6][8]*minus_rot_d_hat[dim][2],
                                                _cov[dim][7] + _cov[6][7]*minus_rot_d_hat[dim][0] + _cov[7][7]*minus_rot_d_hat[dim][1] + _cov[7][8]*minus_rot_d_hat[dim][2],
                                                _cov[dim][8] + _cov[6][8]*minus_rot_d_hat[dim][0] + _cov[7][8]*minus_rot_d_hat[dim][1] + _cov[8][8]*minus_rot_d_hat[dim][2],
                                                _cov[dim][9] + _cov[6][9]*minus_rot_d_hat[dim][0] + _cov[7][9]*minus_rot_d_hat[dim][1] + _cov[8][9]*minus_rot_d_hat[dim][2],
                                                _cov[dim][10] + _cov[6][10]*minus_rot_d_hat[dim][0] + _cov[7][10]*minus_rot_d_hat[dim][1] + _cov[8][10]*minus_rot_d_hat[dim][2],
                                                _cov[dim][11] + _cov[6][11]*minus_rot_d_hat[dim][0] + _cov[7][11]*minus_rot_d_hat[dim][1] + _cov[8][11]*minus_rot_d_hat[dim][2],
                                                _cov[dim][12] + _cov[6][12]*minus_rot_d_hat[dim][0] + _cov[7][12]*minus_rot_d_hat[dim][1] + _cov[8][12]*minus_rot_d_hat[dim][2],
                                                _cov[dim][13] + _cov[6][13]*minus_rot_d_hat[dim][0] + _cov[7][13]*minus_rot_d_hat[dim][1] + _cov[8][13]*minus_rot_d_hat[dim][2],
                                                _cov[dim][14] + _cov[6][14]*minus_rot_d_hat[dim][0] + _cov[7][14]*minus_rot_d_hat[dim][1] + _cov[8][14]*minus_rot_d_hat[dim][2],
                                                _cov[dim][15] + _cov[6][15]*minus_rot_d_hat[dim][0] + _cov[7][15]*minus_rot_d_hat[dim][1] + _cov[8][15]*minus_rot_d_hat[dim][2],
                                                0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};

            if (_control_status.flags.mag) {
                HP[16] = _cov[dim][16] + _cov[6][16]*minus_rot_d_hat[dim][0] + _cov[7][16]*minus_rot_d_hat[dim][1] + _cov[8][16]*minus_rot_d_hat[dim][2];
                HP[17] = _cov[dim][17] + _cov[6][17]*minus_rot_d_hat[dim][0] + _cov[7][17]*minus_rot_d_hat[dim][1] + _cov[8][17]*minus_rot_d_hat[dim][2];
                HP[18] = _cov[dim][18] + _cov[6][18]*minus_rot_d_hat[dim][0] + _cov[7][18]*minus_rot_d_hat[dim][1] + _cov[8][18]*minus_rot_d_hat[dim][2];

                if (_control_status.flags.mag_bias) {
                    HP[19] = _cov[dim][19] + _cov[6][19]*minus_rot_d_hat[dim][0] + _cov[7][19]*minus_rot_d_hat[dim][1] + _cov[8][19]*minus_rot_d_hat[dim][2];
                    HP[20] = _cov[dim][20] + _cov[6][20]*minus_rot_d_hat[dim][0] + _cov[7][20]*minus_rot_d_hat[dim][1] + _cov[8][20]*minus_rot_d_hat[dim][2];
                    HP[21] = _cov[dim][21] + _cov[6][21]*minus_rot_d_hat[dim][0] + _cov[7][21]*minus_rot_d_hat[dim][1] + _cov[8][21]*minus_rot_d_hat[dim][2];                     
                }                                    
            }   

            if (_control_status.flags.wind) {
                HP[22] = _cov[dim][22] + _cov[6][22]*minus_rot_d_hat[dim][0] + _cov[7][22]*minus_rot_d_hat[dim][1] + _cov[8][22]*minus_rot_d_hat[dim][2];
                HP[23] = _cov[dim][23] + _cov[6][23]*minus_rot_d_hat[dim][0] + _cov[7][23]*minus_rot_d_hat[dim][1] + _cov[8][23]*minus_rot_d_hat[dim][2];
            }                                 

            // H * P * H' + R
            const float HPHT_plus_R = HP[dim] + HP[6] * minus_rot_d_hat[dim][0] + HP[7] * minus_rot_d_hat[dim][1] + HP[8] * minus_rot_d_hat[dim][2] + noise_std[dim] * noise_std[dim];

            // h = p + R * dis
            // e = y - h = pos - (p + R * dis)
            const float obs_error = pos[dim] - (_p[dim] + rot_d[dim]);

            /*
            K = P * H' * (H * P * H' + R)^-1
            P = P - K * H * P

            e = y - h = y - (v + R * (w - bg)^ * dis)
            x = x + K * (y - h)
            */
            switch (conservative_posteriori_estimate(HP, HPHT_plus_R, obs_error, gate[dim])) {
                case 0:
                    break;
                case 1:
                    info |= (1 << dim);
                    break;
                case 2:
                    info |= (8 << dim);  
                    reset_covariance_matrix(dim, dim + 1, _q_cov);
                    break;
            }
        }  
        
        regular_covariance_to_symmetric<ESKF::DIM>(0);

        return info;
    }
}