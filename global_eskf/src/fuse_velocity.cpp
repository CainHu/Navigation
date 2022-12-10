//
// Created by Cain on 2022/12/9.
//

#include "geskf.h"
#include <cfloat>
#include <iostream>

namespace eskf {
    using namespace std;

    unsigned char GESKF::fuse_velocity(const Vector3f &vel, const Vector3f &w, const Vector3f &a, 
                                    const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) {
        /*
        s = [p, v, R, bg, ba, g]
        δs = [δp, δv, δθ, δbg, δba, δg]
        s + δs = [p+δp, v+δv, Exp(δθ)*R, bg+δbg, ba+δba, g+δg]

        vel = v + R * (w - bg)^ * dis

        δvel / δv = I
        δvel / δθ = -(R * ((w - bg)^ * dis))^
        δvel / δbg = R * dis^

        H = [O, I, -(R*(w-bg)^*dis)^, R*dis^, O, O]
        */ 

        unsigned char info = 0;

        // R * (w-bg)^ * dis
        const Vector3f rot_w_corr_cross_d = _rot * dis.cross(_bg - w);
        const array<array<float, 3>, 3> rot_d_cross_w_corr_hat {
            {{0.f, rot_w_corr_cross_d[2], -rot_w_corr_cross_d[1]},
            {-rot_w_corr_cross_d[2], 0.f, rot_w_corr_cross_d[0]},
            {rot_w_corr_cross_d[1], -rot_w_corr_cross_d[0], 0.f}}
        };
        
        // Sequential Kalman Filter
        for (unsigned int dim = 0; dim < 3; ++dim) {
            if (!isfinite(noise_std[dim])) {
                continue;
            }

            /*
            K = P * H * (H * P * H' + R)^-1
            x = x + K * (y - h)
            P = (I - K * H) * P

            where, H = [O, I, (R*(dis^*(w-bg)))^, R*dis^, O, O]
                h = v + R * (w-bg)^ * dis
            */

            // R * dis^
            const array<float, 3> rot_d_hat = {_rot(dim, 1)*dis[2] - _rot(dim, 2)*dis[1], 
                                            _rot(dim, 2)*dis[0] - _rot(dim, 0)*dis[2], 
                                            _rot(dim, 0)*dis[1] - _rot(dim, 1)*dis[0]};

            // H * P  or  P * H'
            const unsigned int index = 3 + dim;
            const float cov_4_index = (dim == 2) ? _cov[4][index] : _cov[index][4];
            array<float, ESKF::DIM> HP = {_cov[0][index] + _cov[0][6]*rot_d_cross_w_corr_hat[dim][0] + _cov[0][7]*rot_d_cross_w_corr_hat[dim][1] + _cov[0][8]*rot_d_cross_w_corr_hat[dim][2] + _cov[0][9]*rot_d_hat[0] + _cov[0][10]*rot_d_hat[1] + _cov[0][11]*rot_d_hat[2],
                                            _cov[1][index] + _cov[1][6]*rot_d_cross_w_corr_hat[dim][0] + _cov[1][7]*rot_d_cross_w_corr_hat[dim][1] + _cov[1][8]*rot_d_cross_w_corr_hat[dim][2] + _cov[1][9]*rot_d_hat[0] + _cov[1][10]*rot_d_hat[1] + _cov[1][11]*rot_d_hat[2],
                                            _cov[2][index] + _cov[2][6]*rot_d_cross_w_corr_hat[dim][0] + _cov[2][7]*rot_d_cross_w_corr_hat[dim][1] + _cov[2][8]*rot_d_cross_w_corr_hat[dim][2] + _cov[2][9]*rot_d_hat[0] + _cov[2][10]*rot_d_hat[1] + _cov[2][11]*rot_d_hat[2],
                                            _cov[3][index] + _cov[3][6]*rot_d_cross_w_corr_hat[dim][0] + _cov[3][7]*rot_d_cross_w_corr_hat[dim][1] + _cov[3][8]*rot_d_cross_w_corr_hat[dim][2] + _cov[3][9]*rot_d_hat[0] + _cov[3][10]*rot_d_hat[1] + _cov[3][11]*rot_d_hat[2],
                                            cov_4_index + _cov[4][6]*rot_d_cross_w_corr_hat[dim][0] + _cov[4][7]*rot_d_cross_w_corr_hat[dim][1] + _cov[4][8]*rot_d_cross_w_corr_hat[dim][2] + _cov[4][9]*rot_d_hat[0] + _cov[4][10]*rot_d_hat[1] + _cov[4][11]*rot_d_hat[2],
                                            _cov[index][5] + _cov[5][6]*rot_d_cross_w_corr_hat[dim][0] + _cov[5][7]*rot_d_cross_w_corr_hat[dim][1] + _cov[5][8]*rot_d_cross_w_corr_hat[dim][2] + _cov[5][9]*rot_d_hat[0] + _cov[5][10]*rot_d_hat[1] + _cov[5][11]*rot_d_hat[2],
                                            _cov[index][6] + _cov[6][6]*rot_d_cross_w_corr_hat[dim][0] + _cov[6][7]*rot_d_cross_w_corr_hat[dim][1] + _cov[6][8]*rot_d_cross_w_corr_hat[dim][2] + _cov[6][9]*rot_d_hat[0] + _cov[6][10]*rot_d_hat[1] + _cov[6][11]*rot_d_hat[2],
                                            _cov[index][7] + _cov[6][7]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][7]*rot_d_cross_w_corr_hat[dim][1] + _cov[7][8]*rot_d_cross_w_corr_hat[dim][2] + _cov[7][9]*rot_d_hat[0] + _cov[7][10]*rot_d_hat[1] + _cov[7][11]*rot_d_hat[2],
                                            _cov[index][8] + _cov[6][8]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][8]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][8]*rot_d_cross_w_corr_hat[dim][2] + _cov[8][9]*rot_d_hat[0] + _cov[8][10]*rot_d_hat[1] + _cov[8][11]*rot_d_hat[2],
                                            _cov[index][9] + _cov[6][9]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][9]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][9]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][9]*rot_d_hat[0] + _cov[9][10]*rot_d_hat[1] + _cov[9][11]*rot_d_hat[2],
                                            _cov[index][10] + _cov[6][10]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][10]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][10]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][10]*rot_d_hat[0] + _cov[10][10]*rot_d_hat[1] + _cov[10][11]*rot_d_hat[2],
                                            _cov[index][11] + _cov[6][11]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][11]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][11]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][11]*rot_d_hat[0] + _cov[10][11]*rot_d_hat[1] + _cov[11][11]*rot_d_hat[2],
                                            _cov[index][12] + _cov[6][12]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][12]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][12]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][12]*rot_d_hat[0] + _cov[10][12]*rot_d_hat[1] + _cov[11][12]*rot_d_hat[2],
                                            _cov[index][13] + _cov[6][13]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][13]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][13]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][13]*rot_d_hat[0] + _cov[10][13]*rot_d_hat[1] + _cov[11][13]*rot_d_hat[2],
                                            _cov[index][14] + _cov[6][14]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][14]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][14]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][14]*rot_d_hat[0] + _cov[10][14]*rot_d_hat[1] + _cov[11][14]*rot_d_hat[2],
                                            _cov[index][15] + _cov[6][15]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][15]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][15]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][15]*rot_d_hat[0] + _cov[10][15]*rot_d_hat[1] + _cov[11][15]*rot_d_hat[2],
                                            0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};

            if (_control_status.flags.mag) {
                HP[16] = _cov[index][16] + _cov[6][16]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][16]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][16]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][16]*rot_d_hat[0] + _cov[10][16]*rot_d_hat[1] + _cov[11][16]*rot_d_hat[2];
                HP[17] = _cov[index][17] + _cov[6][17]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][17]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][17]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][17]*rot_d_hat[0] + _cov[10][17]*rot_d_hat[1] + _cov[11][17]*rot_d_hat[2];
                HP[18] = _cov[index][18] + _cov[6][18]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][18]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][18]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][18]*rot_d_hat[0] + _cov[10][18]*rot_d_hat[1] + _cov[11][18]*rot_d_hat[2];                         
                if (_control_status.flags.mag_bias) {
                    HP[19] = _cov[index][19] + _cov[6][19]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][19]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][19]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][19]*rot_d_hat[0] + _cov[10][19]*rot_d_hat[1] + _cov[11][19]*rot_d_hat[2];
                    HP[20] = _cov[index][20] + _cov[6][20]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][20]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][20]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][20]*rot_d_hat[0] + _cov[10][20]*rot_d_hat[1] + _cov[11][20]*rot_d_hat[2];
                    HP[21] = _cov[index][21] + _cov[6][21]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][21]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][21]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][21]*rot_d_hat[0] + _cov[10][21]*rot_d_hat[1] + _cov[11][21]*rot_d_hat[2];
                                            
                }
            }

            if (_control_status.flags.wind) {
                HP[22] = _cov[index][22] + _cov[6][22]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][22]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][22]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][22]*rot_d_hat[0] + _cov[10][22]*rot_d_hat[1] + _cov[11][22]*rot_d_hat[2];
                HP[23] = _cov[index][23] + _cov[6][23]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][23]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][23]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][23]*rot_d_hat[0] + _cov[10][23]*rot_d_hat[1] + _cov[11][23]*rot_d_hat[2];
            }

            // H * P * H' + R
            const float HPHT_plus_R = HP[index] + HP[6] * rot_d_cross_w_corr_hat[dim][0] + HP[7] * rot_d_cross_w_corr_hat[dim][1] + HP[8] * rot_d_cross_w_corr_hat[dim][2] + HP[9] * rot_d_hat[0] + HP[10] * rot_d_hat[1] + HP[11] * rot_d_hat[2] + noise_std[dim] * noise_std[dim];

            // h = v + R * (w - bg)^ * dis
            // e = y - h = y - (v + R * (w - bg)^ * dis)
            const float obs_error = vel[dim] - (_v[dim] + rot_w_corr_cross_d[dim]);

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
                    reset_covariance_matrix(index, index + 1, _q_cov);
                    break;
            }
        }
        
        regular_covariance_to_symmetric<ESKF::DIM>(0);

        return info;
    }
}