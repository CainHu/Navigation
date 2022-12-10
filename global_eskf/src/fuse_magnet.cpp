//
// Created by Cain on 2022/12/9.
//

#include "geskf.h"
#include <cfloat>
#include <iostream>

namespace eskf {
    using namespace std;

    unsigned char GESKF::fuse_magnet(const Vector3f &mag, const Vector3f &w, const Vector3f &a, 
                                    const Vector3f &noise_std, const Vector3f &gate) {
        unsigned char info = 0;                                
        if (!(_control_status.flags.mag || _control_status.flags.dec || _control_status.flags.mag_bias)) {
            return info;
        }        

        // Hb = y - bm
        const Vector3f mag_corr = mag - _bm;

        const float cos_y = cosf(_dec[0]), sin_y = sinf(_dec[0]);
        const float cos_z = cosf(_dec[1]), sin_z = sinf(_dec[1]);

        // Rz * Ry * ex
        const array<float, 3> rz_ry_ex {
            cos_z * cos_y, sin_z * cos_y, -sin_y
        };

        // 1 / H
        const float h_inv = (_h < 1e-6f) ? 1.f / _h_init : 1.f / _h;

        // Sequential Kalman Filter
        for (unsigned int dim = 0; dim < 3; ++dim) {
            if (!isfinite(noise_std[dim])) {
                continue;
            }

            /*
            K = P * H * (H * P * H' + R)^-1
            x = x + K * (y - h)
            P = (I - K * H) * P

            where, H = [O, O, R'*m^, O, O, O, R', I, O]
                h = R' * m
            */

            // R' * RZ * RY * ex
            const float rt_rz_ry_ex = _rot(0, dim) * rz_ry_ex[0] + _rot(1, dim) * rz_ry_ex[1] + _rot(2, dim) * rz_ry_ex[2];

            // R' * RZ * RY * ex / H
            const float rt_rz_ry_ex_h_inv = rt_rz_ry_ex * h_inv;

            // R' * (RZ * RY * ex)^
            const array<float, 3> rt_rz_ry_ex_hat {
                _rot(1, dim) * rz_ry_ex[2] - _rot(2, dim) * rz_ry_ex[1], 
                _rot(2, dim) * rz_ry_ex[0] - _rot(0, dim) * rz_ry_ex[2], 
                _rot(0, dim) * rz_ry_ex[1] - _rot(1, dim) * rz_ry_ex[0]
            };

            const float param_y = rt_rz_ry_ex_hat[1]*cos_z - sin_z*rt_rz_ry_ex_hat[0];

            // H * P  or  P * H'
            const unsigned int index = 19 + dim;
            array<float, ESKF::DIM> HP = {_cov[0][index] * h_inv + _cov[0][6]*rt_rz_ry_ex_hat[0] + _cov[0][7]*rt_rz_ry_ex_hat[1] + _cov[0][8]*rt_rz_ry_ex_hat[2] + _cov[0][16]*rt_rz_ry_ex_h_inv - _cov[0][17]*param_y - _cov[0][18]*rt_rz_ry_ex_hat[2],
                                            _cov[1][index] * h_inv + _cov[1][6]*rt_rz_ry_ex_hat[0] + _cov[1][7]*rt_rz_ry_ex_hat[1] + _cov[1][8]*rt_rz_ry_ex_hat[2] + _cov[1][16]*rt_rz_ry_ex_h_inv - _cov[1][17]*param_y - _cov[1][18]*rt_rz_ry_ex_hat[2],
                                            _cov[2][index] * h_inv + _cov[2][6]*rt_rz_ry_ex_hat[0] + _cov[2][7]*rt_rz_ry_ex_hat[1] + _cov[2][8]*rt_rz_ry_ex_hat[2] + _cov[2][16]*rt_rz_ry_ex_h_inv - _cov[2][17]*param_y - _cov[2][18]*rt_rz_ry_ex_hat[2],
                                            _cov[3][index] * h_inv + _cov[3][6]*rt_rz_ry_ex_hat[0] + _cov[3][7]*rt_rz_ry_ex_hat[1] + _cov[3][8]*rt_rz_ry_ex_hat[2] + _cov[3][16]*rt_rz_ry_ex_h_inv - _cov[3][17]*param_y - _cov[3][18]*rt_rz_ry_ex_hat[2],
                                            _cov[4][index] * h_inv + _cov[4][6]*rt_rz_ry_ex_hat[0] + _cov[4][7]*rt_rz_ry_ex_hat[1] + _cov[4][8]*rt_rz_ry_ex_hat[2] + _cov[4][16]*rt_rz_ry_ex_h_inv - _cov[4][17]*param_y - _cov[4][18]*rt_rz_ry_ex_hat[2],
                                            _cov[5][index] * h_inv + _cov[5][6]*rt_rz_ry_ex_hat[0] + _cov[5][7]*rt_rz_ry_ex_hat[1] + _cov[5][8]*rt_rz_ry_ex_hat[2] + _cov[5][16]*rt_rz_ry_ex_h_inv - _cov[5][17]*param_y - _cov[5][18]*rt_rz_ry_ex_hat[2],
                                            _cov[6][index] * h_inv + _cov[6][6]*rt_rz_ry_ex_hat[0] + _cov[6][7]*rt_rz_ry_ex_hat[1] + _cov[6][8]*rt_rz_ry_ex_hat[2] + _cov[6][16]*rt_rz_ry_ex_h_inv - _cov[6][17]*param_y - _cov[6][18]*rt_rz_ry_ex_hat[2],
                                            _cov[7][index] * h_inv + _cov[6][7]*rt_rz_ry_ex_hat[0] + _cov[7][7]*rt_rz_ry_ex_hat[1] + _cov[7][8]*rt_rz_ry_ex_hat[2] + _cov[7][16]*rt_rz_ry_ex_h_inv - _cov[7][17]*param_y - _cov[7][18]*rt_rz_ry_ex_hat[2],
                                            _cov[8][index] * h_inv + _cov[6][8]*rt_rz_ry_ex_hat[0] + _cov[7][8]*rt_rz_ry_ex_hat[1] + _cov[8][8]*rt_rz_ry_ex_hat[2] + _cov[8][16]*rt_rz_ry_ex_h_inv - _cov[8][17]*param_y - _cov[8][18]*rt_rz_ry_ex_hat[2],
                                            _cov[9][index] * h_inv + _cov[6][9]*rt_rz_ry_ex_hat[0] + _cov[7][9]*rt_rz_ry_ex_hat[1] + _cov[8][9]*rt_rz_ry_ex_hat[2] + _cov[9][16]*rt_rz_ry_ex_h_inv - _cov[9][17]*param_y - _cov[9][18]*rt_rz_ry_ex_hat[2],
                                            _cov[10][index] * h_inv + _cov[6][10]*rt_rz_ry_ex_hat[0] + _cov[7][10]*rt_rz_ry_ex_hat[1] + _cov[8][10]*rt_rz_ry_ex_hat[2] + _cov[10][16]*rt_rz_ry_ex_h_inv - _cov[10][17]*param_y - _cov[10][18]*rt_rz_ry_ex_hat[2],
                                            _cov[11][index] * h_inv + _cov[6][11]*rt_rz_ry_ex_hat[0] + _cov[7][11]*rt_rz_ry_ex_hat[1] + _cov[8][11]*rt_rz_ry_ex_hat[2] + _cov[11][16]*rt_rz_ry_ex_h_inv - _cov[11][17]*param_y - _cov[11][18]*rt_rz_ry_ex_hat[2],
                                            _cov[12][index] * h_inv + _cov[6][12]*rt_rz_ry_ex_hat[0] + _cov[7][12]*rt_rz_ry_ex_hat[1] + _cov[8][12]*rt_rz_ry_ex_hat[2] + _cov[12][16]*rt_rz_ry_ex_h_inv - _cov[12][17]*param_y - _cov[12][18]*rt_rz_ry_ex_hat[2],
                                            _cov[13][index] * h_inv + _cov[6][13]*rt_rz_ry_ex_hat[0] + _cov[7][13]*rt_rz_ry_ex_hat[1] + _cov[8][13]*rt_rz_ry_ex_hat[2] + _cov[13][16]*rt_rz_ry_ex_h_inv - _cov[13][17]*param_y - _cov[13][18]*rt_rz_ry_ex_hat[2],
                                            _cov[14][index] * h_inv + _cov[6][14]*rt_rz_ry_ex_hat[0] + _cov[7][14]*rt_rz_ry_ex_hat[1] + _cov[8][14]*rt_rz_ry_ex_hat[2] + _cov[14][16]*rt_rz_ry_ex_h_inv - _cov[14][17]*param_y - _cov[14][18]*rt_rz_ry_ex_hat[2],
                                            _cov[15][index] * h_inv + _cov[6][15]*rt_rz_ry_ex_hat[0] + _cov[7][15]*rt_rz_ry_ex_hat[1] + _cov[8][15]*rt_rz_ry_ex_hat[2] + _cov[15][16]*rt_rz_ry_ex_h_inv - _cov[15][17]*param_y - _cov[15][18]*rt_rz_ry_ex_hat[2],
                                            0.f, 0.f, 0.f,
                                            0.f, 0.f, 0.f, 0.f, 0.f};

            if (_control_status.flags.mag) {
                HP[16] = _cov[16][index] * h_inv + _cov[6][16]*rt_rz_ry_ex_hat[0] + _cov[7][16]*rt_rz_ry_ex_hat[1] + _cov[8][16]*rt_rz_ry_ex_hat[2] + _cov[16][16]*rt_rz_ry_ex_h_inv - _cov[16][17]*param_y - _cov[16][18]*rt_rz_ry_ex_hat[2];
            }      

            if (_control_status.flags.dec) {
                HP[17] = _cov[17][index] * h_inv + _cov[6][17]*rt_rz_ry_ex_hat[0] + _cov[7][17]*rt_rz_ry_ex_hat[1] + _cov[8][17]*rt_rz_ry_ex_hat[2] + _cov[16][17]*rt_rz_ry_ex_h_inv - _cov[17][17]*param_y - _cov[17][18]*rt_rz_ry_ex_hat[2];
                HP[18] = _cov[18][index] * h_inv + _cov[6][18]*rt_rz_ry_ex_hat[0] + _cov[7][18]*rt_rz_ry_ex_hat[1] + _cov[8][18]*rt_rz_ry_ex_hat[2] + _cov[16][18]*rt_rz_ry_ex_h_inv - _cov[17][18]*param_y - _cov[18][18]*rt_rz_ry_ex_hat[2];
            }                          

            if (_control_status.flags.mag_bias) {
                const float cov_20_index = (dim == 2) ? _cov[20][index] : _cov[index][20];
                HP[19] = _cov[19][index] * h_inv + _cov[6][19]*rt_rz_ry_ex_hat[0] + _cov[7][19]*rt_rz_ry_ex_hat[1] + _cov[8][19]*rt_rz_ry_ex_hat[2] + _cov[16][19]*rt_rz_ry_ex_h_inv - _cov[17][19]*param_y - _cov[18][19]*rt_rz_ry_ex_hat[2];
                HP[20] = cov_20_index * h_inv + _cov[6][20]*rt_rz_ry_ex_hat[0] + _cov[7][20]*rt_rz_ry_ex_hat[1] + _cov[8][20]*rt_rz_ry_ex_hat[2] + _cov[16][20]*rt_rz_ry_ex_h_inv - _cov[17][20]*param_y - _cov[18][20]*rt_rz_ry_ex_hat[2];
                HP[21] = _cov[index][21] * h_inv + _cov[6][21]*rt_rz_ry_ex_hat[0] + _cov[7][21]*rt_rz_ry_ex_hat[1] + _cov[8][21]*rt_rz_ry_ex_hat[2] + _cov[16][21]*rt_rz_ry_ex_h_inv - _cov[17][21]*param_y - _cov[18][21]*rt_rz_ry_ex_hat[2];                           
            }

            if (_control_status.flags.wind) {
                HP[22] = _cov[index][22] * h_inv + _cov[6][22]*rt_rz_ry_ex_hat[0] + _cov[7][22]*rt_rz_ry_ex_hat[1] + _cov[8][22]*rt_rz_ry_ex_hat[2] + _cov[16][22]*rt_rz_ry_ex_h_inv - _cov[17][22]*param_y - _cov[18][22]*rt_rz_ry_ex_hat[2];
                HP[23] = _cov[index][23] * h_inv + _cov[6][23]*rt_rz_ry_ex_hat[0] + _cov[7][23]*rt_rz_ry_ex_hat[1] + _cov[8][23]*rt_rz_ry_ex_hat[2] + _cov[16][23]*rt_rz_ry_ex_h_inv - _cov[17][23]*param_y - _cov[18][23]*rt_rz_ry_ex_hat[2];
            }

            // H * P * H' + R
            const float HPHT_plus_R = HP[index] * h_inv + HP[6] * rt_rz_ry_ex_hat[0] + HP[7] * rt_rz_ry_ex_hat[1] + HP[8] * rt_rz_ry_ex_hat[2] + HP[16] * rt_rz_ry_ex_h_inv - HP[17] * param_y - HP[18] * rt_rz_ry_ex_hat[2] + noise_std[dim] * noise_std[dim];

            // h = m
            // e = (y - bm) - R' * m
            const float obs_error = mag_corr[dim] * h_inv - rt_rz_ry_ex;

            /*
            K = P * H' * (H * P * H' + R)^-1
            P = P - K * H * P

            e = (y - bm) - R'*m
            x = x + K * e
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