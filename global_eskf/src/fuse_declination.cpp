//
// Created by Cain on 2022/12/9.
//

#include "geskf.h"
#include <cfloat>
#include <iostream>

namespace eskf {
    using namespace std;

    unsigned char GESKF::fuse_declination(const Vector2f &dec, const Vector3f &w, const Vector3f &a, const Vector2f &noise_std, const Vector2f &gate) {
        unsigned char info = 0;                                
        if (!_control_status.flags.dec) {
            return info;
        }        

        // Sequential Kalman Filter
        for (unsigned int dim = 0; dim < 2; ++dim) {
            if (!isfinite(noise_std[dim])) {
                continue;
            }
            /*
            K = P * H * (H * P * H' + R)^-1
            x = x + K * (y - h)
            P = (I - K * H) * P

            where, H = [O, O, O, O, O, O, [sinψ, -cosψ, 0], O, O]
                h = [mx, my]
            */

            // H * P  or  P * H'
            const unsigned int index = 17 + dim;
            array<float, ESKF::DIM> HP = {_cov[0][index],
                                            _cov[1][index],
                                            _cov[2][index],
                                            _cov[3][index],
                                            _cov[4][index],
                                            _cov[5][index],
                                            _cov[6][index],
                                            _cov[7][index],
                                            _cov[8][index],
                                            _cov[9][index],
                                            _cov[10][index],
                                            _cov[11][index],
                                            _cov[12][index],
                                            _cov[13][index],
                                            _cov[14][index],
                                            _cov[15][index],
                                            _cov[16][index],
                                            _cov[17][index],
                                            _cov[index][18],
                                            0.f, 0.f, 0.f, 0.f, 0.f};

            if (!_control_status.flags.mag) {
                HP[16] = 0.f;
            }

            if (_control_status.flags.mag_bias) {
                HP[19] = _cov[index][19];
                HP[20] = _cov[index][20];
                HP[21] = _cov[index][21];
            }

            if (_control_status.flags.wind) {
                HP[22] = _cov[index][22];
                HP[23] = _cov[index][23];
            }

            // H * P * H' + R
            const float HPHT_plus_R = HP[index];

            // e = dec_meas - dec
            const float obs_error = normalize_angle(dec[dim] - _dec[dim]);

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