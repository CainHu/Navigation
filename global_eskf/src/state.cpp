//
// Created by Cain on 2022/12/9.
//

#include "geskf.h"
#include <cfloat>
#include <iostream>

namespace eskf {
    using namespace std;

    void GESKF::correct_state() {
        // state: [p, v, bg, ba, g], R
        // error_state : [δp, δv, δθ, δbg, δba, δg]

        /*
        p = p + δp
        v = v + δv
        bg = bg + δbg
        ba = ba + δba
        g = g + δg
        m = m + δm
        bm = bm + δbm
        */
        for (unsigned char i = 0; i < 3; ++i) {
            _p[i] += _error_state[i];
            _v[i] += _error_state[3 + i];
            _bg[i] += _error_state[9 + i];
            _ba[i] += _error_state[12 + i];
            _bm[i] += _error_state[19 + i];
        }
        _g += _error_state[15];
        _h += _error_state[16];
        _dec[0] += _error_state[17];
        _dec[1] += _error_state[18];
        _w[0] += _error_state[22];
        _w[1] += _error_state[23];
        
        // q = Exp(δθ) * q
        Quaternionf delta_q;
        array<float, 3> delta_theta = {_error_state[6], _error_state[7], _error_state[8]};
        quaternion_from_axis_angle(delta_q, delta_theta);
        const Quaternionf q = _q;
        _q = delta_q * q;
        _rot = q;

        // [δp, δv, δθ, δbg, δba, δg] = 0
        for (float &es : _error_state) {
            es = 0.f;
        }
    }
}