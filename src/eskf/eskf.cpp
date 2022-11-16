//
// Created by Cain on 2022/11/11.
//

#include "eskf.h"
#include <cfloat>
#include <iostream>

using namespace std;
using namespace eskf;

void ESKF::initialize() {
    _rot.setIdentity();

    reset_state();
    reset_error_state();
    reset_covariance_matrix(0, 16, _q_cov);
    regular_covariance_to_symmetric();
    reset_accmulator_cov();
}

void ESKF::rotation_from_axis_angle(Matrix3f &r, const Vector3f &a) const {
    /* Rodrigues's Formula:
    * a = n * θ
    * R = cosθ*I + (1 - cosθ)*n*n' + sinθ*n^
    * */
    const float a_norm_square = a.squaredNorm();
    if (a_norm_square < _dt2 * 1e-12f) {
        r.setIdentity();
    } else {
        const float a_norm = sqrtf(a_norm_square);
        const Vector3f a_unit = a / a_norm;
        const float theta = a_norm;
        const float cos_theta = cosf(theta), sin_theta = sinf(theta);
        const float tmp = 1.f - cos_theta;

        const float xx = a_unit(0) * a_unit(0) * tmp;
        const float xy = a_unit(0) * a_unit(1) * tmp;
        const float xz = a_unit(0) * a_unit(2) * tmp;
        const float yy = a_unit(1) * a_unit(1) * tmp;
        const float yz = a_unit(1) * a_unit(2) * tmp;
        const float zz = a_unit(2) * a_unit(2) * tmp;

        const float sx = sin_theta * a_unit(0);
        const float sy = sin_theta * a_unit(1);
        const float sz = sin_theta * a_unit(2);

        r(0, 0) = cos_theta + xx;
        r(0, 1) = xy - sz;
        r(0, 2) = xz + sy;
        r(1, 0) = xy + sz;
        r(1, 1) = cos_theta + yy;
        r(1, 2) = yz - sx;
        r(2, 0) = xz - sy;
        r(2, 1) = yz + sx;
        r(2, 2) = cos_theta + zz;
    }
}

void ESKF::quaternion_from_axis_angle(Quaternionf &q, const Vector3f &a) const {
    /*
    a = n * θ
    q = [cos(θ/2), n*sin(θ/2)]
    */

    const float a_norm_square = a.squaredNorm();
    if (a_norm_square < _dt2 * 1e-12f) {
        q.setIdentity();
    } else {
        const float a_norm = sqrtf(a_norm_square);
        const Vector3f a_unit = a / a_norm;
        const float half_theta = 0.5f * a_norm;
        const float c = cosf(half_theta), s = sinf(half_theta);

        q.w() = c;
        q.x() = s * a_unit(0);
        q.y() = s * a_unit(1);
        q.z() = s * a_unit(2);
    }
}

void ESKF::predict_state(const Vector3f &w, const Vector3f &a) {
    // state: [p, v, bg, ba, g]
    // error_state : [δp, δv, δθ, δbg, δba, δg]

    /* Nominal Model:
    * p = p + v * Δt + 0.5 * (R * (a - ba) + g) * Δt^2
    * v = v + (R * (a - ba) + g) * Δt
    * R = R * Exp((w - bg) * Δt)
    * bg = bg
    * ba = ba
    * g = g
    * */

    // a_body_corr = a - ba
    const Vector3f a_body_corr = a - _ba;

    // a_world = R * (a - ba) + [0;0;g]
    const Vector3f a_world = _rot * a_body_corr + Vector3f(0.f, 0.f, _g);

    // v
    const Vector3f v_last = _v;

    // v' = v + a * Δt
    _v += a_world * _dt;

    // p' = p + 0.5 * (v' + v) * Δt
    _p += 0.5f * (_v + v_last) * _dt;

    // bg' = bg
    // ba' = ba
    // g = g

    // axis_angle = (w - bg) * Δt
    const Vector3f axis_angle = w - _bg;

    // Δq = Exp((w - bg) * Δt)
    Quaternionf delta_q {};
    quaternion_from_axis_angle(delta_q, axis_angle);

    // q' = q * Δq
    const Quaternionf q = _q;
    _q = q * delta_q;
    _q.normalize();

    _rot = _q;
}

void ESKF::predict_covariance(const Vector3f &w, const Vector3f &a) {
    // -R * dt
    const float mrdt00 = -_rot(0, 0) * _dt, mrdt01 = -_rot(0, 1) * _dt, mrdt02 = -_rot(0, 2) * _dt;
    const float mrdt10 = -_rot(1, 0) * _dt, mrdt11 = -_rot(1, 1) * _dt, mrdt12 = -_rot(1, 2) * _dt;
    const float mrdt20 = -_rot(2, 0) * _dt, mrdt21 = -_rot(2, 1) * _dt, mrdt22 = -_rot(2, 2) * _dt;

    // a - ba
    const Vector3f a_corr = a - _ba;
    // X = -R * (a - ba)^ * dt
    const float x00 = mrdt01*a_corr(2) - mrdt02*a_corr(1);
    const float x01 = mrdt02*a_corr(0) - mrdt00*a_corr(2);
    const float x02 = mrdt00*a_corr(1) - mrdt01*a_corr(0);
    const float x10 = mrdt11*a_corr(2) - mrdt12*a_corr(1);
    const float x11 = mrdt12*a_corr(0) - mrdt10*a_corr(2);
    const float x12 = mrdt10*a_corr(1) - mrdt11*a_corr(0);
    const float x20 = mrdt21*a_corr(2) - mrdt22*a_corr(1);
    const float x21 = mrdt22*a_corr(0) - mrdt20*a_corr(2);
    const float x22 = mrdt20*a_corr(1) - mrdt21*a_corr(0);

    // (w - bg) * dt
    const Vector3f axis_angle_corr = (w - _bg) * _dt;
    // Y = Exp(-(w - bg) * Δt)
    const float y00 = 1.f, y01 = axis_angle_corr(2), y02 = -axis_angle_corr(1);
    const float y10 = -axis_angle_corr(2), y11 = 1.f, y12 = axis_angle_corr(0);
    const float y20 = axis_angle_corr(1), y21 = -axis_angle_corr(0), y22 = 1.f;
    // Matrix3f y;
    // rotation_from_axis_angle(y, axis_angle_corr);
    // const float y00 = y(0, 1), y01 = y(0, 1), y02 = y(0, 2);
    // const float y10 = y(1, 1), y11 = y(1, 1), y12 = y(1, 2);
    // const float y20 = y(2, 1), y21 = y(2, 1), y22 = y(2, 2);

    // Equations for covariance matrix prediction, without process noise!
    const float var0 = _cov(3,0) + _cov(3,3)*_dt;
    const float var1 = _cov(4,3)*_dt;
    const float var2 = _cov(4,0) + var1;
    const float var3 = _cov(5,3)*_dt;
    const float var4 = _cov(5,0) + var3;
    const float var5 = _cov(12,0) + _cov(12,3)*_dt;
    const float var6 = _cov(13,0) + _cov(13,3)*_dt;
    const float var7 = _cov(14,0) + _cov(14,3)*_dt;
    const float var8 = _cov(6,0) + _cov(6,3)*_dt;
    const float var9 = _cov(7,0) + _cov(7,3)*_dt;
    const float var10 = _cov(8,0) + _cov(8,3)*_dt;
    const float var11 = _cov(15,0) + _cov(15,3)*_dt;
    const float var12 = _cov(9,0) + _cov(9,3)*_dt;
    const float var13 = _cov(10,0) + _cov(10,3)*_dt;
    const float var14 = _cov(11,0) + _cov(11,3)*_dt;
    const float var15 = _cov(4,1) + _cov(4,4)*_dt;
    const float var16 = _cov(5,4)*_dt;
    const float var17 = _cov(5,1) + var16;
    const float var18 = _cov(12,1) + _cov(12,4)*_dt;
    const float var19 = _cov(13,1) + _cov(13,4)*_dt;
    const float var20 = _cov(14,1) + _cov(14,4)*_dt;
    const float var21 = _cov(6,1) + _cov(6,4)*_dt;
    const float var22 = _cov(7,1) + _cov(7,4)*_dt;
    const float var23 = _cov(8,1) + _cov(8,4)*_dt;
    const float var24 = _cov(15,1) + _cov(15,4)*_dt;
    const float var25 = _cov(9,1) + _cov(9,4)*_dt;
    const float var26 = _cov(10,1) + _cov(10,4)*_dt;
    const float var27 = _cov(11,1) + _cov(11,4)*_dt;
    const float var28 = _cov(5,2) + _cov(5,5)*_dt;
    const float var29 = _cov(12,2) + _cov(12,5)*_dt;
    const float var30 = _cov(13,2) + _cov(13,5)*_dt;
    const float var31 = _cov(14,2) + _cov(14,5)*_dt;
    const float var32 = _cov(6,2) + _cov(6,5)*_dt;
    const float var33 = _cov(7,2) + _cov(7,5)*_dt;
    const float var34 = _cov(8,2) + _cov(8,5)*_dt;
    const float var35 = _cov(15,5)*_dt;
    const float var36 = _cov(15,2) + var35;
    const float var37 = _cov(9,2) + _cov(9,5)*_dt;
    const float var38 = _cov(10,2) + _cov(10,5)*_dt;
    const float var39 = _cov(11,2) + _cov(11,5)*_dt;
    const float var40 = _cov(12,12)*mrdt00 + _cov(13,12)*mrdt01 + _cov(14,12)*mrdt02 + _cov(12,3) + _cov(12,6)*x00 + _cov(12,7)*x01 + _cov(12,8)*x02;
    const float var41 = _cov(13,12)*mrdt00 + _cov(13,13)*mrdt01 + _cov(14,13)*mrdt02 + _cov(13,3) + _cov(13,6)*x00 + _cov(13,7)*x01 + _cov(13,8)*x02;
    const float var42 = _cov(14,12)*mrdt00 + _cov(14,13)*mrdt01 + _cov(14,14)*mrdt02 + _cov(14,3) + _cov(14,6)*x00 + _cov(14,7)*x01 + _cov(14,8)*x02;
    const float var43 = _cov(6,3) + _cov(12,6)*mrdt00 + _cov(13,6)*mrdt01 + _cov(14,6)*mrdt02 + _cov(6,6)*x00 + _cov(7,6)*x01 + _cov(8,6)*x02;
    const float var44 = _cov(7,3) + _cov(7,6)*x00 + _cov(12,7)*mrdt00 + _cov(13,7)*mrdt01 + _cov(14,7)*mrdt02 + _cov(7,7)*x01 + _cov(8,7)*x02;
    const float var45 = _cov(8,3) + _cov(8,6)*x00 + _cov(8,7)*x01 + _cov(12,8)*mrdt00 + _cov(13,8)*mrdt01 + _cov(14,8)*mrdt02 + _cov(8,8)*x02;
    const float var46 = _cov(15,12)*mrdt00 + _cov(15,13)*mrdt01 + _cov(15,14)*mrdt02 + _cov(15,3) + _cov(15,6)*x00 + _cov(15,7)*x01 + _cov(15,8)*x02;
    const float var47 = _cov(9,3) + _cov(9,6)*x00 + _cov(9,7)*x01 + _cov(9,8)*x02 + _cov(12,9)*mrdt00 + _cov(13,9)*mrdt01 + _cov(14,9)*mrdt02;
    const float var48 = _cov(12,10)*mrdt00 + _cov(13,10)*mrdt01 + _cov(14,10)*mrdt02 + _cov(10,3) + _cov(10,6)*x00 + _cov(10,7)*x01 + _cov(10,8)*x02;
    const float var49 = _cov(12,11)*mrdt00 + _cov(13,11)*mrdt01 + _cov(14,11)*mrdt02 + _cov(11,3) + _cov(11,6)*x00 + _cov(11,7)*x01 + _cov(11,8)*x02;
    const float var50 = _cov(12,12)*mrdt10 + _cov(13,12)*mrdt11 + _cov(14,12)*mrdt12 + _cov(12,4) + _cov(12,6)*x10 + _cov(12,7)*x11 + _cov(12,8)*x12;
    const float var51 = _cov(13,12)*mrdt10 + _cov(13,13)*mrdt11 + _cov(14,13)*mrdt12 + _cov(13,4) + _cov(13,6)*x10 + _cov(13,7)*x11 + _cov(13,8)*x12;
    const float var52 = _cov(14,12)*mrdt10 + _cov(14,13)*mrdt11 + _cov(14,14)*mrdt12 + _cov(14,4) + _cov(14,6)*x10 + _cov(14,7)*x11 + _cov(14,8)*x12;
    const float var53 = _cov(6,4) + _cov(12,6)*mrdt10 + _cov(13,6)*mrdt11 + _cov(14,6)*mrdt12 + _cov(6,6)*x10 + _cov(7,6)*x11 + _cov(8,6)*x12;
    const float var54 = _cov(7,4) + _cov(7,6)*x10 + _cov(12,7)*mrdt10 + _cov(13,7)*mrdt11 + _cov(14,7)*mrdt12 + _cov(7,7)*x11 + _cov(8,7)*x12;
    const float var55 = _cov(8,4) + _cov(8,6)*x10 + _cov(8,7)*x11 + _cov(12,8)*mrdt10 + _cov(13,8)*mrdt11 + _cov(14,8)*mrdt12 + _cov(8,8)*x12;
    const float var56 = _cov(15,12)*mrdt10 + _cov(15,13)*mrdt11 + _cov(15,14)*mrdt12 + _cov(15,4) + _cov(15,6)*x10 + _cov(15,7)*x11 + _cov(15,8)*x12;
    const float var57 = _cov(9,4) + _cov(9,6)*x10 + _cov(9,7)*x11 + _cov(9,8)*x12 + _cov(12,9)*mrdt10 + _cov(13,9)*mrdt11 + _cov(14,9)*mrdt12;
    const float var58 = _cov(12,10)*mrdt10 + _cov(13,10)*mrdt11 + _cov(14,10)*mrdt12 + _cov(10,4) + _cov(10,6)*x10 + _cov(10,7)*x11 + _cov(10,8)*x12;
    const float var59 = _cov(12,11)*mrdt10 + _cov(13,11)*mrdt11 + _cov(14,11)*mrdt12 + _cov(11,4) + _cov(11,6)*x10 + _cov(11,7)*x11 + _cov(11,8)*x12;
    const float var60 = _cov(15,12)*mrdt20 + _cov(15,13)*mrdt21 + _cov(15,14)*mrdt22 + _cov(15,15)*_dt + _cov(15,5) + _cov(15,6)*x20 + _cov(15,7)*x21 + _cov(15,8)*x22;
    const float var61 = _cov(12,12)*mrdt20 + _cov(13,12)*mrdt21 + _cov(14,12)*mrdt22 + _cov(15,12)*_dt + _cov(12,5) + _cov(12,6)*x20 + _cov(12,7)*x21 + _cov(12,8)*x22;
    const float var62 = _cov(13,12)*mrdt20 + _cov(13,13)*mrdt21 + _cov(14,13)*mrdt22 + _cov(15,13)*_dt + _cov(13,5) + _cov(13,6)*x20 + _cov(13,7)*x21 + _cov(13,8)*x22;
    const float var63 = _cov(14,12)*mrdt20 + _cov(14,13)*mrdt21 + _cov(14,14)*mrdt22 + _cov(15,14)*_dt + _cov(14,5) + _cov(14,6)*x20 + _cov(14,7)*x21 + _cov(14,8)*x22;
    const float var64 = _cov(6,5) + _cov(12,6)*mrdt20 + _cov(13,6)*mrdt21 + _cov(14,6)*mrdt22 + _cov(15,6)*_dt + _cov(6,6)*x20 + _cov(7,6)*x21 + _cov(8,6)*x22;
    const float var65 = _cov(7,5) + _cov(7,6)*x20 + _cov(12,7)*mrdt20 + _cov(13,7)*mrdt21 + _cov(14,7)*mrdt22 + _cov(15,7)*_dt + _cov(7,7)*x21 + _cov(8,7)*x22;
    const float var66 = _cov(8,5) + _cov(8,6)*x20 + _cov(8,7)*x21 + _cov(12,8)*mrdt20 + _cov(13,8)*mrdt21 + _cov(14,8)*mrdt22 + _cov(15,8)*_dt + _cov(8,8)*x22;
    const float var67 = _cov(15,9)*_dt;
    const float var68 = _cov(9,5) + _cov(9,6)*x20 + _cov(9,7)*x21 + _cov(9,8)*x22 + _cov(12,9)*mrdt20 + _cov(13,9)*mrdt21 + _cov(14,9)*mrdt22 + var67;
    const float var69 = _cov(15,10)*_dt;
    const float var70 = _cov(12,10)*mrdt20 + _cov(13,10)*mrdt21 + _cov(14,10)*mrdt22 + _cov(10,5) + _cov(10,6)*x20 + _cov(10,7)*x21 + _cov(10,8)*x22 + var69;
    const float var71 = _cov(15,11)*_dt;
    const float var72 = _cov(12,11)*mrdt20 + _cov(13,11)*mrdt21 + _cov(14,11)*mrdt22 + _cov(11,5) + _cov(11,6)*x20 + _cov(11,7)*x21 + _cov(11,8)*x22 + var71;
    const float var73 = _cov(6,6)*y00 + _cov(7,6)*y01 + _cov(8,6)*y02 - _cov(9,6)*_dt;
    const float var74 = _cov(7,6)*y00 + _cov(7,7)*y01 + _cov(8,7)*y02 - _cov(9,7)*_dt;
    const float var75 = _cov(8,6)*y00 + _cov(8,7)*y01 + _cov(8,8)*y02 - _cov(9,8)*_dt;
    const float var76 = _cov(9,6)*y00 + _cov(9,7)*y01 + _cov(9,8)*y02 - _cov(9,9)*_dt;
    const float var77 = -_cov(10,9)*_dt;
    const float var78 = _cov(10,6)*y00 + _cov(10,7)*y01 + _cov(10,8)*y02 + var77;
    const float var79 = -_cov(11,9)*_dt;
    const float var80 = _cov(11,6)*y00 + _cov(11,7)*y01 + _cov(11,8)*y02 + var79;
    const float var81 = -_cov(10,6)*_dt + _cov(6,6)*y10 + _cov(7,6)*y11 + _cov(8,6)*y12;
    const float var82 = _cov(7,6)*y10 - _cov(10,7)*_dt + _cov(7,7)*y11 + _cov(8,7)*y12;
    const float var83 = _cov(8,6)*y10 + _cov(8,7)*y11 - _cov(10,8)*_dt + _cov(8,8)*y12;
    const float var84 = -_cov(10,10)*_dt + _cov(10,6)*y10 + _cov(10,7)*y11 + _cov(10,8)*y12;
    const float var85 = -_cov(11,10)*_dt;
    const float var86 = _cov(11,6)*y10 + _cov(11,7)*y11 + _cov(11,8)*y12 + var85;
    const float var87 = -_cov(11,11)*_dt + _cov(11,6)*y20 + _cov(11,7)*y21 + _cov(11,8)*y22;

    _cov(0,0) = _cov(0,0) + _cov(3,0)*_dt + _dt*var0;
    _cov(0,1) = _cov(1,0) + _cov(3,1)*_dt + _dt*var2;
    _cov(1,1) = _cov(1,1) + _cov(4,1)*_dt + _dt*var15;
    _cov(0,2) = _cov(2,0) + _cov(3,2)*_dt + _dt*var4;
    _cov(1,2) = _cov(2,1) + _cov(4,2)*_dt + _dt*var17;
    _cov(2,2) = _cov(2,2) + _cov(5,2)*_dt + _dt*var28;
    _cov(0,3) = mrdt00*var5 + mrdt01*var6 + mrdt02*var7 + var0 + var10*x02 + var8*x00 + var9*x01;
    _cov(1,3) = _cov(3,1) + mrdt00*var18 + mrdt01*var19 + mrdt02*var20 + var1 + var21*x00 + var22*x01 + var23*x02;
    _cov(2,3) = _cov(3,2) + mrdt00*var29 + mrdt01*var30 + mrdt02*var31 + var3 + var32*x00 + var33*x01 + var34*x02;
    _cov(3,3) = _cov(12,3)*mrdt00 + _cov(13,3)*mrdt01 + _cov(14,3)*mrdt02 + _cov(3,3) + _cov(6,3)*x00 + _cov(7,3)*x01 + _cov(8,3)*x02 + mrdt00*var40 + mrdt01*var41 + mrdt02*var42 + var43*x00 + var44*x01 + var45*x02;
    _cov(0,4) = mrdt10*var5 + mrdt11*var6 + mrdt12*var7 + var10*x12 + var2 + var8*x10 + var9*x11;
    _cov(1,4) = mrdt10*var18 + mrdt11*var19 + mrdt12*var20 + var15 + var21*x10 + var22*x11 + var23*x12;
    _cov(2,4) = _cov(4,2) + mrdt10*var29 + mrdt11*var30 + mrdt12*var31 + var16 + var32*x10 + var33*x11 + var34*x12;
    _cov(3,4) = _cov(4,3) + _cov(12,4)*mrdt00 + _cov(13,4)*mrdt01 + _cov(14,4)*mrdt02 + _cov(6,4)*x00 + _cov(7,4)*x01 + _cov(8,4)*x02 + mrdt10*var40 + mrdt11*var41 + mrdt12*var42 + var43*x10 + var44*x11 + var45*x12;
    _cov(4,4) = _cov(12,4)*mrdt10 + _cov(13,4)*mrdt11 + _cov(14,4)*mrdt12 + _cov(4,4) + _cov(6,4)*x10 + _cov(7,4)*x11 + _cov(8,4)*x12 + mrdt10*var50 + mrdt11*var51 + mrdt12*var52 + var53*x10 + var54*x11 + var55*x12;
    _cov(0,5) = _dt*var11 + mrdt20*var5 + mrdt21*var6 + mrdt22*var7 + var10*x22 + var4 + var8*x20 + var9*x21;
    _cov(1,5) = _dt*var24 + mrdt20*var18 + mrdt21*var19 + mrdt22*var20 + var17 + var21*x20 + var22*x21 + var23*x22;
    _cov(2,5) = _dt*var36 + mrdt20*var29 + mrdt21*var30 + mrdt22*var31 + var28 + var32*x20 + var33*x21 + var34*x22;
    _cov(3,5) = _cov(5,3) + _cov(12,5)*mrdt00 + _cov(13,5)*mrdt01 + _cov(14,5)*mrdt02 + _cov(6,5)*x00 + _cov(7,5)*x01 + _cov(8,5)*x02 + _dt*var46 + mrdt20*var40 + mrdt21*var41 + mrdt22*var42 + var43*x20 + var44*x21 + var45*x22;
    _cov(4,5) = _cov(5,4) + _cov(12,5)*mrdt10 + _cov(13,5)*mrdt11 + _cov(14,5)*mrdt12 + _cov(6,5)*x10 + _cov(7,5)*x11 + _cov(8,5)*x12 + _dt*var56 + mrdt20*var50 + mrdt21*var51 + mrdt22*var52 + var53*x20 + var54*x21 + var55*x22;
    _cov(5,5) = _cov(12,5)*mrdt20 + _cov(13,5)*mrdt21 + _cov(14,5)*mrdt22 + _cov(5,5) + _cov(6,5)*x20 + _cov(7,5)*x21 + _cov(8,5)*x22 + _dt*var60 + mrdt20*var61 + mrdt21*var62 + mrdt22*var63 + var35 + var64*x20 + var65*x21 + var66*x22;
    _cov(0,8) = -_dt*var14 + var10*y22 + var8*y20 + var9*y21;
    _cov(1,8) = -_dt*var27 + var21*y20 + var22*y21 + var23*y22;
    _cov(2,8) = -_dt*var39 + var32*y20 + var33*y21 + var34*y22;
    _cov(3,8) = -_dt*var49 + var43*y20 + var44*y21 + var45*y22;
    _cov(4,8) = -_dt*var59 + var53*y20 + var54*y21 + var55*y22;
    _cov(5,8) = -_dt*var72 + var64*y20 + var65*y21 + var66*y22;
    _cov(6,8) = -_dt*var80 + var73*y20 + var74*y21 + var75*y22;
    _cov(7,8) = -_dt*var86 + var81*y20 + var82*y21 + var83*y22;
    _cov(8,8) = -_dt*var87 + y20*(-_cov(11,6)*_dt + _cov(6,6)*y20 + _cov(7,6)*y21 + _cov(8,6)*y22) + y21*(_cov(7,6)*y20 - _cov(11,7)*_dt + _cov(7,7)*y21 + _cov(8,7)*y22) + y22*(_cov(8,6)*y20 + _cov(8,7)*y21 - _cov(11,8)*_dt + _cov(8,8)*y22);
    _cov(0,6) = -_dt*var12 + var10*y02 + var8*y00 + var9*y01;
    _cov(1,6) = -_dt*var25 + var21*y00 + var22*y01 + var23*y02;
    _cov(2,6) = -_dt*var37 + var32*y00 + var33*y01 + var34*y02;
    _cov(3,6) = -_dt*var47 + var43*y00 + var44*y01 + var45*y02;
    _cov(4,6) = -_dt*var57 + var53*y00 + var54*y01 + var55*y02;
    _cov(5,6) = -_dt*var68 + var64*y00 + var65*y01 + var66*y02;
    _cov(6,6) = -_dt*var76 + var73*y00 + var74*y01 + var75*y02;
    _cov(0,7) = -_dt*var13 + var10*y12 + var8*y10 + var9*y11;
    _cov(1,7) = -_dt*var26 + var21*y10 + var22*y11 + var23*y12;
    _cov(2,7) = -_dt*var38 + var32*y10 + var33*y11 + var34*y12;
    _cov(3,7) = -_dt*var48 + var43*y10 + var44*y11 + var45*y12;
    _cov(4,7) = -_dt*var58 + var53*y10 + var54*y11 + var55*y12;
    _cov(5,7) = -_dt*var70 + var64*y10 + var65*y11 + var66*y12;
    _cov(6,7) = -_dt*var78 + var73*y10 + var74*y11 + var75*y12;
    _cov(7,7) = -_dt*var84 + var81*y10 + var82*y11 + var83*y12;
    _cov(0,9) = var12;
    _cov(1,9) = var25;
    _cov(2,9) = var37;
    _cov(3,9) = var47;
    _cov(4,9) = var57;
    _cov(5,9) = var68;
    _cov(6,9) = var76;
    _cov(7,9) = _cov(9,6)*y10 + _cov(9,7)*y11 + _cov(9,8)*y12 + var77;
    _cov(8,9) = _cov(9,6)*y20 + _cov(9,7)*y21 + _cov(9,8)*y22 + var79;
    _cov(9,9) = _cov(9,9);
    _cov(0,10) = var13;
    _cov(1,10) = var26;
    _cov(2,10) = var38;
    _cov(3,10) = var48;
    _cov(4,10) = var58;
    _cov(5,10) = var70;
    _cov(6,10) = var78;
    _cov(7,10) = var84;
    _cov(8,10) = _cov(10,6)*y20 + _cov(10,7)*y21 + _cov(10,8)*y22 + var85;
    _cov(0,11) = var14;
    _cov(1,11) = var27;
    _cov(2,11) = var39;
    _cov(3,11) = var49;
    _cov(4,11) = var59;
    _cov(5,11) = var72;
    _cov(6,11) = var80;
    _cov(7,11) = var86;
    _cov(8,11) = var87;
    _cov(0,12) = var5;
    _cov(1,12) = var18;
    _cov(2,12) = var29;
    _cov(3,12) = var40;
    _cov(4,12) = var50;
    _cov(5,12) = var61;
    _cov(6,12) = _cov(12,6)*y00 + _cov(12,7)*y01 + _cov(12,8)*y02 - _cov(12,9)*_dt;
    _cov(7,12) = -_cov(12,10)*_dt + _cov(12,6)*y10 + _cov(12,7)*y11 + _cov(12,8)*y12;
    _cov(8,12) = -_cov(12,11)*_dt + _cov(12,6)*y20 + _cov(12,7)*y21 + _cov(12,8)*y22;
    _cov(0,13) = var6;
    _cov(1,13) = var19;
    _cov(2,13) = var30;
    _cov(3,13) = var41;
    _cov(4,13) = var51;
    _cov(5,13) = var62;
    _cov(6,13) = _cov(13,6)*y00 + _cov(13,7)*y01 + _cov(13,8)*y02 - _cov(13,9)*_dt;
    _cov(7,13) = -_cov(13,10)*_dt + _cov(13,6)*y10 + _cov(13,7)*y11 + _cov(13,8)*y12;
    _cov(8,13) = -_cov(13,11)*_dt + _cov(13,6)*y20 + _cov(13,7)*y21 + _cov(13,8)*y22;
    _cov(0,14) = var7;
    _cov(1,14) = var20;
    _cov(2,14) = var31;
    _cov(3,14) = var42;
    _cov(4,14) = var52;
    _cov(5,14) = var63;
    _cov(6,14) = _cov(14,6)*y00 + _cov(14,7)*y01 + _cov(14,8)*y02 - _cov(14,9)*_dt;
    _cov(7,14) = -_cov(14,10)*_dt + _cov(14,6)*y10 + _cov(14,7)*y11 + _cov(14,8)*y12;
    _cov(8,14) = -_cov(14,11)*_dt + _cov(14,6)*y20 + _cov(14,7)*y21 + _cov(14,8)*y22;
    _cov(0,15) = var11;
    _cov(1,15) = var24;
    _cov(2,15) = var36;
    _cov(3,15) = var46;
    _cov(4,15) = var56;
    _cov(5,15) = var60;
    _cov(6,15) = _cov(15,6)*y00 + _cov(15,7)*y01 + _cov(15,8)*y02 - var67;
    _cov(7,15) = _cov(15,6)*y10 + _cov(15,7)*y11 + _cov(15,8)*y12 - var69;
    _cov(8,15) = _cov(15,6)*y20 + _cov(15,7)*y21 + _cov(15,8)*y22 - var71;

    // F *P *F' + Q
    _cov(0, 0) = kahan_summation(_cov(0, 0), _q_cov(0) * _dt2, _accumulator_cov(0));
    _cov(1, 1) = kahan_summation(_cov(1, 1), _q_cov(1) * _dt2, _accumulator_cov(1));
    _cov(2, 2) = kahan_summation(_cov(2, 2), _q_cov(2) * _dt2, _accumulator_cov(2));
    if (_q_cov(3) == _q_cov(4)) {
        if (_q_cov(4) == _q_cov(5)) {
            _cov(3, 3) = kahan_summation(_cov(3, 3), _q_cov(3) * _dt2, _accumulator_cov(3));
            _cov(4, 4) = kahan_summation(_cov(4, 4), _q_cov(4) * _dt2, _accumulator_cov(4));
            _cov(5, 5) = kahan_summation(_cov(5, 5), _q_cov(5) * _dt2, _accumulator_cov(5));
        } else {
            float sx = _q_cov(3) * _dt2;
            float sz = _q_cov(5) * _dt2;
            float r00_sx = _rot(0, 0) * sx;
            float r01_sx = _rot(0, 1) * sx;
            float r02_sz = _rot(0, 2) * sz;
            float r00_sx_r10_plus_r01_sx_r11_plus_r02_sz_r12 = r00_sx * _rot(1, 0) + r01_sx * _rot(1, 1) + r02_sz * _rot(1, 2);
            float r00_sx_r20_plus_r01_sx_r21_plus_r02_sz_r22 = r00_sx * _rot(2, 0) + r01_sx * _rot(2, 1) + r02_sz * _rot(2, 2);
            float r10_r20_sx_plus_r11_r21_sx_plus_r12_r22_sz = (_rot(1, 0) * _rot(2, 0) + _rot(1, 1) * _rot(2, 1)) * sx + _rot(1, 2) * _rot(2, 2) * sz;

            _cov(3, 3) += kahan_summation(_cov(3, 3), (_rot(0, 0) * _rot(0, 0) + _rot(0, 1) * _rot(0, 1)) * sx + _rot(0, 2) * _rot(0, 2) * sz, _accumulator_cov(3));
            _cov(3, 4) += r00_sx_r10_plus_r01_sx_r11_plus_r02_sz_r12;
            _cov(3, 5) += r00_sx_r20_plus_r01_sx_r21_plus_r02_sz_r22;
            _cov(4, 4) += kahan_summation(_cov(4, 4), (_rot(1, 0) * _rot(1, 0) + _rot(1, 1) * _rot(1, 1)) * sx + _rot(1, 2) * _rot(1, 2) * sz, _accumulator_cov(4));
            _cov(4, 5) += r10_r20_sx_plus_r11_r21_sx_plus_r12_r22_sz;
            _cov(5, 5) += kahan_summation(_cov(5, 5), (_rot(2, 0) * _rot(2, 0) + _rot(2, 1) * _rot(2, 1)) * sx + _rot(2, 2) * _rot(2, 2) * sz, _accumulator_cov(5));
        }
    } else {
        float sx = _q_cov(3) * _dt2;
        float sy = _q_cov(4) * _dt2;
        float sz = _q_cov(5) * _dt2;
        float r00_sx = _rot(0, 0) * sx;
        float r01_sy = _rot(0, 1) * sy;
        float r02_sz = _rot(0, 2) * sz;
        float r00_sx_r10_plus_r01_sy_r11_plus_r02_sz_r12 = r00_sx * _rot(1, 0) + r01_sy * _rot(1, 1) + r02_sz * _rot(1, 2);
        float r00_sx_r20_plus_r01_sy_r21_plus_r02_sz_r22 = r00_sx * _rot(2, 0) + r01_sy * _rot(2, 1) + r02_sz * _rot(2, 2);
        float r10_r20_sx_plus_r11_r21_sy_plus_r12_r22_sz = _rot(1, 0) * _rot(2, 0) * sx + _rot(1, 1) * _rot(2, 1) * sy + _rot(1, 2) * _rot(2, 2) * sz;
        
        _cov(3, 3) += kahan_summation(_cov(3, 3), _rot(0, 0) * _rot(0, 0) * sx + _rot(0, 1) * _rot(0, 1) * sy + _rot(0, 2) * _rot(0, 2) * sz, _accumulator_cov(3));
        _cov(3, 4) += r00_sx_r10_plus_r01_sy_r11_plus_r02_sz_r12;
        _cov(3, 5) += r00_sx_r20_plus_r01_sy_r21_plus_r02_sz_r22;
        _cov(4, 4) += kahan_summation(_cov(4, 4), _rot(1, 0) * _rot(1, 0) * sx + _rot(1, 1) * _rot(1, 1) * sy + _rot(1, 2) * _rot(1, 2) * sz, _accumulator_cov(4));
        _cov(4, 5) += r10_r20_sx_plus_r11_r21_sy_plus_r12_r22_sz;
        _cov(5, 5) += kahan_summation(_cov(5, 5), _rot(2, 0) * _rot(2, 0) * sx + _rot(2, 1) * _rot(2, 1) * sy + _rot(2, 2) * _rot(2, 2) * sz, _accumulator_cov(5));
    }
    _cov(6, 6) = kahan_summation(_cov(6, 6), _q_cov(6) * _dt2, _accumulator_cov(6));
    _cov(7, 7) = kahan_summation(_cov(7, 7), _q_cov(7) * _dt2, _accumulator_cov(7));
    _cov(8, 8) = kahan_summation(_cov(8, 8), _q_cov(8) * _dt2, _accumulator_cov(8));
    _cov(9, 9) = kahan_summation(_cov(9, 9), _q_cov(9) * _dt2, _accumulator_cov(9));
    _cov(10, 10) = kahan_summation(_cov(10, 10), _q_cov(10) * _dt2, _accumulator_cov(10));
    _cov(11, 11) = kahan_summation(_cov(11, 11), _q_cov(11) * _dt2, _accumulator_cov(11));
    _cov(12, 12) = kahan_summation(_cov(12, 12), _q_cov(12) * _dt2, _accumulator_cov(12));
    _cov(13, 13) = kahan_summation(_cov(13, 13), _q_cov(13) * _dt2, _accumulator_cov(13));
    _cov(14, 14) = kahan_summation(_cov(14, 14), _q_cov(14) * _dt2, _accumulator_cov(14));
    _cov(15, 15) = kahan_summation(_cov(15, 15), _q_cov(15) * _dt2, _accumulator_cov(15));

    regular_covariance_to_symmetric();
}

unsigned char ESKF::fuse_position(const Vector3f &pos, const Vector3f &w, const Vector3f &a, 
                                  const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) {
    /*
    s = (p, v, R, bg, ba, g)
    δs = (δp, δv, δθ, δbg, δba, δg)
    s + δs = (p+δp, v+δv, R*Exp(δθ), bg+δbg, ba+δba, g+δg)

    pos = p + R * dis

    δpos / δp = I
    δpos / δθ = -R * dis^

    H = (I, O, -R*dis^, O, O, O)
    */

    unsigned char info = 0;

    // Sequential Kalman Filter
    for (unsigned int dim = 0; dim < 3; ++dim) {
        if (!isfinite(noise_std(dim))) {
            continue;
        }

        /*
        K = P * H * (H * P * H' + R)^-1
        x = x + K * (y - h)
        P = (I - K * H) * P

        where, H = (I, O, -R*dis^, O, O, O)
               h = p + R * dis
        */

        // -R * dis^
        const Vector3f minus_rot_d_hat = {_rot(dim, 2)*dis(1) - _rot(dim, 1)*dis(2), _rot(dim, 0)*dis(2) - _rot(dim, 2)*dis(0), _rot(dim, 1)*dis(0) - _rot(dim, 0)*dis(1)};

        // H * P  or  P * H'
        const float cov_1_dim = (dim == 2) ? _cov(1, dim) : _cov(dim, 1);
        Matrix<float, 16, 1> HP;
        HP(0) = _cov(0, dim) + _cov(0, 6)*minus_rot_d_hat(0) + _cov(0, 7)*minus_rot_d_hat(1) + _cov(0, 8)*minus_rot_d_hat(2);
        HP(1) = cov_1_dim + _cov(1, 6)*minus_rot_d_hat(0) + _cov(1, 7)*minus_rot_d_hat(1) + _cov(1, 8)*minus_rot_d_hat(2);
        HP(2) = _cov(dim, 2) + _cov(2, 6)*minus_rot_d_hat(0) + _cov(2, 7)*minus_rot_d_hat(1) + _cov(2, 8)*minus_rot_d_hat(2);
        HP(3) = _cov(dim, 3) + _cov(3, 6)*minus_rot_d_hat(0) + _cov(3, 7)*minus_rot_d_hat(1) + _cov(3, 8)*minus_rot_d_hat(2);
        HP(4) = _cov(dim, 4) + _cov(4, 6)*minus_rot_d_hat(0) + _cov(4, 7)*minus_rot_d_hat(1) + _cov(4, 8)*minus_rot_d_hat(2);
        HP(5) = _cov(dim, 5) + _cov(5, 6)*minus_rot_d_hat(0) + _cov(5, 7)*minus_rot_d_hat(1) + _cov(5, 8)*minus_rot_d_hat(2);
        HP(6) = _cov(dim, 6) + _cov(6, 6)*minus_rot_d_hat(0) + _cov(6, 7)*minus_rot_d_hat(1) + _cov(6, 8)*minus_rot_d_hat(2);
        HP(7) = _cov(dim, 7) + _cov(6, 7)*minus_rot_d_hat(0) + _cov(7, 7)*minus_rot_d_hat(1) + _cov(7, 8)*minus_rot_d_hat(2);
        HP(8) = _cov(dim, 8) + _cov(6, 8)*minus_rot_d_hat(0) + _cov(7, 8)*minus_rot_d_hat(1) + _cov(8, 8)*minus_rot_d_hat(2);
        HP(9) = _cov(dim, 9) + _cov(6, 9)*minus_rot_d_hat(0) + _cov(7, 9)*minus_rot_d_hat(1) + _cov(8, 9)*minus_rot_d_hat(2);
        HP(10) = _cov(dim, 10) + _cov(6, 10)*minus_rot_d_hat(0) + _cov(7, 10)*minus_rot_d_hat(1) + _cov(8, 10)*minus_rot_d_hat(2);
        HP(11) = _cov(dim, 11) + _cov(6, 11)*minus_rot_d_hat(0) + _cov(7, 11)*minus_rot_d_hat(1) + _cov(8, 11)*minus_rot_d_hat(2);
        HP(12) = _cov(dim, 12) + _cov(6, 12)*minus_rot_d_hat(0) + _cov(7, 12)*minus_rot_d_hat(1) + _cov(8, 12)*minus_rot_d_hat(2);
        HP(13) = _cov(dim, 13) + _cov(6, 13)*minus_rot_d_hat(0) + _cov(7, 13)*minus_rot_d_hat(1) + _cov(8, 13)*minus_rot_d_hat(2);
        HP(14) = _cov(dim, 14) + _cov(6, 14)*minus_rot_d_hat(0) + _cov(7, 14)*minus_rot_d_hat(1) + _cov(8, 14)*minus_rot_d_hat(2);
        HP(15) = _cov(dim, 15) + _cov(6, 15)*minus_rot_d_hat(0) + _cov(7, 15)*minus_rot_d_hat(1) + _cov(8, 15)*minus_rot_d_hat(2);

        // H * P * H' + R
        const float HPHT_plus_R = HP(dim) + HP(6) * minus_rot_d_hat(0) + HP(7) * minus_rot_d_hat(1) + HP(8) * minus_rot_d_hat(2) + noise_std(dim) * noise_std(dim);

        // h = p + R * dis
        // e = y - h = pos - (p + R * dis)
        const float obs_error = pos(dim) - (_p(dim) + _rot(dim, 0) * dis(0) + _rot(dim, 1) * dis(1) + _rot(dim, 2) * dis(2));

        /*
        K = P * H' * (H * P * H' + R)^-1
        P = P - K * H * P

        e = y - h = y - (v + R * (w - bg)^ * dis)
        x = x + K * (y - h)
        */
        switch (conservative_posteriori_estimate(HP, HPHT_plus_R, obs_error, gate(dim))) {
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
    regular_covariance_to_symmetric();

    return info;
}

unsigned char ESKF::fuse_velocity(const Vector3f &vel, const Vector3f &w, const Vector3f &a, 
                                  const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) {
    /*
    s = (p, v, R, bg, ba, g)
    δs = (δp, δv, δθ, δbg, δba, δg)
    s + δs = (p+δp, v+δv, R*Exp(δθ), bg+δbg, ba+δba, g+δg)

    vel = v + R * (w - bg)^ * dis

    δvel / δv = I
    δvel / δθ = -R * ((w - bg)^ * dis)^
    δvel / δbg = R * dis^

    H = (O, I, -R*((w-bg)^*dis)^, R*dis^, O, O)
    */ 

    unsigned char info = 0;

    // (w-bg)^ * dis
    const Vector3f w_corr = {w(0) - _bg(0), w(1) - _bg(1), w(2) - _bg(2)};
    const Vector3f w_corr_cross_d = {w_corr(1) * dis(2) - w_corr(2) * dis(1), 
                                     w_corr(2) * dis(0) - w_corr(0) * dis(2), 
                                     w_corr(0) * dis(1) - w_corr(1) * dis(0)};
    
    // Sequential Kalman Filter
    for (unsigned int dim = 0; dim < 3; ++dim) {
        if (!isfinite(noise_std(dim))) {
            continue;
        }

        /*
        K = P * H * (H * P * H' + R)^-1
        x = x + K * (y - h)
        P = (I - K * H) * P

        where, H = (O, I, R*(dis^*(w-bg))^, R*dis^, O, O)
               h = v + R * (w-bg)^ * dis
        */

        // R * dis^
        const Vector3f rot_d_hat = {_rot(dim, 1)*dis(2) - _rot(dim, 2)*dis(1), 
                                    _rot(dim, 2)*dis(0) - _rot(dim, 0)*dis(2), 
                                    _rot(dim, 0)*dis(1) - _rot(dim, 1)*dis(0)};

        // R * (dis^*(w-bg))^ = - R * ((w-bg)^*dis)^
        const Vector3f rot_d_cross_w_corr_hat = {_rot(dim, 2)*w_corr_cross_d(1) - _rot(dim, 1)*w_corr_cross_d(2), 
                                                 _rot(dim, 0)*w_corr_cross_d(2) - _rot(dim, 2)*w_corr_cross_d(0), 
                                                 _rot(dim, 1)*w_corr_cross_d(0) - _rot(dim, 0)*w_corr_cross_d(1)};

        // H * P  or  P * H'
        const unsigned int index = 3 + dim;
        const float cov_4_index = (dim == 2) ? _cov(4, index) : _cov(index, 4);
        Matrix<float, 16, 1> HP; 
        HP(0) = _cov(0, index) + _cov(0, 6)*rot_d_cross_w_corr_hat(0) + _cov(0, 7)*rot_d_cross_w_corr_hat(1) + _cov(0, 8)*rot_d_cross_w_corr_hat(2) + _cov(0, 9)*rot_d_hat(0) + _cov(0, 10)*rot_d_hat(1) + _cov(0, 11)*rot_d_hat(2);
        HP(1) = _cov(1, index) + _cov(1, 6)*rot_d_cross_w_corr_hat(0) + _cov(1, 7)*rot_d_cross_w_corr_hat(1) + _cov(1, 8)*rot_d_cross_w_corr_hat(2) + _cov(1, 9)*rot_d_hat(0) + _cov(1, 10)*rot_d_hat(1) + _cov(1, 11)*rot_d_hat(2);
        HP(2) = _cov(2, index) + _cov(2, 6)*rot_d_cross_w_corr_hat(0) + _cov(2, 7)*rot_d_cross_w_corr_hat(1) + _cov(2, 8)*rot_d_cross_w_corr_hat(2) + _cov(2, 9)*rot_d_hat(0) + _cov(2, 10)*rot_d_hat(1) + _cov(2, 11)*rot_d_hat(2);
        HP(3) = _cov(3, index) + _cov(3, 6)*rot_d_cross_w_corr_hat(0) + _cov(3, 7)*rot_d_cross_w_corr_hat(1) + _cov(3, 8)*rot_d_cross_w_corr_hat(2) + _cov(3, 9)*rot_d_hat(0) + _cov(3, 10)*rot_d_hat(1) + _cov(3, 11)*rot_d_hat(2);
        HP(4) = cov_4_index + _cov(4, 6)*rot_d_cross_w_corr_hat(0) + _cov(4, 7)*rot_d_cross_w_corr_hat(1) + _cov(4, 8)*rot_d_cross_w_corr_hat(2) + _cov(4, 9)*rot_d_hat(0) + _cov(4, 10)*rot_d_hat(1) + _cov(4, 11)*rot_d_hat(2);
        HP(5) = _cov(index, 5) + _cov(5, 6)*rot_d_cross_w_corr_hat(0) + _cov(5, 7)*rot_d_cross_w_corr_hat(1) + _cov(5, 8)*rot_d_cross_w_corr_hat(2) + _cov(5, 9)*rot_d_hat(0) + _cov(5, 10)*rot_d_hat(1) + _cov(5, 11)*rot_d_hat(2);
        HP(6) = _cov(index, 6) + _cov(6, 6)*rot_d_cross_w_corr_hat(0) + _cov(6, 7)*rot_d_cross_w_corr_hat(1) + _cov(6, 8)*rot_d_cross_w_corr_hat(2) + _cov(6, 9)*rot_d_hat(0) + _cov(6, 10)*rot_d_hat(1) + _cov(6, 11)*rot_d_hat(2);
        HP(7) = _cov(index, 7) + _cov(6, 7)*rot_d_cross_w_corr_hat(0) + _cov(7, 7)*rot_d_cross_w_corr_hat(1) + _cov(7, 8)*rot_d_cross_w_corr_hat(2) + _cov(7, 9)*rot_d_hat(0) + _cov(7, 10)*rot_d_hat(1) + _cov(7, 11)*rot_d_hat(2);
        HP(8) = _cov(index, 8) + _cov(6, 8)*rot_d_cross_w_corr_hat(0) + _cov(7, 8)*rot_d_cross_w_corr_hat(1) + _cov(8, 8)*rot_d_cross_w_corr_hat(2) + _cov(8, 9)*rot_d_hat(0) + _cov(8, 10)*rot_d_hat(1) + _cov(8, 11)*rot_d_hat(2);
        HP(9) = _cov(index, 9) + _cov(6, 9)*rot_d_cross_w_corr_hat(0) + _cov(7, 9)*rot_d_cross_w_corr_hat(1) + _cov(8, 9)*rot_d_cross_w_corr_hat(2) + _cov(9, 9)*rot_d_hat(0) + _cov(9, 10)*rot_d_hat(1) + _cov(9, 11)*rot_d_hat(2);
        HP(10) = _cov(index, 10) + _cov(6, 10)*rot_d_cross_w_corr_hat(0) + _cov(7, 10)*rot_d_cross_w_corr_hat(1) + _cov(8, 10)*rot_d_cross_w_corr_hat(2) + _cov(9, 10)*rot_d_hat(0) + _cov(10, 10)*rot_d_hat(1) + _cov(10, 11)*rot_d_hat(2);
        HP(11) = _cov(index, 11) + _cov(6, 11)*rot_d_cross_w_corr_hat(0) + _cov(7, 11)*rot_d_cross_w_corr_hat(1) + _cov(8, 11)*rot_d_cross_w_corr_hat(2) + _cov(9, 11)*rot_d_hat(0) + _cov(10, 11)*rot_d_hat(1) + _cov(11, 11)*rot_d_hat(2);
        HP(12) = _cov(index, 12) + _cov(6, 12)*rot_d_cross_w_corr_hat(0) + _cov(7, 12)*rot_d_cross_w_corr_hat(1) + _cov(8, 12)*rot_d_cross_w_corr_hat(2) + _cov(9, 12)*rot_d_hat(0) + _cov(10, 12)*rot_d_hat(1) + _cov(11, 12)*rot_d_hat(2);
        HP(13) = _cov(index, 13) + _cov(6, 13)*rot_d_cross_w_corr_hat(0) + _cov(7, 13)*rot_d_cross_w_corr_hat(1) + _cov(8, 13)*rot_d_cross_w_corr_hat(2) + _cov(9, 13)*rot_d_hat(0) + _cov(10, 13)*rot_d_hat(1) + _cov(11, 13)*rot_d_hat(2);
        HP(14) = _cov(index, 14) + _cov(6, 14)*rot_d_cross_w_corr_hat(0) + _cov(7, 14)*rot_d_cross_w_corr_hat(1) + _cov(8, 14)*rot_d_cross_w_corr_hat(2) + _cov(9, 14)*rot_d_hat(0) + _cov(10, 14)*rot_d_hat(1) + _cov(11, 14)*rot_d_hat(2);
        HP(15) = _cov(index, 15) + _cov(6, 15)*rot_d_cross_w_corr_hat(0) + _cov(7, 15)*rot_d_cross_w_corr_hat(1) + _cov(8, 15)*rot_d_cross_w_corr_hat(2) + _cov(9, 15)*rot_d_hat(0) + _cov(10, 15)*rot_d_hat(1) + _cov(11, 15)*rot_d_hat(2);                                                  
        // H * P * H' + R
        const float HPHT_plus_R = HP(index) + HP(6) * rot_d_cross_w_corr_hat(0) + HP(7) * rot_d_cross_w_corr_hat(1) + HP(8) * rot_d_cross_w_corr_hat(2) + HP(9) * rot_d_hat(0) + HP(10) * rot_d_hat(1) + HP(11) * rot_d_hat(2) + noise_std(dim) * noise_std(dim);

        // h = v + R * (w - bg)^ * dis
        // e = y - h = y - (v + R * (w - bg)^ * dis)
        const float obs_error = vel(dim) - (_v(dim) + _rot(dim, 0) * w_corr_cross_d(0) + _rot(dim, 1) * w_corr_cross_d(1) + _rot(dim, 2) * w_corr_cross_d(2));

        /*
        K = P * H' * (H * P * H' + R)^-1
        P = P - K * H * P

        e = y - h = y - (v + R * (w - bg)^ * dis)
        x = x + K * (y - h)
        */
        switch (conservative_posteriori_estimate(HP, HPHT_plus_R, obs_error, gate(dim))) {
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
    regular_covariance_to_symmetric();

    return info;
}

unsigned char ESKF::conservative_posteriori_estimate(const Matrix<float, 16, 1> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
    /*
    K = P * H' * (H * P * H' + R)^-1
    P = P - K * H * P

    e = y - h = y - (v + R * (w - bg)^ * dis)
    x = x + K * (y - h)
    */
    
    // Don't correct error state unless obs_error / √(H*P*H' + R) >= gate
    if (obs_error * obs_error  >= gate * gate * HPHT_plus_R) {
        return 1;
    }

    // Don't correct error state unless (I - KH) * P > 0 <=> P - PH'(HPH'+R)^-1HP > 0
    for (unsigned char i = 0; i < 16; ++i) {
        if (HP(i) * HP(i) >= HPHT_plus_R * _cov(i, i)) {
            return 2;    
        }
    }
    
    // K = P * H' * (H * P * H' + R)^-1
    Matrix<float, 16, 1> K {};
    for (unsigned char i = 0; i < 16; ++i) {
        K(i) = HP(i) / HPHT_plus_R;
    }

    // x = x + K * e
    for (unsigned char i = 0; i < 16; ++i) {
        _error_state(i) += K(i) * obs_error;
    }

    // P = P - K * H * P
    for (unsigned char i = 0; i < 16; ++i) {
        for (unsigned char j = i; j < 16; ++j) {
            _cov(i, j) -= K(i) * HP(j);
        }
    }
    return 0;
}

void ESKF::correct_state() {
    // for (float &es : _error_state) {
    //     cout << es << endl;
    // }
    // cout << "-------------------------" << endl;
    // state: [p, v, bg, ba, g], R
    // error_state : [δp, δv, δθ, δbg, δba, δg]

    /*
    p = p + δp
    v = v + δv
    bg = bg + δbg
    ba = ba + δba
    g = g + δg
    */
    for (unsigned char i = 0; i < 3; ++i) {
        _p(i) += _error_state(i);
        _v(i) += _error_state(3 + i);
        _bg(i) += _error_state(9 + i);
        _ba(i) += _error_state(12 + i);
    }
    _g += _error_state(15);

    // q = q * Exp(δθ)
    Quaternionf delta_q;
    Vector3f delta_theta = {_error_state(6), _error_state(7), _error_state(8)};
    quaternion_from_axis_angle(delta_q, delta_theta);
    const Quaternionf q = _q;
    _q = q * delta_q;

    _rot = q;

    // (δp, δv, δθ, δbg, δba, δg) = 0
    _error_state.setZero();
}

void ESKF::reset_state() {
    _p.setZero();
    _v.setZero();
    _bg.setZero();
    _ba.setZero();
    _g = _g_init;
    _rot.setIdentity();
    _q.setIdentity();
}

void ESKF::reset_error_state() { 
    _error_state.setZero();
}

void ESKF::reset_priori_covariance_matrix() {
    _q_cov.setZero();
}

void ESKF::reset_covariance_matrix(const unsigned int start_index, const unsigned int end_index, const Matrix<float, 16, 1> &diag_cov) {
    for (unsigned int i = start_index; i < end_index; ++i) {
        // Diaginal
        _cov(i, i) = _q_cov(i) * _dt2;

        // Upper triangular
        for (unsigned int j = start_index; j < i; ++j) {
            _cov(j, i) = 0.f;
        }

        // Columns
        for (unsigned int j = 0; j < start_index; ++j) {
            _cov(j, i) = 0.f;
        }

        // Rows
        for (unsigned int j = end_index; j < 16; ++j) {
            _cov(i, j) = 0.f;
        }
    }
}

void ESKF::reset_accmulator_cov() {
    _accumulator_cov.setZero();
}

void ESKF::regular_covariance_to_symmetric(unsigned int start_index, unsigned int end_index) {
    for (unsigned int i = start_index + 1; i < end_index; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            _cov(i, j) = _cov(j, i);
        }
    }
}

float ESKF::kahan_summation(float sum_previous, float input, float &accumulator) {
    /*
    accumulator中记录了sum_previous + y中, y舍弃掉的部分
    */
    const float y = input - accumulator;
    const float t = sum_previous + y;
    accumulator = (t - sum_previous) - y;
    return t;
}