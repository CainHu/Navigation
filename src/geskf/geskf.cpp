//
// Created by Cain on 2022/11/19.
//

#include "geskf.h"
#include <cfloat>
#include <iostream>

using namespace std;
using namespace geskf;

void GESKF::predict_covariance(const Vector3f &w, const Vector3f &a) {
    // // R * dt               (乘法: 3, 加法: 0)
    // const Matrix3f rot = _rot * _dt;

    // // x = R * (a - ba)^ * dt       (乘法: 9, 加法: 12)
    // const array<float, 3> a_corr = {a[0] - _ba[0], a[1] - _ba[1], a[2] - _ba[2]};
    // const array<array<float, 3>, 3> x = {{{rot(0, 1)*a_corr[2] - rot(0, 2)*a_corr[1], rot(0, 2)*a_corr[0] - rot(0, 0)*a_corr[2], rot(0, 0)*a_corr[1] - rot(0, 1)*a_corr[0]},
    //                                       {rot(1, 1)*a_corr[2] - rot(1, 2)*a_corr[1], rot(1, 2)*a_corr[0] - rot(1, 0)*a_corr[2], rot(1, 0)*a_corr[1] - rot(1, 1)*a_corr[0]},
    //                                       {rot(2, 1)*a_corr[2] - rot(2, 2)*a_corr[1], rot(2, 2)*a_corr[0] - rot(2, 0)*a_corr[2], rot(2, 0)*a_corr[1] - rot(2, 1)*a_corr[0]}}};

    // // y = Exp(-(w - bg) * Δt)      (乘法: 3, 加法: 3)
    // const array<float, 3> axis_angle = {(_bg[0] - w[0]) * _dt, (_bg[1] - w[1]) * _dt, (_bg[2] - w[2]) * _dt};
    // Matrix3f y {};
    // rotation_from_axis_angle(y, axis_angle);

    // // P11_ = P11 + (P21 + P12) * dt + P22 * dt^2
    // _cov[0][0] += _dt*(_cov[0][3] + _cov[0][3] + _cov[3][3]*_dt);
    // _cov[0][1] += _dt*(_cov[0][4] + _cov[1][3] + _cov[3][4]*_dt);
    // _cov[0][2] += _dt*(_cov[0][5] + _cov[2][3] + _cov[3][5]*_dt);
    // _cov[1][1] += _dt*(_cov[1][4] + _cov[1][4] + _cov[4][4]*_dt);
    // _cov[1][2] += _dt*(_cov[1][5] + _cov[2][4] + _cov[4][5]*_dt);
    // _cov[2][2] += _dt*(_cov[2][5] + _cov[2][5] + _cov[5][5]*_dt);

    // // P16_ = P16 + P26 * dt
    // _cov[0][15] += _cov[3][15] * _dt;
    // _cov[1][15] += _cov[4][15] * _dt;
    // _cov[2][15] += _cov[5][15] * _dt;

    // // P15_ = P15 + P25 * dt    (乘法: 9, 加法: 9)
    // _cov[0][12] += _cov[3][12] * _dt;
    // _cov[0][13] += _cov[3][13] * _dt;
    // _cov[0][14] += _cov[3][14] * _dt;
    // _cov[1][12] += _cov[4][12] * _dt;
    // _cov[1][13] += _cov[4][13] * _dt;
    // _cov[1][14] += _cov[4][14] * _dt;
    // _cov[2][12] += _cov[5][12] * _dt;
    // _cov[2][13] += _cov[5][13] * _dt;
    // _cov[2][14] += _cov[5][14] * _dt;

    // // P14_ = P14 + P24 * dt    (乘法: 9, 加法: 9)
    // _cov[0][9] += _cov[3][9] * _dt;
    // _cov[0][10] += _cov[3][10] * _dt;
    // _cov[0][11] += _cov[3][11] * _dt;
    // _cov[1][9] += _cov[4][9] * _dt;
    // _cov[1][10] += _cov[4][10] * _dt;
    // _cov[1][11] += _cov[4][11] * _dt;
    // _cov[2][9] += _cov[5][9] * _dt;
    // _cov[2][10] += _cov[5][10] * _dt;
    // _cov[2][11] += _cov[5][11] * _dt;

    // _cov[6][0] += _cov[3][6] * _dt;
    // _cov[7][0] += _cov[3][7] * _dt;
    // _cov[8][0] += _cov[3][8] * _dt;
    // _cov[6][1] += _cov[4][6] * _dt;
    // _cov[7][1] += _cov[4][7] * _dt;
    // _cov[8][1] += _cov[4][8] * _dt;
    // _cov[6][2] += _cov[5][6] * _dt;
    // _cov[7][2] += _cov[5][7] * _dt;
    // _cov[8][2] += _cov[5][8] * _dt;

    // _cov[0][6] = _cov[6][0] + _cov[0][9]*r[0][0] + _cov[0][10]*r[0][1] + _cov[0][11]*r[0][2];
    // _cov[0][7] = _cov[7][0] + _cov[0][9]*r[1][0] + _cov[0][10]*r[1][1] + _cov[0][11]*r[1][2];
    // _cov[0][8] = _cov[8][0] + _cov[0][9]*r[2][0] + _cov[0][10]*r[2][1] + _cov[0][11]*r[2][2];
    // _cov[1][6] = _cov[6][1] + _cov[1][9]*r[0][0] + _cov[1][10]*r[0][1] + _cov[1][11]*r[0][2];
    // _cov[1][7] = _cov[7][1] + _cov[1][9]*r[1][0] + _cov[1][10]*r[1][1] + _cov[1][11]*r[1][2];
    // _cov[1][8] = _cov[8][1] + _cov[1][9]*r[2][0] + _cov[1][10]*r[2][1] + _cov[1][11]*r[2][2];
    // _cov[2][6] = _cov[6][2] + _cov[2][9]*r[0][0] + _cov[2][10]*r[0][1] + _cov[2][11]*r[0][2];
    // _cov[2][7] = _cov[7][2] + _cov[2][9]*r[1][0] + _cov[2][10]*r[1][1] + _cov[2][11]*r[1][2];
    // _cov[2][8] = _cov[8][2] + _cov[2][9]*r[2][0] + _cov[2][10]*r[2][1] + _cov[2][11]*r[2][2];

    // _cov[3][0] += _cov[3][3] * _dt;
    // _cov[4][0] += _cov[3][4] * _dt;
    // _cov[5][0] += _cov[3][5] * _dt;
    // _cov[3][1] += _cov[4][3] * _dt;
    // _cov[4][1] += _cov[4][4] * _dt;
    // _cov[5][1] += _cov[4][5] * _dt;
    // _cov[3][2] += _cov[5][3] * _dt;
    // _cov[4][2] += _cov[5][4] * _dt;
    // _cov[5][2] += _cov[5][5] * _dt;

    // _cov[0][3] = _cov[3][0] - (_cov[8][0]*v[1]) + _cov[7][0]*v[2] + _cov[0][12]*r[0][0] + _cov[0][13]*r[0][1] + _cov[0][14]*r[0][2];
    // _cov[0][4] = _cov[4][0] + _cov[8][0]*v[0] - _cov[6][0]*v[2] + _cov[0][12]*r[1][0] + _cov[0][13]*r[1][1] + _cov[0][14]*r[1][2];
    // _cov[0][5] = _cov[5][0] - (_cov[7][0]*v[0]) + _cov[6][0]*v[1] + _cov[0][12]*r[2][0] + _cov[0][13]*r[2][1] + _cov[0][14]*r[2][2] + _cov[0][15] * _dt;
    // _cov[1][3] = _cov[3][1] - (_cov[8][1]*v[1]) + _cov[7][1]*v[2] + _cov[1][12]*r[0][0] + _cov[1][13]*r[0][1] + _cov[1][14]*r[0][2];
    // _cov[1][4] = _cov[4][1] + _cov[8][1]*v[0] - _cov[6][1]*v[2] + _cov[1][12]*r[1][0] + _cov[1][13]*r[1][1] + _cov[1][14]*r[1][2];
    // _cov[1][5] = _cov[5][1] - (_cov[7][1]*v[0]) + _cov[6][1]*v[1] + _cov[1][12]*r[2][0] + _cov[1][13]*r[2][1] + _cov[1][14]*r[2][2] + _cov[1][15] * _dt;
    // _cov[2][3] = _cov[3][2] - (_cov[8][2]*v[1]) + _cov[7][2]*v[2] + _cov[2][12]*r[0][0] + _cov[2][13]*r[0][1] + _cov[2][14]*r[0][2];
    // _cov[2][4] = _cov[4][2] + _cov[8][2]*v[0] - _cov[6][2]*v[2] + _cov[2][12]*r[1][0] + _cov[2][13]*r[1][1] + _cov[2][14]*r[1][2];
    // _cov[2][5] = _cov[5][2] -(_cov[7][2]*v[0]) + _cov[6][2]*v[1] + _cov[2][12]*r[2][0] + _cov[2][13]*r[2][1] + _cov[2][14]*r[2][2] + _cov[2][15] * _dt;

    // -(a - ba) * dt
    const Vector3f dv_corr = (_ba - a) * _dt;
    // -R * (a - ba) * dt
    const array<float, 3> dv = {
        _rot(0, 0) * dv_corr[0] + _rot(0, 1) * dv_corr[1] + _rot(0, 2) * dv_corr[2],
        _rot(1, 0) * dv_corr[0] + _rot(1, 1) * dv_corr[1] + _rot(1, 2) * dv_corr[2],
        _rot(2, 0) * dv_corr[0] + _rot(2, 1) * dv_corr[1] + _rot(2, 2) * dv_corr[2],
    };

    // -R * dt
    const array<array<float, 3>, 3> rdt {
        {{-_rot(0, 0) * _dt, -_rot(0, 1) * _dt, -_rot(0, 2) * _dt},
         {-_rot(1, 0) * _dt, -_rot(1, 1) * _dt, -_rot(1, 2) * _dt},
         {-_rot(2, 0) * _dt, -_rot(2, 1) * _dt, -_rot(2, 2) * _dt}}
    };

    // Equations for covariance matrix prediction, without process noise!
    array<float, 90> var;

    var[0] = _cov[3][0] + _cov[3][3]*_dt;
    var[1] = _cov[4][3]*_dt;
    var[2] = _cov[4][0] + var[1];
    var[3] = _cov[5][3]*_dt;
    var[4] = _cov[5][0] + var[3];
    var[5] = _cov[8][0] + _cov[8][3]*_dt;
    var[6] = _cov[12][0] + _cov[12][3]*_dt;
    var[7] = _cov[13][0] + _cov[13][3]*_dt;
    var[8] = _cov[14][0] + _cov[14][3]*_dt;
    var[9] = _cov[7][0] + _cov[7][3]*_dt;
    var[10] = _cov[6][0] + _cov[6][3]*_dt;
    var[11] = _cov[15][0] + _cov[15][3]*_dt;
    var[12] = _cov[9][0] + _cov[9][3]*_dt;
    var[13] = _cov[10][0] + _cov[10][3]*_dt;
    var[14] = _cov[11][0] + _cov[11][3]*_dt;
    var[15] = _cov[4][1] + _cov[4][4]*_dt;
    var[16] = _cov[5][4]*_dt;
    var[17] = _cov[5][1] + var[16];
    var[18] = _cov[8][1] + _cov[8][4]*_dt;
    var[19] = _cov[12][1] + _cov[12][4]*_dt;
    var[20] = _cov[13][1] + _cov[13][4]*_dt;
    var[21] = _cov[14][1] + _cov[14][4]*_dt;
    var[22] = _cov[7][1] + _cov[7][4]*_dt;
    var[23] = _cov[6][1] + _cov[6][4]*_dt;
    var[24] = _cov[15][1] + _cov[15][4]*_dt;
    var[25] = _cov[9][1] + _cov[9][4]*_dt;
    var[26] = _cov[10][1] + _cov[10][4]*_dt;
    var[27] = _cov[11][1] + _cov[11][4]*_dt;
    var[28] = _cov[5][2] + _cov[5][5]*_dt;
    var[29] = _cov[8][2] + _cov[8][5]*_dt;
    var[30] = _cov[12][2] + _cov[12][5]*_dt;
    var[31] = _cov[13][2] + _cov[13][5]*_dt;
    var[32] = _cov[14][2] + _cov[14][5]*_dt;
    var[33] = _cov[7][2] + _cov[7][5]*_dt;
    var[34] = _cov[6][2] + _cov[6][5]*_dt;
    var[35] = _cov[15][5]*_dt;
    var[36] = _cov[15][2] + var[35];
    var[37] = _cov[9][2] + _cov[9][5]*_dt;
    var[38] = _cov[10][2] + _cov[10][5]*_dt;
    var[39] = _cov[11][2] + _cov[11][5]*_dt;
    var[40] = _cov[12][8]*rdt[0][0] + _cov[13][8]*rdt[0][1] + _cov[14][8]*rdt[0][2] + _cov[8][3] - _cov[8][7]*dv[2] + _cov[8][8]*dv[1];
    var[41] = _cov[12][12]*rdt[0][0] + _cov[12][3] - _cov[12][7]*dv[2] + _cov[12][8]*dv[1] + _cov[13][12]*rdt[0][1] + _cov[14][12]*rdt[0][2];
    var[42] = _cov[13][12]*rdt[0][0] + _cov[13][13]*rdt[0][1] + _cov[13][3] - _cov[13][7]*dv[2] + _cov[13][8]*dv[1] + _cov[14][13]*rdt[0][2];
    var[43] = _cov[14][12]*rdt[0][0] + _cov[14][13]*rdt[0][1] + _cov[14][14]*rdt[0][2] + _cov[14][3] - _cov[14][7]*dv[2] + _cov[14][8]*dv[1];
    var[44] = _cov[12][7]*rdt[0][0] + _cov[13][7]*rdt[0][1] + _cov[14][7]*rdt[0][2] + _cov[7][3] - _cov[7][7]*dv[2] + _cov[8][7]*dv[1];
    var[45] = _cov[8][6]*dv[1];
    var[46] = _cov[7][6]*dv[2];
    var[47] = _cov[12][6]*rdt[0][0] + _cov[13][6]*rdt[0][1] + _cov[14][6]*rdt[0][2] + _cov[6][3] + var[45] - var[46];
    var[48] = _cov[15][12]*rdt[0][0] + _cov[15][13]*rdt[0][1] + _cov[15][14]*rdt[0][2] + _cov[15][3] - _cov[15][7]*dv[2] + _cov[15][8]*dv[1];
    var[49] = _cov[12][9]*rdt[0][0];
    var[50] = _cov[13][9]*rdt[0][1] + _cov[14][9]*rdt[0][2] + _cov[9][3] - _cov[9][7]*dv[2] + _cov[9][8]*dv[1] + var[49];
    var[51] = _cov[13][10]*rdt[0][1];
    var[52] = _cov[10][3] - _cov[10][7]*dv[2] + _cov[10][8]*dv[1] + _cov[12][10]*rdt[0][0] + _cov[14][10]*rdt[0][2] + var[51];
    var[53] = _cov[14][11]*rdt[0][2];
    var[54] = _cov[11][3] - _cov[11][7]*dv[2] + _cov[11][8]*dv[1] + _cov[12][11]*rdt[0][0] + _cov[13][11]*rdt[0][1] + var[53];
    var[55] = _cov[12][6]*rdt[1][0] + _cov[13][6]*rdt[1][1] + _cov[14][6]*rdt[1][2] + _cov[6][4] + _cov[6][6]*dv[2] - _cov[8][6]*dv[0];
    var[56] = _cov[12][12]*rdt[1][0] + _cov[12][4] + _cov[12][6]*dv[2] - _cov[12][8]*dv[0] + _cov[13][12]*rdt[1][1] + _cov[14][12]*rdt[1][2];
    var[57] = _cov[13][12]*rdt[1][0] + _cov[13][13]*rdt[1][1] + _cov[13][4] + _cov[13][6]*dv[2] - _cov[13][8]*dv[0] + _cov[14][13]*rdt[1][2];
    var[58] = _cov[14][12]*rdt[1][0] + _cov[14][13]*rdt[1][1] + _cov[14][14]*rdt[1][2] + _cov[14][4] + _cov[14][6]*dv[2] - _cov[14][8]*dv[0];
    var[59] = _cov[12][8]*rdt[1][0] + _cov[13][8]*rdt[1][1] + _cov[14][8]*rdt[1][2] + _cov[8][4] + _cov[8][6]*dv[2] - _cov[8][8]*dv[0];
    var[60] = _cov[15][12]*rdt[1][0] + _cov[15][13]*rdt[1][1] + _cov[15][14]*rdt[1][2] + _cov[15][4] + _cov[15][6]*dv[2] - _cov[15][8]*dv[0];
    var[61] = _cov[8][7]*dv[0];
    var[62] = _cov[12][7]*rdt[1][0] + _cov[13][7]*rdt[1][1] + _cov[14][7]*rdt[1][2] + _cov[7][4] + var[46] - var[61];
    var[63] = _cov[12][9]*rdt[1][0];
    var[64] = _cov[13][9]*rdt[1][1] + _cov[14][9]*rdt[1][2] + _cov[9][4] + _cov[9][6]*dv[2] - _cov[9][8]*dv[0] + var[63];
    var[65] = _cov[13][10]*rdt[1][1];
    var[66] = _cov[10][4] + _cov[10][6]*dv[2] - _cov[10][8]*dv[0] + _cov[12][10]*rdt[1][0] + _cov[14][10]*rdt[1][2] + var[65];
    var[67] = _cov[14][11]*rdt[1][2];
    var[68] = _cov[11][4] + _cov[11][6]*dv[2] - _cov[11][8]*dv[0] + _cov[12][11]*rdt[1][0] + _cov[13][11]*rdt[1][1] + var[67];
    var[69] = _cov[15][12]*rdt[2][0] + _cov[15][13]*rdt[2][1] + _cov[15][14]*rdt[2][2] + _cov[15][15]*_dt + _cov[15][5] - _cov[15][6]*dv[1] + _cov[15][7]*dv[0];
    var[70] = _cov[12][7]*rdt[2][0] + _cov[13][7]*rdt[2][1] + _cov[14][7]*rdt[2][2] + _cov[15][7]*_dt + _cov[7][5] - _cov[7][6]*dv[1] + _cov[7][7]*dv[0];
    var[71] = _cov[12][12]*rdt[2][0] + _cov[12][5] - _cov[12][6]*dv[1] + _cov[12][7]*dv[0] + _cov[13][12]*rdt[2][1] + _cov[14][12]*rdt[2][2] + _cov[15][12]*_dt;
    var[72] = _cov[13][12]*rdt[2][0] + _cov[13][13]*rdt[2][1] + _cov[13][5] - _cov[13][6]*dv[1] + _cov[13][7]*dv[0] + _cov[14][13]*rdt[2][2] + _cov[15][13]*_dt;
    var[73] = _cov[14][12]*rdt[2][0] + _cov[14][13]*rdt[2][1] + _cov[14][14]*rdt[2][2] + _cov[14][5] - _cov[14][6]*dv[1] + _cov[14][7]*dv[0] + _cov[15][14]*_dt;
    var[74] = _cov[12][6]*rdt[2][0] + _cov[13][6]*rdt[2][1] + _cov[14][6]*rdt[2][2] + _cov[15][6]*_dt + _cov[6][5] - _cov[6][6]*dv[1] + _cov[7][6]*dv[0];
    var[75] = _cov[12][9]*rdt[2][0];
    var[76] = _cov[13][9]*rdt[2][1] + _cov[14][9]*rdt[2][2] + _cov[15][9]*_dt + _cov[9][5] - _cov[9][6]*dv[1] + _cov[9][7]*dv[0] + var[75];
    var[77] = _cov[13][10]*rdt[2][1];
    var[78] = _cov[10][5] - _cov[10][6]*dv[1] + _cov[10][7]*dv[0] + _cov[12][10]*rdt[2][0] + _cov[14][10]*rdt[2][2] + _cov[15][10]*_dt + var[77];
    var[79] = _cov[14][11]*rdt[2][2];
    var[80] = _cov[11][5] - _cov[11][6]*dv[1] + _cov[11][7]*dv[0] + _cov[12][11]*rdt[2][0] + _cov[13][11]*rdt[2][1] + _cov[15][11]*_dt + var[79];
    var[81] = _cov[10][9]*rdt[0][1] + _cov[11][9]*rdt[0][2] + _cov[9][6] + _cov[9][9]*rdt[0][0];
    var[82] = _cov[10][10]*rdt[0][1] + _cov[10][6] + _cov[10][9]*rdt[0][0] + _cov[11][10]*rdt[0][2];
    var[83] = _cov[11][10]*rdt[0][1] + _cov[11][11]*rdt[0][2] + _cov[11][6] + _cov[11][9]*rdt[0][0];
    var[84] = _cov[10][9]*rdt[1][1] + _cov[11][9]*rdt[1][2] + _cov[9][7] + _cov[9][9]*rdt[1][0];
    var[85] = _cov[10][10]*rdt[1][1] + _cov[10][7] + _cov[10][9]*rdt[1][0] + _cov[11][10]*rdt[1][2];
    var[86] = _cov[11][10]*rdt[1][1] + _cov[11][11]*rdt[1][2] + _cov[11][7] + _cov[11][9]*rdt[1][0];
    var[87] = _cov[10][9]*rdt[2][1] + _cov[11][9]*rdt[2][2] + _cov[9][8] + _cov[9][9]*rdt[2][0];
    var[88] = _cov[10][10]*rdt[2][1] + _cov[10][8] + _cov[10][9]*rdt[2][0] + _cov[11][10]*rdt[2][2];
    var[89] = _cov[11][10]*rdt[2][1] + _cov[11][11]*rdt[2][2] + _cov[11][8] + _cov[11][9]*rdt[2][0];


    _cov[0][0] = _cov[0][0] + _cov[3][0]*_dt + _dt*var[0];
    _cov[0][1] = _cov[1][0] + _cov[3][1]*_dt + _dt*var[2];
    _cov[1][1] = _cov[1][1] + _cov[4][1]*_dt + _dt*var[15];
    _cov[0][2] = _cov[2][0] + _cov[3][2]*_dt + _dt*var[4];
    _cov[1][2] = _cov[2][1] + _cov[4][2]*_dt + _dt*var[17];
    _cov[2][2] = _cov[2][2] + _cov[5][2]*_dt + _dt*var[28];
    _cov[0][3] = dv[1]*var[5] - dv[2]*var[9] + rdt[0][0]*var[6] + rdt[0][1]*var[7] + rdt[0][2]*var[8] + var[0];
    _cov[1][3] = _cov[3][1] + dv[1]*var[18] - dv[2]*var[22] + rdt[0][0]*var[19] + rdt[0][1]*var[20] + rdt[0][2]*var[21] + var[1];
    _cov[2][3] = _cov[3][2] + dv[1]*var[29] - dv[2]*var[33] + rdt[0][0]*var[30] + rdt[0][1]*var[31] + rdt[0][2]*var[32] + var[3];
    _cov[3][3] = _cov[12][3]*rdt[0][0] + _cov[13][3]*rdt[0][1] + _cov[14][3]*rdt[0][2] + _cov[3][3] - _cov[7][3]*dv[2] + _cov[8][3]*dv[1] + dv[1]*var[40] - dv[2]*var[44] + rdt[0][0]*var[41] + rdt[0][1]*var[42] + rdt[0][2]*var[43];
    _cov[0][4] = -dv[0]*var[5] + dv[2]*var[10] + rdt[1][0]*var[6] + rdt[1][1]*var[7] + rdt[1][2]*var[8] + var[2];
    _cov[1][4] = -dv[0]*var[18] + dv[2]*var[23] + rdt[1][0]*var[19] + rdt[1][1]*var[20] + rdt[1][2]*var[21] + var[15];
    _cov[2][4] = _cov[4][2] - dv[0]*var[29] + dv[2]*var[34] + rdt[1][0]*var[30] + rdt[1][1]*var[31] + rdt[1][2]*var[32] + var[16];
    _cov[3][4] = _cov[12][4]*rdt[0][0] + _cov[13][4]*rdt[0][1] + _cov[14][4]*rdt[0][2] + _cov[4][3] - _cov[7][4]*dv[2] + _cov[8][4]*dv[1] - dv[0]*var[40] + dv[2]*var[47] + rdt[1][0]*var[41] + rdt[1][1]*var[42] + rdt[1][2]*var[43];
    _cov[4][4] = _cov[12][4]*rdt[1][0] + _cov[13][4]*rdt[1][1] + _cov[14][4]*rdt[1][2] + _cov[4][4] + _cov[6][4]*dv[2] - _cov[8][4]*dv[0] - dv[0]*var[59] + dv[2]*var[55] + rdt[1][0]*var[56] + rdt[1][1]*var[57] + rdt[1][2]*var[58];
    _cov[0][5] = _dt*var[11] + dv[0]*var[9] - dv[1]*var[10] + rdt[2][0]*var[6] + rdt[2][1]*var[7] + rdt[2][2]*var[8] + var[4];
    _cov[1][5] = _dt*var[24] + dv[0]*var[22] - dv[1]*var[23] + rdt[2][0]*var[19] + rdt[2][1]*var[20] + rdt[2][2]*var[21] + var[17];
    _cov[2][5] = _dt*var[36] + dv[0]*var[33] - dv[1]*var[34] + rdt[2][0]*var[30] + rdt[2][1]*var[31] + rdt[2][2]*var[32] + var[28];
    _cov[3][5] = _cov[12][5]*rdt[0][0] + _cov[13][5]*rdt[0][1] + _cov[14][5]*rdt[0][2] + _cov[5][3] - _cov[7][5]*dv[2] + _cov[8][5]*dv[1] + _dt*var[48] + dv[0]*var[44] - dv[1]*var[47] + rdt[2][0]*var[41] + rdt[2][1]*var[42] + rdt[2][2]*var[43];
    _cov[4][5] = _cov[12][5]*rdt[1][0] + _cov[13][5]*rdt[1][1] + _cov[14][5]*rdt[1][2] + _cov[5][4] + _cov[6][5]*dv[2] - _cov[8][5]*dv[0] + _dt*var[60] + dv[0]*var[62] - dv[1]*var[55] + rdt[2][0]*var[56] + rdt[2][1]*var[57] + rdt[2][2]*var[58];
    _cov[5][5] = _cov[12][5]*rdt[2][0] + _cov[13][5]*rdt[2][1] + _cov[14][5]*rdt[2][2] + _cov[5][5] - _cov[6][5]*dv[1] + _cov[7][5]*dv[0] + _dt*var[69] + dv[0]*var[70] - dv[1]*var[74] + rdt[2][0]*var[71] + rdt[2][1]*var[72] + rdt[2][2]*var[73] + var[35];
    _cov[0][6] = rdt[0][0]*var[12] + rdt[0][1]*var[13] + rdt[0][2]*var[14] + var[10];
    _cov[1][6] = rdt[0][0]*var[25] + rdt[0][1]*var[26] + rdt[0][2]*var[27] + var[23];
    _cov[2][6] = rdt[0][0]*var[37] + rdt[0][1]*var[38] + rdt[0][2]*var[39] + var[34];
    _cov[3][6] = rdt[0][0]*var[50] + rdt[0][1]*var[52] + rdt[0][2]*var[54] + var[47];
    _cov[4][6] = rdt[0][0]*var[64] + rdt[0][1]*var[66] + rdt[0][2]*var[68] + var[55];
    _cov[5][6] = rdt[0][0]*var[76] + rdt[0][1]*var[78] + rdt[0][2]*var[80] + var[74];
    _cov[6][6] = _cov[10][6]*rdt[0][1] + _cov[11][6]*rdt[0][2] + _cov[6][6] + _cov[9][6]*rdt[0][0] + rdt[0][0]*var[81] + rdt[0][1]*var[82] + rdt[0][2]*var[83];
    _cov[0][7] = rdt[1][0]*var[12] + rdt[1][1]*var[13] + rdt[1][2]*var[14] + var[9];
    _cov[1][7] = rdt[1][0]*var[25] + rdt[1][1]*var[26] + rdt[1][2]*var[27] + var[22];
    _cov[2][7] = rdt[1][0]*var[37] + rdt[1][1]*var[38] + rdt[1][2]*var[39] + var[33];
    _cov[3][7] = rdt[1][0]*var[50] + rdt[1][1]*var[52] + rdt[1][2]*var[54] + var[44];
    _cov[4][7] = rdt[1][0]*var[64] + rdt[1][1]*var[66] + rdt[1][2]*var[68] + var[62];
    _cov[5][7] = rdt[1][0]*var[76] + rdt[1][1]*var[78] + rdt[1][2]*var[80] + var[70];
    _cov[6][7] = _cov[10][7]*rdt[0][1] + _cov[11][7]*rdt[0][2] + _cov[7][6] + _cov[9][7]*rdt[0][0] + rdt[1][0]*var[81] + rdt[1][1]*var[82] + rdt[1][2]*var[83];
    _cov[7][7] = _cov[10][7]*rdt[1][1] + _cov[11][7]*rdt[1][2] + _cov[7][7] + _cov[9][7]*rdt[1][0] + rdt[1][0]*var[84] + rdt[1][1]*var[85] + rdt[1][2]*var[86];
    _cov[0][8] = rdt[2][0]*var[12] + rdt[2][1]*var[13] + rdt[2][2]*var[14] + var[5];
    _cov[1][8] = rdt[2][0]*var[25] + rdt[2][1]*var[26] + rdt[2][2]*var[27] + var[18];
    _cov[2][8] = rdt[2][0]*var[37] + rdt[2][1]*var[38] + rdt[2][2]*var[39] + var[29];
    _cov[3][8] = rdt[2][0]*var[50] + rdt[2][1]*var[52] + rdt[2][2]*var[54] + var[40];
    _cov[4][8] = rdt[2][0]*var[64] + rdt[2][1]*var[66] + rdt[2][2]*var[68] + var[59];
    _cov[5][8] = _cov[12][8]*rdt[2][0] + _cov[13][8]*rdt[2][1] + _cov[14][8]*rdt[2][2] + _cov[15][8]*_dt + _cov[8][5] + rdt[2][0]*var[76] + rdt[2][1]*var[78] + rdt[2][2]*var[80] - var[45] + var[61];
    _cov[6][8] = _cov[10][8]*rdt[0][1] + _cov[11][8]*rdt[0][2] + _cov[8][6] + _cov[9][8]*rdt[0][0] + rdt[2][0]*var[81] + rdt[2][1]*var[82] + rdt[2][2]*var[83];
    _cov[7][8] = _cov[10][8]*rdt[1][1] + _cov[11][8]*rdt[1][2] + _cov[8][7] + _cov[9][8]*rdt[1][0] + rdt[2][0]*var[84] + rdt[2][1]*var[85] + rdt[2][2]*var[86];
    _cov[8][8] = _cov[10][8]*rdt[2][1] + _cov[11][8]*rdt[2][2] + _cov[8][8] + _cov[9][8]*rdt[2][0] + rdt[2][0]*var[87] + rdt[2][1]*var[88] + rdt[2][2]*var[89];
    _cov[0][9] = var[12];
    _cov[1][9] = var[25];
    _cov[2][9] = var[37];
    _cov[3][9] = var[50];
    _cov[4][9] = var[64];
    _cov[5][9] = var[76];
    _cov[6][9] = var[81];
    _cov[7][9] = var[84];
    _cov[8][9] = var[87];
    // _cov[9][9] = _cov[9][9];
    _cov[0][10] = var[13];
    _cov[1][10] = var[26];
    _cov[2][10] = var[38];
    _cov[3][10] = var[52];
    _cov[4][10] = var[66];
    _cov[5][10] = var[78];
    _cov[6][10] = var[82];
    _cov[7][10] = var[85];
    _cov[8][10] = var[88];
    // _cov[9][10] = _cov[10][9];
    // _cov[10][10] = _cov[10][10];
    _cov[0][11] = var[14];
    _cov[1][11] = var[27];
    _cov[2][11] = var[39];
    _cov[3][11] = var[54];
    _cov[4][11] = var[68];
    _cov[5][11] = var[80];
    _cov[6][11] = var[83];
    _cov[7][11] = var[86];
    _cov[8][11] = var[89];
    // _cov[9][11] = _cov[11][9];
    // _cov[10][11] = _cov[11][10];
    // _cov[11][11] = _cov[11][11];
    _cov[0][12] = var[6];
    _cov[1][12] = var[19];
    _cov[2][12] = var[30];
    _cov[3][12] = var[41];
    _cov[4][12] = var[56];
    _cov[5][12] = var[71];
    _cov[6][12] = _cov[12][10]*rdt[0][1] + _cov[12][11]*rdt[0][2] + _cov[12][6] + var[49];
    _cov[7][12] = _cov[12][10]*rdt[1][1] + _cov[12][11]*rdt[1][2] + _cov[12][7] + var[63];
    _cov[8][12] = _cov[12][10]*rdt[2][1] + _cov[12][11]*rdt[2][2] + _cov[12][8] + var[75];
    // _cov[9][12] = _cov[12][9];
    // _cov[10][12] = _cov[12][10];
    // _cov[11][12] = _cov[12][11];
    // _cov[12][12] = _cov[12][12];
    _cov[0][13] = var[7];
    _cov[1][13] = var[20];
    _cov[2][13] = var[31];
    _cov[3][13] = var[42];
    _cov[4][13] = var[57];
    _cov[5][13] = var[72];
    _cov[6][13] = _cov[13][11]*rdt[0][2] + _cov[13][6] + _cov[13][9]*rdt[0][0] + var[51];
    _cov[7][13] = _cov[13][11]*rdt[1][2] + _cov[13][7] + _cov[13][9]*rdt[1][0] + var[65];
    _cov[8][13] = _cov[13][11]*rdt[2][2] + _cov[13][8] + _cov[13][9]*rdt[2][0] + var[77];
    // _cov[9][13] = _cov[13][9];
    // _cov[10][13] = _cov[13][10];
    // _cov[11][13] = _cov[13][11];
    // _cov[12][13] = _cov[13][12];
    // _cov[13][13] = _cov[13][13];
    _cov[0][14] = var[8];
    _cov[1][14] = var[21];
    _cov[2][14] = var[32];
    _cov[3][14] = var[43];
    _cov[4][14] = var[58];
    _cov[5][14] = var[73];
    _cov[6][14] = _cov[14][10]*rdt[0][1] + _cov[14][6] + _cov[14][9]*rdt[0][0] + var[53];
    _cov[7][14] = _cov[14][10]*rdt[1][1] + _cov[14][7] + _cov[14][9]*rdt[1][0] + var[67];
    _cov[8][14] = _cov[14][10]*rdt[2][1] + _cov[14][8] + _cov[14][9]*rdt[2][0] + var[79];
    // _cov[9][14] = _cov[14][9];
    // _cov[10][14] = _cov[14][10];
    // _cov[11][14] = _cov[14][11];
    // _cov[12][14] = _cov[14][12];
    // _cov[13][14] = _cov[14][13];
    // _cov[14][14] = _cov[14][14];
    _cov[0][15] = var[11];
    _cov[1][15] = var[24];
    _cov[2][15] = var[36];
    _cov[3][15] = var[48];
    _cov[4][15] = var[60];
    _cov[5][15] = var[69];
    _cov[6][15] = _cov[15][10]*rdt[0][1] + _cov[15][11]*rdt[0][2] + _cov[15][6] + _cov[15][9]*rdt[0][0];
    _cov[7][15] = _cov[15][10]*rdt[1][1] + _cov[15][11]*rdt[1][2] + _cov[15][7] + _cov[15][9]*rdt[1][0];
    _cov[8][15] = _cov[15][10]*rdt[2][1] + _cov[15][11]*rdt[2][2] + _cov[15][8] + _cov[15][9]*rdt[2][0];
    // _cov[9][15] = _cov[15][9];
    // _cov[10][15] = _cov[15][10];
    // _cov[11][15] = _cov[15][11];
    // _cov[12][15] = _cov[15][12];
    // _cov[13][15] = _cov[15][13];
    // _cov[14][15] = _cov[15][14];
    // _cov[15][15] = _cov[15][15];


    // F *P *F' + Q
    _cov[0][0] = kahan_summation(_cov[0][0], _q_cov[0] * _dt2, _accumulator_cov[0]);
    _cov[1][1] = kahan_summation(_cov[1][1], _q_cov[1] * _dt2, _accumulator_cov[1]);
    _cov[2][2] = kahan_summation(_cov[2][2], _q_cov[2] * _dt2, _accumulator_cov[2]);
    if (_q_cov[3] == _q_cov[4]) {
        if (_q_cov[4] == _q_cov[5]) {
            _cov[3][3] = kahan_summation(_cov[3][3], _q_cov[3] * _dt2, _accumulator_cov[3]);
            _cov[4][4] = kahan_summation(_cov[4][4], _q_cov[4] * _dt2, _accumulator_cov[4]);
            _cov[5][5] = kahan_summation(_cov[5][5], _q_cov[5] * _dt2, _accumulator_cov[5]);
        } else {
            float sx = _q_cov[3] * _dt2;
            float sz = _q_cov[5] * _dt2;
            float r00_sx = _rot(0, 0) * sx;
            float r01_sx = _rot(0, 1) * sx;
            float r02_sz = _rot(0, 2) * sz;
            float r00_sx_r10_plus_r01_sx_r11_plus_r02_sz_r12 = r00_sx * _rot(1, 0) + r01_sx * _rot(1, 1) + r02_sz * _rot(1, 2);
            float r00_sx_r20_plus_r01_sx_r21_plus_r02_sz_r22 = r00_sx * _rot(2, 0) + r01_sx * _rot(2, 1) + r02_sz * _rot(2, 2);
            float r10_r20_sx_plus_r11_r21_sx_plus_r12_r22_sz = (_rot(1, 0) * _rot(2, 0) + _rot(1, 1) * _rot(2, 1)) * sx + _rot(1, 2) * _rot(2, 2) * sz;

            _cov[3][3] += kahan_summation(_cov[3][3], (_rot(0, 0) * _rot(0, 0) + _rot(0, 1) * _rot(0, 1)) * sx + _rot(0, 2) * _rot(0, 2) * sz, _accumulator_cov[3]);
            _cov[3][4] += r00_sx_r10_plus_r01_sx_r11_plus_r02_sz_r12;
            _cov[3][5] += r00_sx_r20_plus_r01_sx_r21_plus_r02_sz_r22;
            _cov[4][4] += kahan_summation(_cov[4][4], (_rot(1, 0) * _rot(1, 0) + _rot(1, 1) * _rot(1, 1)) * sx + _rot(1, 2) * _rot(1, 2) * sz, _accumulator_cov[4]);
            _cov[4][5] += r10_r20_sx_plus_r11_r21_sx_plus_r12_r22_sz;
            _cov[5][5] += kahan_summation(_cov[5][5], (_rot(2, 0) * _rot(2, 0) + _rot(2, 1) * _rot(2, 1)) * sx + _rot(2, 2) * _rot(2, 2) * sz, _accumulator_cov[5]);
        }
    } else {
        float sx = _q_cov[3] * _dt2;
        float sy = _q_cov[4] * _dt2;
        float sz = _q_cov[5] * _dt2;
        float r00_sx = _rot(0, 0) * sx;
        float r01_sy = _rot(0, 1) * sy;
        float r02_sz = _rot(0, 2) * sz;
        float r00_sx_r10_plus_r01_sy_r11_plus_r02_sz_r12 = r00_sx * _rot(1, 0) + r01_sy * _rot(1, 1) + r02_sz * _rot(1, 2);
        float r00_sx_r20_plus_r01_sy_r21_plus_r02_sz_r22 = r00_sx * _rot(2, 0) + r01_sy * _rot(2, 1) + r02_sz * _rot(2, 2);
        float r10_r20_sx_plus_r11_r21_sy_plus_r12_r22_sz = _rot(1, 0) * _rot(2, 0) * sx + _rot(1, 1) * _rot(2, 1) * sy + _rot(1, 2) * _rot(2, 2) * sz;
        
        _cov[3][3] += kahan_summation(_cov[3][3], _rot(0, 0) * _rot(0, 0) * sx + _rot(0, 1) * _rot(0, 1) * sy + _rot(0, 2) * _rot(0, 2) * sz, _accumulator_cov[3]);
        _cov[3][4] += r00_sx_r10_plus_r01_sy_r11_plus_r02_sz_r12;
        _cov[3][5] += r00_sx_r20_plus_r01_sy_r21_plus_r02_sz_r22;
        _cov[4][4] += kahan_summation(_cov[4][4], _rot(1, 0) * _rot(1, 0) * sx + _rot(1, 1) * _rot(1, 1) * sy + _rot(1, 2) * _rot(1, 2) * sz, _accumulator_cov[4]);
        _cov[4][5] += r10_r20_sx_plus_r11_r21_sy_plus_r12_r22_sz;
        _cov[5][5] += kahan_summation(_cov[5][5], _rot(2, 0) * _rot(2, 0) * sx + _rot(2, 1) * _rot(2, 1) * sy + _rot(2, 2) * _rot(2, 2) * sz, _accumulator_cov[5]);
    }
    if (_q_cov[6] == _q_cov[7]) {
        if (_q_cov[7] == _q_cov[8]) {
            _cov[6][6] = kahan_summation(_cov[6][6], _q_cov[6] * _dt2, _accumulator_cov[6]);
            _cov[7][7] = kahan_summation(_cov[7][7], _q_cov[7] * _dt2, _accumulator_cov[7]);
            _cov[8][8] = kahan_summation(_cov[8][8], _q_cov[8] * _dt2, _accumulator_cov[8]);
        } else {
            float sx = _q_cov[6] * _dt2;
            float sz = _q_cov[8] * _dt2;
            float r00_sx = _rot(0, 0) * sx;
            float r01_sx = _rot(0, 1) * sx;
            float r02_sz = _rot(0, 2) * sz;
            float r00_sx_r10_plus_r01_sx_r11_plus_r02_sz_r12 = r00_sx * _rot(1, 0) + r01_sx * _rot(1, 1) + r02_sz * _rot(1, 2);
            float r00_sx_r20_plus_r01_sx_r21_plus_r02_sz_r22 = r00_sx * _rot(2, 0) + r01_sx * _rot(2, 1) + r02_sz * _rot(2, 2);
            float r10_r20_sx_plus_r11_r21_sx_plus_r12_r22_sz = (_rot(1, 0) * _rot(2, 0) + _rot(1, 1) * _rot(2, 1)) * sx + _rot(1, 2) * _rot(2, 2) * sz;

            _cov[6][6] += kahan_summation(_cov[6][6], (_rot(0, 0) * _rot(0, 0) + _rot(0, 1) * _rot(0, 1)) * sx + _rot(0, 2) * _rot(0, 2) * sz, _accumulator_cov[6]);
            _cov[6][7] += r00_sx_r10_plus_r01_sx_r11_plus_r02_sz_r12;
            _cov[6][8] += r00_sx_r20_plus_r01_sx_r21_plus_r02_sz_r22;
            _cov[7][7] += kahan_summation(_cov[7][7], (_rot(1, 0) * _rot(1, 0) + _rot(1, 1) * _rot(1, 1)) * sx + _rot(1, 2) * _rot(1, 2) * sz, _accumulator_cov[7]);
            _cov[7][8] += r10_r20_sx_plus_r11_r21_sx_plus_r12_r22_sz;
            _cov[8][8] += kahan_summation(_cov[8][8], (_rot(2, 0) * _rot(2, 0) + _rot(2, 1) * _rot(2, 1)) * sx + _rot(2, 2) * _rot(2, 2) * sz, _accumulator_cov[8]);
        }
    } else {
        float sx = _q_cov[6] * _dt2;
        float sy = _q_cov[7] * _dt2;
        float sz = _q_cov[8] * _dt2;
        float r00_sx = _rot(0, 0) * sx;
        float r01_sy = _rot(0, 1) * sy;
        float r02_sz = _rot(0, 2) * sz;
        float r00_sx_r10_plus_r01_sy_r11_plus_r02_sz_r12 = r00_sx * _rot(1, 0) + r01_sy * _rot(1, 1) + r02_sz * _rot(1, 2);
        float r00_sx_r20_plus_r01_sy_r21_plus_r02_sz_r22 = r00_sx * _rot(2, 0) + r01_sy * _rot(2, 1) + r02_sz * _rot(2, 2);
        float r10_r20_sx_plus_r11_r21_sy_plus_r12_r22_sz = _rot(1, 0) * _rot(2, 0) * sx + _rot(1, 1) * _rot(2, 1) * sy + _rot(1, 2) * _rot(2, 2) * sz;
        
        _cov[6][6] += kahan_summation(_cov[6][6], _rot(0, 0) * _rot(0, 0) * sx + _rot(0, 1) * _rot(0, 1) * sy + _rot(0, 2) * _rot(0, 2) * sz, _accumulator_cov[6]);
        _cov[6][7] += r00_sx_r10_plus_r01_sy_r11_plus_r02_sz_r12;
        _cov[6][8] += r00_sx_r20_plus_r01_sy_r21_plus_r02_sz_r22;
        _cov[7][7] += kahan_summation(_cov[7][7], _rot(1, 0) * _rot(1, 0) * sx + _rot(1, 1) * _rot(1, 1) * sy + _rot(1, 2) * _rot(1, 2) * sz, _accumulator_cov[7]);
        _cov[7][8] += r10_r20_sx_plus_r11_r21_sy_plus_r12_r22_sz;
        _cov[8][8] += kahan_summation(_cov[8][8], _rot(2, 0) * _rot(2, 0) * sx + _rot(2, 1) * _rot(2, 1) * sy + _rot(2, 2) * _rot(2, 2) * sz, _accumulator_cov[8]);
    }
    _cov[9][9] = kahan_summation(_cov[9][9], _q_cov[9] * _dt2, _accumulator_cov[9]);
    _cov[10][10] = kahan_summation(_cov[10][10], _q_cov[10] * _dt2, _accumulator_cov[10]);
    _cov[11][11] = kahan_summation(_cov[11][11], _q_cov[11] * _dt2, _accumulator_cov[11]);
    _cov[12][12] = kahan_summation(_cov[12][12], _q_cov[12] * _dt2, _accumulator_cov[12]);
    _cov[13][13] = kahan_summation(_cov[13][13], _q_cov[13] * _dt2, _accumulator_cov[13]);
    _cov[14][14] = kahan_summation(_cov[14][14], _q_cov[14] * _dt2, _accumulator_cov[14]);
    _cov[15][15] = kahan_summation(_cov[15][15], _q_cov[15] * _dt2, _accumulator_cov[15]);

    regular_covariance_to_symmetric();
}

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
        const array<float, 16> HP = {_cov[0][dim] + _cov[0][6]*minus_rot_d_hat[dim][0] + _cov[0][7]*minus_rot_d_hat[dim][1] + _cov[0][8]*minus_rot_d_hat[dim][2],
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
                                     _cov[dim][15] + _cov[6][15]*minus_rot_d_hat[dim][0] + _cov[7][15]*minus_rot_d_hat[dim][1] + _cov[8][15]*minus_rot_d_hat[dim][2]};

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
    regular_covariance_to_symmetric();

    return info;
}

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
        const array<float, 16> HP = {_cov[0][index] + _cov[0][6]*rot_d_cross_w_corr_hat[dim][0] + _cov[0][7]*rot_d_cross_w_corr_hat[dim][1] + _cov[0][8]*rot_d_cross_w_corr_hat[dim][2] + _cov[0][9]*rot_d_hat[0] + _cov[0][10]*rot_d_hat[1] + _cov[0][11]*rot_d_hat[2],
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
                                     _cov[index][15] + _cov[6][15]*rot_d_cross_w_corr_hat[dim][0] + _cov[7][15]*rot_d_cross_w_corr_hat[dim][1] + _cov[8][15]*rot_d_cross_w_corr_hat[dim][2] + _cov[9][15]*rot_d_hat[0] + _cov[10][15]*rot_d_hat[1] + _cov[11][15]*rot_d_hat[2]};                                                  
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
    regular_covariance_to_symmetric();

    return info;
}

void GESKF::correct_state() {
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
        _p[i] += _error_state[i];
        _v[i] += _error_state[3 + i];
        _bg[i] += _error_state[9 + i];
        _ba[i] += _error_state[12 + i];
    }
    _g += _error_state[15];
    
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

void GESKF::correct_covariance() {
    array<float, 13> var;

    var[0] = 0.5F*_error_state[8];
    var[1] = 0.5F*_cov[8][0];
    var[2] = 0.5F*_error_state[7];
    var[3] = 0.5F*_error_state[6];
    var[4] = _cov[8][6] - _cov[8][7]*var[0] + _cov[8][8]*var[2];
    var[5] = _cov[7][6] - _cov[7][7]*var[0] + _cov[8][7]*var[2];
    var[6] = _cov[7][6]*var[0];
    var[7] = _cov[8][6]*var[2];
    var[8] = _cov[6][6] - var[6] + var[7];
    var[9] = _cov[6][6]*var[0] + _cov[7][6] - _cov[8][6]*var[3];
    var[10] = _cov[8][6]*var[0] + _cov[8][7] - _cov[8][8]*var[3];
    var[11] = _cov[8][7]*var[3];
    var[12] = _cov[7][7] - var[11] + var[6];

    // _cov[0][0] = _cov[0][0];
    // _cov[0][1] = _cov[1][0];
    // _cov[1][1] = _cov[1][1];
    // _cov[0][2] = _cov[2][0];
    // _cov[1][2] = _cov[2][1];
    // _cov[2][2] = _cov[2][2];
    // _cov[0][3] = _cov[3][0];
    // _cov[1][3] = _cov[3][1];
    // _cov[2][3] = _cov[3][2];
    // _cov[3][3] = _cov[3][3];
    // _cov[0][4] = _cov[4][0];
    // _cov[1][4] = _cov[4][1];
    // _cov[2][4] = _cov[4][2];
    // _cov[3][4] = _cov[4][3];
    // _cov[4][4] = _cov[4][4];
    // _cov[0][5] = _cov[5][0];
    // _cov[1][5] = _cov[5][1];
    // _cov[2][5] = _cov[5][2];
    // _cov[3][5] = _cov[5][3];
    // _cov[4][5] = _cov[5][4];
    // _cov[5][5] = _cov[5][5];
    _cov[0][6] = _cov[6][0] - _cov[7][0]*var[0] + _error_state[7]*var[1];
    _cov[1][6] = _cov[6][1] - _cov[7][1]*var[0] + _cov[8][1]*var[2];
    _cov[2][6] = _cov[6][2] - _cov[7][2]*var[0] + _cov[8][2]*var[2];
    _cov[3][6] = _cov[6][3] - _cov[7][3]*var[0] + _cov[8][3]*var[2];
    _cov[4][6] = _cov[6][4] - _cov[7][4]*var[0] + _cov[8][4]*var[2];
    _cov[5][6] = _cov[6][5] - _cov[7][5]*var[0] + _cov[8][5]*var[2];

    _cov[8][8] = _cov[8][8] + var[11] - var[2]*(-_cov[6][6]*var[2] + _cov[7][6]*var[3] + _cov[8][6]) + var[3]*(-_cov[7][6]*var[2] + _cov[7][7]*var[3] + _cov[8][7]) - var[7];
    _cov[6][6] = -var[0]*var[5] + var[2]*var[4] + var[8];
    _cov[0][7] = _cov[6][0]*var[0] + _cov[7][0] - _error_state[6]*var[1];
    _cov[1][7] = _cov[6][1]*var[0] + _cov[7][1] - _cov[8][1]*var[3];
    _cov[2][7] = _cov[6][2]*var[0] + _cov[7][2] - _cov[8][2]*var[3];
    _cov[3][7] = _cov[6][3]*var[0] + _cov[7][3] - _cov[8][3]*var[3];
    _cov[4][7] = _cov[6][4]*var[0] + _cov[7][4] - _cov[8][4]*var[3];
    _cov[5][7] = _cov[6][5]*var[0] + _cov[7][5] - _cov[8][5]*var[3];
    _cov[6][7] = var[0]*var[8] - var[3]*var[4] + var[5];
    _cov[7][7] = var[0]*var[9] - var[10]*var[3] + var[12];
    _cov[0][8] = -_cov[6][0]*var[2] + _cov[7][0]*var[3] + _cov[8][0];
    _cov[1][8] = -_cov[6][1]*var[2] + _cov[7][1]*var[3] + _cov[8][1];
    _cov[2][8] = -_cov[6][2]*var[2] + _cov[7][2]*var[3] + _cov[8][2];
    _cov[3][8] = -_cov[6][3]*var[2] + _cov[7][3]*var[3] + _cov[8][3];
    _cov[4][8] = -_cov[6][4]*var[2] + _cov[7][4]*var[3] + _cov[8][4];
    _cov[5][8] = -_cov[6][5]*var[2] + _cov[7][5]*var[3] + _cov[8][5];
    _cov[6][8] = -var[2]*var[8] + var[3]*var[5] + var[4];
    _cov[7][8] = var[10] + var[12]*var[3] - var[2]*var[9];
    
    // _cov[0][9] = _cov[9][0];
    // _cov[1][9] = _cov[9][1];
    // _cov[2][9] = _cov[9][2];
    // _cov[3][9] = _cov[9][3];
    // _cov[4][9] = _cov[9][4];
    // _cov[5][9] = _cov[9][5];
    _cov[6][9] = _cov[9][6] - _cov[9][7]*var[0] + _cov[9][8]*var[2];
    _cov[7][9] = _cov[9][6]*var[0] + _cov[9][7] - _cov[9][8]*var[3];
    _cov[8][9] = -_cov[9][6]*var[2] + _cov[9][7]*var[3] + _cov[9][8];
    // _cov[9][9] = _cov[9][9];
    // _cov[0][10] = _cov[10][0];
    // _cov[1][10] = _cov[10][1];
    // _cov[2][10] = _cov[10][2];
    // _cov[3][10] = _cov[10][3];
    // _cov[4][10] = _cov[10][4];
    // _cov[5][10] = _cov[10][5];
    _cov[6][10] = _cov[10][6] - _cov[10][7]*var[0] + _cov[10][8]*var[2];
    _cov[7][10] = _cov[10][6]*var[0] + _cov[10][7] - _cov[10][8]*var[3];
    _cov[8][10] = -_cov[10][6]*var[2] + _cov[10][7]*var[3] + _cov[10][8];
    // _cov[9][10] = _cov[10][9];
    // _cov[10][10] = _cov[10][10];
    // _cov[0][11] = _cov[11][0];
    // _cov[1][11] = _cov[11][1];
    // _cov[2][11] = _cov[11][2];
    // _cov[3][11] = _cov[11][3];
    // _cov[4][11] = _cov[11][4];
    // _cov[5][11] = _cov[11][5];
    _cov[6][11] = _cov[11][6] - _cov[11][7]*var[0] + _cov[11][8]*var[2];
    _cov[7][11] = _cov[11][6]*var[0] + _cov[11][7] - _cov[11][8]*var[3];
    _cov[8][11] = -_cov[11][6]*var[2] + _cov[11][7]*var[3] + _cov[11][8];
    // _cov[9][11] = _cov[11][9];
    // _cov[10][11] = _cov[11][10];
    // _cov[11][11] = _cov[11][11];
    // _cov[0][12] = _cov[12][0];
    // _cov[1][12] = _cov[12][1];
    // _cov[2][12] = _cov[12][2];
    // _cov[3][12] = _cov[12][3];
    // _cov[4][12] = _cov[12][4];
    // _cov[5][12] = _cov[12][5];
    _cov[6][12] = _cov[12][6] - _cov[12][7]*var[0] + _cov[12][8]*var[2];
    _cov[7][12] = _cov[12][6]*var[0] + _cov[12][7] - _cov[12][8]*var[3];
    _cov[8][12] = -_cov[12][6]*var[2] + _cov[12][7]*var[3] + _cov[12][8];
    // _cov[9][12] = _cov[12][9];
    // _cov[10][12] = _cov[12][10];
    // _cov[11][12] = _cov[12][11];
    // _cov[12][12] = _cov[12][12];
    // _cov[0][13] = _cov[13][0];
    // _cov[1][13] = _cov[13][1];
    // _cov[2][13] = _cov[13][2];
    // _cov[3][13] = _cov[13][3];
    // _cov[4][13] = _cov[13][4];
    // _cov[5][13] = _cov[13][5];
    _cov[6][13] = _cov[13][6] - _cov[13][7]*var[0] + _cov[13][8]*var[2];
    _cov[7][13] = _cov[13][6]*var[0] + _cov[13][7] - _cov[13][8]*var[3];
    _cov[8][13] = -_cov[13][6]*var[2] + _cov[13][7]*var[3] + _cov[13][8];
    // _cov[9][13] = _cov[13][9];
    // _cov[10][13] = _cov[13][10];
    // _cov[11][13] = _cov[13][11];
    // _cov[12][13] = _cov[13][12];
    // _cov[13][13] = _cov[13][13];
    // _cov[0][14] = _cov[14][0];
    // _cov[1][14] = _cov[14][1];
    // _cov[2][14] = _cov[14][2];
    // _cov[3][14] = _cov[14][3];
    // _cov[4][14] = _cov[14][4];
    // _cov[5][14] = _cov[14][5];
    _cov[6][14] = _cov[14][6] - _cov[14][7]*var[0] + _cov[14][8]*var[2];
    _cov[7][14] = _cov[14][6]*var[0] + _cov[14][7] - _cov[14][8]*var[3];
    _cov[8][14] = -_cov[14][6]*var[2] + _cov[14][7]*var[3] + _cov[14][8];
    // _cov[9][14] = _cov[14][9];
    // _cov[10][14] = _cov[14][10];
    // _cov[11][14] = _cov[14][11];
    // _cov[12][14] = _cov[14][12];
    // _cov[13][14] = _cov[14][13];
    // _cov[14][14] = _cov[14][14];
    // _cov[0][15] = _cov[15][0];
    // _cov[1][15] = _cov[15][1];
    // _cov[2][15] = _cov[15][2];
    // _cov[3][15] = _cov[15][3];
    // _cov[4][15] = _cov[15][4];
    // _cov[5][15] = _cov[15][5];
    _cov[6][15] = _cov[15][6] - _cov[15][7]*var[0] + _cov[15][8]*var[2];
    _cov[7][15] = _cov[15][6]*var[0] + _cov[15][7] - _cov[15][8]*var[3];
    _cov[8][15] = -_cov[15][6]*var[2] + _cov[15][7]*var[3] + _cov[15][8];
    // _cov[9][15] = _cov[15][9];
    // _cov[10][15] = _cov[15][10];
    // _cov[11][15] = _cov[15][11];
    // _cov[12][15] = _cov[15][12];
    // _cov[13][15] = _cov[15][13];
    // _cov[14][15] = _cov[15][14];
    // _cov[15][15] = _cov[15][15];

    regular_covariance_to_symmetric(6, 9);
}