//
// Created by Cain on 2022/11/19.
//

#include "riekf.h"
#include <cfloat>
#include <iostream>

using namespace std;
using namespace riekf;

void RIEKF::predict_covariance(const Vector3f &w, const Vector3f &a) {
    // g * dt
    const float gdt = _g * _dt;

    // -R * dt
    const array<array<float, 3>, 3> rdt {
        {{-_rot(0, 0) * _dt, -_rot(0, 1) * _dt, -_rot(0, 2) * _dt},
         {-_rot(1, 0) * _dt, -_rot(1, 1) * _dt, -_rot(1, 2) * _dt},
         {-_rot(2, 0) * _dt, -_rot(2, 1) * _dt, -_rot(2, 2) * _dt}}
    };

    // -p^ * R * dt
    const array<array<float, 3>, 3> prdt {
        {{_p[1] * rdt[2][0] - _p[2] * rdt[1][0], _p[1] * rdt[2][1] - _p[2] * rdt[1][1], _p[1] * rdt[2][2] - _p[2] * rdt[1][2]},
         {_p[2] * rdt[0][0] - _p[0] * rdt[2][0], _p[2] * rdt[0][1] - _p[0] * rdt[2][1], _p[2] * rdt[0][2] - _p[0] * rdt[2][2]},
         {_p[0] * rdt[1][0] - _p[1] * rdt[0][0], _p[0] * rdt[1][1] - _p[1] * rdt[0][1], _p[0] * rdt[1][2] - _p[1] * rdt[0][2]}}
    };

    // -v^ * R * dt
    const array<array<float, 3>, 3> vrdt {
        {{_v[1] * rdt[2][0] - _v[2] * rdt[1][0], _v[1] * rdt[2][1] - _v[2] * rdt[1][1], _v[1] * rdt[2][2] - _v[2] * rdt[1][2]},
         {_v[2] * rdt[0][0] - _v[0] * rdt[2][0], _v[2] * rdt[0][1] - _v[0] * rdt[2][1], _v[2] * rdt[0][2] - _v[0] * rdt[2][2]},
         {_v[0] * rdt[1][0] - _v[1] * rdt[0][0], _v[0] * rdt[1][1] - _v[1] * rdt[0][1], _v[0] * rdt[1][2] - _v[1] * rdt[0][2]}}
    };

    // Equations for covariance matrix prediction, without process noise!
    array<float, 116> var;

    var[0] = _cov[10][3]*prdt[0][1] + _cov[11][3]*prdt[0][2] + _cov[3][0] + _cov[3][3]*_dt + _cov[9][3]*prdt[0][0];
    var[1] = _cov[10][9]*prdt[0][1] + _cov[11][9]*prdt[0][2] + _cov[9][0] + _cov[9][3]*_dt + _cov[9][9]*prdt[0][0];
    var[2] = _cov[10][0] + _cov[10][10]*prdt[0][1] + _cov[10][3]*_dt + _cov[10][9]*prdt[0][0] + _cov[11][10]*prdt[0][2];
    var[3] = _cov[11][0] + _cov[11][10]*prdt[0][1] + _cov[11][11]*prdt[0][2] + _cov[11][3]*_dt + _cov[11][9]*prdt[0][0];
    var[4] = _q_cov[6]*prdt[0][0];
    var[5] = _q_cov[7]*prdt[0][1];
    var[6] = _q_cov[8]*prdt[0][2];
    var[7] = _cov[4][3]*_dt;
    var[8] = _cov[10][4]*prdt[0][1] + _cov[11][4]*prdt[0][2] + _cov[4][0] + _cov[9][4]*prdt[0][0] + var[7];
    var[9] = _cov[5][3]*_dt;
    var[10] = _cov[10][5]*prdt[0][1] + _cov[11][5]*prdt[0][2] + _cov[5][0] + _cov[9][5]*prdt[0][0] + var[9];
    var[11] = _cov[12][0] + _cov[12][10]*prdt[0][1] + _cov[12][11]*prdt[0][2] + _cov[12][3]*_dt + _cov[12][9]*prdt[0][0];
    var[12] = _cov[13][0] + _cov[13][10]*prdt[0][1] + _cov[13][11]*prdt[0][2] + _cov[13][3]*_dt + _cov[13][9]*prdt[0][0];
    var[13] = _cov[14][0] + _cov[14][10]*prdt[0][1] + _cov[14][11]*prdt[0][2] + _cov[14][3]*_dt + _cov[14][9]*prdt[0][0];
    var[14] = _cov[10][7]*prdt[0][1] + _cov[11][7]*prdt[0][2] + _cov[7][0] + _cov[7][3]*_dt + _cov[9][7]*prdt[0][0];
    var[15] = _cov[10][6]*prdt[0][1] + _cov[11][6]*prdt[0][2] + _cov[6][0] + _cov[6][3]*_dt + _cov[9][6]*prdt[0][0];
    var[16] = _cov[15][0] + _cov[15][10]*prdt[0][1] + _cov[15][11]*prdt[0][2] + _cov[15][3]*_dt + _cov[15][9]*prdt[0][0];
    var[17] = _cov[10][4]*prdt[1][1] + _cov[11][4]*prdt[1][2] + _cov[4][1] + _cov[4][4]*_dt + _cov[9][4]*prdt[1][0];
    var[18] = _cov[10][9]*prdt[1][1] + _cov[11][9]*prdt[1][2] + _cov[9][1] + _cov[9][4]*_dt + _cov[9][9]*prdt[1][0];
    var[19] = _cov[10][10]*prdt[1][1] + _cov[10][1] + _cov[10][4]*_dt + _cov[10][9]*prdt[1][0] + _cov[11][10]*prdt[1][2];
    var[20] = _cov[11][10]*prdt[1][1] + _cov[11][11]*prdt[1][2] + _cov[11][1] + _cov[11][4]*_dt + _cov[11][9]*prdt[1][0];
    var[21] = _q_cov[6]*prdt[1][0];
    var[22] = _q_cov[7]*prdt[1][1];
    var[23] = _q_cov[8]*prdt[1][2];
    var[24] = _cov[5][4]*_dt;
    var[25] = _cov[10][5]*prdt[1][1] + _cov[11][5]*prdt[1][2] + _cov[5][1] + _cov[9][5]*prdt[1][0] + var[24];
    var[26] = _cov[12][10]*prdt[1][1] + _cov[12][11]*prdt[1][2] + _cov[12][1] + _cov[12][4]*_dt + _cov[12][9]*prdt[1][0];
    var[27] = _cov[13][10]*prdt[1][1] + _cov[13][11]*prdt[1][2] + _cov[13][1] + _cov[13][4]*_dt + _cov[13][9]*prdt[1][0];
    var[28] = _cov[14][10]*prdt[1][1] + _cov[14][11]*prdt[1][2] + _cov[14][1] + _cov[14][4]*_dt + _cov[14][9]*prdt[1][0];
    var[29] = _cov[10][7]*prdt[1][1] + _cov[11][7]*prdt[1][2] + _cov[7][1] + _cov[7][4]*_dt + _cov[9][7]*prdt[1][0];
    var[30] = _cov[10][6]*prdt[1][1] + _cov[11][6]*prdt[1][2] + _cov[6][1] + _cov[6][4]*_dt + _cov[9][6]*prdt[1][0];
    var[31] = _cov[15][10]*prdt[1][1] + _cov[15][11]*prdt[1][2] + _cov[15][1] + _cov[15][4]*_dt + _cov[15][9]*prdt[1][0];
    var[32] = _cov[10][5]*prdt[2][1] + _cov[11][5]*prdt[2][2] + _cov[5][2] + _cov[5][5]*_dt + _cov[9][5]*prdt[2][0];
    var[33] = _cov[10][9]*prdt[2][1] + _cov[11][9]*prdt[2][2] + _cov[9][2] + _cov[9][5]*_dt + _cov[9][9]*prdt[2][0];
    var[34] = _cov[10][10]*prdt[2][1] + _cov[10][2] + _cov[10][5]*_dt + _cov[10][9]*prdt[2][0] + _cov[11][10]*prdt[2][2];
    var[35] = _cov[11][10]*prdt[2][1] + _cov[11][11]*prdt[2][2] + _cov[11][2] + _cov[11][5]*_dt + _cov[11][9]*prdt[2][0];
    var[36] = _q_cov[6]*prdt[2][0];
    var[37] = _q_cov[7]*prdt[2][1];
    var[38] = _q_cov[8]*prdt[2][2];
    var[39] = _cov[12][10]*prdt[2][1] + _cov[12][11]*prdt[2][2] + _cov[12][2] + _cov[12][5]*_dt + _cov[12][9]*prdt[2][0];
    var[40] = _cov[13][10]*prdt[2][1] + _cov[13][11]*prdt[2][2] + _cov[13][2] + _cov[13][5]*_dt + _cov[13][9]*prdt[2][0];
    var[41] = _cov[14][10]*prdt[2][1] + _cov[14][11]*prdt[2][2] + _cov[14][2] + _cov[14][5]*_dt + _cov[14][9]*prdt[2][0];
    var[42] = _cov[10][7]*prdt[2][1] + _cov[11][7]*prdt[2][2] + _cov[7][2] + _cov[7][5]*_dt + _cov[9][7]*prdt[2][0];
    var[43] = _cov[10][6]*prdt[2][1] + _cov[11][6]*prdt[2][2] + _cov[6][2] + _cov[6][5]*_dt + _cov[9][6]*prdt[2][0];
    var[44] = _cov[15][5]*_dt;
    var[45] = _cov[15][10]*prdt[2][1] + _cov[15][11]*prdt[2][2] + _cov[15][2] + _cov[15][9]*prdt[2][0] + var[44];
    var[46] = (rdt[0][0])*(rdt[0][0]);
    var[47] = (rdt[0][1])*(rdt[0][1]);
    var[48] = (rdt[0][2])*(rdt[0][2]);
    var[49] = _cov[12][10]*vrdt[0][1] + _cov[12][11]*vrdt[0][2] + _cov[12][12]*rdt[0][0] + _cov[12][3] - _cov[12][7]*gdt + _cov[12][9]*vrdt[0][0] + _cov[13][12]*rdt[0][1] + _cov[14][12]*rdt[0][2];
    var[50] = _cov[13][10]*vrdt[0][1] + _cov[13][11]*vrdt[0][2] + _cov[13][12]*rdt[0][0] + _cov[13][13]*rdt[0][1] + _cov[13][3] - _cov[13][7]*gdt + _cov[13][9]*vrdt[0][0] + _cov[14][13]*rdt[0][2];
    var[51] = _cov[14][10]*vrdt[0][1] + _cov[14][11]*vrdt[0][2] + _cov[14][12]*rdt[0][0] + _cov[14][13]*rdt[0][1] + _cov[14][14]*rdt[0][2] + _cov[14][3] - _cov[14][7]*gdt + _cov[14][9]*vrdt[0][0];
    var[52] = _cov[12][9]*rdt[0][0];
    var[53] = _cov[10][9]*vrdt[0][1] + _cov[11][9]*vrdt[0][2] + _cov[13][9]*rdt[0][1] + _cov[14][9]*rdt[0][2] + _cov[9][3] - _cov[9][7]*gdt + _cov[9][9]*vrdt[0][0] + var[52];
    var[54] = _cov[13][10]*rdt[0][1];
    var[55] = _cov[10][10]*vrdt[0][1] + _cov[10][3] - _cov[10][7]*gdt + _cov[10][9]*vrdt[0][0] + _cov[11][10]*vrdt[0][2] + _cov[12][10]*rdt[0][0] + _cov[14][10]*rdt[0][2] + var[54];
    var[56] = _cov[14][11]*rdt[0][2];
    var[57] = _cov[11][10]*vrdt[0][1] + _cov[11][11]*vrdt[0][2] + _cov[11][3] - _cov[11][7]*gdt + _cov[11][9]*vrdt[0][0] + _cov[12][11]*rdt[0][0] + _cov[13][11]*rdt[0][1] + var[56];
    var[58] = _cov[10][7]*vrdt[0][1] + _cov[11][7]*vrdt[0][2] + _cov[12][7]*rdt[0][0] + _cov[13][7]*rdt[0][1] + _cov[14][7]*rdt[0][2] + _cov[7][3] - _cov[7][7]*gdt + _cov[9][7]*vrdt[0][0];
    var[59] = _q_cov[3]*rdt[0][0];
    var[60] = _q_cov[4]*rdt[0][1];
    var[61] = _q_cov[5]*rdt[0][2];
    var[62] = _q_cov[6]*vrdt[0][0];
    var[63] = _q_cov[7]*vrdt[0][1];
    var[64] = _q_cov[8]*vrdt[0][2];
    var[65] = _cov[7][6]*gdt;
    var[66] = _cov[10][6]*vrdt[0][1] + _cov[11][6]*vrdt[0][2] + _cov[12][6]*rdt[0][0] + _cov[13][6]*rdt[0][1] + _cov[14][6]*rdt[0][2] + _cov[6][3] + _cov[9][6]*vrdt[0][0] - var[65];
    var[67] = _cov[15][10]*vrdt[0][1] + _cov[15][11]*vrdt[0][2] + _cov[15][12]*rdt[0][0] + _cov[15][13]*rdt[0][1] + _cov[15][14]*rdt[0][2] + _cov[15][3] - _cov[15][7]*gdt + _cov[15][9]*vrdt[0][0];
    var[68] = (rdt[1][0])*(rdt[1][0]);
    var[69] = (rdt[1][1])*(rdt[1][1]);
    var[70] = (rdt[1][2])*(rdt[1][2]);
    var[71] = _cov[10][6]*vrdt[1][1] + _cov[11][6]*vrdt[1][2] + _cov[12][6]*rdt[1][0] + _cov[13][6]*rdt[1][1] + _cov[14][6]*rdt[1][2] + _cov[6][4] + _cov[6][6]*gdt + _cov[9][6]*vrdt[1][0];
    var[72] = _cov[12][10]*vrdt[1][1] + _cov[12][11]*vrdt[1][2] + _cov[12][12]*rdt[1][0] + _cov[12][4] + _cov[12][6]*gdt + _cov[12][9]*vrdt[1][0] + _cov[13][12]*rdt[1][1] + _cov[14][12]*rdt[1][2];
    var[73] = _cov[13][10]*vrdt[1][1] + _cov[13][11]*vrdt[1][2] + _cov[13][12]*rdt[1][0] + _cov[13][13]*rdt[1][1] + _cov[13][4] + _cov[13][6]*gdt + _cov[13][9]*vrdt[1][0] + _cov[14][13]*rdt[1][2];
    var[74] = _cov[14][10]*vrdt[1][1] + _cov[14][11]*vrdt[1][2] + _cov[14][12]*rdt[1][0] + _cov[14][13]*rdt[1][1] + _cov[14][14]*rdt[1][2] + _cov[14][4] + _cov[14][6]*gdt + _cov[14][9]*vrdt[1][0];
    var[75] = _cov[12][9]*rdt[1][0];
    var[76] = _cov[10][9]*vrdt[1][1] + _cov[11][9]*vrdt[1][2] + _cov[13][9]*rdt[1][1] + _cov[14][9]*rdt[1][2] + _cov[9][4] + _cov[9][6]*gdt + _cov[9][9]*vrdt[1][0] + var[75];
    var[77] = _cov[13][10]*rdt[1][1];
    var[78] = _cov[10][10]*vrdt[1][1] + _cov[10][4] + _cov[10][6]*gdt + _cov[10][9]*vrdt[1][0] + _cov[11][10]*vrdt[1][2] + _cov[12][10]*rdt[1][0] + _cov[14][10]*rdt[1][2] + var[77];
    var[79] = _cov[14][11]*rdt[1][2];
    var[80] = _cov[11][10]*vrdt[1][1] + _cov[11][11]*vrdt[1][2] + _cov[11][4] + _cov[11][6]*gdt + _cov[11][9]*vrdt[1][0] + _cov[12][11]*rdt[1][0] + _cov[13][11]*rdt[1][1] + var[79];
    var[81] = rdt[1][0]*rdt[2][0];
    var[82] = rdt[1][1]*rdt[2][1];
    var[83] = rdt[1][2]*rdt[2][2];
    var[84] = _q_cov[6]*vrdt[1][0];
    var[85] = _q_cov[7]*vrdt[1][1];
    var[86] = _q_cov[8]*vrdt[1][2];
    var[87] = _cov[15][10]*vrdt[1][1] + _cov[15][11]*vrdt[1][2] + _cov[15][12]*rdt[1][0] + _cov[15][13]*rdt[1][1] + _cov[15][14]*rdt[1][2] + _cov[15][4] + _cov[15][6]*gdt + _cov[15][9]*vrdt[1][0];
    var[88] = (rdt[2][0])*(rdt[2][0]);
    var[89] = (rdt[2][1])*(rdt[2][1]);
    var[90] = (rdt[2][2])*(rdt[2][2]);
    var[91] = _cov[15][10]*vrdt[2][1] + _cov[15][11]*vrdt[2][2] + _cov[15][12]*rdt[2][0] + _cov[15][13]*rdt[2][1] + _cov[15][14]*rdt[2][2] + _cov[15][15]*_dt + _cov[15][5] + _cov[15][9]*vrdt[2][0];
    var[92] = _cov[12][10]*vrdt[2][1] + _cov[12][11]*vrdt[2][2] + _cov[12][12]*rdt[2][0] + _cov[12][5] + _cov[12][9]*vrdt[2][0] + _cov[13][12]*rdt[2][1] + _cov[14][12]*rdt[2][2] + _cov[15][12]*_dt;
    var[93] = _cov[13][10]*vrdt[2][1] + _cov[13][11]*vrdt[2][2] + _cov[13][12]*rdt[2][0] + _cov[13][13]*rdt[2][1] + _cov[13][5] + _cov[13][9]*vrdt[2][0] + _cov[14][13]*rdt[2][2] + _cov[15][13]*_dt;
    var[94] = _cov[14][10]*vrdt[2][1] + _cov[14][11]*vrdt[2][2] + _cov[14][12]*rdt[2][0] + _cov[14][13]*rdt[2][1] + _cov[14][14]*rdt[2][2] + _cov[14][5] + _cov[14][9]*vrdt[2][0] + _cov[15][14]*_dt;
    var[95] = _cov[12][9]*rdt[2][0];
    var[96] = _cov[10][9]*vrdt[2][1] + _cov[11][9]*vrdt[2][2] + _cov[13][9]*rdt[2][1] + _cov[14][9]*rdt[2][2] + _cov[15][9]*_dt + _cov[9][5] + _cov[9][9]*vrdt[2][0] + var[95];
    var[97] = _cov[13][10]*rdt[2][1];
    var[98] = _cov[10][10]*vrdt[2][1] + _cov[10][5] + _cov[10][9]*vrdt[2][0] + _cov[11][10]*vrdt[2][2] + _cov[12][10]*rdt[2][0] + _cov[14][10]*rdt[2][2] + _cov[15][10]*_dt + var[97];
    var[99] = _cov[14][11]*rdt[2][2];
    var[100] = _cov[11][10]*vrdt[2][1] + _cov[11][11]*vrdt[2][2] + _cov[11][5] + _cov[11][9]*vrdt[2][0] + _cov[12][11]*rdt[2][0] + _cov[13][11]*rdt[2][1] + _cov[15][11]*_dt + var[99];
    var[101] = _q_cov[6]*vrdt[2][0];
    var[102] = _q_cov[7]*vrdt[2][1];
    var[103] = _q_cov[8]*vrdt[2][2];
    var[104] = _cov[10][9]*rdt[0][1] + _cov[11][9]*rdt[0][2] + _cov[9][6] + _cov[9][9]*rdt[0][0];
    var[105] = _cov[10][10]*rdt[0][1] + _cov[10][6] + _cov[10][9]*rdt[0][0] + _cov[11][10]*rdt[0][2];
    var[106] = _cov[11][10]*rdt[0][1] + _cov[11][11]*rdt[0][2] + _cov[11][6] + _cov[11][9]*rdt[0][0];
    var[107] = _q_cov[6]*rdt[0][0];
    var[108] = _q_cov[7]*rdt[0][1];
    var[109] = _q_cov[8]*rdt[0][2];
    var[110] = _cov[10][9]*rdt[1][1] + _cov[11][9]*rdt[1][2] + _cov[9][7] + _cov[9][9]*rdt[1][0];
    var[111] = _cov[10][10]*rdt[1][1] + _cov[10][7] + _cov[10][9]*rdt[1][0] + _cov[11][10]*rdt[1][2];
    var[112] = _cov[11][10]*rdt[1][1] + _cov[11][11]*rdt[1][2] + _cov[11][7] + _cov[11][9]*rdt[1][0];
    var[113] = _cov[10][9]*rdt[2][1] + _cov[11][9]*rdt[2][2] + _cov[9][8] + _cov[9][9]*rdt[2][0];
    var[114] = _cov[10][10]*rdt[2][1] + _cov[10][8] + _cov[10][9]*rdt[2][0] + _cov[11][10]*rdt[2][2];
    var[115] = _cov[11][10]*rdt[2][1] + _cov[11][11]*rdt[2][2] + _cov[11][8] + _cov[11][9]*rdt[2][0];


    _cov[0][0] = _cov[0][0] + _cov[10][0]*prdt[0][1] + _cov[11][0]*prdt[0][2] + _cov[3][0]*_dt + _cov[9][0]*prdt[0][0] + _dt*var[0] + prdt[0][0]*var[1] + prdt[0][1]*var[2] + prdt[0][2]*var[3];    _cov[0][1] = _cov[10][1]*prdt[0][1] + _cov[11][1]*prdt[0][2] + _cov[1][0] + _cov[3][1]*_dt + _cov[9][1]*prdt[0][0] + _dt*var[8] + prdt[1][0]*var[1] + prdt[1][0]*var[4] + prdt[1][1]*var[2] + prdt[1][1]*var[5] + prdt[1][2]*var[3] + prdt[1][2]*var[6];
    _cov[1][1] = _cov[10][1]*prdt[1][1] + _cov[11][1]*prdt[1][2] + _cov[1][1] + _cov[4][1]*_dt + _cov[9][1]*prdt[1][0] + _dt*var[17] + prdt[1][0]*var[18] + prdt[1][1]*var[19] + prdt[1][2]*var[20];    _cov[0][2] = _cov[10][2]*prdt[0][1] + _cov[11][2]*prdt[0][2] + _cov[2][0] + _cov[3][2]*_dt + _cov[9][2]*prdt[0][0] + _dt*var[10] + prdt[2][0]*var[1] + prdt[2][0]*var[4] + prdt[2][1]*var[2] + prdt[2][1]*var[5] + prdt[2][2]*var[3] + prdt[2][2]*var[6];
    _cov[1][2] = _cov[10][2]*prdt[1][1] + _cov[11][2]*prdt[1][2] + _cov[2][1] + _cov[4][2]*_dt + _cov[9][2]*prdt[1][0] + _dt*var[25] + prdt[2][0]*var[18] + prdt[2][0]*var[21] + prdt[2][1]*var[19] + prdt[2][1]*var[22] + prdt[2][2]*var[20] + prdt[2][2]*var[23];
    _cov[2][2] = _cov[10][2]*prdt[2][1] + _cov[11][2]*prdt[2][2] + _cov[2][2] + _cov[5][2]*_dt + _cov[9][2]*prdt[2][0] + _dt*var[32] + prdt[2][0]*var[33] + prdt[2][1]*var[34] + prdt[2][2]*var[35];
    _cov[0][3] = -gdt*var[14] + rdt[0][0]*var[11] + rdt[0][1]*var[12] + rdt[0][2]*var[13] + var[0] + var[1]*vrdt[0][0] + var[2]*vrdt[0][1] + var[3]*vrdt[0][2] + var[4]*vrdt[0][0] + var[5]*vrdt[0][1] + var[6]*vrdt[0][2];
    _cov[1][3] = _cov[10][3]*prdt[1][1] + _cov[11][3]*prdt[1][2] + _cov[3][1] + _cov[9][3]*prdt[1][0] - gdt*var[29] + rdt[0][0]*var[26] + rdt[0][1]*var[27] + rdt[0][2]*var[28] + var[18]*vrdt[0][0] + var[19]*vrdt[0][1] + var[20]*vrdt[0][2] + var[21]*vrdt[0][0] + var[22]*vrdt[0][1] + var[23]*vrdt[0][2] + var[7];
    _cov[2][3] = _cov[10][3]*prdt[2][1] + _cov[11][3]*prdt[2][2] + _cov[3][2] + _cov[9][3]*prdt[2][0] - gdt*var[42] + rdt[0][0]*var[39] + rdt[0][1]*var[40] + rdt[0][2]*var[41] + var[33]*vrdt[0][0] + var[34]*vrdt[0][1] + var[35]*vrdt[0][2] + var[36]*vrdt[0][0] + var[37]*vrdt[0][1] + var[38]*vrdt[0][2] + var[9];
    _cov[3][3] = _cov[10][3]*vrdt[0][1] + _cov[11][3]*vrdt[0][2] + _cov[12][3]*rdt[0][0] + _cov[13][3]*rdt[0][1] + _cov[14][3]*rdt[0][2] + _cov[3][3] - _cov[7][3]*gdt + _cov[9][3]*vrdt[0][0] - gdt*var[58] + rdt[0][0]*var[49] + rdt[0][1]*var[50] + rdt[0][2]*var[51] + var[53]*vrdt[0][0] + var[55]*vrdt[0][1] + var[57]*vrdt[0][2];
    _cov[0][4] = gdt*var[15] + rdt[1][0]*var[11] + rdt[1][1]*var[12] + rdt[1][2]*var[13] + var[1]*vrdt[1][0] + var[2]*vrdt[1][1] + var[3]*vrdt[1][2] + var[4]*vrdt[1][0] + var[5]*vrdt[1][1] + var[6]*vrdt[1][2] + var[8];
    _cov[1][4] = gdt*var[30] + rdt[1][0]*var[26] + rdt[1][1]*var[27] + rdt[1][2]*var[28] + var[17] + var[18]*vrdt[1][0] + var[19]*vrdt[1][1] + var[20]*vrdt[1][2] + var[21]*vrdt[1][0] + var[22]*vrdt[1][1] + var[23]*vrdt[1][2];
    _cov[2][4] = _cov[10][4]*prdt[2][1] + _cov[11][4]*prdt[2][2] + _cov[4][2] + _cov[9][4]*prdt[2][0] + gdt*var[43] + rdt[1][0]*var[39] + rdt[1][1]*var[40] + rdt[1][2]*var[41] + var[24] + var[33]*vrdt[1][0] + var[34]*vrdt[1][1] + var[35]*vrdt[1][2] + var[36]*vrdt[1][0] + var[37]*vrdt[1][1] + var[38]*vrdt[1][2];
    _cov[3][4] = _cov[10][4]*vrdt[0][1] + _cov[11][4]*vrdt[0][2] + _cov[12][4]*rdt[0][0] + _cov[13][4]*rdt[0][1] + _cov[14][4]*rdt[0][2] + _cov[4][3] - _cov[7][4]*gdt + _cov[9][4]*vrdt[0][0] + gdt*var[66] + rdt[1][0]*var[49] + rdt[1][0]*var[59] + rdt[1][1]*var[50] + rdt[1][1]*var[60] + rdt[1][2]*var[51] + rdt[1][2]*var[61] + var[53]*vrdt[1][0] + var[55]*vrdt[1][1] + var[57]*vrdt[1][2] + var[62]*vrdt[1][0] + var[63]*vrdt[1][1] + var[64]*vrdt[1][2];
    _cov[4][4] = _cov[10][4]*vrdt[1][1] + _cov[11][4]*vrdt[1][2] + _cov[12][4]*rdt[1][0] + _cov[13][4]*rdt[1][1] + _cov[14][4]*rdt[1][2] + _cov[4][4] + _cov[6][4]*gdt + _cov[9][4]*vrdt[1][0] + gdt*var[71] + rdt[1][0]*var[72] + rdt[1][1]*var[73] + rdt[1][2]*var[74] + var[76]*vrdt[1][0] + var[78]*vrdt[1][1] + var[80]*vrdt[1][2];
    _cov[0][5] = _dt*var[16] + rdt[2][0]*var[11] + rdt[2][1]*var[12] + rdt[2][2]*var[13] + var[10] + var[1]*vrdt[2][0] + var[2]*vrdt[2][1] + var[3]*vrdt[2][2] + var[4]*vrdt[2][0] + var[5]*vrdt[2][1] + var[6]*vrdt[2][2];
    _cov[1][5] = _dt*var[31] + rdt[2][0]*var[26] + rdt[2][1]*var[27] + rdt[2][2]*var[28] + var[18]*vrdt[2][0] + var[19]*vrdt[2][1] + var[20]*vrdt[2][2] + var[21]*vrdt[2][0] + var[22]*vrdt[2][1] + var[23]*vrdt[2][2] + var[25];
    _cov[2][5] = _dt*var[45] + rdt[2][0]*var[39] + rdt[2][1]*var[40] + rdt[2][2]*var[41] + var[32] + var[33]*vrdt[2][0] + var[34]*vrdt[2][1] + var[35]*vrdt[2][2] + var[36]*vrdt[2][0] + var[37]*vrdt[2][1] + var[38]*vrdt[2][2];
    _cov[3][5] = _cov[10][5]*vrdt[0][1] + _cov[11][5]*vrdt[0][2] + _cov[12][5]*rdt[0][0] + _cov[13][5]*rdt[0][1] + _cov[14][5]*rdt[0][2] + _cov[5][3] - _cov[7][5]*gdt + _cov[9][5]*vrdt[0][0] + _dt*var[67] + rdt[2][0]*var[49] + rdt[2][0]*var[59] + rdt[2][1]*var[50] + rdt[2][1]*var[60] + rdt[2][2]*var[51] + rdt[2][2]*var[61] + var[53]*vrdt[2][0] + var[55]*vrdt[2][1] + var[57]*vrdt[2][2] + var[62]*vrdt[2][0] + var[63]*vrdt[2][1] + var[64]*vrdt[2][2];
    _cov[4][5] = _cov[10][5]*vrdt[1][1] + _cov[11][5]*vrdt[1][2] + _cov[12][5]*rdt[1][0] + _cov[13][5]*rdt[1][1] + _cov[14][5]*rdt[1][2] + _cov[5][4] + _cov[6][5]*gdt + _cov[9][5]*vrdt[1][0] + _dt*var[87] + _q_cov[3]*var[81] + _q_cov[4]*var[82] + _q_cov[5]*var[83] + rdt[2][0]*var[72] + rdt[2][1]*var[73] + rdt[2][2]*var[74] + var[76]*vrdt[2][0] + var[78]*vrdt[2][1] + var[80]*vrdt[2][2] + var[84]*vrdt[2][0] + var[85]*vrdt[2][1] + var[86]*vrdt[2][2];
    _cov[5][5] = _cov[10][5]*vrdt[2][1] + _cov[11][5]*vrdt[2][2] + _cov[12][5]*rdt[2][0] + _cov[13][5]*rdt[2][1] + _cov[14][5]*rdt[2][2] + _cov[5][5] + _cov[9][5]*vrdt[2][0] + _dt*var[91] + rdt[2][0]*var[92] + rdt[2][1]*var[93] + rdt[2][2]*var[94] + var[100]*vrdt[2][2] + var[44] + var[96]*vrdt[2][0] + var[98]*vrdt[2][1];
    _cov[0][6] = rdt[0][0]*var[1] + rdt[0][0]*var[4] + rdt[0][1]*var[2] + rdt[0][1]*var[5] + rdt[0][2]*var[3] + rdt[0][2]*var[6] + var[15];
    _cov[1][6] = rdt[0][0]*var[18] + rdt[0][0]*var[21] + rdt[0][1]*var[19] + rdt[0][1]*var[22] + rdt[0][2]*var[20] + rdt[0][2]*var[23] + var[30];
    _cov[2][6] = rdt[0][0]*var[33] + rdt[0][0]*var[36] + rdt[0][1]*var[34] + rdt[0][1]*var[37] + rdt[0][2]*var[35] + rdt[0][2]*var[38] + var[43];
    _cov[3][6] = rdt[0][0]*var[53] + rdt[0][0]*var[62] + rdt[0][1]*var[55] + rdt[0][1]*var[63] + rdt[0][2]*var[57] + rdt[0][2]*var[64] + var[66];
    _cov[4][6] = rdt[0][0]*var[76] + rdt[0][0]*var[84] + rdt[0][1]*var[78] + rdt[0][1]*var[85] + rdt[0][2]*var[80] + rdt[0][2]*var[86] + var[71];
    _cov[5][6] = _cov[10][6]*vrdt[2][1] + _cov[11][6]*vrdt[2][2] + _cov[12][6]*rdt[2][0] + _cov[13][6]*rdt[2][1] + _cov[14][6]*rdt[2][2] + _cov[15][6]*_dt + _cov[6][5] + _cov[9][6]*vrdt[2][0] + rdt[0][0]*var[101] + rdt[0][0]*var[96] + rdt[0][1]*var[102] + rdt[0][1]*var[98] + rdt[0][2]*var[100] + rdt[0][2]*var[103];
    _cov[6][6] = _cov[10][6]*rdt[0][1] + _cov[11][6]*rdt[0][2] + _cov[6][6] + _cov[9][6]*rdt[0][0] + rdt[0][0]*var[104] + rdt[0][1]*var[105] + rdt[0][2]*var[106];
    _cov[0][7] = rdt[1][0]*var[1] + rdt[1][0]*var[4] + rdt[1][1]*var[2] + rdt[1][1]*var[5] + rdt[1][2]*var[3] + rdt[1][2]*var[6] + var[14];
    _cov[1][7] = rdt[1][0]*var[18] + rdt[1][0]*var[21] + rdt[1][1]*var[19] + rdt[1][1]*var[22] + rdt[1][2]*var[20] + rdt[1][2]*var[23] + var[29];
    _cov[2][7] = rdt[1][0]*var[33] + rdt[1][0]*var[36] + rdt[1][1]*var[34] + rdt[1][1]*var[37] + rdt[1][2]*var[35] + rdt[1][2]*var[38] + var[42];
    _cov[3][7] = rdt[1][0]*var[53] + rdt[1][0]*var[62] + rdt[1][1]*var[55] + rdt[1][1]*var[63] + rdt[1][2]*var[57] + rdt[1][2]*var[64] + var[58];
    _cov[4][7] = _cov[10][7]*vrdt[1][1] + _cov[11][7]*vrdt[1][2] + _cov[12][7]*rdt[1][0] + _cov[13][7]*rdt[1][1] + _cov[14][7]*rdt[1][2] + _cov[7][4] + _cov[9][7]*vrdt[1][0] + rdt[1][0]*var[76] + rdt[1][0]*var[84] + rdt[1][1]*var[78] + rdt[1][1]*var[85] + rdt[1][2]*var[80] + rdt[1][2]*var[86] + var[65];
    _cov[5][7] = _cov[10][7]*vrdt[2][1] + _cov[11][7]*vrdt[2][2] + _cov[12][7]*rdt[2][0] + _cov[13][7]*rdt[2][1] + _cov[14][7]*rdt[2][2] + _cov[15][7]*_dt + _cov[7][5] + _cov[9][7]*vrdt[2][0] + rdt[1][0]*var[101] + rdt[1][0]*var[96] + rdt[1][1]*var[102] + rdt[1][1]*var[98] + rdt[1][2]*var[100] + rdt[1][2]*var[103];
    _cov[6][7] = _cov[10][7]*rdt[0][1] + _cov[11][7]*rdt[0][2] + _cov[7][6] + _cov[9][7]*rdt[0][0] + rdt[1][0]*var[104] + rdt[1][0]*var[107] + rdt[1][1]*var[105] + rdt[1][1]*var[108] + rdt[1][2]*var[106] + rdt[1][2]*var[109];
    _cov[7][7] = _cov[10][7]*rdt[1][1] + _cov[11][7]*rdt[1][2] + _cov[7][7] + _cov[9][7]*rdt[1][0] + rdt[1][0]*var[110] + rdt[1][1]*var[111] + rdt[1][2]*var[112];
    _cov[0][8] = _cov[10][8]*prdt[0][1] + _cov[11][8]*prdt[0][2] + _cov[8][0] + _cov[8][3]*_dt + _cov[9][8]*prdt[0][0] + rdt[2][0]*var[1] + rdt[2][0]*var[4] + rdt[2][1]*var[2] + rdt[2][1]*var[5] + rdt[2][2]*var[3] + rdt[2][2]*var[6];
    _cov[1][8] = _cov[10][8]*prdt[1][1] + _cov[11][8]*prdt[1][2] + _cov[8][1] + _cov[8][4]*_dt + _cov[9][8]*prdt[1][0] + rdt[2][0]*var[18] + rdt[2][0]*var[21] + rdt[2][1]*var[19] + rdt[2][1]*var[22] + rdt[2][2]*var[20] + rdt[2][2]*var[23];
    _cov[2][8] = _cov[10][8]*prdt[2][1] + _cov[11][8]*prdt[2][2] + _cov[8][2] + _cov[8][5]*_dt + _cov[9][8]*prdt[2][0] + rdt[2][0]*var[33] + rdt[2][0]*var[36] + rdt[2][1]*var[34] + rdt[2][1]*var[37] + rdt[2][2]*var[35] + rdt[2][2]*var[38];
    _cov[3][8] = _cov[10][8]*vrdt[0][1] + _cov[11][8]*vrdt[0][2] + _cov[12][8]*rdt[0][0] + _cov[13][8]*rdt[0][1] + _cov[14][8]*rdt[0][2] + _cov[8][3] - _cov[8][7]*gdt + _cov[9][8]*vrdt[0][0] + rdt[2][0]*var[53] + rdt[2][0]*var[62] + rdt[2][1]*var[55] + rdt[2][1]*var[63] + rdt[2][2]*var[57] + rdt[2][2]*var[64];
    _cov[4][8] = _cov[10][8]*vrdt[1][1] + _cov[11][8]*vrdt[1][2] + _cov[12][8]*rdt[1][0] + _cov[13][8]*rdt[1][1] + _cov[14][8]*rdt[1][2] + _cov[8][4] + _cov[8][6]*gdt + _cov[9][8]*vrdt[1][0] + rdt[2][0]*var[76] + rdt[2][0]*var[84] + rdt[2][1]*var[78] + rdt[2][1]*var[85] + rdt[2][2]*var[80] + rdt[2][2]*var[86];
    _cov[5][8] = _cov[10][8]*vrdt[2][1] + _cov[11][8]*vrdt[2][2] + _cov[12][8]*rdt[2][0] + _cov[13][8]*rdt[2][1] + _cov[14][8]*rdt[2][2] + _cov[15][8]*_dt + _cov[8][5] + _cov[9][8]*vrdt[2][0] + rdt[2][0]*var[101] + rdt[2][0]*var[96] + rdt[2][1]*var[102] + rdt[2][1]*var[98] + rdt[2][2]*var[100] + rdt[2][2]*var[103];
    _cov[6][8] = _cov[10][8]*rdt[0][1] + _cov[11][8]*rdt[0][2] + _cov[8][6] + _cov[9][8]*rdt[0][0] + rdt[2][0]*var[104] + rdt[2][0]*var[107] + rdt[2][1]*var[105] + rdt[2][1]*var[108] + rdt[2][2]*var[106] + rdt[2][2]*var[109];
    _cov[7][8] = _cov[10][8]*rdt[1][1] + _cov[11][8]*rdt[1][2] + _cov[8][7] + _cov[9][8]*rdt[1][0] + _q_cov[6]*var[81] + _q_cov[7]*var[82] + _q_cov[8]*var[83] + rdt[2][0]*var[110] + rdt[2][1]*var[111] + rdt[2][2]*var[112];
    _cov[8][8] = _cov[10][8]*rdt[2][1] + _cov[11][8]*rdt[2][2] + _cov[8][8] + _cov[9][8]*rdt[2][0] + _q_cov[6]*var[88] + _q_cov[7]*var[89] + _q_cov[8]*var[90] + rdt[2][0]*var[113] + rdt[2][1]*var[114] + rdt[2][2]*var[115];
    _cov[0][9] = var[1];
    _cov[1][9] = var[18];
    _cov[2][9] = var[33];
    _cov[3][9] = var[53];
    _cov[4][9] = var[76];
    _cov[5][9] = var[96];
    _cov[6][9] = var[104];
    _cov[7][9] = var[110];
    _cov[8][9] = var[113];
    // _cov[9][9] = _cov[9][9];
    _cov[0][10] = var[2];
    _cov[1][10] = var[19];
    _cov[2][10] = var[34];
    _cov[3][10] = var[55];
    _cov[4][10] = var[78];
    _cov[5][10] = var[98];
    _cov[6][10] = var[105];
    _cov[7][10] = var[111];
    _cov[8][10] = var[114];
    // _cov[9][10] = _cov[10][9];
    // _cov[10][10] = _cov[10][10];
    _cov[0][11] = var[3];
    _cov[1][11] = var[20];
    _cov[2][11] = var[35];
    _cov[3][11] = var[57];
    _cov[4][11] = var[80];
    _cov[5][11] = var[100];
    _cov[6][11] = var[106];
    _cov[7][11] = var[112];
    _cov[8][11] = var[115];
    // _cov[9][11] = _cov[11][9];
    // _cov[10][11] = _cov[11][10];
    // _cov[11][11] = _cov[11][11];
    _cov[0][12] = var[11];
    _cov[1][12] = var[26];
    _cov[2][12] = var[39];
    _cov[3][12] = var[49];
    _cov[4][12] = var[72];
    _cov[5][12] = var[92];
    _cov[6][12] = _cov[12][10]*rdt[0][1] + _cov[12][11]*rdt[0][2] + _cov[12][6] + var[52];
    _cov[7][12] = _cov[12][10]*rdt[1][1] + _cov[12][11]*rdt[1][2] + _cov[12][7] + var[75];
    _cov[8][12] = _cov[12][10]*rdt[2][1] + _cov[12][11]*rdt[2][2] + _cov[12][8] + var[95];
    // _cov[9][12] = _cov[12][9];
    // _cov[10][12] = _cov[12][10];
    // _cov[11][12] = _cov[12][11];
    // _cov[12][12] = _cov[12][12];
    _cov[0][13] = var[12];
    _cov[1][13] = var[27];
    _cov[2][13] = var[40];
    _cov[3][13] = var[50];
    _cov[4][13] = var[73];
    _cov[5][13] = var[93];
    _cov[6][13] = _cov[13][11]*rdt[0][2] + _cov[13][6] + _cov[13][9]*rdt[0][0] + var[54];
    _cov[7][13] = _cov[13][11]*rdt[1][2] + _cov[13][7] + _cov[13][9]*rdt[1][0] + var[77];
    _cov[8][13] = _cov[13][11]*rdt[2][2] + _cov[13][8] + _cov[13][9]*rdt[2][0] + var[97];
    // _cov[9][13] = _cov[13][9];
    // _cov[10][13] = _cov[13][10];
    // _cov[11][13] = _cov[13][11];
    // _cov[12][13] = _cov[13][12];
    // _cov[13][13] = _cov[13][13];
    _cov[0][14] = var[13];
    _cov[1][14] = var[28];
    _cov[2][14] = var[41];
    _cov[3][14] = var[51];
    _cov[4][14] = var[74];
    _cov[5][14] = var[94];
    _cov[6][14] = _cov[14][10]*rdt[0][1] + _cov[14][6] + _cov[14][9]*rdt[0][0] + var[56];
    _cov[7][14] = _cov[14][10]*rdt[1][1] + _cov[14][7] + _cov[14][9]*rdt[1][0] + var[79];
    _cov[8][14] = _cov[14][10]*rdt[2][1] + _cov[14][8] + _cov[14][9]*rdt[2][0] + var[99];
    // _cov[9][14] = _cov[14][9];
    // _cov[10][14] = _cov[14][10];
    // _cov[11][14] = _cov[14][11];
    // _cov[12][14] = _cov[14][12];
    // _cov[13][14] = _cov[14][13];
    // _cov[14][14] = _cov[14][14];
    _cov[0][15] = var[16];
    _cov[1][15] = var[31];
    _cov[2][15] = var[45];
    _cov[3][15] = var[67];
    _cov[4][15] = var[87];
    _cov[5][15] = var[91];
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

    _cov[0][0] = kahan_summation(_cov[0][0], _q_cov[6]*(prdt[0][0])*(prdt[0][0]) + _q_cov[7]*(prdt[0][1])*(prdt[0][1]) + _q_cov[8]*(prdt[0][2])*(prdt[0][2]), _accumulator_cov[0]);
    _cov[1][1] = kahan_summation(_cov[1][1], _q_cov[6]*(prdt[1][0])*(prdt[1][0]) + _q_cov[7]*(prdt[1][1])*(prdt[1][1]) + _q_cov[8]*(prdt[1][2])*(prdt[1][2]), _accumulator_cov[1]);
    _cov[2][2] = kahan_summation(_cov[2][2], _q_cov[6]*(prdt[2][0])*(prdt[2][0]) + _q_cov[7]*(prdt[2][1])*(prdt[2][1]) + _q_cov[8]*(prdt[2][2])*(prdt[2][2]), _accumulator_cov[2]);
    _cov[3][3] = kahan_summation(_cov[3][3], _q_cov[3]*var[46] + _q_cov[4]*var[47] + _q_cov[5]*var[48] + _q_cov[6]*(vrdt[0][0])*(vrdt[0][0]) + _q_cov[7]*(vrdt[0][1])*(vrdt[0][1]) + _q_cov[8]*(vrdt[0][2])*(vrdt[0][2]), _accumulator_cov[3]);
    _cov[4][4] = kahan_summation(_cov[4][4], _q_cov[3]*var[68] + _q_cov[4]*var[69] + _q_cov[5]*var[70] + _q_cov[6]*(vrdt[1][0])*(vrdt[1][0]) + _q_cov[7]*(vrdt[1][1])*(vrdt[1][1]) + _q_cov[8]*(vrdt[1][2])*(vrdt[1][2]), _accumulator_cov[4]);
    _cov[5][5] = kahan_summation(_cov[5][5], _q_cov[3]*var[88] + _q_cov[4]*var[89] + _q_cov[5]*var[90] + _q_cov[6]*(vrdt[2][0])*(vrdt[2][0]) + _q_cov[7]*(vrdt[2][1])*(vrdt[2][1]) + _q_cov[8]*(vrdt[2][2])*(vrdt[2][2]), _accumulator_cov[5]);
    _cov[6][6] = kahan_summation(_cov[6][6], _q_cov[6]*var[46] + _q_cov[7]*var[47] + _q_cov[8]*var[48], _accumulator_cov[6]);
    _cov[7][7] = kahan_summation(_cov[7][7], _q_cov[6]*var[68] + _q_cov[7]*var[69] + _q_cov[8]*var[70], _accumulator_cov[7]);
    _cov[8][8] = kahan_summation(_cov[8][8], _q_cov[6]*var[88] + _q_cov[7]*var[89] + _q_cov[8]*var[90], _accumulator_cov[8]);
    _cov[9][9] = kahan_summation(_cov[9][9], _q_cov[9] * _dt2, _accumulator_cov[9]);
    _cov[10][10] = kahan_summation(_cov[10][10], _q_cov[10] * _dt2, _accumulator_cov[10]);
    _cov[11][11] = kahan_summation(_cov[11][11], _q_cov[11] * _dt2, _accumulator_cov[11]);
    _cov[12][12] = kahan_summation(_cov[12][12], _q_cov[12] * _dt2, _accumulator_cov[12]);
    _cov[13][13] = kahan_summation(_cov[13][13], _q_cov[13] * _dt2, _accumulator_cov[13]);
    _cov[14][14] = kahan_summation(_cov[14][14], _q_cov[14] * _dt2, _accumulator_cov[14]);
    _cov[15][15] = kahan_summation(_cov[15][15], _q_cov[15] * _dt2, _accumulator_cov[15]);

    regular_covariance_to_symmetric();
}

unsigned char RIEKF::fuse_position(const Vector3f &pos, const Vector3f &w, const Vector3f &a, 
                                  const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) {
    /*
    s = [p, v, R, bg, ba, g]
    δs = [δp, δv, δθ, δbg, δba, δg]
    s + δs = [Exp(δθ)*(p+δp), Exp(δθ)*(v+δv), Exp(δθ)*R, bg+δbg, ba+δba, g+δg]

    pos = p + R * dis

    δpos / δp = I
    δpos / δθ = -(p + R * dis)^

    H = [I, O, -(p + R*dis)^, O, O, O]
    */

    const Vector3f h = _p + _rot * dis;
    const array<array<float, 3>, 3> minus_h_hat {
        {{0.f, h[2], -h[1]},
         {-h[2], 0.f, h[0]},
         {h[1], -h[0], 0.f}}
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

        where, H = [I, O, -(p + R*dis)^, O, O, O]
               h = p + R * dis
        */

        // H * P  or  P * H'
        const float cov_1_dim = (dim == 2) ? _cov[1][dim] : _cov[dim][1];
        const array<float, 16> HP = {_cov[0][dim] + _cov[0][6]*minus_h_hat[dim][0] + _cov[0][7]*minus_h_hat[dim][1] + _cov[0][8]*minus_h_hat[dim][2],
                                     cov_1_dim + _cov[1][6]*minus_h_hat[dim][0] + _cov[1][7]*minus_h_hat[dim][1] + _cov[1][8]*minus_h_hat[dim][2],
                                     _cov[dim][2] + _cov[2][6]*minus_h_hat[dim][0] + _cov[2][7]*minus_h_hat[dim][1] + _cov[2][8]*minus_h_hat[dim][2],
                                     _cov[dim][3] + _cov[3][6]*minus_h_hat[dim][0] + _cov[3][7]*minus_h_hat[dim][1] + _cov[3][8]*minus_h_hat[dim][2],
                                     _cov[dim][4] + _cov[4][6]*minus_h_hat[dim][0] + _cov[4][7]*minus_h_hat[dim][1] + _cov[4][8]*minus_h_hat[dim][2],
                                     _cov[dim][5] + _cov[5][6]*minus_h_hat[dim][0] + _cov[5][7]*minus_h_hat[dim][1] + _cov[5][8]*minus_h_hat[dim][2],
                                     _cov[dim][6] + _cov[6][6]*minus_h_hat[dim][0] + _cov[6][7]*minus_h_hat[dim][1] + _cov[6][8]*minus_h_hat[dim][2],
                                     _cov[dim][7] + _cov[6][7]*minus_h_hat[dim][0] + _cov[7][7]*minus_h_hat[dim][1] + _cov[7][8]*minus_h_hat[dim][2],
                                     _cov[dim][8] + _cov[6][8]*minus_h_hat[dim][0] + _cov[7][8]*minus_h_hat[dim][1] + _cov[8][8]*minus_h_hat[dim][2],
                                     _cov[dim][9] + _cov[6][9]*minus_h_hat[dim][0] + _cov[7][9]*minus_h_hat[dim][1] + _cov[8][9]*minus_h_hat[dim][2],
                                     _cov[dim][10] + _cov[6][10]*minus_h_hat[dim][0] + _cov[7][10]*minus_h_hat[dim][1] + _cov[8][10]*minus_h_hat[dim][2],
                                     _cov[dim][11] + _cov[6][11]*minus_h_hat[dim][0] + _cov[7][11]*minus_h_hat[dim][1] + _cov[8][11]*minus_h_hat[dim][2],
                                     _cov[dim][12] + _cov[6][12]*minus_h_hat[dim][0] + _cov[7][12]*minus_h_hat[dim][1] + _cov[8][12]*minus_h_hat[dim][2],
                                     _cov[dim][13] + _cov[6][13]*minus_h_hat[dim][0] + _cov[7][13]*minus_h_hat[dim][1] + _cov[8][13]*minus_h_hat[dim][2],
                                     _cov[dim][14] + _cov[6][14]*minus_h_hat[dim][0] + _cov[7][14]*minus_h_hat[dim][1] + _cov[8][14]*minus_h_hat[dim][2],
                                     _cov[dim][15] + _cov[6][15]*minus_h_hat[dim][0] + _cov[7][15]*minus_h_hat[dim][1] + _cov[8][15]*minus_h_hat[dim][2]};

        // H * P * H' + R
        const float HPHT_plus_R = HP[dim] + HP[6] * minus_h_hat[dim][0] + HP[7] * minus_h_hat[dim][1] + HP[8] * minus_h_hat[dim][2] + noise_std[dim] * noise_std[dim];

        // h = p + R * dis
        // e = y - h = pos - (p + R * dis)
        const float obs_error = pos[dim] - h[dim];

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

unsigned char RIEKF::fuse_velocity(const Vector3f &vel, const Vector3f &w, const Vector3f &a, 
                                  const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) {
    /*
    s = [p, v, R, bg, ba, g]
    δs = [δp, δv, δθ, δbg, δba, δg]
    s + δs = [Exp(δθ)*(p+δp), Exp(δθ)*(v+δv), Exp(δθ)*R, bg+δbg, ba+δba, g+δg]

    vel = v + R * (w - bg)^ * dis

    δvel / δv = I
    δvel / δθ = -(v + R*((w - bg)^*dis))^
    δvel / δbg = R * dis^

    H = [O, I, -(v + R*(w-bg)^*dis)^, R*dis^, O, O]
    */ 

    unsigned char info = 0;

    // v + R * (w-bg)^ * dis
    const Vector3f h = _v + _rot * dis.cross(_bg - w);
    const array<array<float, 3>, 3> minus_h_hat {
        {{0.f, h[2], -h[1]},
         {-h[2], 0.f, h[0]},
         {h[1], -h[0], 0.f}}
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

        where, H = [O, I, (v + R*(dis^*(w-bg)))^, R*dis^, O, O]
               h = v + R * (w-bg)^ * dis
        */

        // R * dis^
        const array<float, 3> rot_d_hat = {_rot(dim, 1)*dis[2] - _rot(dim, 2)*dis[1], 
                                           _rot(dim, 2)*dis[0] - _rot(dim, 0)*dis[2], 
                                           _rot(dim, 0)*dis[1] - _rot(dim, 1)*dis[0]};

        // H * P  or  P * H'
        const unsigned int index = 3 + dim;
        const float cov_4_index = (dim == 2) ? _cov[4][index] : _cov[index][4];
        const array<float, 16> HP = {_cov[0][index] + _cov[0][6]*minus_h_hat[dim][0] + _cov[0][7]*minus_h_hat[dim][1] + _cov[0][8]*minus_h_hat[dim][2] + _cov[0][9]*rot_d_hat[0] + _cov[0][10]*rot_d_hat[1] + _cov[0][11]*rot_d_hat[2],
                                     _cov[1][index] + _cov[1][6]*minus_h_hat[dim][0] + _cov[1][7]*minus_h_hat[dim][1] + _cov[1][8]*minus_h_hat[dim][2] + _cov[1][9]*rot_d_hat[0] + _cov[1][10]*rot_d_hat[1] + _cov[1][11]*rot_d_hat[2],
                                     _cov[2][index] + _cov[2][6]*minus_h_hat[dim][0] + _cov[2][7]*minus_h_hat[dim][1] + _cov[2][8]*minus_h_hat[dim][2] + _cov[2][9]*rot_d_hat[0] + _cov[2][10]*rot_d_hat[1] + _cov[2][11]*rot_d_hat[2],
                                     _cov[3][index] + _cov[3][6]*minus_h_hat[dim][0] + _cov[3][7]*minus_h_hat[dim][1] + _cov[3][8]*minus_h_hat[dim][2] + _cov[3][9]*rot_d_hat[0] + _cov[3][10]*rot_d_hat[1] + _cov[3][11]*rot_d_hat[2],
                                     cov_4_index + _cov[4][6]*minus_h_hat[dim][0] + _cov[4][7]*minus_h_hat[dim][1] + _cov[4][8]*minus_h_hat[dim][2] + _cov[4][9]*rot_d_hat[0] + _cov[4][10]*rot_d_hat[1] + _cov[4][11]*rot_d_hat[2],
                                     _cov[index][5] + _cov[5][6]*minus_h_hat[dim][0] + _cov[5][7]*minus_h_hat[dim][1] + _cov[5][8]*minus_h_hat[dim][2] + _cov[5][9]*rot_d_hat[0] + _cov[5][10]*rot_d_hat[1] + _cov[5][11]*rot_d_hat[2],
                                     _cov[index][6] + _cov[6][6]*minus_h_hat[dim][0] + _cov[6][7]*minus_h_hat[dim][1] + _cov[6][8]*minus_h_hat[dim][2] + _cov[6][9]*rot_d_hat[0] + _cov[6][10]*rot_d_hat[1] + _cov[6][11]*rot_d_hat[2],
                                     _cov[index][7] + _cov[6][7]*minus_h_hat[dim][0] + _cov[7][7]*minus_h_hat[dim][1] + _cov[7][8]*minus_h_hat[dim][2] + _cov[7][9]*rot_d_hat[0] + _cov[7][10]*rot_d_hat[1] + _cov[7][11]*rot_d_hat[2],
                                     _cov[index][8] + _cov[6][8]*minus_h_hat[dim][0] + _cov[7][8]*minus_h_hat[dim][1] + _cov[8][8]*minus_h_hat[dim][2] + _cov[8][9]*rot_d_hat[0] + _cov[8][10]*rot_d_hat[1] + _cov[8][11]*rot_d_hat[2],
                                     _cov[index][9] + _cov[6][9]*minus_h_hat[dim][0] + _cov[7][9]*minus_h_hat[dim][1] + _cov[8][9]*minus_h_hat[dim][2] + _cov[9][9]*rot_d_hat[0] + _cov[9][10]*rot_d_hat[1] + _cov[9][11]*rot_d_hat[2],
                                     _cov[index][10] + _cov[6][10]*minus_h_hat[dim][0] + _cov[7][10]*minus_h_hat[dim][1] + _cov[8][10]*minus_h_hat[dim][2] + _cov[9][10]*rot_d_hat[0] + _cov[10][10]*rot_d_hat[1] + _cov[10][11]*rot_d_hat[2],
                                     _cov[index][11] + _cov[6][11]*minus_h_hat[dim][0] + _cov[7][11]*minus_h_hat[dim][1] + _cov[8][11]*minus_h_hat[dim][2] + _cov[9][11]*rot_d_hat[0] + _cov[10][11]*rot_d_hat[1] + _cov[11][11]*rot_d_hat[2],
                                     _cov[index][12] + _cov[6][12]*minus_h_hat[dim][0] + _cov[7][12]*minus_h_hat[dim][1] + _cov[8][12]*minus_h_hat[dim][2] + _cov[9][12]*rot_d_hat[0] + _cov[10][12]*rot_d_hat[1] + _cov[11][12]*rot_d_hat[2],
                                     _cov[index][13] + _cov[6][13]*minus_h_hat[dim][0] + _cov[7][13]*minus_h_hat[dim][1] + _cov[8][13]*minus_h_hat[dim][2] + _cov[9][13]*rot_d_hat[0] + _cov[10][13]*rot_d_hat[1] + _cov[11][13]*rot_d_hat[2],
                                     _cov[index][14] + _cov[6][14]*minus_h_hat[dim][0] + _cov[7][14]*minus_h_hat[dim][1] + _cov[8][14]*minus_h_hat[dim][2] + _cov[9][14]*rot_d_hat[0] + _cov[10][14]*rot_d_hat[1] + _cov[11][14]*rot_d_hat[2],
                                     _cov[index][15] + _cov[6][15]*minus_h_hat[dim][0] + _cov[7][15]*minus_h_hat[dim][1] + _cov[8][15]*minus_h_hat[dim][2] + _cov[9][15]*rot_d_hat[0] + _cov[10][15]*rot_d_hat[1] + _cov[11][15]*rot_d_hat[2]};                                                  
        // H * P * H' + R
        const float HPHT_plus_R = HP[index] + HP[6] * minus_h_hat[dim][0] + HP[7] * minus_h_hat[dim][1] + HP[8] * minus_h_hat[dim][2] + HP[9] * rot_d_hat[0] + HP[10] * rot_d_hat[1] + HP[11] * rot_d_hat[2] + noise_std[dim] * noise_std[dim];

        // h = v + R * (w - bg)^ * dis
        // e = y - h = y - (v + R * (w - bg)^ * dis)
        const float obs_error = vel[dim] - h[dim];

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

void RIEKF::correct_state() {
    // for (float &es : _error_state) {
    //     cout << es << endl;
    // }
    // cout << "-------------------------" << endl;
    // state: [p, v, bg, ba, g], R
    // error_state : [δp, δv, δθ, δbg, δba, δg]

    // q = Exp(δθ) * q
    Quaternionf delta_q;
    array<float, 3> delta_theta = {_error_state[6], _error_state[7], _error_state[8]};
    quaternion_from_axis_angle(delta_q, delta_theta);
    const Quaternionf q = _q;
    _q = delta_q * q;
    _rot = q;

    /*
    p = p + Exp(δθ) * (p + δp)
    v = v + Exp(δθ) * (v + δv)
    */
    const Matrix3f delta_r = delta_q.toRotationMatrix();
    const Vector3f p_corr(_p[0] + _error_state[0], _p[1] + _error_state[1], _p[2] + _error_state[2]);
    const Vector3f v_corr(_v[0] + _error_state[3], _v[1] + _error_state[4], _v[2] + _error_state[5]);
    _p = delta_r * p_corr;
    _v = delta_r * v_corr;

    /*
    bg = bg + δbg
    ba = ba + δba
    g = g + δg
    */
    for (unsigned char i = 0; i < 3; ++i) {
        _bg[i] += _error_state[9 + i];
        _ba[i] += _error_state[12 + i];
    }
    _g += _error_state[15];

    // [δp, δv, δθ, δbg, δba, δg] = 0
    for (float &es : _error_state) {
        es = 0.f;
    }
}

void RIEKF::correct_covariance() {
    array<float, 58> var;

    var[0] = _cov[2][0] - _cov[2][1]*_error_state[8] + _cov[2][2]*_error_state[7];
    var[1] = _cov[1][0] - _cov[1][1]*_error_state[8] + _cov[2][1]*_error_state[7];
    var[2] = _cov[2][0]*_error_state[7];
    var[3] = _cov[1][0]*_error_state[8];
    var[4] = _cov[0][0] + var[2] - var[3];
    var[5] = _cov[5][0] - _cov[5][1]*_error_state[8] + _cov[5][2]*_error_state[7];
    var[6] = _cov[4][0] - _cov[4][1]*_error_state[8] + _cov[4][2]*_error_state[7];
    var[7] = _cov[3][0] - _cov[3][1]*_error_state[8] + _cov[3][2]*_error_state[7];
    var[8] = _cov[8][0] - _cov[8][1]*_error_state[8] + _cov[8][2]*_error_state[7];
    var[9] = 0.5F*var[8];
    var[10] = _cov[7][0] - _cov[7][1]*_error_state[8] + _cov[7][2]*_error_state[7];
    var[11] = 0.5F*_error_state[8];
    var[12] = _cov[6][0] - _cov[6][1]*_error_state[8] + _cov[6][2]*_error_state[7];
    var[13] = 0.5F*_error_state[6];
    var[14] = 0.5F*_error_state[7];
    var[15] = _cov[0][0]*_error_state[8] + _cov[1][0] - _cov[2][0]*_error_state[6];
    var[16] = _cov[2][0]*_error_state[8] + _cov[2][1] - _cov[2][2]*_error_state[6];
    var[17] = _cov[2][1]*_error_state[6];
    var[18] = _cov[1][1] - var[17] + var[3];
    var[19] = _cov[5][0]*_error_state[8] + _cov[5][1] - _cov[5][2]*_error_state[6];
    var[20] = _cov[4][0]*_error_state[8] + _cov[4][1] - _cov[4][2]*_error_state[6];
    var[21] = _cov[3][0]*_error_state[8] + _cov[3][1] - _cov[3][2]*_error_state[6];
    var[22] = _cov[8][0]*_error_state[8] + _cov[8][1] - _cov[8][2]*_error_state[6];
    var[23] = _cov[7][0]*_error_state[8] + _cov[7][1] - _cov[7][2]*_error_state[6];
    var[24] = _cov[6][0]*_error_state[8] + _cov[6][1] - _cov[6][2]*_error_state[6];
    var[25] = -_cov[5][0]*_error_state[7] + _cov[5][1]*_error_state[6] + _cov[5][2];
    var[26] = -_cov[4][0]*_error_state[7] + _cov[4][1]*_error_state[6] + _cov[4][2];
    var[27] = -_cov[3][0]*_error_state[7] + _cov[3][1]*_error_state[6] + _cov[3][2];
    var[28] = -_cov[8][0]*_error_state[7] + _cov[8][1]*_error_state[6] + _cov[8][2];
    var[29] = -_cov[7][0]*_error_state[7] + _cov[7][1]*_error_state[6] + _cov[7][2];
    var[30] = -_cov[6][0]*_error_state[7] + _cov[6][1]*_error_state[6] + _cov[6][2];
    var[31] = _cov[5][3] - _cov[5][4]*_error_state[8] + _cov[5][5]*_error_state[7];
    var[32] = _cov[4][3] - _cov[4][4]*_error_state[8] + _cov[5][4]*_error_state[7];
    var[33] = _cov[5][3]*_error_state[7];
    var[34] = _cov[4][3]*_error_state[8];
    var[35] = _cov[3][3] + var[33] - var[34];
    var[36] = _cov[8][3] - _cov[8][4]*_error_state[8] + _cov[8][5]*_error_state[7];
    var[37] = _cov[7][3] - _cov[7][4]*_error_state[8] + _cov[7][5]*_error_state[7];
    var[38] = _cov[6][3] - _cov[6][4]*_error_state[8] + _cov[6][5]*_error_state[7];
    var[39] = _cov[3][3]*_error_state[8] + _cov[4][3] - _cov[5][3]*_error_state[6];
    var[40] = _cov[5][3]*_error_state[8] + _cov[5][4] - _cov[5][5]*_error_state[6];
    var[41] = _cov[5][4]*_error_state[6];
    var[42] = _cov[4][4] + var[34] - var[41];
    var[43] = _cov[8][3]*_error_state[8] + _cov[8][4] - _cov[8][5]*_error_state[6];
    var[44] = _cov[7][3]*_error_state[8] + _cov[7][4] - _cov[7][5]*_error_state[6];
    var[45] = _cov[6][3]*_error_state[8] + _cov[6][4] - _cov[6][5]*_error_state[6];
    var[46] = -_cov[8][3]*_error_state[7] + _cov[8][4]*_error_state[6] + _cov[8][5];
    var[47] = -_cov[7][3]*_error_state[7] + _cov[7][4]*_error_state[6] + _cov[7][5];
    var[48] = -_cov[6][3]*_error_state[7] + _cov[6][4]*_error_state[6] + _cov[6][5];
    var[49] = _cov[8][6] - _cov[8][7]*var[11] + _cov[8][8]*var[14];
    var[50] = _cov[7][6] - _cov[7][7]*var[11] + _cov[8][7]*var[14];
    var[51] = _cov[7][6]*var[11];
    var[52] = _cov[8][6]*var[14];
    var[53] = _cov[6][6] - var[51] + var[52];
    var[54] = _cov[6][6]*var[11] + _cov[7][6] - _cov[8][6]*var[13];
    var[55] = _cov[8][6]*var[11] + _cov[8][7] - _cov[8][8]*var[13];
    var[56] = _cov[8][7]*var[13];
    var[57] = _cov[7][7] + var[51] - var[56];

    _cov[2][2] = _cov[2][2] + _error_state[6]*(-_cov[1][0]*_error_state[7] + _cov[1][1]*_error_state[6] + _cov[2][1]) - _error_state[7]*(-_cov[0][0]*_error_state[7] + _cov[1][0]*_error_state[6] + _cov[2][0]) + var[17] - var[2];
    _cov[0][0] = _error_state[7]*var[0] - _error_state[8]*var[1] + var[4];
    _cov[0][1] = -_error_state[6]*var[0] + _error_state[8]*var[4] + var[1];
    _cov[1][1] = -_error_state[6]*var[16] + _error_state[8]*var[15] + var[18];
    _cov[0][2] = _error_state[6]*var[1] - _error_state[7]*var[4] + var[0];
    _cov[1][2] = _error_state[6]*var[18] - _error_state[7]*var[15] + var[16];
    _cov[0][3] = _error_state[7]*var[5] - _error_state[8]*var[6] + var[7];
    _cov[1][3] = _error_state[7]*var[19] - _error_state[8]*var[20] + var[21];
    _cov[2][3] = _error_state[7]*var[25] - _error_state[8]*var[26] + var[27];

    _cov[5][5] = _cov[5][5] + _error_state[6]*(-_cov[4][3]*_error_state[7] + _cov[4][4]*_error_state[6] + _cov[5][4]) - _error_state[7]*(-_cov[3][3]*_error_state[7] + _cov[4][3]*_error_state[6] + _cov[5][3]) - var[33] + var[41];
    _cov[3][3] = _error_state[7]*var[31] - _error_state[8]*var[32] + var[35];
    _cov[0][4] = -_error_state[6]*var[5] + _error_state[8]*var[7] + var[6];
    _cov[1][4] = -_error_state[6]*var[19] + _error_state[8]*var[21] + var[20];
    _cov[2][4] = -_error_state[6]*var[25] + _error_state[8]*var[27] + var[26];
    _cov[3][4] = -_error_state[6]*var[31] + _error_state[8]*var[35] + var[32];
    _cov[4][4] = -_error_state[6]*var[40] + _error_state[8]*var[39] + var[42];
    _cov[0][5] = _error_state[6]*var[6] - _error_state[7]*var[7] + var[5];
    _cov[1][5] = _error_state[6]*var[20] - _error_state[7]*var[21] + var[19];
    _cov[2][5] = _error_state[6]*var[26] - _error_state[7]*var[27] + var[25];
    _cov[3][5] = _error_state[6]*var[32] - _error_state[7]*var[35] + var[31];
    _cov[4][5] = _error_state[6]*var[42] - _error_state[7]*var[39] + var[40];
    _cov[0][6] = _error_state[7]*var[9] - var[10]*var[11] + var[12];
    _cov[1][6] = -var[11]*var[23] + var[14]*var[22] + var[24];
    _cov[2][6] = -var[11]*var[29] + var[14]*var[28] + var[30];
    _cov[3][6] = -var[11]*var[37] + var[14]*var[36] + var[38];
    _cov[4][6] = -var[11]*var[44] + var[14]*var[43] + var[45];
    _cov[5][6] = -var[11]*var[47] + var[14]*var[46] + var[48];

    _cov[8][8] = _cov[8][8] + var[13]*(-_cov[7][6]*var[14] + _cov[7][7]*var[13] + _cov[8][7]) - var[14]*(-_cov[6][6]*var[14] + _cov[7][6]*var[13] + _cov[8][6]) - var[52] + var[56];
    _cov[6][6] = -var[11]*var[50] + var[14]*var[49] + var[53];
    _cov[0][7] = -_error_state[6]*var[9] + var[10] + var[11]*var[12];
    _cov[1][7] = var[11]*var[24] - var[13]*var[22] + var[23];
    _cov[2][7] = var[11]*var[30] - var[13]*var[28] + var[29];
    _cov[3][7] = var[11]*var[38] - var[13]*var[36] + var[37];
    _cov[4][7] = var[11]*var[45] - var[13]*var[43] + var[44];
    _cov[5][7] = var[11]*var[48] - var[13]*var[46] + var[47];
    _cov[6][7] = var[11]*var[53] - var[13]*var[49] + var[50];
    _cov[7][7] = var[11]*var[54] - var[13]*var[55] + var[57];
    _cov[0][8] = var[10]*var[13] - var[12]*var[14] + var[8];
    _cov[1][8] = var[13]*var[23] - var[14]*var[24] + var[22];
    _cov[2][8] = var[13]*var[29] - var[14]*var[30] + var[28];
    _cov[3][8] = var[13]*var[37] - var[14]*var[38] + var[36];
    _cov[4][8] = var[13]*var[44] - var[14]*var[45] + var[43];
    _cov[5][8] = var[13]*var[47] - var[14]*var[48] + var[46];
    _cov[6][8] = var[13]*var[50] - var[14]*var[53] + var[49];
    _cov[7][8] = var[13]*var[57] - var[14]*var[54] + var[55];
    _cov[0][9] = _cov[9][0] - _cov[9][1]*_error_state[8] + _cov[9][2]*_error_state[7];
    _cov[1][9] = _cov[9][0]*_error_state[8] + _cov[9][1] - _cov[9][2]*_error_state[6];
    _cov[2][9] = -_cov[9][0]*_error_state[7] + _cov[9][1]*_error_state[6] + _cov[9][2];
    _cov[3][9] = _cov[9][3] - _cov[9][4]*_error_state[8] + _cov[9][5]*_error_state[7];
    _cov[4][9] = _cov[9][3]*_error_state[8] + _cov[9][4] - _cov[9][5]*_error_state[6];
    _cov[5][9] = -_cov[9][3]*_error_state[7] + _cov[9][4]*_error_state[6] + _cov[9][5];
    _cov[6][9] = _cov[9][6] - _cov[9][7]*var[11] + _cov[9][8]*var[14];
    _cov[7][9] = _cov[9][6]*var[11] + _cov[9][7] - _cov[9][8]*var[13];
    _cov[8][9] = -_cov[9][6]*var[14] + _cov[9][7]*var[13] + _cov[9][8];
    // _cov[9][9] = _cov[9][9];
    _cov[0][10] = _cov[10][0] - _cov[10][1]*_error_state[8] + _cov[10][2]*_error_state[7];
    _cov[1][10] = _cov[10][0]*_error_state[8] + _cov[10][1] - _cov[10][2]*_error_state[6];
    _cov[2][10] = -_cov[10][0]*_error_state[7] + _cov[10][1]*_error_state[6] + _cov[10][2];
    _cov[3][10] = _cov[10][3] - _cov[10][4]*_error_state[8] + _cov[10][5]*_error_state[7];
    _cov[4][10] = _cov[10][3]*_error_state[8] + _cov[10][4] - _cov[10][5]*_error_state[6];
    _cov[5][10] = -_cov[10][3]*_error_state[7] + _cov[10][4]*_error_state[6] + _cov[10][5];
    _cov[6][10] = _cov[10][6] - _cov[10][7]*var[11] + _cov[10][8]*var[14];
    _cov[7][10] = _cov[10][6]*var[11] + _cov[10][7] - _cov[10][8]*var[13];
    _cov[8][10] = -_cov[10][6]*var[14] + _cov[10][7]*var[13] + _cov[10][8];
    // _cov[9][10] = _cov[10][9];
    // _cov[10][10] = _cov[10][10];
    _cov[0][11] = _cov[11][0] - _cov[11][1]*_error_state[8] + _cov[11][2]*_error_state[7];
    _cov[1][11] = _cov[11][0]*_error_state[8] + _cov[11][1] - _cov[11][2]*_error_state[6];
    _cov[2][11] = -_cov[11][0]*_error_state[7] + _cov[11][1]*_error_state[6] + _cov[11][2];
    _cov[3][11] = _cov[11][3] - _cov[11][4]*_error_state[8] + _cov[11][5]*_error_state[7];
    _cov[4][11] = _cov[11][3]*_error_state[8] + _cov[11][4] - _cov[11][5]*_error_state[6];
    _cov[5][11] = -_cov[11][3]*_error_state[7] + _cov[11][4]*_error_state[6] + _cov[11][5];
    _cov[6][11] = _cov[11][6] - _cov[11][7]*var[11] + _cov[11][8]*var[14];
    _cov[7][11] = _cov[11][6]*var[11] + _cov[11][7] - _cov[11][8]*var[13];
    _cov[8][11] = -_cov[11][6]*var[14] + _cov[11][7]*var[13] + _cov[11][8];
    // _cov[9][11] = _cov[11][9];
    // _cov[10][11] = _cov[11][10];
    // _cov[11][11] = _cov[11][11];
    _cov[0][12] = _cov[12][0] - _cov[12][1]*_error_state[8] + _cov[12][2]*_error_state[7];
    _cov[1][12] = _cov[12][0]*_error_state[8] + _cov[12][1] - _cov[12][2]*_error_state[6];
    _cov[2][12] = -_cov[12][0]*_error_state[7] + _cov[12][1]*_error_state[6] + _cov[12][2];
    _cov[3][12] = _cov[12][3] - _cov[12][4]*_error_state[8] + _cov[12][5]*_error_state[7];
    _cov[4][12] = _cov[12][3]*_error_state[8] + _cov[12][4] - _cov[12][5]*_error_state[6];
    _cov[5][12] = -_cov[12][3]*_error_state[7] + _cov[12][4]*_error_state[6] + _cov[12][5];
    _cov[6][12] = _cov[12][6] - _cov[12][7]*var[11] + _cov[12][8]*var[14];
    _cov[7][12] = _cov[12][6]*var[11] + _cov[12][7] - _cov[12][8]*var[13];
    _cov[8][12] = -_cov[12][6]*var[14] + _cov[12][7]*var[13] + _cov[12][8];
    // _cov[9][12] = _cov[12][9];
    // _cov[10][12] = _cov[12][10];
    // _cov[11][12] = _cov[12][11];
    // _cov[12][12] = _cov[12][12];
    _cov[0][13] = _cov[13][0] - _cov[13][1]*_error_state[8] + _cov[13][2]*_error_state[7];
    _cov[1][13] = _cov[13][0]*_error_state[8] + _cov[13][1] - _cov[13][2]*_error_state[6];
    _cov[2][13] = -_cov[13][0]*_error_state[7] + _cov[13][1]*_error_state[6] + _cov[13][2];
    _cov[3][13] = _cov[13][3] - _cov[13][4]*_error_state[8] + _cov[13][5]*_error_state[7];
    _cov[4][13] = _cov[13][3]*_error_state[8] + _cov[13][4] - _cov[13][5]*_error_state[6];
    _cov[5][13] = -_cov[13][3]*_error_state[7] + _cov[13][4]*_error_state[6] + _cov[13][5];
    _cov[6][13] = _cov[13][6] - _cov[13][7]*var[11] + _cov[13][8]*var[14];
    _cov[7][13] = _cov[13][6]*var[11] + _cov[13][7] - _cov[13][8]*var[13];
    _cov[8][13] = -_cov[13][6]*var[14] + _cov[13][7]*var[13] + _cov[13][8];
    // _cov[9][13] = _cov[13][9];
    // _cov[10][13] = _cov[13][10];
    // _cov[11][13] = _cov[13][11];
    // _cov[12][13] = _cov[13][12];
    // _cov[13][13] = _cov[13][13];
    _cov[0][14] = _cov[14][0] - _cov[14][1]*_error_state[8] + _cov[14][2]*_error_state[7];
    _cov[1][14] = _cov[14][0]*_error_state[8] + _cov[14][1] - _cov[14][2]*_error_state[6];
    _cov[2][14] = -_cov[14][0]*_error_state[7] + _cov[14][1]*_error_state[6] + _cov[14][2];
    _cov[3][14] = _cov[14][3] - _cov[14][4]*_error_state[8] + _cov[14][5]*_error_state[7];
    _cov[4][14] = _cov[14][3]*_error_state[8] + _cov[14][4] - _cov[14][5]*_error_state[6];
    _cov[5][14] = -_cov[14][3]*_error_state[7] + _cov[14][4]*_error_state[6] + _cov[14][5];
    _cov[6][14] = _cov[14][6] - _cov[14][7]*var[11] + _cov[14][8]*var[14];
    _cov[7][14] = _cov[14][6]*var[11] + _cov[14][7] - _cov[14][8]*var[13];
    _cov[8][14] = -_cov[14][6]*var[14] + _cov[14][7]*var[13] + _cov[14][8];
    // _cov[9][14] = _cov[14][9];
    // _cov[10][14] = _cov[14][10];
    // _cov[11][14] = _cov[14][11];
    // _cov[12][14] = _cov[14][12];
    // _cov[13][14] = _cov[14][13];
    // _cov[14][14] = _cov[14][14];
    _cov[0][15] = _cov[15][0] - _cov[15][1]*_error_state[8] + _cov[15][2]*_error_state[7];
    _cov[1][15] = _cov[15][0]*_error_state[8] + _cov[15][1] - _cov[15][2]*_error_state[6];
    _cov[2][15] = -_cov[15][0]*_error_state[7] + _cov[15][1]*_error_state[6] + _cov[15][2];
    _cov[3][15] = _cov[15][3] - _cov[15][4]*_error_state[8] + _cov[15][5]*_error_state[7];
    _cov[4][15] = _cov[15][3]*_error_state[8] + _cov[15][4] - _cov[15][5]*_error_state[6];
    _cov[5][15] = -_cov[15][3]*_error_state[7] + _cov[15][4]*_error_state[6] + _cov[15][5];
    _cov[6][15] = _cov[15][6] - _cov[15][7]*var[11] + _cov[15][8]*var[14];
    _cov[7][15] = _cov[15][6]*var[11] + _cov[15][7] - _cov[15][8]*var[13];
    _cov[8][15] = -_cov[15][6]*var[14] + _cov[15][7]*var[13] + _cov[15][8];
    // _cov[9][15] = _cov[15][9];
    // _cov[10][15] = _cov[15][10];
    // _cov[11][15] = _cov[15][11];
    // _cov[12][15] = _cov[15][12];
    // _cov[13][15] = _cov[15][13];
    // _cov[14][15] = _cov[15][14];
    // _cov[15][15] = _cov[15][15];

    for (unsigned char i = 0; i < 9; ++i) {
        for (unsigned char j = i + 1; j < 16; ++j) {
            _cov[j][i] = _cov[i][j];
        }
    }
}