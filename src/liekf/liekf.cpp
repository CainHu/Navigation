//
// Created by Cain on 2022/11/11.
//

#include "liekf.h"
#include <cfloat>
#include <iostream>

using namespace std;
using namespace liekf;

void LIEKF::predict_covariance(const Vector3f &w, const Vector3f &a) {
    const Vector3f dv_corr = (a - _ba) * _dt;
    const Vector3f da_corr = (w - _bg) * _dt;

    const array<array<float, 3>, 3> x {
        {{0.f, dv_corr[2], -dv_corr[1]},
         {-dv_corr[2], 0.f, dv_corr[0]},
         {dv_corr[1], -dv_corr[0], 0.f}}
    };

    const array<array<float, 3>, 3> y {
        {{1.f, da_corr[2], -da_corr[1]},
         {-da_corr[2], 1.f, da_corr[0]},
         {da_corr[1], -da_corr[0], 1.f}}
    };

    const array<float, 3> r {
        _rot(2, 0) * _dt, _rot(2, 1) * _dt, _rot(2, 2) * _dt
    };

    array<float, 117> var;
    // Equations for covariance matrix prediction, without process noise!
    var[0] = _cov[3][0]*y[0][0] + _cov[3][1]*y[0][1] + _cov[3][2]*y[0][2] + _cov[3][3]*_dt;
    var[1] = _cov[0][0]*y[0][0] + _cov[1][0]*y[0][1] + _cov[2][0]*y[0][2] + _cov[3][0]*_dt;
    var[2] = _cov[1][0]*y[0][0] + _cov[1][1]*y[0][1] + _cov[2][1]*y[0][2] + _cov[3][1]*_dt;
    var[3] = _cov[2][0]*y[0][0] + _cov[2][1]*y[0][1] + _cov[2][2]*y[0][2] + _cov[3][2]*_dt;
    var[4] = _cov[4][3]*_dt;
    var[5] = _cov[4][0]*y[0][0] + _cov[4][1]*y[0][1] + _cov[4][2]*y[0][2] + var[4];
    var[6] = _cov[5][3]*_dt;
    var[7] = _cov[5][0]*y[0][0] + _cov[5][1]*y[0][1] + _cov[5][2]*y[0][2] + var[6];
    var[8] = _cov[15][0]*y[0][0] + _cov[15][1]*y[0][1] + _cov[15][2]*y[0][2] + _cov[15][3]*_dt;
    var[9] = _cov[6][0]*y[0][0] + _cov[6][1]*y[0][1] + _cov[6][2]*y[0][2] + _cov[6][3]*_dt;
    var[10] = _cov[7][0]*y[0][0] + _cov[7][1]*y[0][1] + _cov[7][2]*y[0][2] + _cov[7][3]*_dt;
    var[11] = _cov[8][0]*y[0][0] + _cov[8][1]*y[0][1] + _cov[8][2]*y[0][2] + _cov[8][3]*_dt;
    var[12] = _cov[12][3]*_dt;
    var[13] = _cov[12][0]*y[0][0] + _cov[12][1]*y[0][1] + _cov[12][2]*y[0][2] + var[12];
    var[14] = _cov[13][3]*_dt;
    var[15] = _cov[13][0]*y[0][0] + _cov[13][1]*y[0][1] + _cov[13][2]*y[0][2] + var[14];
    var[16] = _cov[14][3]*_dt;
    var[17] = _cov[14][0]*y[0][0] + _cov[14][1]*y[0][1] + _cov[14][2]*y[0][2] + var[16];
    var[18] = _cov[9][0]*y[0][0] + _cov[9][1]*y[0][1] + _cov[9][2]*y[0][2] + _cov[9][3]*_dt;
    var[19] = _cov[10][0]*y[0][0] + _cov[10][1]*y[0][1] + _cov[10][2]*y[0][2] + _cov[10][3]*_dt;
    var[20] = _cov[11][0]*y[0][0] + _cov[11][1]*y[0][1] + _cov[11][2]*y[0][2] + _cov[11][3]*_dt;
    var[21] = _cov[4][0]*y[1][0] + _cov[4][1]*y[1][1] + _cov[4][2]*y[1][2] + _cov[4][4]*_dt;
    var[22] = _cov[0][0]*y[1][0] + _cov[1][0]*y[1][1] + _cov[2][0]*y[1][2] + _cov[4][0]*_dt;
    var[23] = _cov[1][0]*y[1][0] + _cov[1][1]*y[1][1] + _cov[2][1]*y[1][2] + _cov[4][1]*_dt;
    var[24] = _cov[2][0]*y[1][0] + _cov[2][1]*y[1][1] + _cov[2][2]*y[1][2] + _cov[4][2]*_dt;
    var[25] = _cov[5][4]*_dt;
    var[26] = _cov[5][0]*y[1][0] + _cov[5][1]*y[1][1] + _cov[5][2]*y[1][2] + var[25];
    var[27] = _cov[15][0]*y[1][0] + _cov[15][1]*y[1][1] + _cov[15][2]*y[1][2] + _cov[15][4]*_dt;
    var[28] = _cov[6][0]*y[1][0] + _cov[6][1]*y[1][1] + _cov[6][2]*y[1][2] + _cov[6][4]*_dt;
    var[29] = _cov[7][0]*y[1][0] + _cov[7][1]*y[1][1] + _cov[7][2]*y[1][2] + _cov[7][4]*_dt;
    var[30] = _cov[8][0]*y[1][0] + _cov[8][1]*y[1][1] + _cov[8][2]*y[1][2] + _cov[8][4]*_dt;
    var[31] = _cov[3][0]*y[1][0] + _cov[3][1]*y[1][1] + _cov[3][2]*y[1][2] + var[4];
    var[32] = _cov[12][4]*_dt;
    var[33] = _cov[12][0]*y[1][0] + _cov[12][1]*y[1][1] + _cov[12][2]*y[1][2] + var[32];
    var[34] = _cov[13][4]*_dt;
    var[35] = _cov[13][0]*y[1][0] + _cov[13][1]*y[1][1] + _cov[13][2]*y[1][2] + var[34];
    var[36] = _cov[14][4]*_dt;
    var[37] = _cov[14][0]*y[1][0] + _cov[14][1]*y[1][1] + _cov[14][2]*y[1][2] + var[36];
    var[38] = _cov[9][0]*y[1][0] + _cov[9][1]*y[1][1] + _cov[9][2]*y[1][2] + _cov[9][4]*_dt;
    var[39] = _cov[10][0]*y[1][0] + _cov[10][1]*y[1][1] + _cov[10][2]*y[1][2] + _cov[10][4]*_dt;
    var[40] = _cov[11][0]*y[1][0] + _cov[11][1]*y[1][1] + _cov[11][2]*y[1][2] + _cov[11][4]*_dt;
    var[41] = _cov[5][0]*y[2][0] + _cov[5][1]*y[2][1] + _cov[5][2]*y[2][2] + _cov[5][5]*_dt;
    var[42] = _cov[15][0]*y[2][0] + _cov[15][1]*y[2][1] + _cov[15][2]*y[2][2] + _cov[15][5]*_dt;
    var[43] = _cov[6][0]*y[2][0] + _cov[6][1]*y[2][1] + _cov[6][2]*y[2][2] + _cov[6][5]*_dt;
    var[44] = _cov[7][0]*y[2][0] + _cov[7][1]*y[2][1] + _cov[7][2]*y[2][2] + _cov[7][5]*_dt;
    var[45] = _cov[8][0]*y[2][0] + _cov[8][1]*y[2][1] + _cov[8][2]*y[2][2] + _cov[8][5]*_dt;
    var[46] = _cov[3][0]*y[2][0] + _cov[3][1]*y[2][1] + _cov[3][2]*y[2][2] + var[6];
    var[47] = _cov[4][0]*y[2][0] + _cov[4][1]*y[2][1] + _cov[4][2]*y[2][2] + var[25];
    var[48] = _cov[12][5]*_dt;
    var[49] = _cov[12][0]*y[2][0] + _cov[12][1]*y[2][1] + _cov[12][2]*y[2][2] + var[48];
    var[50] = _cov[13][5]*_dt;
    var[51] = _cov[13][0]*y[2][0] + _cov[13][1]*y[2][1] + _cov[13][2]*y[2][2] + var[50];
    var[52] = _cov[14][5]*_dt;
    var[53] = _cov[14][0]*y[2][0] + _cov[14][1]*y[2][1] + _cov[14][2]*y[2][2] + var[52];
    var[54] = _cov[9][0]*y[2][0] + _cov[9][1]*y[2][1] + _cov[9][2]*y[2][2] + _cov[9][5]*_dt;
    var[55] = _cov[10][0]*y[2][0] + _cov[10][1]*y[2][1] + _cov[10][2]*y[2][2] + _cov[10][5]*_dt;
    var[56] = _cov[11][0]*y[2][0] + _cov[11][1]*y[2][1] + _cov[11][2]*y[2][2] + _cov[11][5]*_dt;
    var[57] = -_cov[15][12]*_dt + _cov[15][15]*r[0] + _cov[15][3]*y[0][0] + _cov[15][4]*y[0][1] + _cov[15][5]*y[0][2] + _cov[15][6]*x[0][0] + _cov[15][7]*x[0][1] + _cov[15][8]*x[0][2];
    var[58] = -_cov[12][6]*_dt + _cov[15][6]*r[0] + _cov[6][3]*y[0][0] + _cov[6][4]*y[0][1] + _cov[6][5]*y[0][2] + _cov[6][6]*x[0][0] + _cov[7][6]*x[0][1] + _cov[8][6]*x[0][2];
    var[59] = -_cov[12][7]*_dt + _cov[15][7]*r[0] + _cov[7][3]*y[0][0] + _cov[7][4]*y[0][1] + _cov[7][5]*y[0][2] + _cov[7][6]*x[0][0] + _cov[7][7]*x[0][1] + _cov[8][7]*x[0][2];
    var[60] = -_cov[12][8]*_dt + _cov[15][8]*r[0] + _cov[8][3]*y[0][0] + _cov[8][4]*y[0][1] + _cov[8][5]*y[0][2] + _cov[8][6]*x[0][0] + _cov[8][7]*x[0][1] + _cov[8][8]*x[0][2];
    var[61] = _cov[15][3]*r[0] + _cov[3][3]*y[0][0] + _cov[4][3]*y[0][1] + _cov[5][3]*y[0][2] + _cov[6][3]*x[0][0] + _cov[7][3]*x[0][1] + _cov[8][3]*x[0][2] - var[12];
    var[62] = _cov[15][4]*r[0] + _cov[4][3]*y[0][0] + _cov[4][4]*y[0][1] + _cov[5][4]*y[0][2] + _cov[6][4]*x[0][0] + _cov[7][4]*x[0][1] + _cov[8][4]*x[0][2] - var[32];
    var[63] = _cov[15][5]*r[0] + _cov[5][3]*y[0][0] + _cov[5][4]*y[0][1] + _cov[5][5]*y[0][2] + _cov[6][5]*x[0][0] + _cov[7][5]*x[0][1] + _cov[8][5]*x[0][2] - var[48];
    var[64] = -_cov[12][12]*_dt + _cov[12][3]*y[0][0] + _cov[12][4]*y[0][1] + _cov[12][5]*y[0][2] + _cov[12][6]*x[0][0] + _cov[12][7]*x[0][1] + _cov[12][8]*x[0][2] + _cov[15][12]*r[0];
    var[65] = -_cov[13][12]*_dt;
    var[66] = _cov[13][3]*y[0][0] + _cov[13][4]*y[0][1] + _cov[13][5]*y[0][2] + _cov[13][6]*x[0][0] + _cov[13][7]*x[0][1] + _cov[13][8]*x[0][2] + _cov[15][13]*r[0] + var[65];
    var[67] = -_cov[14][12]*_dt;
    var[68] = _cov[14][3]*y[0][0] + _cov[14][4]*y[0][1] + _cov[14][5]*y[0][2] + _cov[14][6]*x[0][0] + _cov[14][7]*x[0][1] + _cov[14][8]*x[0][2] + _cov[15][14]*r[0] + var[67];
    var[69] = -_cov[12][9]*_dt;
    var[70] = _cov[15][9]*r[0] + _cov[9][3]*y[0][0] + _cov[9][4]*y[0][1] + _cov[9][5]*y[0][2] + _cov[9][6]*x[0][0] + _cov[9][7]*x[0][1] + _cov[9][8]*x[0][2] + var[69];
    var[71] = -_cov[12][10]*_dt;
    var[72] = _cov[10][3]*y[0][0] + _cov[10][4]*y[0][1] + _cov[10][5]*y[0][2] + _cov[10][6]*x[0][0] + _cov[10][7]*x[0][1] + _cov[10][8]*x[0][2] + _cov[15][10]*r[0] + var[71];
    var[73] = -_cov[12][11]*_dt;
    var[74] = _cov[11][3]*y[0][0] + _cov[11][4]*y[0][1] + _cov[11][5]*y[0][2] + _cov[11][6]*x[0][0] + _cov[11][7]*x[0][1] + _cov[11][8]*x[0][2] + _cov[15][11]*r[0] + var[73];
    var[75] = -_cov[15][13]*_dt + _cov[15][15]*r[1] + _cov[15][3]*y[1][0] + _cov[15][4]*y[1][1] + _cov[15][5]*y[1][2] + _cov[15][6]*x[1][0] + _cov[15][7]*x[1][1] + _cov[15][8]*x[1][2];
    var[76] = -_cov[13][6]*_dt + _cov[15][6]*r[1] + _cov[6][3]*y[1][0] + _cov[6][4]*y[1][1] + _cov[6][5]*y[1][2] + _cov[6][6]*x[1][0] + _cov[7][6]*x[1][1] + _cov[8][6]*x[1][2];
    var[77] = -_cov[13][7]*_dt + _cov[15][7]*r[1] + _cov[7][3]*y[1][0] + _cov[7][4]*y[1][1] + _cov[7][5]*y[1][2] + _cov[7][6]*x[1][0] + _cov[7][7]*x[1][1] + _cov[8][7]*x[1][2];
    var[78] = -_cov[13][8]*_dt + _cov[15][8]*r[1] + _cov[8][3]*y[1][0] + _cov[8][4]*y[1][1] + _cov[8][5]*y[1][2] + _cov[8][6]*x[1][0] + _cov[8][7]*x[1][1] + _cov[8][8]*x[1][2];
    var[79] = _cov[15][3]*r[1] + _cov[3][3]*y[1][0] + _cov[4][3]*y[1][1] + _cov[5][3]*y[1][2] + _cov[6][3]*x[1][0] + _cov[7][3]*x[1][1] + _cov[8][3]*x[1][2] - var[14];
    var[80] = _cov[15][4]*r[1] + _cov[4][3]*y[1][0] + _cov[4][4]*y[1][1] + _cov[5][4]*y[1][2] + _cov[6][4]*x[1][0] + _cov[7][4]*x[1][1] + _cov[8][4]*x[1][2] - var[34];
    var[81] = _cov[15][5]*r[1] + _cov[5][3]*y[1][0] + _cov[5][4]*y[1][1] + _cov[5][5]*y[1][2] + _cov[6][5]*x[1][0] + _cov[7][5]*x[1][1] + _cov[8][5]*x[1][2] - var[50];
    var[82] = -_cov[13][13]*_dt + _cov[13][3]*y[1][0] + _cov[13][4]*y[1][1] + _cov[13][5]*y[1][2] + _cov[13][6]*x[1][0] + _cov[13][7]*x[1][1] + _cov[13][8]*x[1][2] + _cov[15][13]*r[1];
    var[83] = -_cov[14][13]*_dt;
    var[84] = _cov[14][3]*y[1][0] + _cov[14][4]*y[1][1] + _cov[14][5]*y[1][2] + _cov[14][6]*x[1][0] + _cov[14][7]*x[1][1] + _cov[14][8]*x[1][2] + _cov[15][14]*r[1] + var[83];
    var[85] = -_cov[13][9]*_dt;
    var[86] = _cov[15][9]*r[1] + _cov[9][3]*y[1][0] + _cov[9][4]*y[1][1] + _cov[9][5]*y[1][2] + _cov[9][6]*x[1][0] + _cov[9][7]*x[1][1] + _cov[9][8]*x[1][2] + var[85];
    var[87] = -_cov[13][10]*_dt;
    var[88] = _cov[10][3]*y[1][0] + _cov[10][4]*y[1][1] + _cov[10][5]*y[1][2] + _cov[10][6]*x[1][0] + _cov[10][7]*x[1][1] + _cov[10][8]*x[1][2] + _cov[15][10]*r[1] + var[87];
    var[89] = -_cov[13][11]*_dt;
    var[90] = _cov[11][3]*y[1][0] + _cov[11][4]*y[1][1] + _cov[11][5]*y[1][2] + _cov[11][6]*x[1][0] + _cov[11][7]*x[1][1] + _cov[11][8]*x[1][2] + _cov[15][11]*r[1] + var[89];
    var[91] = -_cov[15][14]*_dt + _cov[15][15]*r[2] + _cov[15][3]*y[2][0] + _cov[15][4]*y[2][1] + _cov[15][5]*y[2][2] + _cov[15][6]*x[2][0] + _cov[15][7]*x[2][1] + _cov[15][8]*x[2][2];
    var[92] = -_cov[14][6]*_dt + _cov[15][6]*r[2] + _cov[6][3]*y[2][0] + _cov[6][4]*y[2][1] + _cov[6][5]*y[2][2] + _cov[6][6]*x[2][0] + _cov[7][6]*x[2][1] + _cov[8][6]*x[2][2];
    var[93] = -_cov[14][7]*_dt + _cov[15][7]*r[2] + _cov[7][3]*y[2][0] + _cov[7][4]*y[2][1] + _cov[7][5]*y[2][2] + _cov[7][6]*x[2][0] + _cov[7][7]*x[2][1] + _cov[8][7]*x[2][2];
    var[94] = -_cov[14][8]*_dt + _cov[15][8]*r[2] + _cov[8][3]*y[2][0] + _cov[8][4]*y[2][1] + _cov[8][5]*y[2][2] + _cov[8][6]*x[2][0] + _cov[8][7]*x[2][1] + _cov[8][8]*x[2][2];
    var[95] = -_cov[14][14]*_dt + _cov[14][3]*y[2][0] + _cov[14][4]*y[2][1] + _cov[14][5]*y[2][2] + _cov[14][6]*x[2][0] + _cov[14][7]*x[2][1] + _cov[14][8]*x[2][2] + _cov[15][14]*r[2];
    var[96] = -_cov[14][9]*_dt;
    var[97] = _cov[15][9]*r[2] + _cov[9][3]*y[2][0] + _cov[9][4]*y[2][1] + _cov[9][5]*y[2][2] + _cov[9][6]*x[2][0] + _cov[9][7]*x[2][1] + _cov[9][8]*x[2][2] + var[96];
    var[98] = -_cov[14][10]*_dt;
    var[99] = _cov[10][3]*y[2][0] + _cov[10][4]*y[2][1] + _cov[10][5]*y[2][2] + _cov[10][6]*x[2][0] + _cov[10][7]*x[2][1] + _cov[10][8]*x[2][2] + _cov[15][10]*r[2] + var[98];
    var[100] = -_cov[14][11]*_dt;
    var[101] = _cov[11][3]*y[2][0] + _cov[11][4]*y[2][1] + _cov[11][5]*y[2][2] + _cov[11][6]*x[2][0] + _cov[11][7]*x[2][1] + _cov[11][8]*x[2][2] + _cov[15][11]*r[2] + var[100];
    var[102] = _cov[6][6]*y[0][0] + _cov[7][6]*y[0][1] + _cov[8][6]*y[0][2] - _cov[9][6]*_dt;
    var[103] = _cov[7][6]*y[0][0] + _cov[7][7]*y[0][1] + _cov[8][7]*y[0][2] - _cov[9][7]*_dt;
    var[104] = _cov[8][6]*y[0][0] + _cov[8][7]*y[0][1] + _cov[8][8]*y[0][2] - _cov[9][8]*_dt;
    var[105] = _cov[9][6]*y[0][0] + _cov[9][7]*y[0][1] + _cov[9][8]*y[0][2] - _cov[9][9]*_dt;
    var[106] = -_cov[10][9]*_dt;
    var[107] = _cov[10][6]*y[0][0] + _cov[10][7]*y[0][1] + _cov[10][8]*y[0][2] + var[106];
    var[108] = -_cov[11][9]*_dt;
    var[109] = _cov[11][6]*y[0][0] + _cov[11][7]*y[0][1] + _cov[11][8]*y[0][2] + var[108];
    var[110] = -_cov[10][6]*_dt + _cov[6][6]*y[1][0] + _cov[7][6]*y[1][1] + _cov[8][6]*y[1][2];
    var[111] = -_cov[10][7]*_dt + _cov[7][6]*y[1][0] + _cov[7][7]*y[1][1] + _cov[8][7]*y[1][2];
    var[112] = -_cov[10][8]*_dt + _cov[8][6]*y[1][0] + _cov[8][7]*y[1][1] + _cov[8][8]*y[1][2];
    var[113] = -_cov[10][10]*_dt + _cov[10][6]*y[1][0] + _cov[10][7]*y[1][1] + _cov[10][8]*y[1][2];
    var[114] = -_cov[11][10]*_dt;
    var[115] = _cov[11][6]*y[1][0] + _cov[11][7]*y[1][1] + _cov[11][8]*y[1][2] + var[114];
    var[116] = -_cov[11][11]*_dt + _cov[11][6]*y[2][0] + _cov[11][7]*y[2][1] + _cov[11][8]*y[2][2];

    _cov[2][2] = _dt*var[41] + y[2][0]*(_cov[0][0]*y[2][0] + _cov[1][0]*y[2][1] + _cov[2][0]*y[2][2] + _cov[5][0]*_dt) + y[2][1]*(_cov[1][0]*y[2][0] + _cov[1][1]*y[2][1] + _cov[2][1]*y[2][2] + _cov[5][1]*_dt) + y[2][2]*(_cov[2][0]*y[2][0] + _cov[2][1]*y[2][1] + _cov[2][2]*y[2][2] + _cov[5][2]*_dt);
    _cov[0][0] = _dt*var[0] + var[1]*y[0][0] + var[2]*y[0][1] + var[3]*y[0][2];
    _cov[0][1] = _dt*var[5] + var[1]*y[1][0] + var[2]*y[1][1] + var[3]*y[1][2];
    _cov[1][1] = _dt*var[21] + var[22]*y[1][0] + var[23]*y[1][1] + var[24]*y[1][2];
    _cov[0][2] = _dt*var[7] + var[1]*y[2][0] + var[2]*y[2][1] + var[3]*y[2][2];
    _cov[1][2] = _dt*var[26] + var[22]*y[2][0] + var[23]*y[2][1] + var[24]*y[2][2];

    _cov[5][5] = -_dt*var[95] + r[2]*var[91] + var[92]*x[2][0] + var[93]*x[2][1] + var[94]*x[2][2] + y[2][0]*(_cov[15][3]*r[2] + _cov[3][3]*y[2][0] + _cov[4][3]*y[2][1] + _cov[5][3]*y[2][2] + _cov[6][3]*x[2][0] + _cov[7][3]*x[2][1] + _cov[8][3]*x[2][2] - var[16]) + y[2][1]*(_cov[15][4]*r[2] + _cov[4][3]*y[2][0] + _cov[4][4]*y[2][1] + _cov[5][4]*y[2][2] + _cov[6][4]*x[2][0] + _cov[7][4]*x[2][1] + _cov[8][4]*x[2][2] - var[36]) + y[2][2]*(_cov[15][5]*r[2] + _cov[5][3]*y[2][0] + _cov[5][4]*y[2][1] + _cov[5][5]*y[2][2] + _cov[6][5]*x[2][0] + _cov[7][5]*x[2][1] + _cov[8][5]*x[2][2] - var[52]);
    _cov[0][3] = -_dt*var[13] + r[0]*var[8] + var[0]*y[0][0] + var[10]*x[0][1] + var[11]*x[0][2] + var[5]*y[0][1] + var[7]*y[0][2] + var[9]*x[0][0];
    _cov[1][3] = -_dt*var[33] + r[0]*var[27] + var[21]*y[0][1] + var[26]*y[0][2] + var[28]*x[0][0] + var[29]*x[0][1] + var[30]*x[0][2] + var[31]*y[0][0];
    _cov[2][3] = -_dt*var[49] + r[0]*var[42] + var[41]*y[0][2] + var[43]*x[0][0] + var[44]*x[0][1] + var[45]*x[0][2] + var[46]*y[0][0] + var[47]*y[0][1];
    _cov[3][3] = -_dt*var[64] + r[0]*var[57] + var[58]*x[0][0] + var[59]*x[0][1] + var[60]*x[0][2] + var[61]*y[0][0] + var[62]*y[0][1] + var[63]*y[0][2];
    _cov[0][4] = -_dt*var[15] + r[1]*var[8] + var[0]*y[1][0] + var[10]*x[1][1] + var[11]*x[1][2] + var[5]*y[1][1] + var[7]*y[1][2] + var[9]*x[1][0];
    _cov[1][4] = -_dt*var[35] + r[1]*var[27] + var[21]*y[1][1] + var[26]*y[1][2] + var[28]*x[1][0] + var[29]*x[1][1] + var[30]*x[1][2] + var[31]*y[1][0];
    _cov[2][4] = -_dt*var[51] + r[1]*var[42] + var[41]*y[1][2] + var[43]*x[1][0] + var[44]*x[1][1] + var[45]*x[1][2] + var[46]*y[1][0] + var[47]*y[1][1];
    _cov[3][4] = -_dt*var[66] + r[1]*var[57] + var[58]*x[1][0] + var[59]*x[1][1] + var[60]*x[1][2] + var[61]*y[1][0] + var[62]*y[1][1] + var[63]*y[1][2];
    _cov[4][4] = -_dt*var[82] + r[1]*var[75] + var[76]*x[1][0] + var[77]*x[1][1] + var[78]*x[1][2] + var[79]*y[1][0] + var[80]*y[1][1] + var[81]*y[1][2];
    _cov[0][5] = -_dt*var[17] + r[2]*var[8] + var[0]*y[2][0] + var[10]*x[2][1] + var[11]*x[2][2] + var[5]*y[2][1] + var[7]*y[2][2] + var[9]*x[2][0];
    _cov[1][5] = -_dt*var[37] + r[2]*var[27] + var[21]*y[2][1] + var[26]*y[2][2] + var[28]*x[2][0] + var[29]*x[2][1] + var[30]*x[2][2] + var[31]*y[2][0];
    _cov[2][5] = -_dt*var[53] + r[2]*var[42] + var[41]*y[2][2] + var[43]*x[2][0] + var[44]*x[2][1] + var[45]*x[2][2] + var[46]*y[2][0] + var[47]*y[2][1];
    _cov[3][5] = -_dt*var[68] + r[2]*var[57] + var[58]*x[2][0] + var[59]*x[2][1] + var[60]*x[2][2] + var[61]*y[2][0] + var[62]*y[2][1] + var[63]*y[2][2];
    _cov[4][5] = -_dt*var[84] + r[2]*var[75] + var[76]*x[2][0] + var[77]*x[2][1] + var[78]*x[2][2] + var[79]*y[2][0] + var[80]*y[2][1] + var[81]*y[2][2];

    _cov[8][8] = -_dt*var[116] + y[2][0]*(-_cov[11][6]*_dt + _cov[6][6]*y[2][0] + _cov[7][6]*y[2][1] + _cov[8][6]*y[2][2]) + y[2][1]*(-_cov[11][7]*_dt + _cov[7][6]*y[2][0] + _cov[7][7]*y[2][1] + _cov[8][7]*y[2][2]) + y[2][2]*(-_cov[11][8]*_dt + _cov[8][6]*y[2][0] + _cov[8][7]*y[2][1] + _cov[8][8]*y[2][2]);
    _cov[0][6] = -_dt*var[18] + var[10]*y[0][1] + var[11]*y[0][2] + var[9]*y[0][0];
    _cov[1][6] = -_dt*var[38] + var[28]*y[0][0] + var[29]*y[0][1] + var[30]*y[0][2];
    _cov[2][6] = -_dt*var[54] + var[43]*y[0][0] + var[44]*y[0][1] + var[45]*y[0][2];
    _cov[3][6] = -_dt*var[70] + var[58]*y[0][0] + var[59]*y[0][1] + var[60]*y[0][2];
    _cov[4][6] = -_dt*var[86] + var[76]*y[0][0] + var[77]*y[0][1] + var[78]*y[0][2];
    _cov[5][6] = -_dt*var[97] + var[92]*y[0][0] + var[93]*y[0][1] + var[94]*y[0][2];
    _cov[6][6] = -_dt*var[105] + var[102]*y[0][0] + var[103]*y[0][1] + var[104]*y[0][2];
    _cov[0][7] = -_dt*var[19] + var[10]*y[1][1] + var[11]*y[1][2] + var[9]*y[1][0];
    _cov[1][7] = -_dt*var[39] + var[28]*y[1][0] + var[29]*y[1][1] + var[30]*y[1][2];
    _cov[2][7] = -_dt*var[55] + var[43]*y[1][0] + var[44]*y[1][1] + var[45]*y[1][2];
    _cov[3][7] = -_dt*var[72] + var[58]*y[1][0] + var[59]*y[1][1] + var[60]*y[1][2];
    _cov[4][7] = -_dt*var[88] + var[76]*y[1][0] + var[77]*y[1][1] + var[78]*y[1][2];
    _cov[5][7] = -_dt*var[99] + var[92]*y[1][0] + var[93]*y[1][1] + var[94]*y[1][2];
    _cov[6][7] = -_dt*var[107] + var[102]*y[1][0] + var[103]*y[1][1] + var[104]*y[1][2];
    _cov[7][7] = -_dt*var[113] + var[110]*y[1][0] + var[111]*y[1][1] + var[112]*y[1][2];
    _cov[0][8] = -_dt*var[20] + var[10]*y[2][1] + var[11]*y[2][2] + var[9]*y[2][0];
    _cov[1][8] = -_dt*var[40] + var[28]*y[2][0] + var[29]*y[2][1] + var[30]*y[2][2];
    _cov[2][8] = -_dt*var[56] + var[43]*y[2][0] + var[44]*y[2][1] + var[45]*y[2][2];
    _cov[3][8] = -_dt*var[74] + var[58]*y[2][0] + var[59]*y[2][1] + var[60]*y[2][2];
    _cov[4][8] = -_dt*var[90] + var[76]*y[2][0] + var[77]*y[2][1] + var[78]*y[2][2];
    _cov[5][8] = -_dt*var[101] + var[92]*y[2][0] + var[93]*y[2][1] + var[94]*y[2][2];
    _cov[6][8] = -_dt*var[109] + var[102]*y[2][0] + var[103]*y[2][1] + var[104]*y[2][2];
    _cov[7][8] = -_dt*var[115] + var[110]*y[2][0] + var[111]*y[2][1] + var[112]*y[2][2];

    _cov[0][9] = var[18];
    _cov[1][9] = var[38];
    _cov[2][9] = var[54];
    _cov[3][9] = var[70];
    _cov[4][9] = var[86];
    _cov[5][9] = var[97];
    _cov[6][9] = var[105];
    _cov[7][9] = _cov[9][6]*y[1][0] + _cov[9][7]*y[1][1] + _cov[9][8]*y[1][2] + var[106];
    _cov[8][9] = _cov[9][6]*y[2][0] + _cov[9][7]*y[2][1] + _cov[9][8]*y[2][2] + var[108];

    _cov[0][10] = var[19];
    _cov[1][10] = var[39];
    _cov[2][10] = var[55];
    _cov[3][10] = var[72];
    _cov[4][10] = var[88];
    _cov[5][10] = var[99];
    _cov[6][10] = var[107];
    _cov[7][10] = var[113];
    _cov[8][10] = _cov[10][6]*y[2][0] + _cov[10][7]*y[2][1] + _cov[10][8]*y[2][2] + var[114];

    _cov[0][11] = var[20];
    _cov[1][11] = var[40];
    _cov[2][11] = var[56];
    _cov[3][11] = var[74];
    _cov[4][11] = var[90];
    _cov[5][11] = var[101];
    _cov[6][11] = var[109];
    _cov[7][11] = var[115];
    _cov[8][11] = var[116];

    // F *P *F' + Q
    add_processing_covariance<12>(0);

    if (_control_status.flags.acc_x_bias) {
        _cov[0][12] = var[13];
        _cov[1][12] = var[33];
        _cov[2][12] = var[49];
        _cov[3][12] = var[64];
        _cov[4][12] = _cov[12][3]*y[1][0] + _cov[12][4]*y[1][1] + _cov[12][5]*y[1][2] + _cov[12][6]*x[1][0] + _cov[12][7]*x[1][1] + _cov[12][8]*x[1][2] + _cov[15][12]*r[1] + var[65];
        _cov[5][12] = _cov[12][3]*y[2][0] + _cov[12][4]*y[2][1] + _cov[12][5]*y[2][2] + _cov[12][6]*x[2][0] + _cov[12][7]*x[2][1] + _cov[12][8]*x[2][2] + _cov[15][12]*r[2] + var[67];
        _cov[6][12] = _cov[12][6]*y[0][0] + _cov[12][7]*y[0][1] + _cov[12][8]*y[0][2] + var[69];
        _cov[7][12] = _cov[12][6]*y[1][0] + _cov[12][7]*y[1][1] + _cov[12][8]*y[1][2] + var[71];
        _cov[8][12] = _cov[12][6]*y[2][0] + _cov[12][7]*y[2][1] + _cov[12][8]*y[2][2] + var[73];

        _cov[12][12] = kahan_summation(_cov[12][12], _q_cov[12] * _dt2, _accumulator_cov[12]);
    }

    if (_control_status.flags.acc_y_bias) {
        _cov[0][13] = var[15];
        _cov[1][13] = var[35];
        _cov[2][13] = var[51];
        _cov[3][13] = var[66];
        _cov[4][13] = var[82];
        _cov[5][13] = _cov[13][3]*y[2][0] + _cov[13][4]*y[2][1] + _cov[13][5]*y[2][2] + _cov[13][6]*x[2][0] + _cov[13][7]*x[2][1] + _cov[13][8]*x[2][2] + _cov[15][13]*r[2] + var[83];
        _cov[6][13] = _cov[13][6]*y[0][0] + _cov[13][7]*y[0][1] + _cov[13][8]*y[0][2] + var[85];
        _cov[7][13] = _cov[13][6]*y[1][0] + _cov[13][7]*y[1][1] + _cov[13][8]*y[1][2] + var[87];
        _cov[8][13] = _cov[13][6]*y[2][0] + _cov[13][7]*y[2][1] + _cov[13][8]*y[2][2] + var[89];

        _cov[13][13] = kahan_summation(_cov[13][13], _q_cov[13] * _dt2, _accumulator_cov[13]);
    }

    if (_control_status.flags.acc_z_bias) {
        _cov[0][14] = var[17];
        _cov[1][14] = var[37];
        _cov[2][14] = var[53];
        _cov[3][14] = var[68];
        _cov[4][14] = var[84];
        _cov[5][14] = var[95];
        _cov[6][14] = _cov[14][6]*y[0][0] + _cov[14][7]*y[0][1] + _cov[14][8]*y[0][2] + var[96];
        _cov[7][14] = _cov[14][6]*y[1][0] + _cov[14][7]*y[1][1] + _cov[14][8]*y[1][2] + var[98];
        _cov[8][14] = _cov[14][6]*y[2][0] + _cov[14][7]*y[2][1] + _cov[14][8]*y[2][2] + var[100];

        _cov[14][14] = kahan_summation(_cov[14][14], _q_cov[14] * _dt2, _accumulator_cov[14]);
    }

    if (_control_status.flags.grav) {
        _cov[0][15] = var[8];
        _cov[1][15] = var[27];
        _cov[2][15] = var[42];
        _cov[3][15] = var[57];
        _cov[4][15] = var[75];
        _cov[5][15] = var[91];
        _cov[6][15] = _cov[15][6]*y[0][0] + _cov[15][7]*y[0][1] + _cov[15][8]*y[0][2] - _cov[15][9]*_dt;
        _cov[7][15] = -_cov[15][10]*_dt + _cov[15][6]*y[1][0] + _cov[15][7]*y[1][1] + _cov[15][8]*y[1][2];
        _cov[8][15] = -_cov[15][11]*_dt + _cov[15][6]*y[2][0] + _cov[15][7]*y[2][1] + _cov[15][8]*y[2][2];

        _cov[15][15] = kahan_summation(_cov[15][15], _q_cov[15] * _dt2, _accumulator_cov[15]);
    }

    if (_control_status.flags.mag) {
        _cov[0][16] = _cov[16][0]*y[0][0] + _cov[16][1]*y[0][1] + _cov[16][2]*y[0][2] + _cov[16][3]*_dt;
        _cov[1][16] = _cov[16][0]*y[1][0] + _cov[16][1]*y[1][1] + _cov[16][2]*y[1][2] + _cov[16][4]*_dt;
        _cov[2][16] = _cov[16][0]*y[2][0] + _cov[16][1]*y[2][1] + _cov[16][2]*y[2][2] + _cov[16][5]*_dt;
        _cov[3][16] = -_cov[16][12]*_dt + _cov[16][15]*r[0] + _cov[16][3]*y[0][0] + _cov[16][4]*y[0][1] + _cov[16][5]*y[0][2] + _cov[16][6]*x[0][0] + _cov[16][7]*x[0][1] + _cov[16][8]*x[0][2];
        _cov[4][16] = -_cov[16][13]*_dt + _cov[16][15]*r[1] + _cov[16][3]*y[1][0] + _cov[16][4]*y[1][1] + _cov[16][5]*y[1][2] + _cov[16][6]*x[1][0] + _cov[16][7]*x[1][1] + _cov[16][8]*x[1][2];
        _cov[5][16] = -_cov[16][14]*_dt + _cov[16][15]*r[2] + _cov[16][3]*y[2][0] + _cov[16][4]*y[2][1] + _cov[16][5]*y[2][2] + _cov[16][6]*x[2][0] + _cov[16][7]*x[2][1] + _cov[16][8]*x[2][2];
        _cov[6][16] = _cov[16][6]*y[0][0] + _cov[16][7]*y[0][1] + _cov[16][8]*y[0][2] - _cov[16][9]*_dt;
        _cov[7][16] = -_cov[16][10]*_dt + _cov[16][6]*y[1][0] + _cov[16][7]*y[1][1] + _cov[16][8]*y[1][2];
        _cov[8][16] = -_cov[16][11]*_dt + _cov[16][6]*y[2][0] + _cov[16][7]*y[2][1] + _cov[16][8]*y[2][2];

        _cov[0][17] = _cov[17][0]*y[0][0] + _cov[17][1]*y[0][1] + _cov[17][2]*y[0][2] + _cov[17][3]*_dt;
        _cov[1][17] = _cov[17][0]*y[1][0] + _cov[17][1]*y[1][1] + _cov[17][2]*y[1][2] + _cov[17][4]*_dt;
        _cov[2][17] = _cov[17][0]*y[2][0] + _cov[17][1]*y[2][1] + _cov[17][2]*y[2][2] + _cov[17][5]*_dt;
        _cov[3][17] = -_cov[17][12]*_dt + _cov[17][15]*r[0] + _cov[17][3]*y[0][0] + _cov[17][4]*y[0][1] + _cov[17][5]*y[0][2] + _cov[17][6]*x[0][0] + _cov[17][7]*x[0][1] + _cov[17][8]*x[0][2];
        _cov[4][17] = -_cov[17][13]*_dt + _cov[17][15]*r[1] + _cov[17][3]*y[1][0] + _cov[17][4]*y[1][1] + _cov[17][5]*y[1][2] + _cov[17][6]*x[1][0] + _cov[17][7]*x[1][1] + _cov[17][8]*x[1][2];
        _cov[5][17] = -_cov[17][14]*_dt + _cov[17][15]*r[2] + _cov[17][3]*y[2][0] + _cov[17][4]*y[2][1] + _cov[17][5]*y[2][2] + _cov[17][6]*x[2][0] + _cov[17][7]*x[2][1] + _cov[17][8]*x[2][2];
        _cov[6][17] = _cov[17][6]*y[0][0] + _cov[17][7]*y[0][1] + _cov[17][8]*y[0][2] - _cov[17][9]*_dt;
        _cov[7][17] = -_cov[17][10]*_dt + _cov[17][6]*y[1][0] + _cov[17][7]*y[1][1] + _cov[17][8]*y[1][2];
        _cov[8][17] = -_cov[17][11]*_dt + _cov[17][6]*y[2][0] + _cov[17][7]*y[2][1] + _cov[17][8]*y[2][2];

        _cov[0][18] = _cov[18][0]*y[0][0] + _cov[18][1]*y[0][1] + _cov[18][2]*y[0][2] + _cov[18][3]*_dt;
        _cov[1][18] = _cov[18][0]*y[1][0] + _cov[18][1]*y[1][1] + _cov[18][2]*y[1][2] + _cov[18][4]*_dt;
        _cov[2][18] = _cov[18][0]*y[2][0] + _cov[18][1]*y[2][1] + _cov[18][2]*y[2][2] + _cov[18][5]*_dt;
        _cov[3][18] = -_cov[18][12]*_dt + _cov[18][15]*r[0] + _cov[18][3]*y[0][0] + _cov[18][4]*y[0][1] + _cov[18][5]*y[0][2] + _cov[18][6]*x[0][0] + _cov[18][7]*x[0][1] + _cov[18][8]*x[0][2];
        _cov[4][18] = -_cov[18][13]*_dt + _cov[18][15]*r[1] + _cov[18][3]*y[1][0] + _cov[18][4]*y[1][1] + _cov[18][5]*y[1][2] + _cov[18][6]*x[1][0] + _cov[18][7]*x[1][1] + _cov[18][8]*x[1][2];
        _cov[5][18] = -_cov[18][14]*_dt + _cov[18][15]*r[2] + _cov[18][3]*y[2][0] + _cov[18][4]*y[2][1] + _cov[18][5]*y[2][2] + _cov[18][6]*x[2][0] + _cov[18][7]*x[2][1] + _cov[18][8]*x[2][2];
        _cov[6][18] = _cov[18][6]*y[0][0] + _cov[18][7]*y[0][1] + _cov[18][8]*y[0][2] - _cov[18][9]*_dt;
        _cov[7][18] = -_cov[18][10]*_dt + _cov[18][6]*y[1][0] + _cov[18][7]*y[1][1] + _cov[18][8]*y[1][2];
        _cov[8][18] = -_cov[18][11]*_dt + _cov[18][6]*y[2][0] + _cov[18][7]*y[2][1] + _cov[18][8]*y[2][2];

        add_processing_covariance<3>(16);

        if (_control_status.flags.mag_bias) {
            _cov[0][19] = _cov[19][0]*y[0][0] + _cov[19][1]*y[0][1] + _cov[19][2]*y[0][2] + _cov[19][3]*_dt;
            _cov[1][19] = _cov[19][0]*y[1][0] + _cov[19][1]*y[1][1] + _cov[19][2]*y[1][2] + _cov[19][4]*_dt;
            _cov[2][19] = _cov[19][0]*y[2][0] + _cov[19][1]*y[2][1] + _cov[19][2]*y[2][2] + _cov[19][5]*_dt;
            _cov[3][19] = -_cov[19][12]*_dt + _cov[19][15]*r[0] + _cov[19][3]*y[0][0] + _cov[19][4]*y[0][1] + _cov[19][5]*y[0][2] + _cov[19][6]*x[0][0] + _cov[19][7]*x[0][1] + _cov[19][8]*x[0][2];
            _cov[4][19] = -_cov[19][13]*_dt + _cov[19][15]*r[1] + _cov[19][3]*y[1][0] + _cov[19][4]*y[1][1] + _cov[19][5]*y[1][2] + _cov[19][6]*x[1][0] + _cov[19][7]*x[1][1] + _cov[19][8]*x[1][2];
            _cov[5][19] = -_cov[19][14]*_dt + _cov[19][15]*r[2] + _cov[19][3]*y[2][0] + _cov[19][4]*y[2][1] + _cov[19][5]*y[2][2] + _cov[19][6]*x[2][0] + _cov[19][7]*x[2][1] + _cov[19][8]*x[2][2];
            _cov[6][19] = _cov[19][6]*y[0][0] + _cov[19][7]*y[0][1] + _cov[19][8]*y[0][2] - _cov[19][9]*_dt;
            _cov[7][19] = -_cov[19][10]*_dt + _cov[19][6]*y[1][0] + _cov[19][7]*y[1][1] + _cov[19][8]*y[1][2];
            _cov[8][19] = -_cov[19][11]*_dt + _cov[19][6]*y[2][0] + _cov[19][7]*y[2][1] + _cov[19][8]*y[2][2];

            _cov[0][20] = _cov[20][0]*y[0][0] + _cov[20][1]*y[0][1] + _cov[20][2]*y[0][2] + _cov[20][3]*_dt;
            _cov[1][20] = _cov[20][0]*y[1][0] + _cov[20][1]*y[1][1] + _cov[20][2]*y[1][2] + _cov[20][4]*_dt;
            _cov[2][20] = _cov[20][0]*y[2][0] + _cov[20][1]*y[2][1] + _cov[20][2]*y[2][2] + _cov[20][5]*_dt;
            _cov[3][20] = -_cov[20][12]*_dt + _cov[20][15]*r[0] + _cov[20][3]*y[0][0] + _cov[20][4]*y[0][1] + _cov[20][5]*y[0][2] + _cov[20][6]*x[0][0] + _cov[20][7]*x[0][1] + _cov[20][8]*x[0][2];
            _cov[4][20] = -_cov[20][13]*_dt + _cov[20][15]*r[1] + _cov[20][3]*y[1][0] + _cov[20][4]*y[1][1] + _cov[20][5]*y[1][2] + _cov[20][6]*x[1][0] + _cov[20][7]*x[1][1] + _cov[20][8]*x[1][2];
            _cov[5][20] = -_cov[20][14]*_dt + _cov[20][15]*r[2] + _cov[20][3]*y[2][0] + _cov[20][4]*y[2][1] + _cov[20][5]*y[2][2] + _cov[20][6]*x[2][0] + _cov[20][7]*x[2][1] + _cov[20][8]*x[2][2];
            _cov[6][20] = _cov[20][6]*y[0][0] + _cov[20][7]*y[0][1] + _cov[20][8]*y[0][2] - _cov[20][9]*_dt;
            _cov[7][20] = -_cov[20][10]*_dt + _cov[20][6]*y[1][0] + _cov[20][7]*y[1][1] + _cov[20][8]*y[1][2];
            _cov[8][20] = -_cov[20][11]*_dt + _cov[20][6]*y[2][0] + _cov[20][7]*y[2][1] + _cov[20][8]*y[2][2];

            _cov[0][21] = _cov[21][0]*y[0][0] + _cov[21][1]*y[0][1] + _cov[21][2]*y[0][2] + _cov[21][3]*_dt;
            _cov[1][21] = _cov[21][0]*y[1][0] + _cov[21][1]*y[1][1] + _cov[21][2]*y[1][2] + _cov[21][4]*_dt;
            _cov[2][21] = _cov[21][0]*y[2][0] + _cov[21][1]*y[2][1] + _cov[21][2]*y[2][2] + _cov[21][5]*_dt;
            _cov[3][21] = -_cov[21][12]*_dt + _cov[21][15]*r[0] + _cov[21][3]*y[0][0] + _cov[21][4]*y[0][1] + _cov[21][5]*y[0][2] + _cov[21][6]*x[0][0] + _cov[21][7]*x[0][1] + _cov[21][8]*x[0][2];
            _cov[4][21] = -_cov[21][13]*_dt + _cov[21][15]*r[1] + _cov[21][3]*y[1][0] + _cov[21][4]*y[1][1] + _cov[21][5]*y[1][2] + _cov[21][6]*x[1][0] + _cov[21][7]*x[1][1] + _cov[21][8]*x[1][2];
            _cov[5][21] = -_cov[21][14]*_dt + _cov[21][15]*r[2] + _cov[21][3]*y[2][0] + _cov[21][4]*y[2][1] + _cov[21][5]*y[2][2] + _cov[21][6]*x[2][0] + _cov[21][7]*x[2][1] + _cov[21][8]*x[2][2];
            _cov[6][21] = _cov[21][6]*y[0][0] + _cov[21][7]*y[0][1] + _cov[21][8]*y[0][2] - _cov[21][9]*_dt;
            _cov[7][21] = -_cov[21][10]*_dt + _cov[21][6]*y[1][0] + _cov[21][7]*y[1][1] + _cov[21][8]*y[1][2];
            _cov[8][21] = -_cov[21][11]*_dt + _cov[21][6]*y[2][0] + _cov[21][7]*y[2][1] + _cov[21][8]*y[2][2];

            add_processing_covariance<3>(19);
        }
    }

    if (_control_status.flags.wind) {
        _cov[0][22] = _cov[22][0]*y[0][0] + _cov[22][1]*y[0][1] + _cov[22][2]*y[0][2] + _cov[22][3]*_dt;
        _cov[1][22] = _cov[22][0]*y[1][0] + _cov[22][1]*y[1][1] + _cov[22][2]*y[1][2] + _cov[22][4]*_dt;
        _cov[2][22] = _cov[22][0]*y[2][0] + _cov[22][1]*y[2][1] + _cov[22][2]*y[2][2] + _cov[22][5]*_dt;
        _cov[3][22] = -_cov[22][12]*_dt + _cov[22][15]*r[0] + _cov[22][3]*y[0][0] + _cov[22][4]*y[0][1] + _cov[22][5]*y[0][2] + _cov[22][6]*x[0][0] + _cov[22][7]*x[0][1] + _cov[22][8]*x[0][2];
        _cov[4][22] = -_cov[22][13]*_dt + _cov[22][15]*r[1] + _cov[22][3]*y[1][0] + _cov[22][4]*y[1][1] + _cov[22][5]*y[1][2] + _cov[22][6]*x[1][0] + _cov[22][7]*x[1][1] + _cov[22][8]*x[1][2];
        _cov[5][22] = -_cov[22][14]*_dt + _cov[22][15]*r[2] + _cov[22][3]*y[2][0] + _cov[22][4]*y[2][1] + _cov[22][5]*y[2][2] + _cov[22][6]*x[2][0] + _cov[22][7]*x[2][1] + _cov[22][8]*x[2][2];
        _cov[6][22] = _cov[22][6]*y[0][0] + _cov[22][7]*y[0][1] + _cov[22][8]*y[0][2] - _cov[22][9]*_dt;
        _cov[7][22] = -_cov[22][10]*_dt + _cov[22][6]*y[1][0] + _cov[22][7]*y[1][1] + _cov[22][8]*y[1][2];
        _cov[8][22] = -_cov[22][11]*_dt + _cov[22][6]*y[2][0] + _cov[22][7]*y[2][1] + _cov[22][8]*y[2][2];

        _cov[0][23] = _cov[23][0]*y[0][0] + _cov[23][1]*y[0][1] + _cov[23][2]*y[0][2] + _cov[23][3]*_dt;
        _cov[1][23] = _cov[23][0]*y[1][0] + _cov[23][1]*y[1][1] + _cov[23][2]*y[1][2] + _cov[23][4]*_dt;
        _cov[2][23] = _cov[23][0]*y[2][0] + _cov[23][1]*y[2][1] + _cov[23][2]*y[2][2] + _cov[23][5]*_dt;
        _cov[3][23] = -_cov[23][12]*_dt + _cov[23][15]*r[0] + _cov[23][3]*y[0][0] + _cov[23][4]*y[0][1] + _cov[23][5]*y[0][2] + _cov[23][6]*x[0][0] + _cov[23][7]*x[0][1] + _cov[23][8]*x[0][2];
        _cov[4][23] = -_cov[23][13]*_dt + _cov[23][15]*r[1] + _cov[23][3]*y[1][0] + _cov[23][4]*y[1][1] + _cov[23][5]*y[1][2] + _cov[23][6]*x[1][0] + _cov[23][7]*x[1][1] + _cov[23][8]*x[1][2];
        _cov[5][23] = -_cov[23][14]*_dt + _cov[23][15]*r[2] + _cov[23][3]*y[2][0] + _cov[23][4]*y[2][1] + _cov[23][5]*y[2][2] + _cov[23][6]*x[2][0] + _cov[23][7]*x[2][1] + _cov[23][8]*x[2][2];
        _cov[6][23] = _cov[23][6]*y[0][0] + _cov[23][7]*y[0][1] + _cov[23][8]*y[0][2] - _cov[23][9]*_dt;
        _cov[7][23] = -_cov[23][10]*_dt + _cov[23][6]*y[1][0] + _cov[23][7]*y[1][1] + _cov[23][8]*y[1][2];
        _cov[8][23] = -_cov[23][11]*_dt + _cov[23][6]*y[2][0] + _cov[23][7]*y[2][1] + _cov[23][8]*y[2][2];

        add_processing_covariance<2>(22);
    }

    regular_covariance_to_symmetric<ESKF::dim>(0);
}

unsigned char LIEKF::fuse_position(const Vector3f &pos, const Vector3f &w, const Vector3f &a, 
                                  const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) {
    /*
    s = [p, v, R, bg, ba, g]
    δs = [δp, δv, δθ, δbg, δba, δg]
    s + δs = [p+R*Exp(δθ)*δp, v+R*Exp(δθ)*δv, R*Exp(δθ), bg+δbg, ba+δba, g+δg]

    pos = p + R * dis

    δpos / δp = I
    δpos / δθ = -dis^

    H = [I, O, -dis^, O, O, O]
    */

    unsigned char info = 0;

    // -dis^
    const array<array<float, 3>, 3> minus_d_hat = {
        {{0.f, dis[2], -dis[1]},
         {-dis[2], 0.f, dis[0]},
         {dis[1], -dis[0], 0.f}}
    };

    // pos - p
    const array<float, 3> e_p {pos[0] - _p[0], pos[1] - _p[1], pos[2] - _p[2]};

    // Sequential Kalman Filter
    for (unsigned int dim = 0; dim < 3; ++dim) {
        if (!isfinite(noise_std[dim])) {
            continue;
        }

        /*
        K = P * H * (H * P * H' + R)^-1
        x = x + K * (y - h)
        P = (I - K * H) * P

        where, H = [I, O, -dis^, O, O, O]
               h = p + R * dis
        */
        
        // H * P  or  P * H'
        const float cov_1_dim = (dim == 2) ? _cov[1][dim] : _cov[dim][1];
        array<float, ESKF::dim> HP = {_cov[0][dim] + _cov[0][6]*minus_d_hat[dim][0] + _cov[0][7]*minus_d_hat[dim][1] + _cov[0][8]*minus_d_hat[dim][2],
                                        cov_1_dim + _cov[1][6]*minus_d_hat[dim][0] + _cov[1][7]*minus_d_hat[dim][1] + _cov[1][8]*minus_d_hat[dim][2],
                                        _cov[dim][2] + _cov[2][6]*minus_d_hat[dim][0] + _cov[2][7]*minus_d_hat[dim][1] + _cov[2][8]*minus_d_hat[dim][2],
                                        _cov[dim][3] + _cov[3][6]*minus_d_hat[dim][0] + _cov[3][7]*minus_d_hat[dim][1] + _cov[3][8]*minus_d_hat[dim][2],
                                        _cov[dim][4] + _cov[4][6]*minus_d_hat[dim][0] + _cov[4][7]*minus_d_hat[dim][1] + _cov[4][8]*minus_d_hat[dim][2],
                                        _cov[dim][5] + _cov[5][6]*minus_d_hat[dim][0] + _cov[5][7]*minus_d_hat[dim][1] + _cov[5][8]*minus_d_hat[dim][2],
                                        _cov[dim][6] + _cov[6][6]*minus_d_hat[dim][0] + _cov[6][7]*minus_d_hat[dim][1] + _cov[6][8]*minus_d_hat[dim][2],
                                        _cov[dim][7] + _cov[6][7]*minus_d_hat[dim][0] + _cov[7][7]*minus_d_hat[dim][1] + _cov[7][8]*minus_d_hat[dim][2],
                                        _cov[dim][8] + _cov[6][8]*minus_d_hat[dim][0] + _cov[7][8]*minus_d_hat[dim][1] + _cov[8][8]*minus_d_hat[dim][2],
                                        _cov[dim][9] + _cov[6][9]*minus_d_hat[dim][0] + _cov[7][9]*minus_d_hat[dim][1] + _cov[8][9]*minus_d_hat[dim][2],
                                        _cov[dim][10] + _cov[6][10]*minus_d_hat[dim][0] + _cov[7][10]*minus_d_hat[dim][1] + _cov[8][10]*minus_d_hat[dim][2],
                                        _cov[dim][11] + _cov[6][11]*minus_d_hat[dim][0] + _cov[7][11]*minus_d_hat[dim][1] + _cov[8][11]*minus_d_hat[dim][2],
                                        _cov[dim][12] + _cov[6][12]*minus_d_hat[dim][0] + _cov[7][12]*minus_d_hat[dim][1] + _cov[8][12]*minus_d_hat[dim][2],
                                        _cov[dim][13] + _cov[6][13]*minus_d_hat[dim][0] + _cov[7][13]*minus_d_hat[dim][1] + _cov[8][13]*minus_d_hat[dim][2],
                                        _cov[dim][14] + _cov[6][14]*minus_d_hat[dim][0] + _cov[7][14]*minus_d_hat[dim][1] + _cov[8][14]*minus_d_hat[dim][2],
                                        _cov[dim][15] + _cov[6][15]*minus_d_hat[dim][0] + _cov[7][15]*minus_d_hat[dim][1] + _cov[8][15]*minus_d_hat[dim][2],
                                        0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};

        if (_control_status.flags.mag) {
            HP[16] = _cov[dim][16] + _cov[6][16]*minus_d_hat[dim][0] + _cov[7][16]*minus_d_hat[dim][1] + _cov[8][16]*minus_d_hat[dim][2];
            HP[17] = _cov[dim][17] + _cov[6][17]*minus_d_hat[dim][0] + _cov[7][17]*minus_d_hat[dim][1] + _cov[8][17]*minus_d_hat[dim][2];
            HP[18] = _cov[dim][18] + _cov[6][18]*minus_d_hat[dim][0] + _cov[7][18]*minus_d_hat[dim][1] + _cov[8][18]*minus_d_hat[dim][2];

            if (_control_status.flags.mag_bias) {
                HP[19] = _cov[dim][19] + _cov[6][19]*minus_d_hat[dim][0] + _cov[7][19]*minus_d_hat[dim][1] + _cov[8][19]*minus_d_hat[dim][2];
                HP[20] = _cov[dim][20] + _cov[6][20]*minus_d_hat[dim][0] + _cov[7][20]*minus_d_hat[dim][1] + _cov[8][20]*minus_d_hat[dim][2];
                HP[21] = _cov[dim][21] + _cov[6][21]*minus_d_hat[dim][0] + _cov[7][21]*minus_d_hat[dim][1] + _cov[8][21]*minus_d_hat[dim][2];                    
            }                                    
        }   

        if (_control_status.flags.wind) {
            HP[22] = _cov[dim][22] + _cov[6][22]*minus_d_hat[dim][0] + _cov[7][22]*minus_d_hat[dim][1] + _cov[8][22]*minus_d_hat[dim][2];
            HP[23] = _cov[dim][23] + _cov[6][23]*minus_d_hat[dim][0] + _cov[7][23]*minus_d_hat[dim][1] + _cov[8][23]*minus_d_hat[dim][2];  
        }       

        // H * P * H' + R
        const float HPHT_plus_R = HP[dim] + HP[6] * minus_d_hat[dim][0] + HP[7] * minus_d_hat[dim][1] + HP[8] * minus_d_hat[dim][2] + noise_std[dim] * noise_std[dim];

        // e = R' * (y - p) - dis
        const float obs_error = _rot(0, dim) * e_p[0] + _rot(1, dim) * e_p[1] + _rot(2, dim) * e_p[2] - dis[dim];

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

    regular_covariance_to_symmetric<ESKF::dim>(0);

    return info;
}

unsigned char LIEKF::fuse_velocity(const Vector3f &vel, const Vector3f &w, const Vector3f &a, 
                                  const Vector3f &dis, const Vector3f &noise_std, const Vector3f &gate) {
    /*
    s = [p, v, R, bg, ba, g]
    δs = [δp, δv, δθ, δbg, δba, δg]
    s + δs = [p+R*Exp(δθ)*δp, v+R*Exp(δθ)*δv, R*Exp(δθ), bg+δbg, ba+δba, g+δg]

    vel = v + R * (w - bg)^ * dis

    δvel / δv = I
    δvel / δθ = -((w - bg)^ * dis)^
    δvel / δbg = dis^

    H = [O, I, -((w-bg)^*dis)^, dis^, O, O]
    */ 

    unsigned char info = 0;

    // dis^
    const array<array<float, 3>, 3> d_hat = {
        {{0.f, -dis[2], dis[1]},
         {dis[2], 0.f, -dis[0]},
         {-dis[1], dis[0], 0.f}}
    };

    // -(w-bg)^ * dis = dis^ * (w-bg)
    const array<float, 3> w_corr = {w[0] - _bg[0], w[1] - _bg[1], w[2] - _bg[2]};
    const array<float, 3> d_cross_w_corr = {w_corr[2] * dis[1] - w_corr[1] * dis[2], 
                                            w_corr[0] * dis[2] - w_corr[2] * dis[0], 
                                            w_corr[1] * dis[0] - w_corr[0] * dis[1]};

    // (dis^ * (w-bg))^
    const array<array<float, 3>, 3> d_cross_w_corr_hat {
        {{0.f, -d_cross_w_corr[2], d_cross_w_corr[1]},
         {d_cross_w_corr[2], 0.f, -d_cross_w_corr[0]},
         {-d_cross_w_corr[1], d_cross_w_corr[0], 0.f}}
    };

    // vel - v
    const array<float, 3> e_v {vel[0] - _v[0], vel[1] - _v[1], vel[2] - _v[2]};
    
    // Sequential Kalman Filter
    for (unsigned int dim = 0; dim < 3; ++dim) {
        if (!isfinite(noise_std[dim])) {
            continue;
        }

        /*
        K = P * H * (H * P * H' + R)^-1
        x = x + K * (y - h)
        P = (I - K * H) * P

        where, H = [O, I, (dis^*(w-bg))^, dis^, O, O]
               h = v + R * (w-bg)^ * dis
        */

        // H * P  or  P * H'
        const unsigned int index = 3 + dim;
        const float cov_4_index = (dim == 2) ? _cov[4][index] : _cov[index][4];
        array<float, ESKF::dim> HP = {_cov[0][index] + _cov[0][6]*d_cross_w_corr_hat[dim][0] + _cov[0][7]*d_cross_w_corr_hat[dim][1] + _cov[0][8]*d_cross_w_corr_hat[dim][2] + _cov[0][9]*d_hat[dim][0] + _cov[0][10]*d_hat[dim][1] + _cov[0][11]*d_hat[dim][2],
                                        _cov[1][index] + _cov[1][6]*d_cross_w_corr_hat[dim][0] + _cov[1][7]*d_cross_w_corr_hat[dim][1] + _cov[1][8]*d_cross_w_corr_hat[dim][2] + _cov[1][9]*d_hat[dim][0] + _cov[1][10]*d_hat[dim][1] + _cov[1][11]*d_hat[dim][2],
                                        _cov[2][index] + _cov[2][6]*d_cross_w_corr_hat[dim][0] + _cov[2][7]*d_cross_w_corr_hat[dim][1] + _cov[2][8]*d_cross_w_corr_hat[dim][2] + _cov[2][9]*d_hat[dim][0] + _cov[2][10]*d_hat[dim][1] + _cov[2][11]*d_hat[dim][2],
                                        _cov[3][index] + _cov[3][6]*d_cross_w_corr_hat[dim][0] + _cov[3][7]*d_cross_w_corr_hat[dim][1] + _cov[3][8]*d_cross_w_corr_hat[dim][2] + _cov[3][9]*d_hat[dim][0] + _cov[3][10]*d_hat[dim][1] + _cov[3][11]*d_hat[dim][2],
                                        cov_4_index + _cov[4][6]*d_cross_w_corr_hat[dim][0] + _cov[4][7]*d_cross_w_corr_hat[dim][1] + _cov[4][8]*d_cross_w_corr_hat[dim][2] + _cov[4][9]*d_hat[dim][0] + _cov[4][10]*d_hat[dim][1] + _cov[4][11]*d_hat[dim][2],
                                        _cov[index][5] + _cov[5][6]*d_cross_w_corr_hat[dim][0] + _cov[5][7]*d_cross_w_corr_hat[dim][1] + _cov[5][8]*d_cross_w_corr_hat[dim][2] + _cov[5][9]*d_hat[dim][0] + _cov[5][10]*d_hat[dim][1] + _cov[5][11]*d_hat[dim][2],
                                        _cov[index][6] + _cov[6][6]*d_cross_w_corr_hat[dim][0] + _cov[6][7]*d_cross_w_corr_hat[dim][1] + _cov[6][8]*d_cross_w_corr_hat[dim][2] + _cov[6][9]*d_hat[dim][0] + _cov[6][10]*d_hat[dim][1] + _cov[6][11]*d_hat[dim][2],
                                        _cov[index][7] + _cov[6][7]*d_cross_w_corr_hat[dim][0] + _cov[7][7]*d_cross_w_corr_hat[dim][1] + _cov[7][8]*d_cross_w_corr_hat[dim][2] + _cov[7][9]*d_hat[dim][0] + _cov[7][10]*d_hat[dim][1] + _cov[7][11]*d_hat[dim][2],
                                        _cov[index][8] + _cov[6][8]*d_cross_w_corr_hat[dim][0] + _cov[7][8]*d_cross_w_corr_hat[dim][1] + _cov[8][8]*d_cross_w_corr_hat[dim][2] + _cov[8][9]*d_hat[dim][0] + _cov[8][10]*d_hat[dim][1] + _cov[8][11]*d_hat[dim][2],
                                        _cov[index][9] + _cov[6][9]*d_cross_w_corr_hat[dim][0] + _cov[7][9]*d_cross_w_corr_hat[dim][1] + _cov[8][9]*d_cross_w_corr_hat[dim][2] + _cov[9][9]*d_hat[dim][0] + _cov[9][10]*d_hat[dim][1] + _cov[9][11]*d_hat[dim][2],
                                        _cov[index][10] + _cov[6][10]*d_cross_w_corr_hat[dim][0] + _cov[7][10]*d_cross_w_corr_hat[dim][1] + _cov[8][10]*d_cross_w_corr_hat[dim][2] + _cov[9][10]*d_hat[dim][0] + _cov[10][10]*d_hat[dim][1] + _cov[10][11]*d_hat[dim][2],
                                        _cov[index][11] + _cov[6][11]*d_cross_w_corr_hat[dim][0] + _cov[7][11]*d_cross_w_corr_hat[dim][1] + _cov[8][11]*d_cross_w_corr_hat[dim][2] + _cov[9][11]*d_hat[dim][0] + _cov[10][11]*d_hat[dim][1] + _cov[11][11]*d_hat[dim][2],
                                        _cov[index][12] + _cov[6][12]*d_cross_w_corr_hat[dim][0] + _cov[7][12]*d_cross_w_corr_hat[dim][1] + _cov[8][12]*d_cross_w_corr_hat[dim][2] + _cov[9][12]*d_hat[dim][0] + _cov[10][12]*d_hat[dim][1] + _cov[11][12]*d_hat[dim][2],
                                        _cov[index][13] + _cov[6][13]*d_cross_w_corr_hat[dim][0] + _cov[7][13]*d_cross_w_corr_hat[dim][1] + _cov[8][13]*d_cross_w_corr_hat[dim][2] + _cov[9][13]*d_hat[dim][0] + _cov[10][13]*d_hat[dim][1] + _cov[11][13]*d_hat[dim][2],
                                        _cov[index][14] + _cov[6][14]*d_cross_w_corr_hat[dim][0] + _cov[7][14]*d_cross_w_corr_hat[dim][1] + _cov[8][14]*d_cross_w_corr_hat[dim][2] + _cov[9][14]*d_hat[dim][0] + _cov[10][14]*d_hat[dim][1] + _cov[11][14]*d_hat[dim][2],
                                        _cov[index][15] + _cov[6][15]*d_cross_w_corr_hat[dim][0] + _cov[7][15]*d_cross_w_corr_hat[dim][1] + _cov[8][15]*d_cross_w_corr_hat[dim][2] + _cov[9][15]*d_hat[dim][0] + _cov[10][15]*d_hat[dim][1] + _cov[11][15]*d_hat[dim][2],
                                        0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};         

        if (_control_status.flags.mag) {
            HP[16] = _cov[index][16] + _cov[6][16]*d_cross_w_corr_hat[dim][0] + _cov[7][16]*d_cross_w_corr_hat[dim][1] + _cov[8][16]*d_cross_w_corr_hat[dim][2] + _cov[9][16]*d_hat[dim][0] + _cov[10][16]*d_hat[dim][1] + _cov[11][16]*d_hat[dim][2];   
            HP[17] = _cov[index][17] + _cov[6][17]*d_cross_w_corr_hat[dim][0] + _cov[7][17]*d_cross_w_corr_hat[dim][1] + _cov[8][17]*d_cross_w_corr_hat[dim][2] + _cov[9][17]*d_hat[dim][0] + _cov[10][17]*d_hat[dim][1] + _cov[11][17]*d_hat[dim][2];  
            HP[18] = _cov[index][18] + _cov[6][18]*d_cross_w_corr_hat[dim][0] + _cov[7][18]*d_cross_w_corr_hat[dim][1] + _cov[8][18]*d_cross_w_corr_hat[dim][2] + _cov[9][18]*d_hat[dim][0] + _cov[10][18]*d_hat[dim][1] + _cov[11][18]*d_hat[dim][2];                       
            if (_control_status.flags.mag_bias) {
                HP[19] = _cov[index][19] + _cov[6][19]*d_cross_w_corr_hat[dim][0] + _cov[7][19]*d_cross_w_corr_hat[dim][1] + _cov[8][19]*d_cross_w_corr_hat[dim][2] + _cov[9][19]*d_hat[dim][0] + _cov[10][19]*d_hat[dim][1] + _cov[11][19]*d_hat[dim][2];   
                HP[20] = _cov[index][20] + _cov[6][20]*d_cross_w_corr_hat[dim][0] + _cov[7][20]*d_cross_w_corr_hat[dim][1] + _cov[8][20]*d_cross_w_corr_hat[dim][2] + _cov[9][20]*d_hat[dim][0] + _cov[10][20]*d_hat[dim][1] + _cov[11][20]*d_hat[dim][2];  
                HP[21] = _cov[index][21] + _cov[6][21]*d_cross_w_corr_hat[dim][0] + _cov[7][21]*d_cross_w_corr_hat[dim][1] + _cov[8][21]*d_cross_w_corr_hat[dim][2] + _cov[9][21]*d_hat[dim][0] + _cov[10][21]*d_hat[dim][1] + _cov[11][21]*d_hat[dim][2];
                                        
            }
        }

        if (_control_status.flags.wind) {
            HP[22] = _cov[index][22] + _cov[6][22]*d_cross_w_corr_hat[dim][0] + _cov[7][22]*d_cross_w_corr_hat[dim][1] + _cov[8][22]*d_cross_w_corr_hat[dim][2] + _cov[9][22]*d_hat[dim][0] + _cov[10][22]*d_hat[dim][1] + _cov[11][22]*d_hat[dim][2];  
            HP[23] = _cov[index][23] + _cov[6][23]*d_cross_w_corr_hat[dim][0] + _cov[7][23]*d_cross_w_corr_hat[dim][1] + _cov[8][23]*d_cross_w_corr_hat[dim][2] + _cov[9][23]*d_hat[dim][0] + _cov[10][23]*d_hat[dim][1] + _cov[11][23]*d_hat[dim][2];
        }

        // H * P * H' + R
        const float HPHT_plus_R = HP[index] + HP[6] * d_cross_w_corr_hat[dim][0] + HP[7] * d_cross_w_corr_hat[dim][1] + HP[8] * d_cross_w_corr_hat[dim][2] + HP[9] * d_hat[dim][0] + HP[10] * d_hat[dim][1] + HP[11] * d_hat[dim][2] + noise_std[dim] * noise_std[dim];

        // h = v + R * (w - bg)^ * dis
        // e = y - h = y - (v + R * (w - bg)^ * dis)
        const float obs_error = _rot(0, dim) * e_v[0] + _rot(1, dim) * e_v[1] + _rot(2, dim) * e_v[2] + d_cross_w_corr[dim];

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

    regular_covariance_to_symmetric<ESKF::dim>(0);

    return info;
}

unsigned char LIEKF::fuse_magnet(const Vector3f &mag, const Vector3f &w, const Vector3f &a, 
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

    // R' * RZ * RY * ex
    const array<float, 3> rt_rz_ry_ex {
        _rot(0, 0) * rz_ry_ex[0] + _rot(1, 0) * rz_ry_ex[1] + _rot(2, 0) * rz_ry_ex[2],
        _rot(0, 1) * rz_ry_ex[0] + _rot(1, 1) * rz_ry_ex[1] + _rot(2, 1) * rz_ry_ex[2],
        _rot(0, 2) * rz_ry_ex[0] + _rot(1, 2) * rz_ry_ex[1] + _rot(2, 2) * rz_ry_ex[2]
    };

    // (R' * RZ * RY * ex)^
    const array<array<float, 3>, 3> rt_rz_ry_ex_hat {
        {{0.f, -rt_rz_ry_ex[2], rt_rz_ry_ex[1]},
         {rt_rz_ry_ex[2], 0.f, -rt_rz_ry_ex[0]},
         {-rt_rz_ry_ex[1], rt_rz_ry_ex[0], 0.f}}
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

        // R' * RZ * RY * ex / H
        const float rt_rz_ry_ex_h_inv = rt_rz_ry_ex[dim] * h_inv;

        // R' * (RZ * RY * ex)^
        const array<float, 3> rt__rz_ry_ex_hat {
            _rot(1, dim) * rz_ry_ex[2] - _rot(2, dim) * rz_ry_ex[1], 
            _rot(2, dim) * rz_ry_ex[0] - _rot(0, dim) * rz_ry_ex[2], 
            _rot(0, dim) * rz_ry_ex[1] - _rot(1, dim) * rz_ry_ex[0]
        };

        const float param_y = rt__rz_ry_ex_hat[1]*cos_z - sin_z*rt__rz_ry_ex_hat[0];

        // H * P  or  P * H'
        const unsigned int index = 19 + dim;
        array<float, ESKF::dim> HP = {_cov[0][index] * h_inv + _cov[0][6]*rt_rz_ry_ex_hat[dim][0] + _cov[0][7]*rt_rz_ry_ex_hat[dim][1] + _cov[0][8]*rt_rz_ry_ex_hat[dim][2] + _cov[0][16]*rt_rz_ry_ex_h_inv - _cov[0][17]*param_y - _cov[0][18]*rt__rz_ry_ex_hat[2],
                                        _cov[1][index] * h_inv + _cov[1][6]*rt_rz_ry_ex_hat[dim][0] + _cov[1][7]*rt_rz_ry_ex_hat[dim][1] + _cov[1][8]*rt_rz_ry_ex_hat[dim][2] + _cov[1][16]*rt_rz_ry_ex_h_inv - _cov[1][17]*param_y - _cov[1][18]*rt__rz_ry_ex_hat[2],
                                        _cov[2][index] * h_inv + _cov[2][6]*rt_rz_ry_ex_hat[dim][0] + _cov[2][7]*rt_rz_ry_ex_hat[dim][1] + _cov[2][8]*rt_rz_ry_ex_hat[dim][2] + _cov[2][16]*rt_rz_ry_ex_h_inv - _cov[2][17]*param_y - _cov[2][18]*rt__rz_ry_ex_hat[2],
                                        _cov[3][index] * h_inv + _cov[3][6]*rt_rz_ry_ex_hat[dim][0] + _cov[3][7]*rt_rz_ry_ex_hat[dim][1] + _cov[3][8]*rt_rz_ry_ex_hat[dim][2] + _cov[3][16]*rt_rz_ry_ex_h_inv - _cov[3][17]*param_y - _cov[3][18]*rt__rz_ry_ex_hat[2],
                                        _cov[4][index] * h_inv + _cov[4][6]*rt_rz_ry_ex_hat[dim][0] + _cov[4][7]*rt_rz_ry_ex_hat[dim][1] + _cov[4][8]*rt_rz_ry_ex_hat[dim][2] + _cov[4][16]*rt_rz_ry_ex_h_inv - _cov[4][17]*param_y - _cov[4][18]*rt__rz_ry_ex_hat[2],
                                        _cov[5][index] * h_inv + _cov[5][6]*rt_rz_ry_ex_hat[dim][0] + _cov[5][7]*rt_rz_ry_ex_hat[dim][1] + _cov[5][8]*rt_rz_ry_ex_hat[dim][2] + _cov[5][16]*rt_rz_ry_ex_h_inv - _cov[5][17]*param_y - _cov[5][18]*rt__rz_ry_ex_hat[2],
                                        _cov[6][index] * h_inv + _cov[6][6]*rt_rz_ry_ex_hat[dim][0] + _cov[6][7]*rt_rz_ry_ex_hat[dim][1] + _cov[6][8]*rt_rz_ry_ex_hat[dim][2] + _cov[6][16]*rt_rz_ry_ex_h_inv - _cov[6][17]*param_y - _cov[6][18]*rt__rz_ry_ex_hat[2],
                                        _cov[7][index] * h_inv + _cov[6][7]*rt_rz_ry_ex_hat[dim][0] + _cov[7][7]*rt_rz_ry_ex_hat[dim][1] + _cov[7][8]*rt_rz_ry_ex_hat[dim][2] + _cov[7][16]*rt_rz_ry_ex_h_inv - _cov[7][17]*param_y - _cov[7][18]*rt__rz_ry_ex_hat[2],
                                        _cov[8][index] * h_inv + _cov[6][8]*rt_rz_ry_ex_hat[dim][0] + _cov[7][8]*rt_rz_ry_ex_hat[dim][1] + _cov[8][8]*rt_rz_ry_ex_hat[dim][2] + _cov[8][16]*rt_rz_ry_ex_h_inv - _cov[8][17]*param_y - _cov[8][18]*rt__rz_ry_ex_hat[2],
                                        _cov[9][index] * h_inv + _cov[6][9]*rt_rz_ry_ex_hat[dim][0] + _cov[7][9]*rt_rz_ry_ex_hat[dim][1] + _cov[8][9]*rt_rz_ry_ex_hat[dim][2] + _cov[9][16]*rt_rz_ry_ex_h_inv - _cov[9][17]*param_y - _cov[9][18]*rt__rz_ry_ex_hat[2],
                                        _cov[10][index] * h_inv + _cov[6][10]*rt_rz_ry_ex_hat[dim][0] + _cov[7][10]*rt_rz_ry_ex_hat[dim][1] + _cov[8][10]*rt_rz_ry_ex_hat[dim][2] + _cov[10][16]*rt_rz_ry_ex_h_inv - _cov[10][17]*param_y - _cov[10][18]*rt__rz_ry_ex_hat[2],
                                        _cov[11][index] * h_inv + _cov[6][11]*rt_rz_ry_ex_hat[dim][0] + _cov[7][11]*rt_rz_ry_ex_hat[dim][1] + _cov[8][11]*rt_rz_ry_ex_hat[dim][2] + _cov[11][16]*rt_rz_ry_ex_h_inv - _cov[11][17]*param_y - _cov[11][18]*rt__rz_ry_ex_hat[2],
                                        _cov[12][index] * h_inv + _cov[6][12]*rt_rz_ry_ex_hat[dim][0] + _cov[7][12]*rt_rz_ry_ex_hat[dim][1] + _cov[8][12]*rt_rz_ry_ex_hat[dim][2] + _cov[12][16]*rt_rz_ry_ex_h_inv - _cov[12][17]*param_y - _cov[12][18]*rt__rz_ry_ex_hat[2],
                                        _cov[13][index] * h_inv + _cov[6][13]*rt_rz_ry_ex_hat[dim][0] + _cov[7][13]*rt_rz_ry_ex_hat[dim][1] + _cov[8][13]*rt_rz_ry_ex_hat[dim][2] + _cov[13][16]*rt_rz_ry_ex_h_inv - _cov[13][17]*param_y - _cov[13][18]*rt__rz_ry_ex_hat[2],
                                        _cov[14][index] * h_inv + _cov[6][14]*rt_rz_ry_ex_hat[dim][0] + _cov[7][14]*rt_rz_ry_ex_hat[dim][1] + _cov[8][14]*rt_rz_ry_ex_hat[dim][2] + _cov[14][16]*rt_rz_ry_ex_h_inv - _cov[14][17]*param_y - _cov[14][18]*rt__rz_ry_ex_hat[2],
                                        _cov[15][index] * h_inv + _cov[6][15]*rt_rz_ry_ex_hat[dim][0] + _cov[7][15]*rt_rz_ry_ex_hat[dim][1] + _cov[8][15]*rt_rz_ry_ex_hat[dim][2] + _cov[15][16]*rt_rz_ry_ex_h_inv - _cov[15][17]*param_y - _cov[15][18]*rt__rz_ry_ex_hat[2],
                                        0.f, 0.f, 0.f,
                                        0.f, 0.f, 0.f, 0.f, 0.f};

        if (_control_status.flags.mag) {
            HP[16] = _cov[16][index] * h_inv + _cov[6][16]*rt_rz_ry_ex_hat[dim][0] + _cov[7][16]*rt_rz_ry_ex_hat[dim][1] + _cov[8][16]*rt_rz_ry_ex_hat[dim][2] + _cov[16][16]*rt_rz_ry_ex_h_inv - _cov[16][17]*param_y - _cov[16][18]*rt__rz_ry_ex_hat[2];
        }      

        if (_control_status.flags.dec) {
            HP[17] = _cov[17][index] * h_inv + _cov[6][17]*rt_rz_ry_ex_hat[dim][0] + _cov[7][17]*rt_rz_ry_ex_hat[dim][1] + _cov[8][17]*rt_rz_ry_ex_hat[dim][2] + _cov[16][17]*rt_rz_ry_ex_h_inv - _cov[17][17]*param_y - _cov[17][18]*rt__rz_ry_ex_hat[2];
            HP[18] = _cov[18][index] * h_inv + _cov[6][18]*rt_rz_ry_ex_hat[dim][0] + _cov[7][18]*rt_rz_ry_ex_hat[dim][1] + _cov[8][18]*rt_rz_ry_ex_hat[dim][2] + _cov[16][18]*rt_rz_ry_ex_h_inv - _cov[17][18]*param_y - _cov[18][18]*rt__rz_ry_ex_hat[2];
        }                          

        if (_control_status.flags.mag_bias) {
            const float cov_20_index = (dim == 2) ? _cov[20][index] : _cov[index][20];
            HP[19] = _cov[19][index] * h_inv + _cov[6][19]*rt_rz_ry_ex_hat[dim][0] + _cov[7][19]*rt_rz_ry_ex_hat[dim][1] + _cov[8][19]*rt_rz_ry_ex_hat[dim][2] + _cov[16][19]*rt_rz_ry_ex_h_inv - _cov[17][19]*param_y - _cov[18][19]*rt__rz_ry_ex_hat[2];
            HP[20] = cov_20_index * h_inv + _cov[6][20]*rt_rz_ry_ex_hat[dim][0] + _cov[7][20]*rt_rz_ry_ex_hat[dim][1] + _cov[8][20]*rt_rz_ry_ex_hat[dim][2] + _cov[16][20]*rt_rz_ry_ex_h_inv - _cov[17][20]*param_y - _cov[18][20]*rt__rz_ry_ex_hat[2];
            HP[21] = _cov[index][21] * h_inv + _cov[6][21]*rt_rz_ry_ex_hat[dim][0] + _cov[7][21]*rt_rz_ry_ex_hat[dim][1] + _cov[8][21]*rt_rz_ry_ex_hat[dim][2] + _cov[16][21]*rt_rz_ry_ex_h_inv - _cov[17][21]*param_y - _cov[18][21]*rt__rz_ry_ex_hat[2];                           
        }

        if (_control_status.flags.wind) {
            HP[22] = _cov[index][22] * h_inv + _cov[6][22]*rt_rz_ry_ex_hat[dim][0] + _cov[7][22]*rt_rz_ry_ex_hat[dim][1] + _cov[8][22]*rt_rz_ry_ex_hat[dim][2] + _cov[16][22]*rt_rz_ry_ex_h_inv - _cov[17][22]*param_y - _cov[18][22]*rt__rz_ry_ex_hat[2];
            HP[23] = _cov[index][23] * h_inv + _cov[6][23]*rt_rz_ry_ex_hat[dim][0] + _cov[7][23]*rt_rz_ry_ex_hat[dim][1] + _cov[8][23]*rt_rz_ry_ex_hat[dim][2] + _cov[16][23]*rt_rz_ry_ex_h_inv - _cov[17][23]*param_y - _cov[18][23]*rt__rz_ry_ex_hat[2];
        }

        // H * P * H' + R
        const float HPHT_plus_R = HP[index] * h_inv + HP[6] * rt_rz_ry_ex_hat[dim][0] + HP[7] * rt_rz_ry_ex_hat[dim][1] + HP[8] * rt_rz_ry_ex_hat[dim][2] + HP[16] * rt_rz_ry_ex_h_inv - HP[17] * param_y - HP[18] * rt__rz_ry_ex_hat[2] + noise_std[dim] * noise_std[dim];
        // h = m
        // e = (y - bm) - R' * m
        const float obs_error = mag_corr[dim] * h_inv - rt_rz_ry_ex[dim];

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
    
    regular_covariance_to_symmetric<ESKF::dim>(0);

    return info;                 
}

unsigned char LIEKF::fuse_declination(const Vector2f &dec, const Vector3f &w, const Vector3f &a, const Vector2f &noise_std, const Vector2f &gate) {
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
        array<float, ESKF::dim> HP = {_cov[0][index],
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
    
    regular_covariance_to_symmetric<ESKF::dim>(0);

    return info;
}

void LIEKF::correct_state() {
    // for (float &es : _error_state) {
    //     cout << es << endl;
    // }
    // cout << "-------------------------" << endl;
    // state: [p, v, bg, ba, g], R
    // error_state : [δp, δv, δθ, δbg, δba, δg]

    // q = q * Exp(δθ)
    Quaternionf delta_q;
    array<float, 3> delta_theta = {_error_state[6], _error_state[7], _error_state[8]};
    quaternion_from_axis_angle(delta_q, delta_theta);
    const Quaternionf q = _q;
    _q = q * delta_q;
    _rot = q;

    /*
    p = p + R*Exp(δθ)*δp
    v = v + R*Exp(δθ)*δv
    bg = bg + δbg
    ba = ba + δba
    g = g + δg
    m = m + δm
    bm = bm + δbm
    */
    for (unsigned char i = 0; i < 3; ++i) {
        _p[i] += _rot(i, 0) * _error_state[0] + _rot(i, 1) * _error_state[1] + _rot(i, 2) * _error_state[2];
        _v[i] += _rot(i, 0) * _error_state[3] + _rot(i, 1) * _error_state[4] + _rot(i, 2) * _error_state[5];
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

    // [δp, δv, δθ, δbg, δba, δg] = 0
    for (float &es : _error_state) {
        es = 0.f;
    }
}

void LIEKF::correct_covariance() {
    array<float, 58> var;

    var[0] = _cov[1][0] + _cov[1][1]*_error_state[8] - _cov[2][1]*_error_state[7];
    var[1] = _cov[2][0] + _cov[2][1]*_error_state[8] - _cov[2][2]*_error_state[7];
    var[2] = _cov[1][0]*_error_state[8];
    var[3] = _cov[2][0]*_error_state[7];
    var[4] = _cov[0][0] + var[2] - var[3];
    var[5] = _cov[4][0] + _cov[4][1]*_error_state[8] - _cov[4][2]*_error_state[7];
    var[6] = _cov[5][0] + _cov[5][1]*_error_state[8] - _cov[5][2]*_error_state[7];
    var[7] = _cov[3][0] + _cov[3][1]*_error_state[8] - _cov[3][2]*_error_state[7];
    var[8] = _cov[7][0] + _cov[7][1]*_error_state[8] - _cov[7][2]*_error_state[7];
    var[9] = 0.5F*_error_state[8];
    var[10] = _cov[8][0] + _cov[8][1]*_error_state[8] - _cov[8][2]*_error_state[7];
    var[11] = 0.5F*var[10];
    var[12] = _cov[6][0] + _cov[6][1]*_error_state[8] - _cov[6][2]*_error_state[7];
    var[13] = 0.5F*_error_state[7];
    var[14] = 0.5F*_error_state[6];
    var[15] = -_cov[2][0]*_error_state[8] + _cov[2][1] + _cov[2][2]*_error_state[6];
    var[16] = -_cov[0][0]*_error_state[8] + _cov[1][0] + _cov[2][0]*_error_state[6];
    var[17] = _cov[2][1]*_error_state[6];
    var[18] = _cov[1][1] + var[17] - var[2];
    var[19] = -_cov[4][0]*_error_state[8] + _cov[4][1] + _cov[4][2]*_error_state[6];
    var[20] = -_cov[5][0]*_error_state[8] + _cov[5][1] + _cov[5][2]*_error_state[6];
    var[21] = -_cov[3][0]*_error_state[8] + _cov[3][1] + _cov[3][2]*_error_state[6];
    var[22] = -_cov[7][0]*_error_state[8] + _cov[7][1] + _cov[7][2]*_error_state[6];
    var[23] = -_cov[8][0]*_error_state[8] + _cov[8][1] + _cov[8][2]*_error_state[6];
    var[24] = -_cov[6][0]*_error_state[8] + _cov[6][1] + _cov[6][2]*_error_state[6];
    var[25] = _cov[4][0]*_error_state[7] - _cov[4][1]*_error_state[6] + _cov[4][2];
    var[26] = _cov[5][0]*_error_state[7] - _cov[5][1]*_error_state[6] + _cov[5][2];
    var[27] = _cov[3][0]*_error_state[7] - _cov[3][1]*_error_state[6] + _cov[3][2];
    var[28] = _cov[7][0]*_error_state[7] - _cov[7][1]*_error_state[6] + _cov[7][2];
    var[29] = _cov[8][0]*_error_state[7] - _cov[8][1]*_error_state[6] + _cov[8][2];
    var[30] = _cov[6][0]*_error_state[7] - _cov[6][1]*_error_state[6] + _cov[6][2];
    var[31] = _cov[4][3] + _cov[4][4]*_error_state[8] - _cov[5][4]*_error_state[7];
    var[32] = _cov[5][3] + _cov[5][4]*_error_state[8] - _cov[5][5]*_error_state[7];
    var[33] = _cov[4][3]*_error_state[8];
    var[34] = _cov[5][3]*_error_state[7];
    var[35] = _cov[3][3] + var[33] - var[34];
    var[36] = _cov[7][3] + _cov[7][4]*_error_state[8] - _cov[7][5]*_error_state[7];
    var[37] = _cov[8][3] + _cov[8][4]*_error_state[8] - _cov[8][5]*_error_state[7];
    var[38] = _cov[6][3] + _cov[6][4]*_error_state[8] - _cov[6][5]*_error_state[7];
    var[39] = -_cov[5][3]*_error_state[8] + _cov[5][4] + _cov[5][5]*_error_state[6];
    var[40] = -_cov[3][3]*_error_state[8] + _cov[4][3] + _cov[5][3]*_error_state[6];
    var[41] = _cov[5][4]*_error_state[6];
    var[42] = _cov[4][4] - var[33] + var[41];
    var[43] = -_cov[7][3]*_error_state[8] + _cov[7][4] + _cov[7][5]*_error_state[6];
    var[44] = -_cov[8][3]*_error_state[8] + _cov[8][4] + _cov[8][5]*_error_state[6];
    var[45] = -_cov[6][3]*_error_state[8] + _cov[6][4] + _cov[6][5]*_error_state[6];
    var[46] = _cov[7][3]*_error_state[7] - _cov[7][4]*_error_state[6] + _cov[7][5];
    var[47] = _cov[8][3]*_error_state[7] - _cov[8][4]*_error_state[6] + _cov[8][5];
    var[48] = _cov[6][3]*_error_state[7] - _cov[6][4]*_error_state[6] + _cov[6][5];
    var[49] = _cov[7][6] + _cov[7][7]*var[9] - _cov[8][7]*var[13];
    var[50] = _cov[8][6] + _cov[8][7]*var[9] - _cov[8][8]*var[13];
    var[51] = _cov[7][6]*var[9];
    var[52] = _cov[8][6]*var[13];
    var[53] = _cov[6][6] + var[51] - var[52];
    var[54] = -_cov[8][6]*var[9] + _cov[8][7] + _cov[8][8]*var[14];
    var[55] = -_cov[6][6]*var[9] + _cov[7][6] + _cov[8][6]*var[14];
    var[56] = _cov[8][7]*var[14];
    var[57] = _cov[7][7] - var[51] + var[56];


    _cov[2][2] = _cov[2][2] - _error_state[6]*(_cov[1][0]*_error_state[7] - _cov[1][1]*_error_state[6] + _cov[2][1]) + _error_state[7]*(_cov[0][0]*_error_state[7] - _cov[1][0]*_error_state[6] + _cov[2][0]) - var[17] + var[3];
    _cov[0][0] = -_error_state[7]*var[1] + _error_state[8]*var[0] + var[4];
    _cov[0][1] = _error_state[6]*var[1] - _error_state[8]*var[4] + var[0];
    _cov[1][1] = _error_state[6]*var[15] - _error_state[8]*var[16] + var[18];
    _cov[0][2] = -_error_state[6]*var[0] + _error_state[7]*var[4] + var[1];
    _cov[1][2] = -_error_state[6]*var[18] + _error_state[7]*var[16] + var[15];
    
    _cov[0][3] = -_error_state[7]*var[6] + _error_state[8]*var[5] + var[7];
    _cov[1][3] = -_error_state[7]*var[20] + _error_state[8]*var[19] + var[21];
    _cov[2][3] = -_error_state[7]*var[26] + _error_state[8]*var[25] + var[27];

    _cov[5][5] = _cov[5][5] - _error_state[6]*(_cov[4][3]*_error_state[7] - _cov[4][4]*_error_state[6] + _cov[5][4]) + _error_state[7]*(_cov[3][3]*_error_state[7] - _cov[4][3]*_error_state[6] + _cov[5][3]) + var[34] - var[41];
    _cov[3][3] = -_error_state[7]*var[32] + _error_state[8]*var[31] + var[35];
    _cov[0][4] = _error_state[6]*var[6] - _error_state[8]*var[7] + var[5];
    _cov[1][4] = _error_state[6]*var[20] - _error_state[8]*var[21] + var[19];
    _cov[2][4] = _error_state[6]*var[26] - _error_state[8]*var[27] + var[25];
    _cov[3][4] = _error_state[6]*var[32] - _error_state[8]*var[35] + var[31];
    _cov[4][4] = _error_state[6]*var[39] - _error_state[8]*var[40] + var[42];
    _cov[0][5] = -_error_state[6]*var[5] + _error_state[7]*var[7] + var[6];
    _cov[1][5] = -_error_state[6]*var[19] + _error_state[7]*var[21] + var[20];
    _cov[2][5] = -_error_state[6]*var[25] + _error_state[7]*var[27] + var[26];
    _cov[3][5] = -_error_state[6]*var[31] + _error_state[7]*var[35] + var[32];
    _cov[4][5] = -_error_state[6]*var[42] + _error_state[7]*var[40] + var[39];
    
    _cov[0][6] = -_error_state[7]*var[11] + var[12] + var[8]*var[9];
    _cov[1][6] = -var[13]*var[23] + var[22]*var[9] + var[24];
    _cov[2][6] = -var[13]*var[29] + var[28]*var[9] + var[30];
    _cov[3][6] = -var[13]*var[37] + var[36]*var[9] + var[38];
    _cov[4][6] = -var[13]*var[44] + var[43]*var[9] + var[45];
    _cov[5][6] = -var[13]*var[47] + var[46]*var[9] + var[48];

    _cov[8][8] = _cov[8][8] + var[13]*(_cov[6][6]*var[13] - _cov[7][6]*var[14] + _cov[8][6]) - var[14]*(_cov[7][6]*var[13] - _cov[7][7]*var[14] + _cov[8][7]) + var[52] - var[56];
    _cov[6][6] = -var[13]*var[50] + var[49]*var[9] + var[53];
    _cov[0][7] = _error_state[6]*var[11] - var[12]*var[9] + var[8];
    _cov[1][7] = var[14]*var[23] + var[22] - var[24]*var[9];
    _cov[2][7] = var[14]*var[29] + var[28] - var[30]*var[9];
    _cov[3][7] = var[14]*var[37] + var[36] - var[38]*var[9];
    _cov[4][7] = var[14]*var[44] + var[43] - var[45]*var[9];
    _cov[5][7] = var[14]*var[47] + var[46] - var[48]*var[9];
    _cov[6][7] = var[14]*var[50] + var[49] - var[53]*var[9];
    _cov[7][7] = var[14]*var[54] - var[55]*var[9] + var[57];
    _cov[0][8] = var[10] + var[12]*var[13] - var[14]*var[8];
    _cov[1][8] = var[13]*var[24] - var[14]*var[22] + var[23];
    _cov[2][8] = var[13]*var[30] - var[14]*var[28] + var[29];
    _cov[3][8] = var[13]*var[38] - var[14]*var[36] + var[37];
    _cov[4][8] = var[13]*var[45] - var[14]*var[43] + var[44];
    _cov[5][8] = var[13]*var[48] - var[14]*var[46] + var[47];
    _cov[6][8] = var[13]*var[53] - var[14]*var[49] + var[50];
    _cov[7][8] = var[13]*var[55] - var[14]*var[57] + var[54];
    
    _cov[0][9] = _cov[9][0] + _cov[9][1]*_error_state[8] - _cov[9][2]*_error_state[7];
    _cov[1][9] = -_cov[9][0]*_error_state[8] + _cov[9][1] + _cov[9][2]*_error_state[6];
    _cov[2][9] = _cov[9][0]*_error_state[7] - _cov[9][1]*_error_state[6] + _cov[9][2];
    _cov[3][9] = _cov[9][3] + _cov[9][4]*_error_state[8] - _cov[9][5]*_error_state[7];
    _cov[4][9] = -_cov[9][3]*_error_state[8] + _cov[9][4] + _cov[9][5]*_error_state[6];
    _cov[5][9] = _cov[9][3]*_error_state[7] - _cov[9][4]*_error_state[6] + _cov[9][5];
    _cov[6][9] = _cov[9][6] + _cov[9][7]*var[9] - _cov[9][8]*var[13];
    _cov[7][9] = -_cov[9][6]*var[9] + _cov[9][7] + _cov[9][8]*var[14];
    _cov[8][9] = _cov[9][6]*var[13] - _cov[9][7]*var[14] + _cov[9][8];
    // _cov[9][9] = _cov[9][9];
    _cov[0][10] = _cov[10][0] + _cov[10][1]*_error_state[8] - _cov[10][2]*_error_state[7];
    _cov[1][10] = -_cov[10][0]*_error_state[8] + _cov[10][1] + _cov[10][2]*_error_state[6];
    _cov[2][10] = _cov[10][0]*_error_state[7] - _cov[10][1]*_error_state[6] + _cov[10][2];
    _cov[3][10] = _cov[10][3] + _cov[10][4]*_error_state[8] - _cov[10][5]*_error_state[7];
    _cov[4][10] = -_cov[10][3]*_error_state[8] + _cov[10][4] + _cov[10][5]*_error_state[6];
    _cov[5][10] = _cov[10][3]*_error_state[7] - _cov[10][4]*_error_state[6] + _cov[10][5];
    _cov[6][10] = _cov[10][6] + _cov[10][7]*var[9] - _cov[10][8]*var[13];
    _cov[7][10] = -_cov[10][6]*var[9] + _cov[10][7] + _cov[10][8]*var[14];
    _cov[8][10] = _cov[10][6]*var[13] - _cov[10][7]*var[14] + _cov[10][8];
    // _cov[9][10] = _cov[10][9];
    // _cov[10][10] = _cov[10][10];
    _cov[0][11] = _cov[11][0] + _cov[11][1]*_error_state[8] - _cov[11][2]*_error_state[7];
    _cov[1][11] = -_cov[11][0]*_error_state[8] + _cov[11][1] + _cov[11][2]*_error_state[6];
    _cov[2][11] = _cov[11][0]*_error_state[7] - _cov[11][1]*_error_state[6] + _cov[11][2];
    _cov[3][11] = _cov[11][3] + _cov[11][4]*_error_state[8] - _cov[11][5]*_error_state[7];
    _cov[4][11] = -_cov[11][3]*_error_state[8] + _cov[11][4] + _cov[11][5]*_error_state[6];
    _cov[5][11] = _cov[11][3]*_error_state[7] - _cov[11][4]*_error_state[6] + _cov[11][5];
    _cov[6][11] = _cov[11][6] + _cov[11][7]*var[9] - _cov[11][8]*var[13];
    _cov[7][11] = -_cov[11][6]*var[9] + _cov[11][7] + _cov[11][8]*var[14];
    _cov[8][11] = _cov[11][6]*var[13] - _cov[11][7]*var[14] + _cov[11][8];
    // _cov[9][11] = _cov[11][9];
    // _cov[10][11] = _cov[11][10];
    // _cov[11][11] = _cov[11][11];
    _cov[0][12] = _cov[12][0] + _cov[12][1]*_error_state[8] - _cov[12][2]*_error_state[7];
    _cov[1][12] = -_cov[12][0]*_error_state[8] + _cov[12][1] + _cov[12][2]*_error_state[6];
    _cov[2][12] = _cov[12][0]*_error_state[7] - _cov[12][1]*_error_state[6] + _cov[12][2];
    _cov[3][12] = _cov[12][3] + _cov[12][4]*_error_state[8] - _cov[12][5]*_error_state[7];
    _cov[4][12] = -_cov[12][3]*_error_state[8] + _cov[12][4] + _cov[12][5]*_error_state[6];
    _cov[5][12] = _cov[12][3]*_error_state[7] - _cov[12][4]*_error_state[6] + _cov[12][5];
    _cov[6][12] = _cov[12][6] + _cov[12][7]*var[9] - _cov[12][8]*var[13];
    _cov[7][12] = -_cov[12][6]*var[9] + _cov[12][7] + _cov[12][8]*var[14];
    _cov[8][12] = _cov[12][6]*var[13] - _cov[12][7]*var[14] + _cov[12][8];
    // _cov[9][12] = _cov[12][9];
    // _cov[10][12] = _cov[12][10];
    // _cov[11][12] = _cov[12][11];
    // _cov[12][12] = _cov[12][12];
    _cov[0][13] = _cov[13][0] + _cov[13][1]*_error_state[8] - _cov[13][2]*_error_state[7];
    _cov[1][13] = -_cov[13][0]*_error_state[8] + _cov[13][1] + _cov[13][2]*_error_state[6];
    _cov[2][13] = _cov[13][0]*_error_state[7] - _cov[13][1]*_error_state[6] + _cov[13][2];
    _cov[3][13] = _cov[13][3] + _cov[13][4]*_error_state[8] - _cov[13][5]*_error_state[7];
    _cov[4][13] = -_cov[13][3]*_error_state[8] + _cov[13][4] + _cov[13][5]*_error_state[6];
    _cov[5][13] = _cov[13][3]*_error_state[7] - _cov[13][4]*_error_state[6] + _cov[13][5];
    _cov[6][13] = _cov[13][6] + _cov[13][7]*var[9] - _cov[13][8]*var[13];
    _cov[7][13] = -_cov[13][6]*var[9] + _cov[13][7] + _cov[13][8]*var[14];
    _cov[8][13] = _cov[13][6]*var[13] - _cov[13][7]*var[14] + _cov[13][8];
    // _cov[9][13] = _cov[13][9];
    // _cov[10][13] = _cov[13][10];
    // _cov[11][13] = _cov[13][11];
    // _cov[12][13] = _cov[13][12];
    // _cov[13][13] = _cov[13][13];
    _cov[0][14] = _cov[14][0] + _cov[14][1]*_error_state[8] - _cov[14][2]*_error_state[7];
    _cov[1][14] = -_cov[14][0]*_error_state[8] + _cov[14][1] + _cov[14][2]*_error_state[6];
    _cov[2][14] = _cov[14][0]*_error_state[7] - _cov[14][1]*_error_state[6] + _cov[14][2];
    _cov[3][14] = _cov[14][3] + _cov[14][4]*_error_state[8] - _cov[14][5]*_error_state[7];
    _cov[4][14] = -_cov[14][3]*_error_state[8] + _cov[14][4] + _cov[14][5]*_error_state[6];
    _cov[5][14] = _cov[14][3]*_error_state[7] - _cov[14][4]*_error_state[6] + _cov[14][5];
    _cov[6][14] = _cov[14][6] + _cov[14][7]*var[9] - _cov[14][8]*var[13];
    _cov[7][14] = -_cov[14][6]*var[9] + _cov[14][7] + _cov[14][8]*var[14];
    _cov[8][14] = _cov[14][6]*var[13] - _cov[14][7]*var[14] + _cov[14][8];
    // _cov[9][14] = _cov[14][9];
    // _cov[10][14] = _cov[14][10];
    // _cov[11][14] = _cov[14][11];
    // _cov[12][14] = _cov[14][12];
    // _cov[13][14] = _cov[14][13];
    // _cov[14][14] = _cov[14][14];
    _cov[0][15] = _cov[15][0] + _cov[15][1]*_error_state[8] - _cov[15][2]*_error_state[7];
    _cov[1][15] = -_cov[15][0]*_error_state[8] + _cov[15][1] + _cov[15][2]*_error_state[6];
    _cov[2][15] = _cov[15][0]*_error_state[7] - _cov[15][1]*_error_state[6] + _cov[15][2];
    _cov[3][15] = _cov[15][3] + _cov[15][4]*_error_state[8] - _cov[15][5]*_error_state[7];
    _cov[4][15] = -_cov[15][3]*_error_state[8] + _cov[15][4] + _cov[15][5]*_error_state[6];
    _cov[5][15] = _cov[15][3]*_error_state[7] - _cov[15][4]*_error_state[6] + _cov[15][5];
    _cov[6][15] = _cov[15][6] + _cov[15][7]*var[9] - _cov[15][8]*var[13];
    _cov[7][15] = -_cov[15][6]*var[9] + _cov[15][7] + _cov[15][8]*var[14];
    _cov[8][15] = _cov[15][6]*var[13] - _cov[15][7]*var[14] + _cov[15][8];
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