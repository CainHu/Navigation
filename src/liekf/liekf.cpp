//
// Created by Cain on 2022/11/11.
//

#include "liekf.h"
#include <cfloat>
#include <iostream>

using namespace std;
using namespace liekf;

void LIEKF::predict_covariance(const Vector3f &w, const Vector3f &a) {
    const float dang0 = (w[0] - _bg[0]) * _dt, dang1 = (w[1] - _bg[1]) * _dt, dang2 = (w[2] - _bg[2]) * _dt;
    const float dv0 = (a[0] - _ba[0]) * _dt, dv1 = (a[1] - _ba[1]) * _dt, dv2 = (a[2] - _ba[2]) * _dt;
    const float r0 = _rot(2, 0) * _dt,  r1 = _rot(2, 1) * _dt, r2 = _rot(2, 2) * _dt;

    array<array<float, 3>, 3> XP11 {
        {{_cov[0][0] - _cov[2][0]*dang1 + _cov[1][0]*dang2,
          _cov[0][1] - _cov[2][1]*dang1 + _cov[1][1]*dang2,
          _cov[0][2] - _cov[2][2]*dang1 + _cov[1][2]*dang2},
         {_cov[1][0] + _cov[2][0]*dang0 - _cov[0][0]*dang2,
          _cov[1][1] + _cov[2][1]*dang0 - _cov[0][1]*dang2,
          _cov[1][2] + _cov[2][2]*dang0 - _cov[0][2]*dang2},
         {_cov[2][0] - _cov[1][0]*dang0 + _cov[0][0]*dang1,
          _cov[2][1] - _cov[1][1]*dang0 + _cov[0][1]*dang1,
          _cov[2][2] - _cov[1][2]*dang0 + _cov[0][2]*dang1}}
    };
    array<array<float, 3>, 3> XP12 {
        {{_cov[0][3] - _cov[2][3]*dang1 + _cov[1][3]*dang2,
          _cov[0][4] - _cov[2][4]*dang1 + _cov[1][4]*dang2,
          _cov[0][5] - _cov[2][5]*dang1 + _cov[1][5]*dang2},
         {_cov[1][3] + _cov[2][3]*dang0 - _cov[0][3]*dang2,
          _cov[1][4] + _cov[2][4]*dang0 - _cov[0][4]*dang2,
          _cov[1][5] + _cov[2][5]*dang0 - _cov[0][5]*dang2},
         {_cov[2][3] - _cov[1][3]*dang0 + _cov[0][3]*dang1,
          _cov[2][4] - _cov[1][4]*dang0 + _cov[0][4]*dang1,
          _cov[2][5] - _cov[1][5]*dang0 + _cov[0][5]*dang1}}
    };
    array<array<float, 3>, 3> XP13 {
        {{_cov[0][6] - _cov[2][6]*dang1 + _cov[1][6]*dang2,
          _cov[0][7] - _cov[2][7]*dang1 + _cov[1][7]*dang2,
          _cov[0][8] - _cov[2][8]*dang1 + _cov[1][8]*dang2},
         {_cov[1][6] + _cov[2][6]*dang0 - _cov[0][6]*dang2,
          _cov[1][7] + _cov[2][7]*dang0 - _cov[0][7]*dang2,
          _cov[1][8] + _cov[2][8]*dang0 - _cov[0][8]*dang2},
         {_cov[2][6] - _cov[1][6]*dang0 + _cov[0][6]*dang1,
          _cov[2][7] - _cov[1][7]*dang0 + _cov[0][7]*dang1,
          _cov[2][8] - _cov[1][8]*dang0 + _cov[0][8]*dang1}}
    };
    array<array<float, 3>, 3> XP14 {
        {{_cov[0][9] - _cov[2][9]*dang1 + _cov[1][9]*dang2,
          _cov[0][10] - _cov[2][10]*dang1 + _cov[1][10]*dang2,
          _cov[0][11] - _cov[2][11]*dang1 + _cov[1][11]*dang2},
         {_cov[1][9] + _cov[2][9]*dang0 - _cov[0][9]*dang2,
          _cov[1][10] + _cov[2][10]*dang0 - _cov[0][10]*dang2,
          _cov[1][11] + _cov[2][11]*dang0 - _cov[0][11]*dang2},
         {_cov[2][9] - _cov[1][9]*dang0 + _cov[0][9]*dang1,
          _cov[2][10] - _cov[1][10]*dang0 + _cov[0][10]*dang1,
          _cov[2][11] - _cov[1][11]*dang0 + _cov[0][11]*dang1}}
    };
    array<array<float, 3>, 3> XP15 {
        {{_cov[0][12] - _cov[2][12]*dang1 + _cov[1][12]*dang2,
          _cov[0][13] - _cov[2][13]*dang1 + _cov[1][13]*dang2,
          _cov[0][14] - _cov[2][14]*dang1 + _cov[1][14]*dang2},
         {_cov[1][12] + _cov[2][12]*dang0 - _cov[0][12]*dang2,
          _cov[1][13] + _cov[2][13]*dang0 - _cov[0][13]*dang2,
          _cov[1][14] + _cov[2][14]*dang0 - _cov[0][14]*dang2},
         {_cov[2][12] - _cov[1][12]*dang0 + _cov[0][12]*dang1,
          _cov[2][13] - _cov[1][13]*dang0 + _cov[0][13]*dang1,
          _cov[2][14] - _cov[1][14]*dang0 + _cov[0][14]*dang1}}
    };
    array<float, 3> XP16 {
        _cov[0][15] - _cov[2][15]*dang1 + _cov[1][15]*dang2,
        _cov[1][15] + _cov[2][15]*dang0 - _cov[0][15]*dang2,
        _cov[2][15] - _cov[1][15]*dang0 + _cov[0][15]*dang1
    };
    float XP11XT00 = XP11[0][0] - XP11[0][2]*dang1 + XP11[0][1]*dang2;
    float XP11XT01 = XP11[0][1] + XP11[0][2]*dang0 - XP11[0][0]*dang2;
    float XP11XT02 = XP11[0][2] - XP11[0][1]*dang0 + XP11[0][0]*dang1;
    float XP11XT11 = XP11[1][1] + XP11[1][2]*dang0 - XP11[1][0]*dang2;
    float XP11XT12 = XP11[1][2] - XP11[1][1]*dang0 + XP11[1][0]*dang1;
    float XP11XT22 = XP11[2][2] - XP11[2][1]*dang0 + XP11[2][0]*dang1;

    _cov[0][0] = XP11XT00 + (XP12[0][0] + XP12[0][0] + _cov[3][3] * _dt) * _dt;
    _cov[0][1] = XP11XT01 + (XP12[0][1] + XP12[1][0] + _cov[3][4] * _dt) * _dt;
    _cov[0][2] = XP11XT02 + (XP12[0][2] + XP12[2][0] + _cov[3][5] * _dt) * _dt;
    _cov[1][1] = XP11XT11 + (XP12[1][1] + XP12[1][1] + _cov[4][4] * _dt) * _dt;
    _cov[1][2] = XP11XT12 + (XP12[1][2] + XP12[2][1] + _cov[4][5] * _dt) * _dt;
    _cov[2][2] = XP11XT22 + (XP12[2][2] + XP12[2][2] + _cov[5][5] * _dt) * _dt;

    _cov[0][15] = XP16[0] + _cov[3][15] * _dt;
    _cov[1][15] = XP16[1] + _cov[4][15] * _dt;
    _cov[2][15] = XP16[2] + _cov[5][15] * _dt;

    _cov[0][12] = XP15[0][0] + _cov[3][12] * _dt;
    _cov[0][13] = XP15[0][1] + _cov[3][13] * _dt;
    _cov[0][14] = XP15[0][2] + _cov[3][14] * _dt;
    _cov[1][12] = XP15[1][0] + _cov[4][12] * _dt;
    _cov[1][13] = XP15[1][1] + _cov[4][13] * _dt;
    _cov[1][14] = XP15[1][2] + _cov[4][14] * _dt;
    _cov[2][12] = XP15[2][0] + _cov[5][12] * _dt;
    _cov[2][13] = XP15[2][1] + _cov[5][13] * _dt;
    _cov[2][14] = XP15[2][2] + _cov[5][14] * _dt;

    _cov[0][9] = XP14[0][0] + _cov[3][9] * _dt;
    _cov[0][10] = XP14[0][1] + _cov[3][10] * _dt;
    _cov[0][11] = XP14[0][2] + _cov[3][11] * _dt;
    _cov[1][9] = XP14[1][0] + _cov[4][9] * _dt;
    _cov[1][10] = XP14[1][1] + _cov[4][10] * _dt;
    _cov[1][11] = XP14[1][2] + _cov[4][11] * _dt;
    _cov[2][9] = XP14[2][0] + _cov[5][9] * _dt;
    _cov[2][10] = XP14[2][1] + _cov[5][10] * _dt;
    _cov[2][11] = XP14[2][2] + _cov[5][11] * _dt;

    const array<array<float, 3>, 3> c13 {
        {{XP13[0][0] + _cov[3][6] * _dt, XP13[0][1] + _cov[3][7] * _dt, XP13[0][2] + _cov[3][8] * _dt},
         {XP13[1][0] + _cov[4][6] * _dt, XP13[1][1] + _cov[4][7] * _dt, XP13[1][2] + _cov[4][8] * _dt},
         {XP13[2][0] + _cov[5][6] * _dt, XP13[2][1] + _cov[5][7] * _dt, XP13[2][2] + _cov[5][8] * _dt}}
    };

    const array<array<float, 3>, 3> c12 {
        {{XP12[0][0] + _cov[3][3] * _dt, XP12[0][1] + _cov[3][4] * _dt, XP12[0][2] + _cov[3][5] * _dt},
         {XP12[1][0] + _cov[4][3] * _dt, XP12[1][1] + _cov[4][4] * _dt, XP12[1][2] + _cov[4][5] * _dt},
         {XP12[2][0] + _cov[5][3] * _dt, XP12[2][1] + _cov[5][4] * _dt, XP12[2][2] + _cov[5][5] * _dt}}
    };

    _cov[0][6] = c13[0][0] - c13[0][2]*dang1 + c13[0][1]*dang2 - _cov[0][9] * _dt;
    _cov[0][7] = c13[0][1] + c13[0][2]*dang0 - c13[0][0]*dang2 - _cov[0][10] * _dt;
    _cov[0][8] = c13[0][2] - c13[0][1]*dang0 + c13[0][0]*dang1 - _cov[0][11] * _dt;
    _cov[1][6] = c13[1][0] - c13[1][2]*dang1 + c13[1][1]*dang2 - _cov[1][9] * _dt;
    _cov[1][7] = c13[1][1] + c13[1][2]*dang0 - c13[1][0]*dang2 - _cov[1][10] * _dt;
    _cov[1][8] = c13[1][2] - c13[1][1]*dang0 + c13[1][0]*dang1 - _cov[1][11] * _dt;
    _cov[2][6] = c13[2][0] - c13[2][2]*dang1 + c13[2][1]*dang2 - _cov[2][9] * _dt;
    _cov[2][7] = c13[2][1] + c13[2][2]*dang0 - c13[2][0]*dang2 - _cov[2][10] * _dt;
    _cov[2][8] = c13[2][2] - c13[2][1]*dang0 + c13[2][0]*dang1 - _cov[2][11] * _dt;

    _cov[0][3] = c12[0][0] - c12[0][2]*dang1 + c12[0][1]*dang2 - c13[0][2]*dv1 + c13[0][1]*dv2 - _cov[0][12] * _dt + _cov[0][15] * r0;
    _cov[0][4] = c12[0][1] + c12[0][2]*dang0 - c12[0][0]*dang2 + c13[0][2]*dv0 - c13[0][0]*dv2 - _cov[0][13] * _dt + _cov[0][15] * r1;
    _cov[0][5] = c12[0][2] - c12[0][1]*dang0 + c12[0][0]*dang1 - c13[0][1]*dv0 + c13[0][0]*dv1 - _cov[0][14] * _dt + _cov[0][15] * r2;
    _cov[1][3] = c12[1][0] - c12[1][2]*dang1 + c12[1][1]*dang2 - c13[1][2]*dv1 + c13[1][1]*dv2 - _cov[1][12] * _dt + _cov[1][15] * r0;
    _cov[1][4] = c12[1][1] + c12[1][2]*dang0 - c12[1][0]*dang2 + c13[1][2]*dv0 - c13[1][0]*dv2 - _cov[1][13] * _dt + _cov[1][15] * r1;
    _cov[1][5] = c12[1][2] - c12[1][1]*dang0 + c12[1][0]*dang1 - c13[1][1]*dv0 + c13[1][0]*dv1 - _cov[1][14] * _dt + _cov[1][15] * r2;
    _cov[2][3] = c12[2][0] - c12[2][2]*dang1 + c12[2][1]*dang2 - c13[2][2]*dv1 + c13[2][1]*dv2 - _cov[2][12] * _dt + _cov[2][15] * r0;
    _cov[2][4] = c12[2][1] + c12[2][2]*dang0 - c12[2][0]*dang2 + c13[2][2]*dv0 - c13[2][0]*dv2 - _cov[2][13] * _dt + _cov[2][15] * r1;
    _cov[2][5] = c12[2][2] - c12[2][1]*dang0 + c12[2][0]*dang1 - c13[2][1]*dv0 + c13[2][0]*dv1 - _cov[2][14] * _dt + _cov[2][15] * r2;

    // XP26
    const array<float, 3> XP26 {
        _cov[3][15] - _cov[5][15]*dang1 + _cov[4][15]*dang2,
        _cov[4][15] + _cov[5][15]*dang0 - _cov[3][15]*dang2,
        _cov[5][15] - _cov[4][15]*dang0 + _cov[3][15]*dang1
    };

    // XP25
    const array<array<float, 3>, 3> XP25 {
        {{_cov[3][12] - _cov[5][12]*dang1 + _cov[4][12]*dang2,
          _cov[3][13] - _cov[5][13]*dang1 + _cov[4][13]*dang2,
          _cov[3][14] - _cov[5][14]*dang1 + _cov[4][14]*dang2},
         {_cov[4][12] + _cov[5][12]*dang0 - _cov[3][12]*dang2,
          _cov[4][13] + _cov[5][13]*dang0 - _cov[3][13]*dang2,
          _cov[4][14] + _cov[5][14]*dang0 - _cov[3][14]*dang2},
         {_cov[5][12] - _cov[4][12]*dang0 + _cov[3][12]*dang1,
          _cov[5][13] - _cov[4][13]*dang0 + _cov[3][13]*dang1,
          _cov[5][14] - _cov[4][14]*dang0 + _cov[3][14]*dang1}}
    };

    // XP24
    const array<array<float, 3>, 3> XP24 {
        {{_cov[3][9] - _cov[5][9]*dang1 + _cov[4][9]*dang2,
          _cov[3][10] - _cov[5][10]*dang1 + _cov[4][10]*dang2,
          _cov[3][11] - _cov[5][11]*dang1 + _cov[4][11]*dang2},
         {_cov[4][9] + _cov[5][9]*dang0 - _cov[3][9]*dang2,
          _cov[4][10] + _cov[5][10]*dang0 - _cov[3][10]*dang2,
          _cov[4][11] + _cov[5][11]*dang0 - _cov[3][11]*dang2},
         {_cov[5][9] - _cov[4][9]*dang0 + _cov[3][9]*dang1,
          _cov[5][10] - _cov[4][10]*dang0 + _cov[3][10]*dang1,
          _cov[5][11] - _cov[4][11]*dang0 + _cov[3][11]*dang1}}
    };

    // XP23
    const array<array<float, 3>, 3> XP23 {
        {{_cov[3][6] - _cov[5][6]*dang1 + _cov[4][6]*dang2,
          _cov[3][7] - _cov[5][7]*dang1 + _cov[4][7]*dang2,
          _cov[3][8] - _cov[5][8]*dang1 + _cov[4][8]*dang2},
         {_cov[4][6] + _cov[5][6]*dang0 - _cov[3][6]*dang2,
          _cov[4][7] + _cov[5][7]*dang0 - _cov[3][7]*dang2,
          _cov[4][8] + _cov[5][8]*dang0 - _cov[3][8]*dang2},
         {_cov[5][6] - _cov[4][6]*dang0 + _cov[3][6]*dang1,
          _cov[5][7] - _cov[4][7]*dang0 + _cov[3][7]*dang1,
          _cov[5][8] - _cov[4][8]*dang0 + _cov[3][8]*dang1}}
    };

    // XP22
    const array<array<float, 3>, 3> XP22 {
        {{_cov[3][3] - _cov[5][3]*dang1 + _cov[4][3]*dang2,
          _cov[3][4] - _cov[5][4]*dang1 + _cov[4][4]*dang2,
          _cov[3][5] - _cov[5][5]*dang1 + _cov[4][5]*dang2},
         {_cov[4][3] + _cov[5][3]*dang0 - _cov[3][3]*dang2,
          _cov[4][4] + _cov[5][4]*dang0 - _cov[3][4]*dang2,
          _cov[4][5] + _cov[5][5]*dang0 - _cov[3][5]*dang2},
         {_cov[5][3] - _cov[4][3]*dang0 + _cov[3][3]*dang1,
          _cov[5][4] - _cov[4][4]*dang0 + _cov[3][4]*dang1,
          _cov[5][5] - _cov[4][5]*dang0 + _cov[3][5]*dang1}}
    };

    // YP36
    const array<float, 3> YP36 {
        -(_cov[8][15]*dv1) + _cov[7][15]*dv2,
        _cov[8][15]*dv0 - _cov[6][15]*dv2
        -(_cov[7][15]*dv0) + _cov[6][15]*dv1
    };

    // YP35
    const array<array<float, 3>, 3> YP35 {
        {{-(_cov[8][12]*dv1) + _cov[7][12]*dv2,
          -(_cov[8][13]*dv1) + _cov[7][13]*dv2,
          -(_cov[8][14]*dv1) + _cov[7][14]*dv2},
         {_cov[8][12]*dv0 - _cov[6][12]*dv2,
          _cov[8][13]*dv0 - _cov[6][13]*dv2,
          _cov[8][14]*dv0 - _cov[6][14]*dv2},
         {-(_cov[7][12]*dv0) + _cov[6][12]*dv1,
          -(_cov[7][13]*dv0) + _cov[6][13]*dv1,
          -(_cov[7][14]*dv0) + _cov[6][14]*dv1}}
    };

    // YP34
    const array<array<float, 3>, 3> YP34 {
        {{-(_cov[8][9]*dv1) + _cov[7][9]*dv2,
          -(_cov[8][10]*dv1) + _cov[7][10]*dv2,
          -(_cov[8][11]*dv1) + _cov[7][11]*dv2},
         {_cov[8][9]*dv0 - _cov[6][9]*dv2,
          _cov[8][10]*dv0 - _cov[6][10]*dv2,
          _cov[8][11]*dv0 - _cov[6][11]*dv2},
         {-(_cov[7][9]*dv0) + _cov[6][9]*dv1,
          -(_cov[7][10]*dv0) + _cov[6][10]*dv1,
          -(_cov[7][11]*dv0) + _cov[6][11]*dv1}}
    };

    // YP33
    const array<array<float, 3>, 3> YP33 {
        {{-(_cov[8][6]*dv1) + _cov[7][6]*dv2,
          -(_cov[8][7]*dv1) + _cov[7][7]*dv2,
          -(_cov[8][8]*dv1) + _cov[7][8]*dv2},
         {_cov[8][6]*dv0 - _cov[6][6]*dv2,
          _cov[8][7]*dv0 - _cov[6][7]*dv2,
          _cov[8][8]*dv0 - _cov[6][8]*dv2},
         {-(_cov[7][6]*dv0) + _cov[6][6]*dv1,
          -(_cov[7][7]*dv0) + _cov[6][7]*dv1,
          -(_cov[7][8]*dv0) + _cov[6][8]*dv1}}
    };

    const array<array<float, 3>, 3> YP32 {
        {{-(_cov[3][8]*dv1) + _cov[3][7]*dv2,
          -(_cov[4][8]*dv1) + _cov[4][7]*dv2,
          -(_cov[5][8]*dv1) + _cov[5][7]*dv2},
         {_cov[3][8]*dv0 - _cov[3][6]*dv2,
          _cov[4][8]*dv0 - _cov[4][6]*dv2,
          _cov[5][8]*dv0 - _cov[5][6]*dv2},
         {-(_cov[3][7]*dv0) + _cov[3][6]*dv1,
          -(_cov[4][7]*dv0) + _cov[4][6]*dv1,
          -(_cov[5][7]*dv0) + _cov[5][6]*dv1}}
    };

    // XP22XT
    const float XP22XT00 = XP22[0][0] - XP22[0][2]*dang1 + XP22[0][1]*dang2;
    const float XP22XT01 = XP22[0][1] + XP22[0][2]*dang0 - XP22[0][0]*dang2;
    const float XP22XT02 = XP22[0][2] - XP22[0][1]*dang0 + XP22[0][0]*dang1;
    const float XP22XT11 = XP22[1][1] + XP22[1][2]*dang0 - XP22[1][0]*dang2;
    const float XP22XT12 = XP22[1][2] - XP22[1][1]*dang0 + XP22[1][0]*dang1;
    const float XP22XT22 = XP22[2][2] - XP22[2][1]*dang0 + XP22[2][0]*dang1;

    // XP23YT
    const array<array<float, 3>, 3> XP23YT {
        {{-(XP23[0][2]*dv1) + XP23[0][1]*dv2,
          XP23[0][2]*dv0 - XP23[0][0]*dv2,
          -(XP23[0][1]*dv0) + XP23[0][0]*dv1},
         {-(XP23[1][2]*dv1) + XP23[1][1]*dv2,
          XP23[1][2]*dv0 - XP23[1][0]*dv2,
          -(XP23[1][1]*dv0) + XP23[1][0]*dv1},
         {-(XP23[2][2]*dv1) + XP23[2][1]*dv2,
          XP23[2][2]*dv0 - XP23[2][0]*dv2,
          -(XP23[2][1]*dv0) + XP23[2][0]*dv1}}
    };

    // YP33YT
    const float YP33YT00 = -(YP33[0][2]*dv1) + YP33[0][1]*dv2;
    const float YP33YT01 = YP33[0][2]*dv0 - YP33[0][0]*dv2;
    const float YP33YT02 = -(YP33[0][1]*dv0) + YP33[0][0]*dv1;
    const float YP33YT11 = YP33[1][2]*dv0 - YP33[1][0]*dv2;
    const float YP33YT12 = -(YP33[1][1]*dv0) + YP33[1][0]*dv1;
    const float YP33YT22 = -(YP33[2][1]*dv0) + YP33[2][0]*dv1;

    _cov[3][15] = XP26[0] + YP36[0] - _cov[12][15] * _dt + r0 * _cov[15][15];
    _cov[4][15] = XP26[1] + YP36[1] - _cov[13][15] * _dt + r1 * _cov[15][15];
    _cov[5][15] = XP26[2] + YP36[2] - _cov[14][15] * _dt + r2 * _cov[15][15];

    _cov[3][12] = XP25[0][0] + YP35[0][0] - _cov[12][12] * _dt + r0 * _cov[12][15];
    _cov[3][13] = XP25[0][1] + YP35[0][1] - _cov[12][13] * _dt + r0 * _cov[13][15];
    _cov[3][14] = XP25[0][2] + YP35[0][2] - _cov[12][14] * _dt + r0 * _cov[14][15];
    _cov[4][12] = XP25[1][0] + YP35[1][0] - _cov[12][13] * _dt + r1 * _cov[12][15];
    _cov[4][13] = XP25[1][1] + YP35[1][1] - _cov[13][13] * _dt + r1 * _cov[13][15];
    _cov[4][14] = XP25[1][2] + YP35[1][2] - _cov[13][14] * _dt + r1 * _cov[14][15];
    _cov[5][12] = XP25[2][0] + YP35[2][0] - _cov[12][14] * _dt + r2 * _cov[12][15];
    _cov[5][13] = XP25[2][1] + YP35[2][1] - _cov[13][14] * _dt + r2 * _cov[13][15];
    _cov[5][14] = XP25[2][2] + YP35[2][2] - _cov[14][14] * _dt + r2 * _cov[14][15];

    _cov[3][9] = XP24[0][0] + YP34[0][0] - _cov[12][9] * _dt + r0 * _cov[15][9];
    _cov[3][10] = XP24[0][1] + YP34[0][1] - _cov[12][10] * _dt + r0 * _cov[15][10];
    _cov[3][11] = XP24[0][2] + YP34[0][2] - _cov[12][11] * _dt + r0 * _cov[15][11];
    _cov[4][9] = XP24[1][0] + YP34[1][0] - _cov[13][9] * _dt + r1 * _cov[15][9];
    _cov[4][10] = XP24[1][1] + YP34[1][1] - _cov[13][10] * _dt + r1 * _cov[15][10];
    _cov[4][11] = XP24[1][2] + YP34[1][2] - _cov[13][11] * _dt + r1 * _cov[15][11];
    _cov[5][9] = XP24[2][0] + YP34[2][0] - _cov[14][9] * _dt + r2 * _cov[15][9];
    _cov[5][10] = XP24[2][1] + YP34[2][1] - _cov[14][10] * _dt + r2 * _cov[15][10];
    _cov[5][11] = XP24[2][2] + YP34[2][2] - _cov[14][11] * _dt + r2 * _cov[15][11];

    const array<array<float, 3>, 3> c23 {
        {{XP23[0][0] + YP33[0][0] - _cov[12][6] * _dt + r0 * _cov[15][6],
          XP23[0][1] + YP33[0][1] - _cov[12][7] * _dt + r0 * _cov[15][7],
          XP23[0][2] + YP33[0][2] - _cov[12][8] * _dt + r0 * _cov[15][8]},
         {XP23[1][0] + YP33[1][0] - _cov[13][6] * _dt + r1 * _cov[15][6],
          XP23[1][1] + YP33[1][1] - _cov[13][7] * _dt + r1 * _cov[15][7],
          XP23[1][2] + YP33[1][2] - _cov[13][8] * _dt + r1 * _cov[15][8]},
         {XP23[2][0] + YP33[2][0] - _cov[14][6] * _dt + r2 * _cov[15][6],
          XP23[2][1] + YP33[2][1] - _cov[14][7] * _dt + r2 * _cov[15][7],
          XP23[2][2] + YP33[2][2] - _cov[14][8] * _dt + r2 * _cov[15][8]}}
    };

    _cov[3][6] = c23[0][0] - c23[0][2]*dang1 + c23[0][1]*dang2;
    _cov[3][7] = c23[0][1] + c23[0][2]*dang0 - c23[0][0]*dang2;
    _cov[3][8] = c23[0][2] - c23[0][1]*dang0 + c23[0][0]*dang1;
    _cov[4][6] = c23[1][0] - c23[1][2]*dang1 + c23[1][1]*dang2;
    _cov[4][7] = c23[1][1] + c23[1][2]*dang0 - c23[1][0]*dang2;
    _cov[4][8] = c23[1][2] - c23[1][1]*dang0 + c23[1][0]*dang1;
    _cov[5][6] = c23[2][0] - c23[2][2]*dang1 + c23[2][1]*dang2;
    _cov[5][7] = c23[2][1] + c23[2][2]*dang0 - c23[2][0]*dang2;
    _cov[5][8] = c23[2][2] - c23[2][1]*dang0 + c23[2][0]*dang1;

    _cov[3][3] = XP22XT00 + ((XP23YT[0][0] + XP23YT[0][0]) - (XP25[0][0] + XP25[0][0]) * _dt + (XP26[0] * r0 + XP26[0] * r0))
                 + (YP33YT00 + _cov[12][12] * _dt2 - (YP35[0][0] + YP35[0][0]) * _dt + (_cov[6][15] * r0 + _cov[6][15] * r0) + r0 * _cov[15][15] * r0); 
    _cov[3][4] = XP22XT01 + ((XP23YT[0][1] + XP23YT[1][0]) - (XP25[0][1] + XP25[1][0]) * _dt + (XP26[0] * r1 + XP26[1] * r0))
                 + (YP33YT01 + _cov[12][13] * _dt2 - (YP35[0][1] + YP35[1][0]) * _dt + (_cov[7][15] * r0 + _cov[6][15] * r1) + r0 * _cov[15][15] * r1); 
    _cov[3][5] = XP22XT02 + ((XP23YT[0][2] + XP23YT[2][0]) - (XP25[0][2] + XP25[2][0]) * _dt + (XP26[0] * r2 + XP26[2] * r0))
                 + (YP33YT02 + _cov[12][14] * _dt2 - (YP35[0][2] + YP35[2][0]) * _dt + (_cov[8][15] * r0 + _cov[6][15] * r2) + r0 * _cov[15][15] * r2); 
    _cov[4][4] = XP22XT11 + ((XP23YT[1][1] + XP23YT[1][1]) - (XP25[1][1] + XP25[1][1]) * _dt + (XP26[1] * r1 + XP26[1] * r1))
                 + (YP33YT11 + _cov[13][13] * _dt2 - (YP35[1][1] + YP35[1][1]) * _dt + (_cov[7][15] * r1 + _cov[7][15] * r1) + r1 * _cov[15][15] * r1); 
    _cov[4][5] = XP22XT12 + ((XP23YT[1][2] + XP23YT[2][1]) - (XP25[1][2] + XP25[2][1]) * _dt + (XP26[1] * r2 + XP26[2] * r1))
                 + (YP33YT12 + _cov[13][14] * _dt2 - (YP35[1][2] + YP35[2][1]) * _dt + (_cov[8][15] * r1 + _cov[7][15] * r2) + r1 * _cov[15][15] * r2); 
    _cov[5][5] = XP22XT22 + ((XP23YT[2][2] + XP23YT[2][2]) - (XP25[2][2] + XP25[2][2]) * _dt + (XP26[2] * r2 + XP26[2] * r2))
                 + (YP33YT22 + _cov[14][14] * _dt2 - (YP35[2][2] + YP35[2][2]) * _dt + (_cov[8][15] * r2 + _cov[8][15] * r2) + r2 * _cov[15][15] * r2); 

    _cov[6][15] = _cov[6][15] - _cov[8][15]*dang1 + _cov[7][15]*dang2 - _cov[9][15] * _dt;
    _cov[7][15] = _cov[7][15] + _cov[8][15]*dang0 - _cov[6][15]*dang2 - _cov[10][15] * _dt;
    _cov[8][15] = _cov[8][15] - _cov[7][15]*dang0 + _cov[6][15]*dang1 - _cov[11][15] * _dt;

    _cov[6][12] = _cov[6][12] - _cov[8][12]*dang1 + _cov[7][12]*dang2 - _cov[9][12] * _dt;
    _cov[6][13] = _cov[6][13] - _cov[8][13]*dang1 + _cov[7][13]*dang2 - _cov[9][13] * _dt;
    _cov[6][14] = _cov[6][14] - _cov[8][14]*dang1 + _cov[7][14]*dang2 - _cov[9][14] * _dt;
    _cov[7][12] = _cov[7][12] + _cov[8][12]*dang0 - _cov[6][12]*dang2 - _cov[10][12] * _dt;
    _cov[7][13] = _cov[7][13] + _cov[8][13]*dang0 - _cov[6][13]*dang2 - _cov[10][13] * _dt;
    _cov[7][14] = _cov[7][14] + _cov[8][14]*dang0 - _cov[6][14]*dang2 - _cov[10][14] * _dt;
    _cov[8][12] = _cov[8][12] - _cov[7][12]*dang0 + _cov[6][12]*dang1 - _cov[11][12] * _dt;
    _cov[8][13] = _cov[8][13] - _cov[7][13]*dang0 + _cov[6][13]*dang1 - _cov[11][13] * _dt;
    _cov[8][14] = _cov[8][14] - _cov[7][14]*dang0 + _cov[6][14]*dang1 - _cov[11][14] * _dt;

    const array<array<float, 3>, 3> XP34 {
        {{_cov[6][9] - _cov[8][9]*dang1 + _cov[7][9]*dang2,
          _cov[6][10] - _cov[8][10]*dang1 + _cov[7][10]*dang2,
          _cov[6][11] - _cov[8][11]*dang1 + _cov[7][11]*dang2},
         {_cov[7][9] + _cov[8][9]*dang0 - _cov[6][9]*dang2,
          _cov[7][10] + _cov[8][10]*dang0 - _cov[6][10]*dang2,
          _cov[7][11] + _cov[8][11]*dang0 - _cov[6][11]*dang2},
         {_cov[8][9] - _cov[7][9]*dang0 + _cov[6][9]*dang1,
          _cov[8][10] - _cov[7][10]*dang0 + _cov[6][10]*dang1,
          _cov[8][11] - _cov[7][11]*dang0 + _cov[6][11]*dang1}}
    };

    const array<array<float, 3>, 3> XP33 {
        {{_cov[6][6] - _cov[8][6]*dang1 + _cov[7][6]*dang2,
          _cov[6][7] - _cov[8][7]*dang1 + _cov[7][7]*dang2,
          _cov[6][8] - _cov[8][8]*dang1 + _cov[7][8]*dang2},
         {_cov[7][6] + _cov[8][6]*dang0 - _cov[6][6]*dang2,
          _cov[7][7] + _cov[8][7]*dang0 - _cov[6][7]*dang2,
          _cov[7][8] + _cov[8][8]*dang0 - _cov[6][8]*dang2},
         {_cov[8][6] - _cov[7][6]*dang0 + _cov[6][6]*dang1,
          _cov[8][7] - _cov[7][7]*dang0 + _cov[6][7]*dang1,
          _cov[8][8] - _cov[7][8]*dang0 + _cov[6][8]*dang1}}
    };

    _cov[6][9] = XP34[0][0] - _cov[9][9] * _dt;
    _cov[6][10] = XP34[0][1] - _cov[9][10] * _dt;
    _cov[6][11] = XP34[0][2] - _cov[9][11] * _dt;
    _cov[7][9] = XP34[1][0] - _cov[10][9] * _dt;
    _cov[7][10] = XP34[1][1] - _cov[10][10] * _dt;
    _cov[7][11] = XP34[1][2] - _cov[10][11] * _dt;
    _cov[8][9] = XP34[2][0] - _cov[11][9] * _dt;
    _cov[8][10] = XP34[2][1] - _cov[11][10] * _dt;
    _cov[8][11] = XP34[2][2] - _cov[11][11] * _dt;

    const float XP33XT00 = XP33[0][0] - XP33[0][2]*dang1 + XP33[0][1]*dang2;
    const float XP33XT01 = XP33[0][1] + XP33[0][2]*dang0 - XP33[0][0]*dang2;
    const float XP33XT02 = XP33[0][2] - XP33[0][1]*dang0 + XP33[0][0]*dang1;
    const float XP33XT11 = XP33[1][1] + XP33[1][2]*dang0 - XP33[1][0]*dang2;
    const float XP33XT12 = XP33[1][2] - XP33[1][1]*dang0 + XP33[1][0]*dang1;
    const float XP33XT22 = XP33[2][2] - XP33[2][1]*dang0 + XP33[2][0]*dang1;

    _cov[6][6] = XP33XT00 - (XP34[0][0] + XP34[0][0]) * _dt + _cov[9][9] * _dt2;
    _cov[6][7] = XP33XT01 - (XP34[0][1] + XP34[1][0]) * _dt + _cov[9][10] * _dt2;
    _cov[6][8] = XP33XT02 - (XP34[0][2] + XP34[2][0]) * _dt + _cov[9][11] * _dt2;
    _cov[7][7] = XP33XT11 - (XP34[1][1] + XP34[1][1]) * _dt + _cov[10][10] * _dt2;
    _cov[7][8] = XP33XT12 - (XP34[1][2] + XP34[2][1]) * _dt + _cov[10][11] * _dt2;
    _cov[8][8] = XP33XT22 - (XP34[2][2] + XP34[2][2]) * _dt + _cov[11][11] * _dt2;

    // // Equations for covariance matrix prediction, without process noise!
    // const float var0 = _cov[1][0] + _cov[1][1]*dang2 - _cov[2][1]*dang1 + _cov[3][1]*_dt;
    // const float var1 = _cov[3][0] + _cov[3][1]*dang2 - _cov[3][2]*dang1 + _cov[3][3]*_dt;
    // const float var2 = _cov[2][0] + _cov[2][1]*dang2 - _cov[2][2]*dang1 + _cov[3][2]*_dt;
    // const float var3 = _cov[1][0]*dang2;
    // const float var4 = _cov[2][0]*dang1;
    // const float var5 = _cov[0][0] + _cov[3][0]*_dt + var3 - var4;
    // const float var6 = _cov[4][3]*_dt;
    // const float var7 = _cov[4][0] + _cov[4][1]*dang2 - _cov[4][2]*dang1 + var6;
    // const float var8 = _cov[5][3]*_dt;
    // const float var9 = _cov[5][0] + _cov[5][1]*dang2 - _cov[5][2]*dang1 + var8;
    // const float var10 = _cov[7][0] + _cov[7][1]*dang2 - _cov[7][2]*dang1 + _cov[7][3]*_dt;
    // const float var11 = _cov[15][0] + _cov[15][1]*dang2 - _cov[15][2]*dang1 + _cov[15][3]*_dt;
    // const float var12 = _cov[12][3]*_dt;
    // const float var13 = _cov[12][0] + _cov[12][1]*dang2 - _cov[12][2]*dang1 + var12;
    // const float var14 = _cov[8][0] + _cov[8][1]*dang2 - _cov[8][2]*dang1 + _cov[8][3]*_dt;
    // const float var15 = _cov[13][3]*_dt;
    // const float var16 = _cov[13][0] + _cov[13][1]*dang2 - _cov[13][2]*dang1 + var15;
    // const float var17 = _cov[6][0] + _cov[6][1]*dang2 - _cov[6][2]*dang1 + _cov[6][3]*_dt;
    // const float var18 = _cov[14][3]*_dt;
    // const float var19 = _cov[14][0] + _cov[14][1]*dang2 - _cov[14][2]*dang1 + var18;
    // const float var20 = _cov[9][0] + _cov[9][1]*dang2 - _cov[9][2]*dang1 + _cov[9][3]*_dt;
    // const float var21 = _cov[10][0] + _cov[10][1]*dang2 - _cov[10][2]*dang1 + _cov[10][3]*_dt;
    // const float var22 = _cov[11][0] + _cov[11][1]*dang2 - _cov[11][2]*dang1 + _cov[11][3]*_dt;
    // const float var23 = -_cov[2][0]*dang2 + _cov[2][1] + _cov[2][2]*dang0 + _cov[4][2]*_dt;
    // const float var24 = -_cov[4][0]*dang2 + _cov[4][1] + _cov[4][2]*dang0 + _cov[4][4]*_dt;
    // const float var25 = -_cov[0][0]*dang2 + _cov[1][0] + _cov[2][0]*dang0 + _cov[4][0]*_dt;
    // const float var26 = _cov[2][1]*dang0;
    // const float var27 = _cov[1][1] + _cov[4][1]*_dt + var26 - var3;
    // const float var28 = _cov[5][4]*_dt;
    // const float var29 = -_cov[5][0]*dang2 + _cov[5][1] + _cov[5][2]*dang0 + var28;
    // const float var30 = -_cov[7][0]*dang2 + _cov[7][1] + _cov[7][2]*dang0 + _cov[7][4]*_dt;
    // const float var31 = -_cov[15][0]*dang2 + _cov[15][1] + _cov[15][2]*dang0 + _cov[15][4]*_dt;
    // const float var32 = _cov[12][4]*_dt;
    // const float var33 = -_cov[12][0]*dang2 + _cov[12][1] + _cov[12][2]*dang0 + var32;
    // const float var34 = -_cov[8][0]*dang2 + _cov[8][1] + _cov[8][2]*dang0 + _cov[8][4]*_dt;
    // const float var35 = -_cov[3][0]*dang2 + _cov[3][1] + _cov[3][2]*dang0 + var6;
    // const float var36 = _cov[13][4]*_dt;
    // const float var37 = -_cov[13][0]*dang2 + _cov[13][1] + _cov[13][2]*dang0 + var36;
    // const float var38 = -_cov[6][0]*dang2 + _cov[6][1] + _cov[6][2]*dang0 + _cov[6][4]*_dt;
    // const float var39 = _cov[14][4]*_dt;
    // const float var40 = -_cov[14][0]*dang2 + _cov[14][1] + _cov[14][2]*dang0 + var39;
    // const float var41 = -_cov[9][0]*dang2 + _cov[9][1] + _cov[9][2]*dang0 + _cov[9][4]*_dt;
    // const float var42 = -_cov[10][0]*dang2 + _cov[10][1] + _cov[10][2]*dang0 + _cov[10][4]*_dt;
    // const float var43 = -_cov[11][0]*dang2 + _cov[11][1] + _cov[11][2]*dang0 + _cov[11][4]*_dt;
    // const float var44 = _cov[5][0]*dang1 - _cov[5][1]*dang0 + _cov[5][2] + _cov[5][5]*_dt;
    // const float var45 = _cov[4][0]*dang1 - _cov[4][1]*dang0 + _cov[4][2] + var28;
    // const float var46 = _cov[7][0]*dang1 - _cov[7][1]*dang0 + _cov[7][2] + _cov[7][5]*_dt;
    // const float var47 = _cov[15][0]*dang1 - _cov[15][1]*dang0 + _cov[15][2] + _cov[15][5]*_dt;
    // const float var48 = _cov[12][5]*_dt;
    // const float var49 = _cov[12][0]*dang1 - _cov[12][1]*dang0 + _cov[12][2] + var48;
    // const float var50 = _cov[8][0]*dang1 - _cov[8][1]*dang0 + _cov[8][2] + _cov[8][5]*_dt;
    // const float var51 = _cov[3][0]*dang1 - _cov[3][1]*dang0 + _cov[3][2] + var8;
    // const float var52 = _cov[13][5]*_dt;
    // const float var53 = _cov[13][0]*dang1 - _cov[13][1]*dang0 + _cov[13][2] + var52;
    // const float var54 = _cov[6][0]*dang1 - _cov[6][1]*dang0 + _cov[6][2] + _cov[6][5]*_dt;
    // const float var55 = _cov[14][5]*_dt;
    // const float var56 = _cov[14][0]*dang1 - _cov[14][1]*dang0 + _cov[14][2] + var55;
    // const float var57 = _cov[9][0]*dang1 - _cov[9][1]*dang0 + _cov[9][2] + _cov[9][5]*_dt;
    // const float var58 = _cov[10][0]*dang1 - _cov[10][1]*dang0 + _cov[10][2] + _cov[10][5]*_dt;
    // const float var59 = _cov[11][0]*dang1 - _cov[11][1]*dang0 + _cov[11][2] + _cov[11][5]*_dt;
    // const float var60 = _cov[4][3] + _cov[15][4]*r0 + _cov[4][4]*dang2 - _cov[5][4]*dang1 + _cov[7][4]*dv2 - _cov[8][4]*dv1 - var32;
    // const float var61 = _cov[7][3] + _cov[7][4]*dang2 - _cov[7][5]*dang1 - _cov[12][7]*_dt + _cov[15][7]*r0 + _cov[7][7]*dv2 - _cov[8][7]*dv1;
    // const float var62 = -_cov[15][12]*_dt + _cov[15][15]*r0 + _cov[15][3] + _cov[15][4]*dang2 - _cov[15][5]*dang1 + _cov[15][7]*dv2 - _cov[15][8]*dv1;
    // const float var63 = _cov[5][3] + _cov[5][4]*dang2 + _cov[15][5]*r0 - _cov[5][5]*dang1 + _cov[7][5]*dv2 - _cov[8][5]*dv1 - var48;
    // const float var64 = -_cov[12][12]*_dt + _cov[15][12]*r0 + _cov[12][3] + _cov[12][4]*dang2 - _cov[12][5]*dang1 + _cov[12][7]*dv2 - _cov[12][8]*dv1;
    // const float var65 = _cov[8][3] + _cov[8][4]*dang2 - _cov[8][5]*dang1 + _cov[8][7]*dv2 - _cov[12][8]*_dt + _cov[15][8]*r0 - _cov[8][8]*dv1;
    // const float var66 = _cov[4][3]*dang2;
    // const float var67 = _cov[5][3]*dang1;
    // const float var68 = _cov[15][3]*r0 + _cov[3][3] + _cov[7][3]*dv2 - _cov[8][3]*dv1 - var12 + var66 - var67;
    // const float var69 = -_cov[13][12]*_dt;
    // const float var70 = _cov[15][13]*r0 + _cov[13][3] + _cov[13][4]*dang2 - _cov[13][5]*dang1 + _cov[13][7]*dv2 - _cov[13][8]*dv1 + var69;
    // const float var71 = _cov[7][6]*dv2;
    // const float var72 = _cov[8][6]*dv1;
    // const float var73 = _cov[6][3] + _cov[6][4]*dang2 - _cov[6][5]*dang1 - _cov[12][6]*_dt + _cov[15][6]*r0 + var71 - var72;
    // const float var74 = -_cov[14][12]*_dt;
    // const float var75 = _cov[15][14]*r0 + _cov[14][3] + _cov[14][4]*dang2 - _cov[14][5]*dang1 + _cov[14][7]*dv2 - _cov[14][8]*dv1 + var74;
    // const float var76 = -_cov[12][9]*_dt;
    // const float var77 = _cov[9][3] + _cov[9][4]*dang2 - _cov[9][5]*dang1 + _cov[9][7]*dv2 - _cov[9][8]*dv1 + _cov[15][9]*r0 + var76;
    // const float var78 = -_cov[12][10]*_dt;
    // const float var79 = _cov[15][10]*r0 + _cov[10][3] + _cov[10][4]*dang2 - _cov[10][5]*dang1 + _cov[10][7]*dv2 - _cov[10][8]*dv1 + var78;
    // const float var80 = -_cov[12][11]*_dt;
    // const float var81 = _cov[15][11]*r0 + _cov[11][3] + _cov[11][4]*dang2 - _cov[11][5]*dang1 + _cov[11][7]*dv2 - _cov[11][8]*dv1 + var80;
    // const float var82 = -_cov[5][3]*dang2 + _cov[5][4] + _cov[15][5]*r1 + _cov[5][5]*dang0 - _cov[6][5]*dv2 + _cov[8][5]*dv0 - var52;
    // const float var83 = -_cov[8][3]*dang2 + _cov[8][4] + _cov[8][5]*dang0 - _cov[8][6]*dv2 - _cov[13][8]*_dt + _cov[15][8]*r1 + _cov[8][8]*dv0;
    // const float var84 = -_cov[15][13]*_dt + _cov[15][15]*r1 - _cov[15][3]*dang2 + _cov[15][4] + _cov[15][5]*dang0 - _cov[15][6]*dv2 + _cov[15][8]*dv0;
    // const float var85 = _cov[15][3]*r1 - _cov[3][3]*dang2 + _cov[4][3] + _cov[5][3]*dang0 - _cov[6][3]*dv2 + _cov[8][3]*dv0 - var15;
    // const float var86 = -_cov[13][13]*_dt + _cov[15][13]*r1 - _cov[13][3]*dang2 + _cov[13][4] + _cov[13][5]*dang0 - _cov[13][6]*dv2 + _cov[13][8]*dv0;
    // const float var87 = -_cov[6][3]*dang2 + _cov[6][4] + _cov[6][5]*dang0 - _cov[13][6]*_dt + _cov[15][6]*r1 - _cov[6][6]*dv2 + _cov[8][6]*dv0;
    // const float var88 = _cov[5][4]*dang0;
    // const float var89 = _cov[15][4]*r1 + _cov[4][4] - _cov[6][4]*dv2 + _cov[8][4]*dv0 - var36 - var66 + var88;
    // const float var90 = -_cov[14][13]*_dt;
    // const float var91 = _cov[15][14]*r1 - _cov[14][3]*dang2 + _cov[14][4] + _cov[14][5]*dang0 - _cov[14][6]*dv2 + _cov[14][8]*dv0 + var90;
    // const float var92 = _cov[8][7]*dv0;
    // const float var93 = -_cov[7][3]*dang2 + _cov[7][4] + _cov[7][5]*dang0 - _cov[13][7]*_dt + _cov[15][7]*r1 - var71 + var92;
    // const float var94 = -_cov[13][9]*_dt;
    // const float var95 = -_cov[9][3]*dang2 + _cov[9][4] + _cov[9][5]*dang0 - _cov[9][6]*dv2 + _cov[9][8]*dv0 + _cov[15][9]*r1 + var94;
    // const float var96 = -_cov[13][10]*_dt;
    // const float var97 = _cov[15][10]*r1 - _cov[10][3]*dang2 + _cov[10][4] + _cov[10][5]*dang0 - _cov[10][6]*dv2 + _cov[10][8]*dv0 + var96;
    // const float var98 = -_cov[13][11]*_dt;
    // const float var99 = _cov[15][11]*r1 - _cov[11][3]*dang2 + _cov[11][4] + _cov[11][5]*dang0 - _cov[11][6]*dv2 + _cov[11][8]*dv0 + var98;
    // const float var100 = _cov[6][3]*dang1 - _cov[6][4]*dang0 + _cov[6][5] - _cov[14][6]*_dt + _cov[15][6]*r2 + _cov[6][6]*dv1 - _cov[7][6]*dv0;
    // const float var101 = -_cov[15][14]*_dt + _cov[15][15]*r2 + _cov[15][3]*dang1 - _cov[15][4]*dang0 + _cov[15][5] + _cov[15][6]*dv1 - _cov[15][7]*dv0;
    // const float var102 = -_cov[14][14]*_dt + _cov[15][14]*r2 + _cov[14][3]*dang1 - _cov[14][4]*dang0 + _cov[14][5] + _cov[14][6]*dv1 - _cov[14][7]*dv0;
    // const float var103 = _cov[7][3]*dang1 - _cov[7][4]*dang0 + _cov[7][5] + _cov[7][6]*dv1 - _cov[14][7]*_dt + _cov[15][7]*r2 - _cov[7][7]*dv0;
    // const float var104 = _cov[8][3]*dang1 - _cov[8][4]*dang0 + _cov[8][5] - _cov[14][8]*_dt + _cov[15][8]*r2 + var72 - var92;
    // const float var105 = -_cov[14][9]*_dt;
    // const float var106 = _cov[9][3]*dang1 - _cov[9][4]*dang0 + _cov[9][5] + _cov[9][6]*dv1 - _cov[9][7]*dv0 + _cov[15][9]*r2 + var105;
    // const float var107 = -_cov[14][10]*_dt;
    // const float var108 = _cov[15][10]*r2 + _cov[10][3]*dang1 - _cov[10][4]*dang0 + _cov[10][5] + _cov[10][6]*dv1 - _cov[10][7]*dv0 + var107;
    // const float var109 = -_cov[14][11]*_dt;
    // const float var110 = _cov[15][11]*r2 + _cov[11][3]*dang1 - _cov[11][4]*dang0 + _cov[11][5] + _cov[11][6]*dv1 - _cov[11][7]*dv0 + var109;
    // const float var111 = _cov[7][6] + _cov[7][7]*dang2 - _cov[8][7]*dang1 - _cov[9][7]*_dt;
    // const float var112 = _cov[8][6] + _cov[8][7]*dang2 - _cov[8][8]*dang1 - _cov[9][8]*_dt;
    // const float var113 = _cov[9][6] + _cov[9][7]*dang2 - _cov[9][8]*dang1 - _cov[9][9]*_dt;
    // const float var114 = _cov[7][6]*dang2;
    // const float var115 = _cov[8][6]*dang1;
    // const float var116 = _cov[6][6] - _cov[9][6]*_dt + var114 - var115;
    // const float var117 = -_cov[10][9]*_dt;
    // const float var118 = _cov[10][6] + _cov[10][7]*dang2 - _cov[10][8]*dang1 + var117;
    // const float var119 = -_cov[11][9]*_dt;
    // const float var120 = _cov[11][6] + _cov[11][7]*dang2 - _cov[11][8]*dang1 + var119;
    // const float var121 = _cov[8][7]*dang0;
    // const float var122 = _cov[10][7]*_dt;
    // const float var123 = _cov[10][6]*_dt + _cov[6][6]*dang2 - _cov[7][6] - _cov[8][6]*dang0;
    // const float var124 = _cov[10][10]*_dt;
    // const float var125 = _cov[10][6]*dang2;
    // const float var126 = _cov[10][8]*dang0;
    // const float var127 = _cov[8][6]*dang2;
    // const float var128 = _cov[10][8]*_dt;
    // const float var129 = _cov[8][8]*dang0;
    // const float var130 = -var121;
    // const float var131 = _cov[11][10]*_dt;
    // const float var132 = _cov[11][6]*dang2;
    // const float var133 = _cov[11][8]*dang0;
    // const float var134 = -var131;
    // const float var135 = _cov[11][11]*_dt;
    // const float var136 = _cov[11][7]*dang0;
    // const float var137 = _cov[11][6]*dang1;

    // _cov[2][2] = _cov[2][2] + _cov[5][2]*_dt - dang0*(_cov[1][0]*dang1 - _cov[1][1]*dang0 + _cov[2][1] + _cov[5][1]*_dt) + dang1*(_cov[0][0]*dang1 - _cov[1][0]*dang0 + _cov[2][0] + _cov[5][0]*_dt) + _dt*var44 - var26 + var4;
    // _cov[0][0] = -dang1*var2 + dang2*var0 + _dt*var1 + var5;
    // _cov[0][1] = dang0*var2 - dang2*var5 + _dt*var7 + var0;
    // _cov[1][1] = dang0*var23 - dang2*var25 + _dt*var24 + var27;
    // _cov[0][2] = -dang0*var0 + dang1*var5 + _dt*var9 + var2;
    // _cov[1][2] = -dang0*var27 + dang1*var25 + _dt*var29 + var23;

    // _cov[5][5] = _cov[15][5]*r2 + _cov[5][5] + _cov[6][5]*dv1 - _cov[7][5]*dv0 - dang0*(_cov[4][3]*dang1 + _cov[15][4]*r2 - _cov[4][4]*dang0 + _cov[5][4] + _cov[6][4]*dv1 - _cov[7][4]*dv0 - var39) + dang1*(_cov[15][3]*r2 + _cov[3][3]*dang1 - _cov[4][3]*dang0 + _cov[5][3] + _cov[6][3]*dv1 - _cov[7][3]*dv0 - var18) - _dt*var102 - dv0*var103 + dv1*var100 + r2*var101 - var55 + var67 - var88;
    // _cov[0][3] = -dang1*var9 + dang2*var7 - _dt*var13 - dv1*var14 + dv2*var10 + r0*var11 + var1;
    // _cov[1][3] = -dang1*var29 + dang2*var24 - _dt*var33 - dv1*var34 + dv2*var30 + r0*var31 + var35;
    // _cov[2][3] = -dang1*var44 + dang2*var45 - _dt*var49 - dv1*var50 + dv2*var46 + r0*var47 + var51;
    // _cov[3][3] = -dang1*var63 + dang2*var60 - _dt*var64 - dv1*var65 + dv2*var61 + r0*var62 + var68;
    // _cov[0][4] = dang0*var9 - dang2*var1 - _dt*var16 + dv0*var14 - dv2*var17 + r1*var11 + var7;
    // _cov[1][4] = dang0*var29 - dang2*var35 - _dt*var37 + dv0*var34 - dv2*var38 + r1*var31 + var24;
    // _cov[2][4] = dang0*var44 - dang2*var51 - _dt*var53 + dv0*var50 - dv2*var54 + r1*var47 + var45;
    // _cov[3][4] = dang0*var63 - dang2*var68 - _dt*var70 + dv0*var65 - dv2*var73 + r1*var62 + var60;
    // _cov[4][4] = dang0*var82 - dang2*var85 - _dt*var86 + dv0*var83 - dv2*var87 + r1*var84 + var89;
    // _cov[0][5] = -dang0*var7 + dang1*var1 - _dt*var19 - dv0*var10 + dv1*var17 + r2*var11 + var9;
    // _cov[1][5] = -dang0*var24 + dang1*var35 - _dt*var40 - dv0*var30 + dv1*var38 + r2*var31 + var29;
    // _cov[2][5] = -dang0*var45 + dang1*var51 - _dt*var56 - dv0*var46 + dv1*var54 + r2*var47 + var44;
    // _cov[3][5] = -dang0*var60 + dang1*var68 - _dt*var75 - dv0*var61 + dv1*var73 + r2*var62 + var63;
    // _cov[4][5] = -dang0*var89 + dang1*var85 - _dt*var91 - dv0*var93 + dv1*var87 + r2*var84 + var82;

    // _cov[7][8] = _cov[8][7] + dang0*(-_cov[7][7] + var114 + var122 + var130) - dang1*var123 + _dt*(-_cov[11][7] + var131 + var132 - var133) - var127 - var128 + var129;
    // _cov[8][8] = -_cov[11][8]*_dt + _cov[8][8] - dang0*(_cov[7][6]*dang1 - _cov[11][7]*_dt - _cov[7][7]*dang0 + _cov[8][7]) - dang1*(_cov[11][6]*_dt - _cov[6][6]*dang1 + _cov[7][6]*dang0 - _cov[8][6]) + _dt*(-_cov[11][8] + var135 + var136 - var137) + var115 + var130;
    // _cov[7][7] = _cov[7][7] - dang0*(-_cov[8][7] + var127 + var128 - var129) + dang2*var123 + _dt*(-_cov[10][7] + var124 + var125 - var126) - var114 + var121 - var122;
    // _cov[0][6] = -dang1*var14 + dang2*var10 - _dt*var20 + var17;
    // _cov[1][6] = -dang1*var34 + dang2*var30 - _dt*var41 + var38;
    // _cov[2][6] = -dang1*var50 + dang2*var46 - _dt*var57 + var54;
    // _cov[3][6] = -dang1*var65 + dang2*var61 - _dt*var77 + var73;
    // _cov[4][6] = -dang1*var83 + dang2*var93 - _dt*var95 + var87;
    // _cov[5][6] = -dang1*var104 + dang2*var103 - _dt*var106 + var100;
    // _cov[6][6] = -dang1*var112 + dang2*var111 - _dt*var113 + var116;
    // _cov[0][7] = dang0*var14 - dang2*var17 - _dt*var21 + var10;
    // _cov[1][7] = dang0*var34 - dang2*var38 - _dt*var42 + var30;
    // _cov[2][7] = dang0*var50 - dang2*var54 - _dt*var58 + var46;
    // _cov[3][7] = dang0*var65 - dang2*var73 - _dt*var79 + var61;
    // _cov[4][7] = dang0*var83 - dang2*var87 - _dt*var97 + var93;
    // _cov[5][7] = dang0*var104 - dang2*var100 - _dt*var108 + var103;
    // _cov[6][7] = dang0*var112 - dang2*var116 - _dt*var118 + var111;
    // _cov[0][8] = -dang0*var10 + dang1*var17 - _dt*var22 + var14;
    // _cov[1][8] = -dang0*var30 + dang1*var38 - _dt*var43 + var34;
    // _cov[2][8] = -dang0*var46 + dang1*var54 - _dt*var59 + var50;
    // _cov[3][8] = -dang0*var61 + dang1*var73 - _dt*var81 + var65;
    // _cov[4][8] = -dang0*var93 + dang1*var87 - _dt*var99 + var83;
    // _cov[5][8] = -dang0*var103 + dang1*var100 - _dt*var110 + var104;
    // _cov[6][8] = -dang0*var111 + dang1*var116 - _dt*var120 + var112;

    // _cov[0][9] = var20;
    // _cov[1][9] = var41;
    // _cov[2][9] = var57;
    // _cov[3][9] = var77;
    // _cov[4][9] = var95;
    // _cov[5][9] = var106;
    // _cov[6][9] = var113;
    // _cov[7][9] = -_cov[9][6]*dang2 + _cov[9][7] + _cov[9][8]*dang0 + var117;
    // _cov[8][9] = _cov[9][6]*dang1 - _cov[9][7]*dang0 + _cov[9][8] + var119;
    // _cov[0][10] = var21;
    // _cov[1][10] = var42;
    // _cov[2][10] = var58;
    // _cov[3][10] = var79;
    // _cov[4][10] = var97;
    // _cov[5][10] = var108;
    // _cov[6][10] = var118;
    // _cov[7][10] = _cov[10][7] - var124 - var125 + var126;
    // _cov[8][10] = _cov[10][6]*dang1 - _cov[10][7]*dang0 + _cov[10][8] + var134;
    // _cov[0][11] = var22;
    // _cov[1][11] = var43;
    // _cov[2][11] = var59;
    // _cov[3][11] = var81;
    // _cov[4][11] = var99;
    // _cov[5][11] = var110;
    // _cov[6][11] = var120;
    // _cov[7][11] = _cov[11][7] - var132 + var133 + var134;
    // _cov[8][11] = _cov[11][8] - var135 - var136 + var137;

    // _cov[0][12] = var13;
    // _cov[1][12] = var33;
    // _cov[2][12] = var49;
    // _cov[3][12] = var64;
    // _cov[4][12] = _cov[15][12]*r1 - _cov[12][3]*dang2 + _cov[12][4] + _cov[12][5]*dang0 - _cov[12][6]*dv2 + _cov[12][8]*dv0 + var69;
    // _cov[5][12] = _cov[15][12]*r2 + _cov[12][3]*dang1 - _cov[12][4]*dang0 + _cov[12][5] + _cov[12][6]*dv1 - _cov[12][7]*dv0 + var74;
    // _cov[6][12] = _cov[12][6] + _cov[12][7]*dang2 - _cov[12][8]*dang1 + var76;
    // _cov[7][12] = -_cov[12][6]*dang2 + _cov[12][7] + _cov[12][8]*dang0 + var78;
    // _cov[8][12] = _cov[12][6]*dang1 - _cov[12][7]*dang0 + _cov[12][8] + var80;
    // _cov[0][13] = var16;
    // _cov[1][13] = var37;
    // _cov[2][13] = var53;
    // _cov[3][13] = var70;
    // _cov[4][13] = var86;
    // _cov[5][13] = _cov[15][13]*r2 + _cov[13][3]*dang1 - _cov[13][4]*dang0 + _cov[13][5] + _cov[13][6]*dv1 - _cov[13][7]*dv0 + var90;
    // _cov[6][13] = _cov[13][6] + _cov[13][7]*dang2 - _cov[13][8]*dang1 + var94;
    // _cov[7][13] = -_cov[13][6]*dang2 + _cov[13][7] + _cov[13][8]*dang0 + var96;
    // _cov[8][13] = _cov[13][6]*dang1 - _cov[13][7]*dang0 + _cov[13][8] + var98;
    // _cov[0][14] = var19;
    // _cov[1][14] = var40;
    // _cov[2][14] = var56;
    // _cov[3][14] = var75;
    // _cov[4][14] = var91;
    // _cov[5][14] = var102;
    // _cov[6][14] = _cov[14][6] + _cov[14][7]*dang2 - _cov[14][8]*dang1 + var105;
    // _cov[7][14] = -_cov[14][6]*dang2 + _cov[14][7] + _cov[14][8]*dang0 + var107;
    // _cov[8][14] = _cov[14][6]*dang1 - _cov[14][7]*dang0 + _cov[14][8] + var109;

    // _cov[0][15] = var11;
    // _cov[1][15] = var31;
    // _cov[2][15] = var47;
    // _cov[3][15] = var62;
    // _cov[4][15] = var84;
    // _cov[5][15] = var101;
    // _cov[6][15] = _cov[15][6] + _cov[15][7]*dang2 - _cov[15][8]*dang1 - _cov[15][9]*_dt;
    // _cov[7][15] = -_cov[15][10]*_dt - _cov[15][6]*dang2 + _cov[15][7] + _cov[15][8]*dang0;
    // _cov[8][15] = -_cov[15][11]*_dt + _cov[15][6]*dang1 - _cov[15][7]*dang0 + _cov[15][8];    

    // F *P *F' + Q
    _cov[0][0] = kahan_summation(_cov[0][0], _q_cov[0] * _dt2, _accumulator_cov[0]);
    _cov[1][1] = kahan_summation(_cov[1][1], _q_cov[1] * _dt2, _accumulator_cov[1]);
    _cov[2][2] = kahan_summation(_cov[2][2], _q_cov[2] * _dt2, _accumulator_cov[2]);
    _cov[3][3] = kahan_summation(_cov[3][3], _q_cov[3] * _dt2, _accumulator_cov[3]);
    _cov[4][4] = kahan_summation(_cov[4][4], _q_cov[4] * _dt2, _accumulator_cov[4]);
    _cov[5][5] = kahan_summation(_cov[5][5], _q_cov[5] * _dt2, _accumulator_cov[5]);
    _cov[6][6] = kahan_summation(_cov[6][6], _q_cov[6] * _dt2, _accumulator_cov[6]);
    _cov[7][7] = kahan_summation(_cov[7][7], _q_cov[7] * _dt2, _accumulator_cov[7]);
    _cov[8][8] = kahan_summation(_cov[8][8], _q_cov[8] * _dt2, _accumulator_cov[8]);
    _cov[9][9] = kahan_summation(_cov[9][9], _q_cov[9] * _dt2, _accumulator_cov[9]);
    _cov[10][10] = kahan_summation(_cov[10][10], _q_cov[10] * _dt2, _accumulator_cov[10]);
    _cov[11][11] = kahan_summation(_cov[11][11], _q_cov[11] * _dt2, _accumulator_cov[11]);
    _cov[12][12] = kahan_summation(_cov[12][12], _q_cov[12] * _dt2, _accumulator_cov[12]);
    _cov[13][13] = kahan_summation(_cov[13][13], _q_cov[13] * _dt2, _accumulator_cov[13]);
    _cov[14][14] = kahan_summation(_cov[14][14], _q_cov[14] * _dt2, _accumulator_cov[14]);
    _cov[15][15] = kahan_summation(_cov[15][15], _q_cov[15] * _dt2, _accumulator_cov[15]);

    regular_covariance_to_symmetric();
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
        const array<float, 16> HP = {_cov[0][dim] + _cov[0][6]*minus_d_hat[dim][0] + _cov[0][7]*minus_d_hat[dim][1] + _cov[0][8]*minus_d_hat[dim][2],
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
                                     _cov[dim][15] + _cov[6][15]*minus_d_hat[dim][0] + _cov[7][15]*minus_d_hat[dim][1] + _cov[8][15]*minus_d_hat[dim][2]};

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
    regular_covariance_to_symmetric();

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
        const array<float, 16> HP = {_cov[0][index] + _cov[0][6]*d_cross_w_corr_hat[dim][0] + _cov[0][7]*d_cross_w_corr_hat[dim][1] + _cov[0][8]*d_cross_w_corr_hat[dim][2] + _cov[0][9]*d_hat[dim][0] + _cov[0][10]*d_hat[dim][1] + _cov[0][11]*d_hat[dim][2],
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
                                     _cov[index][15] + _cov[6][15]*d_cross_w_corr_hat[dim][0] + _cov[7][15]*d_cross_w_corr_hat[dim][1] + _cov[8][15]*d_cross_w_corr_hat[dim][2] + _cov[9][15]*d_hat[dim][0] + _cov[10][15]*d_hat[dim][1] + _cov[11][15]*d_hat[dim][2]};                                                  
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
    regular_covariance_to_symmetric();

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
    */
    for (unsigned char i = 0; i < 3; ++i) {
        _p[i] += _rot(i, 0) * _error_state[0] + _rot(i, 1) * _error_state[1] + _rot(i, 2) * _error_state[2];
        _v[i] += _rot(i, 0) * _error_state[3] + _rot(i, 1) * _error_state[4] + _rot(i, 2) * _error_state[5];
        _bg[i] += _error_state[9 + i];
        _ba[i] += _error_state[12 + i];
    }
    _g += _error_state[15];

    // [δp, δv, δθ, δbg, δba, δg] = 0
    for (float &es : _error_state) {
        es = 0.f;
    }
}