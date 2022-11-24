// #include "eskf.h"
// #include <cfloat>
// #include <iostream>

// using namespace std;
// using namespace eskf;

// unsigned char conservative_posteriori_estimate_000000(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 12; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }

//     conservative_posteriori_estimate<12>(0, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_000001(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 13; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }

//     conservative_posteriori_estimate<13>(0, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_000010(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 12; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }
//     if (HP[13] * HP[13] > HPHT_plus_R * _cov[13][13]) {
//         return 2;    
//     }

//     conservative_posteriori_estimate<12>(0, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(13, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_000011(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 14; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }

//     conservative_posteriori_estimate<14>(0, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_000100(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 12; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }
//     if (HP[14] * HP[14] > HPHT_plus_R * _cov[14][14]) {
//         return 2;    
//     }

//     conservative_posteriori_estimate<12>(0, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(14, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_000101(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 12; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }
//     if (HP[12] * HP[12] > HPHT_plus_R * _cov[12][12]) {
//         return 2;    
//     }
//     if (HP[14] * HP[14] > HPHT_plus_R * _cov[14][14]) {
//         return 2;    
//     }

//     conservative_posteriori_estimate<12>(0, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(12, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(14, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_000110(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 12; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }
//     if (HP[13] * HP[13] > HPHT_plus_R * _cov[13][13]) {
//         return 2;    
//     }
//     if (HP[14] * HP[14] > HPHT_plus_R * _cov[14][14]) {
//         return 2;    
//     }

//     conservative_posteriori_estimate<12>(0, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(13, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(14, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_000111(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 15; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }

//     conservative_posteriori_estimate<15>(0, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_001000(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 12; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }
//     if (HP[15] * HP[15] > HPHT_plus_R * _cov[15][15]) {
//         return 2;    
//     }

//     conservative_posteriori_estimate<12>(0, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(15, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_001001(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 13; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }
//     if (HP[15] * HP[15] > HPHT_plus_R * _cov[15][15]) {
//         return 2;    
//     }

//     conservative_posteriori_estimate<13>(0, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(15, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_001010(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 12; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }
//     if (HP[13] * HP[13] > HPHT_plus_R * _cov[13][13]) {
//         return 2;    
//     }
//     if (HP[15] * HP[15] > HPHT_plus_R * _cov[15][15]) {
//         return 2;    
//     }

//     conservative_posteriori_estimate<12>(0, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(13, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(15, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_001011(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 14; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }
//     if (HP[15] * HP[15] > HPHT_plus_R * _cov[15][15]) {
//         return 2;    
//     }

//     conservative_posteriori_estimate<14>(0, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(15, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_001100(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 12; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }
//     if (HP[14] * HP[14] > HPHT_plus_R * _cov[14][14]) {
//         return 2;    
//     }
//     if (HP[15] * HP[15] > HPHT_plus_R * _cov[15][15]) {
//         return 2;    
//     }

//     conservative_posteriori_estimate<12>(0, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(14, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(15, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_001101(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 12; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }
//     if (HP[12] * HP[12] > HPHT_plus_R * _cov[12][12]) {
//         return 2;    
//     }
//     if (HP[14] * HP[14] > HPHT_plus_R * _cov[14][14]) {
//         return 2;    
//     }
//     if (HP[15] * HP[15] > HPHT_plus_R * _cov[15][15]) {
//         return 2;    
//     }

//     conservative_posteriori_estimate<13>(0, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(14, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(15, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_001110(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 12; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }
//     if (HP[13] * HP[13] > HPHT_plus_R * _cov[13][13]) {
//         return 2;    
//     }
//     if (HP[14] * HP[14] > HPHT_plus_R * _cov[14][14]) {
//         return 2;    
//     }
//     if (HP[15] * HP[15] > HPHT_plus_R * _cov[15][15]) {
//         return 2;    
//     }

//     conservative_posteriori_estimate<12>(0, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(13, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(14, HP, HPHT_plus_R, obs_error);
//     conservative_posteriori_estimate<1>(15, HP, HPHT_plus_R, obs_error);
// }

// unsigned char conservative_posteriori_estimate_001111(const array<float, ESKF::dim> &HP, const float &HPHT_plus_R, const float &obs_error, const float &gate) {
//     /*
//     K = P * H' * (H * P * H' + R)^-1
//     P = P - K * H * P

//     e = y - h = y - (v + R * (w - bg)^ * dis)
//     x = x + K * (y - h)
//     */
    
//     // Don't correct error state unless obs_error / √(H*P*H' + R) > gate
//     if (obs_error * obs_error  > gate * gate * HPHT_plus_R) {
//         return 1;
//     }

//     // Don't correct error state unless (I - KH) * P >= 0 <=> P - PH'(HPH'+R)^-1HP >= 0
//     for (unsigned char i = 0; i < 16; ++i) {
//         if (HP[i] * HP[i] > HPHT_plus_R * _cov[i][i]) {
//             return 2;    
//         }
//     }

//     conservative_posteriori_estimate<16>(0, HP, HPHT_plus_R, obs_error);
// }