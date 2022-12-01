#include "gvar.h"

// WGS-84 model
const double GM = 3.986004415e14;
const double Re = 6.378136998405e6;
const double wie = 7.2921151467e-5;

const double ff = 1. / 298.257223563;
const double ee = sqrt(2. * ff - ff * ff);
const double e2 = ee * ee;
const double Rp = (1. - ff) * Re;

// gravity, ug
const double ge = 9.780325333434361;
const double gp = 9.832184935381024;
const double g0 = ge;
const double ug = g0 * 1e-6;

// angle unit
const double arcdeg = M_PI / 180.;
const double arcmin = arcdeg / 60.;
const double arcsec = arcmin / 60.;

// hour, deg/hour, deg / sqrt(hour)
const double hur = 3600.;
const double dph = arcdeg / hur;
const double dpsh = arcdeg / sqrt(hur);

// ug/sqrt(Hz)
const double ugpsHz = ug / sqrt(1.);



