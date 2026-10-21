// Independent calls to pinned native single-frequency slip helpers.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/algorithms/rtk_slip_detection.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
    using namespace libgnss::rtk_slip_detection;
    const double wavelength = 0.19029367279836487;
    std::cout << std::setprecision(17);
    std::cout << "# name,kind,inputs...,native_result; finite valid inputs only\n";
    auto rate = [&](const char* name, double r, double b, double w) {
        std::cout << name << ",rate," << r << ',' << b << ',' << w << ','
                  << singleDifferenceRangeRateMps(r,b,w) << '\n';
    };
    auto cmc = [&](const char* name, double r, double b, double rp, double bp, double w) {
        std::cout << name << ",cmc," << r << ',' << b << ',' << rp << ',' << bp << ',' << w << ','
                  << singleDifferenceCodeMinusPhaseM(r,b,rp,bp,w) << '\n';
    };
    auto code = [&](const char* name, double p, double c, double t) {
        std::cout << name << ",code," << p << ',' << c << ',' << t << ','
                  << detectCodeSlip(p,c,t) << '\n';
    };
    auto doppler = [&](const char* name, double p, double c, double r, double dt, double t) {
        std::cout << name << ",doppler," << p << ',' << c << ',' << r << ',' << dt << ',' << t << ','
                  << detectDopplerSlip(p,c,r,dt,t) << '\n';
    };
    rate("negative_doppler",-1000,-998,wavelength);
    rate("opposite_rate",1000,998,wavelength);
    rate("zero_rate",0,0,wavelength);
    rate("equal_nonzero",-3421.23,-3421.23,wavelength);
    rate("other_wavelength",-502.5,-500.25,0.24421021342456825);
    cmc("ordinary",20200000.40,20199998.10,106151720,106151708,wavelength);
    cmc("zero_difference",0,0,0,0,wavelength);
    cmc("negative_bias",22000000,22000005,110000000,110000040,wavelength);
    cmc("fractional",24000001.125,24000000.375,120000003.25,120000000.5,wavelength);
    cmc("clock_common_mode",24001001.125,24001000.375,120000003.25,120000000.5,wavelength);
    code("small",1.5,3.9,3);
    code("jump",1.5,6.75,3);
    code("negative_jump",6.75,1.5,3);
    code("exact_5m",0,5,5);
    code("above_5m",0,std::nextafter(5.0,INFINITY),5);
    code("below_5m",0,std::nextafter(5.0,0.0),5);
    code("unchanged_nonzero",-1000,-1000,5);
    code("exact_negative",0,-5,5);
    doppler("consistent",10,10.15,0.15,1,0.2);
    doppler("phase_jump",10,10.4,0.05,1,0.2);
    doppler("negative_phase_jump",10,9.6,-0.05,1,0.2);
    doppler("exact_threshold",0,0.2,0,1,0.2);
    doppler("above_threshold",0,std::nextafter(0.2,INFINITY),0,1,0.2);
    doppler("below_threshold",0,std::nextafter(0.2,0.0),0,1,0.2);
    doppler("half_second",-100,-99.5,1,0.5,0.2);
    doppler("five_seconds",100,90,-2,5,0.2);
    doppler("zero_doppler",0,1,0,1,0.2);
    doppler("zero_values",0,0,0,1,0.2);
}
