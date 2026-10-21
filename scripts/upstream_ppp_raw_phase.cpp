// Synthetic relative-code gauge consistent with the independent IF PPP fixture.
// Expected raw FLOAT ambiguities use pinned native IF coefficients/wavelengths.
// No claim of real product clock compatibility or phase OSB availability.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/core/constants.hpp>
#include <fstream>
#include <iomanip>
#include <string>

int main(int argc,char**argv) {
    if(argc!=2)return 2;
    const std::string root=argv[1];
    std::ofstream bia(root+"/synthetic_ppp_raw_phase_dsb.bsx"),truth(root+"/upstream_ppp_raw_phase.csv");
    bia << "%=BIA 1.00 SYN 2024:039:00000 SYN 2024:039:00000 2024:040:00000 R 00000018\n"
        << "+BIAS/DESCRIPTION\n TIME_SYSTEM G\n BIAS_MODE RELATIVE\n"
        << " CLOCK_DATUM SYNTHETIC_GPS_IF_C1W_C2W\n-BIAS/DESCRIPTION\n+BIAS/SOLUTION\n" << std::setprecision(17);
    truth << "satellite,calibrated_IF_ambiguity_m,unmodelled_phase_IF_bias_m,raw_FLOAT_ambiguity_m\n" << std::setprecision(17);
    const double a=libgnss::constants::IFLC_F1,b=libgnss::constants::IFLC_F2;
    for(int prn:{1,2,7,8,9,18,19,25,26}) {
        const double code1c=3.2*(prn%3)-2.0,code2w=-1.1*(prn%4)+.7;
        const double code1w=-b/a*code2w;
        for(int edge=0;edge<2;++edge) {
            const double d=edge==0?code1w-code2w:code1c-code1w;
            bia << " DSB G" << std::setfill('0') << std::setw(3) << 80+prn << " G" << std::setw(2) << prn << std::setfill(' ')
                << (edge==0?"           C1W C2W":"           C1C C1W")
                << " 2024:039:00600 2024:039:02520 ns " << d/libgnss::constants::SPEED_OF_LIGHT*1e9 << " 0.001\n";
        }
        const double phase1=.015*(prn%3)-.007,phase2=.012*(prn%4)-.020;
        const double ambiguity=a*libgnss::constants::GPS_L1_WAVELENGTH*(12+prn)+b*libgnss::constants::GPS_L2_WAVELENGTH*(-7+2*prn);
        const double bias=a*phase1+b*phase2;
        truth << 'G' << std::setfill('0') << std::setw(2) << prn << std::setfill(' ') << ',' << ambiguity << ',' << bias << ',' << ambiguity+bias << '\n';
    }
    bia << "-BIAS/SOLUTION\n%=ENDBIA\n";
    return bia&&truth?0:3;
}
