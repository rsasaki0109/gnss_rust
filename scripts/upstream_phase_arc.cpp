// Native GF threshold helper plus independent ordinary-PPP predicate reference.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include "ppp_internal.hpp"
#include <libgnss++/core/constants.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace libgnss;
int main(int argc,char** argv) {
    if (argc!=3) return 2;
    std::ofstream out(argv[1]),obs(argv[2]);
    if (!out || !obs) return 3;
    out << "configured,gf_threshold,mw_threshold,previous_gf,current_gf,previous_mw,current_mw,have_previous_mw,have_current_mw,gf_slip,mw_slip\n";
    out << std::setprecision(17);
    for (double configured : {0.05,0.2,0.5,0.75}) {
        const double gf=ppp_internal::geometryFreeSlipThresholdMeters(false,configured);
        const double mw=std::max(10.0,configured*100.0);
        for (int i=0;i<12;++i) {
            double cg=0.0,cm=0.0;
            bool pg=true,pc=true;
            if (i==1) cg=gf;
            if (i==2) cg=std::nextafter(gf,INFINITY);
            if (i==3) cg=-gf;
            if (i==4) cg=std::nextafter(-gf,-INFINITY);
            if (i==5) cm=mw;
            if (i==6) cm=std::nextafter(mw,INFINITY);
            if (i==7) cm=-mw;
            if (i==8) cm=std::nextafter(-mw,-INFINITY);
            if (i==9) {cg=gf*2;cm=mw*2;}
            if (i==10) {pg=false;cm=mw*2;}
            if (i==11) {pc=false;cm=mw*2;}
            // These strict comparisons are the literal ordinary-PPP conditions
            // from ppp_state.cpp. PPPProcessor itself is not executed here.
            const bool gs=std::abs(cg)>gf;
            const bool ms=pg && pc && std::abs(cm)>mw;
            out << configured << ',' << gf << ',' << mw << ",0," << cg << ",0," << cm << ',' << pg << ',' << pc << ',' << gs << ',' << ms << '\n';
        }
    }
    auto header=[&](const std::string& value,const std::string& label) {
        obs << std::left << std::setw(60) << value << label << '\n' << std::right;
    };
    header("     3.04           O                   G","RINEX VERSION / TYPE");
    header("Native-component continuous phase-arc synthetic fixture","PGM / RUN BY / DATE");
    header("G    4 C1C L1C C2W L2W","SYS / # / OBS TYPES");
    header("","END OF HEADER");
    const double f1=constants::GPS_L1_FREQ,f2=constants::GPS_L2_FREQ;
    const double lambda1=constants::SPEED_OF_LIGHT/f1,lambda2=constants::SPEED_OF_LIGHT/f2;
    // 16 epochs, with slow ionosphere and constant integer arcs except explicit
    // slips. G02 remains continuous until the receiver-wide power/event epochs.
    for (int i=0;i<16;++i) {
        const int seconds=i*30+(i>=14?150:0);
        const int flag=i==12?6:(i==10?1:0);
        const int count=i==12?1:(i==8?1:2);
        obs << "> 2024 02 08 00 " << std::setfill('0') << std::setw(2) << seconds/60
            << ' ' << std::setfill(' ') << std::fixed << std::setprecision(7)
            << std::setw(10) << double(seconds%60) << "  " << flag << ' ' << count << '\n';
        for (int sat=1;sat<=2;++sat) {
            if ((i==8 && sat==1) || (i==12 && sat==2)) continue;
            const double range=23000000.125+sat*1000+seconds*5.0;
            const double ion=0.001*seconds;
            const double ion2=ion*f1*f1/(f2*f2);
            // LLI slip at epoch 2; unannounced phase jump at epoch 6.
            const double n1=12+(sat==1 && i>=2?4:0)+(sat==1 && i>=6?5:0);
            const double n2=-7;
            const double p1=range+ion+(sat==1 && i==15?25.0:0.0),p2=range+ion2;
            const double l1=(range-ion)/lambda1+n1,l2=(range-ion2)/lambda2+n2;
            auto code=[&](double value) {obs << std::setprecision(5) << std::setw(14) << value << "  ";};
            auto phase=[&](double value,int lli,bool missing) {
                if (missing) obs << std::string(16,' ');
                else obs << std::setprecision(4) << std::setw(14) << value << lli << '7';
            };
            obs << (sat==1?"G01":"G02");
            code(p1);phase(l1,sat==1 && i==2?1:0,false);
            code(p2);phase(l2,sat==1 && i==11?2:0,sat==1 && i==4);
            obs << '\n';
        }
    }
    return out && obs ? 0 : 4;
}
