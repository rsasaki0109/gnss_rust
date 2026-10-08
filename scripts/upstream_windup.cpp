// Native wind-up/Sun components and independent synthetic PPP phase injection.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#define GNSS_RUST_PRECISE_TRANSMIT_MAIN precise_transmit_fixture_main
#include "upstream_precise_transmit.cpp"
#include <libgnss++/algorithms/ppp.hpp>
#include <filesystem>
Vector3d nativeWindupSun(const GNSSTime&);
int main(int argc,char**argv) {
    if(argc!=3) return 2;
    const std::string root=argv[1],fixtures=argv[2];
    std::ofstream solar(root+"/upstream_windup_sun.csv"),model(root+"/upstream_windup.csv"),trace(root+"/upstream_ppp_windup.csv");
    solar << "week,tow,x,y,z\n" << std::setprecision(17);
    model << "rx_x,rx_y,rx_z,sat_x,sat_y,sat_z,sun_x,sun_y,sun_z,previous,cycles\n" << std::setprecision(17);
    // Real helper and native wind-up across hemispheres, fractional time,
    // GPS week boundaries and independent previous cycle gauges.
    for(int week:{0,1042,2300,2301}) for(double tow:{0.0,346200.25,604799.75}) {
        const auto sun=nativeWindupSun(GNSSTime(week,tow));
        solar << week << ',' << tow << ',' << sun.x() << ',' << sun.y() << ',' << sun.z() << '\n';
        for(double latitude:{-70.0,-35.0,0.0,35.0,70.0}) {
            const auto rx=geodetic2ecef(latitude*M_PI/180,139*M_PI/180,45);
            for(int prn:{1,2}) for(double previous:{-2.25,0.0,3.75}) {
                auto e=make_eph(GNSSSystem::GPS,prn);e.omega0=0.31*prn;e.m0=0.77*prn;
                Vector3d sat;double unused;state(e,GNSSTime(2300,346200),sat,unused);
                for(const auto& v:{rx,sat,sun}) for(int i=0;i<3;++i) model << v(i) << ',';
                model << previous << ',' << ppp_utils::calculatePhaseWindup(rx,sat,sun,previous) << '\n';
            }
        }
    }
    const Vector3d receiver=geodetic2ecef(35*M_PI/180,139*M_PI/180,45.0);
    PreciseProducts products;
    if(!products.loadSP3File(fixtures+"/synthetic_ppp.sp3") || !products.loadClockFile(fixtures+"/synthetic_ppp.clk")) return 3;
    trace << "epoch,sat,week,tow,sun_x,sun_y,sun_z,previous,cycles,if_add_m\n" << std::setprecision(17);
    // Add one common native wind-up cycle count to each raw carrier frequency.
    // Code, bias and products are inherited byte-for-byte; no antenna/tide.
    for(bool noisy:{false,true}) {
        std::ifstream input(fixtures+(noisy?"/synthetic_ppp_noisy.obs":"/synthetic_ppp.obs"));
        std::ofstream output(root+(noisy?"/synthetic_ppp_windup_noisy.obs":"/synthetic_ppp_windup.obs"));
        std::map<std::string,double> previous;
        int index=-1; std::string line;
        while(std::getline(input,line)) {
            if(!line.empty() && line[0]=='>') ++index;
            if(index<0 || line.empty() || line[0]!='G') {output << line << '\n';continue;}
            const auto id=line.substr(0,3);
            const SatelliteId satellite(GNSSSystem::GPS,std::stoi(id.substr(1)));
            const GNSSTime time(2300,346200.0+30*index);
            const auto tx=transmit(products,satellite,time,receiver,false);
            const auto sun=nativeWindupSun(time);
            const double prior=previous[id],cycles=ppp_utils::calculatePhaseWindup(receiver,tx.p,sun,prior);
            previous[id]=cycles;
            if(!noisy) {
                trace << index << ',' << id << ",2300," << time.tow << ',' << sun.x() << ',' << sun.y() << ',' << sun.z()
                    << ',' << prior << ',' << cycles << ',' << -cycles*constants::SPEED_OF_LIGHT/(constants::GPS_L1_FREQ+constants::GPS_L2_FREQ) << '\n';
            }
            for(int field:{1,3}) {
                const int start=3+16*field;
                const double raw=std::stod(line.substr(start,14));
                std::ostringstream value;value << std::fixed << std::setprecision(4) << std::setw(14) << raw+cycles;
                line.replace(start,14,value.str());
            }
            output << line << '\n';
        }
        if(index!=63 || !input.eof() || !output) return 4;
    }
    return solar && model && trace ? 0 : 5;
}
