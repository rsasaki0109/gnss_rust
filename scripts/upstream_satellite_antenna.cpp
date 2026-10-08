// Actual pinned satellite ANTEX loader and literal PPP nominal-yaw IF-PCO
// rotation, precise geometry/Niell/wind-up. Synthetic, not PPPProcessor.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#define GNSS_RUST_PRECISE_TRANSMIT_MAIN precise_transmit_fixture_main
#include "upstream_precise_transmit.cpp"
#include <libgnss++/algorithms/ppp.hpp>
#include <libgnss++/algorithms/ppp_osr.hpp>
#include <libgnss++/models/troposphere.hpp>
bool nativeSatelliteCalibrations(const std::string&,std::vector<libgnss::SatelliteAntexEntry>&);
bool nativeReceiverCalibration(const std::string&,std::array<Vector3d,2>&,std::array<libgnss::ReceiverPcvGrid,2>&);
Vector3d nativeWindupSun(const GNSSTime&);
void record(std::ostream& out,const std::string& value,const std::string& label) {
    out << std::left << std::setw(60) << value << label << '\n';
}
// Literal rotation from pinned ppp_corrections.cpp, independent of Rust.
Vector3d pco_ecef(const Vector3d& com,const Vector3d& sun,const Vector3d& body) {
    const Vector3d z=-com.normalized(),y=z.cross(sun-com).normalized(),x=y.cross(z);
    return body.x()*x+body.y()*y+body.z()*z;
}
Vector3d apc(const Vector3d& com,const Vector3d& sun,const Vector3d& body) {
    return com+pco_ecef(com,sun,body);
}
double pcv(const libgnss::ReceiverPcvGrid& g,double elevation) {
    const double node=(90-elevation*180/M_PI-g.zen1_deg)/g.dzen_deg;
    const int last=static_cast<int>(g.noazi_m.size())-1;
    if(node<=0) return g.noazi_m.front(); if(node>=last) return g.noazi_m.back();
    const int i=static_cast<int>(std::floor(node)); const double f=node-i;
    return g.noazi_m[i]*(1-f)+g.noazi_m[i+1]*f;
}
Vector3d body_if(const libgnss::SatelliteAntexEntry& c,libgnss::SignalType s1,libgnss::SignalType s2,double f1,double f2) {
    const auto coefficients=libgnss::ppp_utils::getIonosphereFreeCoefficients(f1,f2);
    return coefficients.first*c.body_offsets_m.at(s1)+coefficients.second*c.body_offsets_m.at(s2);
}
#ifndef GNSS_RUST_SATELLITE_ANTENNA_MAIN
#define GNSS_RUST_SATELLITE_ANTENNA_MAIN main
#endif
int GNSS_RUST_SATELLITE_ANTENNA_MAIN(int argc,char**argv) {
    if(argc!=3) return 2;
    const std::string root=argv[1],fixtures=argv[2];
    const Vector3d rx=geodetic2ecef(35*M_PI/180,139*M_PI/180,45.0);
    const auto satellites=visible(rx);
    {
        std::ofstream out(root+"/synthetic_satellite.atx");
        record(out,"     1.4            G","ANTEX VERSION / SYST");record(out,"A","PCV TYPE / REFANT");record(out,"","END OF HEADER");
        for(const auto& e:satellites) {
            const int prn=e.satellite.prn;
            record(out,"","START OF ANTENNA");std::ostringstream type;
            type << std::left << std::setw(20) << "BLOCK TEST" << std::setw(20) << e.satellite.toString()
                 << 'G' << std::right << std::setfill('0') << std::setw(3) << (100+prn);
            record(out,type.str(),"TYPE / SERIAL NO");record(out,"     0.0","DAZI");record(out,"     0.0    14.0     2.0","ZEN1 / ZEN2 / DZEN");
            record(out,"     3","# OF FREQUENCIES");
            record(out,"  2020     1     1     0     0    0.0000000","VALID FROM");record(out,"  2030     1     1     0     0    0.0000000","VALID UNTIL");
            for(int f=0;f<3;++f) {
                const std::string code=f==0?"G01":f==1?"G02":"G05";
                const Vector3d p=f==0?Vector3d(30*prn,-17*prn,2100+4*prn):f==1?Vector3d(-22*prn,14*prn,2400+7*prn):Vector3d(8*prn,31*prn,2200+5*prn);
                record(out,"   "+code,"START OF FREQUENCY");std::ostringstream pco;
                pco << std::fixed << std::setprecision(2);for(int i=0;i<3;++i) pco << std::setw(10) << p(i);
                record(out,pco.str(),"NORTH / EAST / UP");out << "   NOAZI";
                // Nonzero satellite PCV is retained but not applied by this
                // initial native-equivalent PCO path.
                for(int i=0;i<8;++i) out << std::right << std::fixed << std::setprecision(2) << std::setw(8) << (1.0+f+0.5*i);
                out << '\n';record(out,"   "+code,"END OF FREQUENCY");
            }
            record(out,"","END OF ANTENNA");
        }
    }
    // The same fictional orbit samples now declare COM. No real APC product
    // is converted by editing a comment: synthetic observations below are
    // generated from their new declared COM and independently rotated PCO.
    {
        std::ifstream in(fixtures+"/synthetic_ppp.sp3");std::ofstream out(root+"/synthetic_ppp_com.sp3");std::string line;
        while(std::getline(in,line)) {
            if(line.find("ANTENNA PHASE CENTER")!=std::string::npos) line="/* Synthetic satellite positions: CENTRE OF MASS";
            out << line << '\n';
        }
    }
    std::vector<libgnss::SatelliteAntexEntry> entries;
    if(!nativeSatelliteCalibrations(root+"/synthetic_satellite.atx",entries) || entries.size()!=9) return 3;
    std::map<std::string,libgnss::SatelliteAntexEntry> calibrations;for(const auto& e:entries) calibrations.emplace(e.satellite.toString(),e);
    std::array<Vector3d,2> rx_enu;std::array<libgnss::ReceiverPcvGrid,2> rx_grids;
    if(!nativeReceiverCalibration(fixtures+"/synthetic_receiver.atx",rx_enu,rx_grids)) return 4;
    PreciseProducts products;
    if(!products.loadSP3File(root+"/synthetic_ppp_com.sp3") || !products.loadClockFile(fixtures+"/synthetic_ppp.clk") || products.orbits_at_antenna_phase_center) return 5;
    using libgnss::SignalType;
    const double f1=constants::GPS_L1_FREQ,f2=constants::GPS_L2_FREQ,f5=constants::GPS_L5_FREQ;
    const auto coeff=libgnss::ppp_utils::getIonosphereFreeCoefficients(f1,f2);
    std::ofstream model(root+"/upstream_satellite_antenna.csv");model << std::setprecision(17)
      << "sat,week,tow,signal1,signal2,com_x,com_y,com_z,sun_x,sun_y,sun_z,body_x,body_y,body_z,offset_x,offset_y,offset_z,apc_x,apc_y,apc_z,range,elevation,mapping,predicted,hx,hy,hz,clock_s\n";
    for(int i:{0,31,63}) for(const auto& e:satellites) for(int pair=0;pair<3;++pair) {
        const GNSSTime t(2300,346200+30*i);const auto tx=transmit(products,e.satellite,t,rx,false);const auto sun=nativeWindupSun(t);
        const auto body=pair==0?body_if(calibrations.at(e.satellite.toString()),SignalType::GPS_L1CA,SignalType::GPS_L2C,f1,f2)
           :pair==1?body_if(calibrations.at(e.satellite.toString()),SignalType::GPS_L2C,SignalType::GPS_L1CA,f2,f1)
           :body_if(calibrations.at(e.satellite.toString()),SignalType::GPS_L1CA,SignalType::GPS_L5,f1,f5);
        const auto pos=apc(tx.p,sun,body);const Vector3d offset=pco_ecef(tx.p,sun,body);
        const auto enu=ecef2enu(pos-rx,35*M_PI/180,139*M_PI/180);const double el=std::atan2(enu.z(),std::hypot(enu.x(),enu.y()));
        const double mapping=libgnss::models::niellHydrostaticMapping(35*M_PI/180,45.0,el,39),range=geodist(pos,rx);
        const auto los=(pos-rx).normalized();const double common=range+(2e-6+i*5e-10)*constants::SPEED_OF_LIGHT-tx.total*constants::SPEED_OF_LIGHT+mapping*2.4;
        model << e.satellite.toString() << ",2300," << t.tow << ',' << (pair==1?"2W":"1C") << ',' << (pair==0?"2W":pair==1?"1C":"5Q");
        for(const auto& v:{tx.p,sun,body,offset,pos}) for(int j=0;j<3;++j) model << ',' << v(j);
        model << ',' << range << ',' << el << ',' << mapping << ',' << common
              << ',' << -los.x()-constants::OMEGA_E*pos.y()/constants::SPEED_OF_LIGHT
              << ',' << -los.y()+constants::OMEGA_E*pos.x()/constants::SPEED_OF_LIGHT << ',' << -los.z() << ',' << tx.total << '\n';
    }
    std::ofstream trace(root+"/upstream_ppp_satellite_antenna.csv");trace << std::setprecision(17)
        << "epoch,sat,apc_x,apc_y,apc_z,range,elevation,mapping,predicted,windup,receiver_if_add,admitted\n";
    for(bool noisy:{false,true}) {
        std::ifstream in(fixtures+(noisy?"/synthetic_ppp_noisy.obs":"/synthetic_ppp.obs"));
        std::ofstream out(root+(noisy?"/synthetic_ppp_satellite_antenna_noisy.obs":"/synthetic_ppp_satellite_antenna.obs"));
        std::map<std::string,double> previous;int index=-1;std::string line;
        while(std::getline(in,line)) {
            if(!line.empty() && line[0]=='>') ++index;
            if(index<0 || line.empty() || line[0]!='G') {out << line << '\n';continue;}
            const auto id=line.substr(0,3);const SatelliteId sat(GNSSSystem::GPS,std::stoi(id.substr(1)));
            const GNSSTime t(2300,346200+30*index);const auto tx=transmit(products,sat,t,rx,false);const auto sun=nativeWindupSun(t);
            const Vector3d body=body_if(calibrations.at(id),SignalType::GPS_L1CA,SignalType::GPS_L2C,f1,f2),pos=apc(tx.p,sun,body);
            const auto enu=ecef2enu(pos-rx,35*M_PI/180,139*M_PI/180);
            const double az=std::atan2(enu.x(),enu.y()),el=std::atan2(enu.z(),std::hypot(enu.x(),enu.y()));
            const double mapping=libgnss::models::niellHydrostaticMapping(35*M_PI/180,45.0,el,39);
            const double range=geodist(pos,rx),common=range+(2e-6+index*5e-10)*constants::SPEED_OF_LIGHT-tx.total*constants::SPEED_OF_LIGHT+mapping*2.4;
            const double cycles=libgnss::ppp_utils::calculatePhaseWindup(rx,pos,sun,previous[id]);previous[id]=cycles;
            std::array<double,2> rx_add;
            for(int f=0;f<2;++f) {
                const Vector3d neu(rx_enu[f].y(),rx_enu[f].x(),rx_enu[f].z());
                rx_add[f]=-libgnss::clasReceiverAntennaCorrectionMeters(Vector3d(0.025,-0.035,0.8),neu,pcv(rx_grids[f],el),az,el);
                const double lambda=constants::SPEED_OF_LIGHT/(f?f2:f1);
                const int prn=sat.prn;
                const double ion1=4+0.01*prn+0.002*(30*index),ion=f?ion1*f1*f1/(f2*f2):ion1;
                const double code_bias=f?-1.1*(prn%4)+0.7:3.2*(prn%3)-2.0;
                const double phase_bias=f?0.012*(prn%4)-0.020:0.015*(prn%3)-0.007;
                const double code_noise=noisy?0.05*std::sin(0.37*index+0.53*prn+f):0.0;
                const double phase_noise=noisy?0.001*std::sin(0.29*index+0.31*prn+f):0.0;
                for(int j=0;j<2;++j) {
                    // Generate the complete new native model before a single
                    // RINEX rounding, rather than re-rounding old carriers.
                    const int start=3+16*(2*f+j);
                    const double raw=j?(common-rx_add[f]-ion+phase_bias+phase_noise)/lambda+(f?-7+2*prn:12+prn)+cycles
                        :common-rx_add[f]+ion+code_bias+code_noise;
                    std::ostringstream value;value << std::fixed << std::setprecision(j?4:5) << std::setw(14)
                        << raw;
                    line.replace(start,14,value.str());
                }
            }
            if(!noisy) trace << index << ',' << id << ',' << pos.x() << ',' << pos.y() << ',' << pos.z() << ',' << range << ',' << el << ',' << mapping
                << ',' << range+(2e-6+index*5e-10)*constants::SPEED_OF_LIGHT-tx.total*constants::SPEED_OF_LIGHT+mapping*2.4
                << ',' << cycles << ',' << coeff.first*rx_add[0]+coeff.second*rx_add[1] << ',' << (el>=10*M_PI/180) << '\n';
            out << line << '\n';
        }
        if(index!=63 || !in.eof() || !out) return 6;
    }
    return model && trace?0:7;
}
