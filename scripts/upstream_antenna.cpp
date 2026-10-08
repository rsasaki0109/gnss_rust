// Native ANTEX loader/range correction, literal native NOAZI interpolation,
// and independent receiver-only injection into synthetic APC PPP products.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#define GNSS_RUST_PRECISE_TRANSMIT_MAIN precise_transmit_fixture_main
#include "upstream_precise_transmit.cpp"
#include <libgnss++/algorithms/ppp.hpp>
#include <libgnss++/algorithms/ppp_osr.hpp>
bool nativeReceiverCalibration(const std::string&,std::array<Vector3d,2>&,
                              std::array<libgnss::ReceiverPcvGrid,2>&);
// Literal receiverAntennaPcvMeters lookup; independent of the Rust code.
double pcv(const libgnss::ReceiverPcvGrid& g,double elevation) {
    const double node=(90.0-elevation*180/M_PI-g.zen1_deg)/g.dzen_deg;
    const int last=static_cast<int>(g.noazi_m.size())-1;
    if(node<=0) return g.noazi_m.front();
    if(node>=last) return g.noazi_m.back();
    const int i=static_cast<int>(std::floor(node));const double f=node-i;
    return g.noazi_m[i]*(1-f)+g.noazi_m[i+1]*f;
}
void record(std::ostream& out,std::string value,const std::string& label) {
    out << std::left << std::setw(60) << value << label << '\n';
}
int main(int argc,char**argv) {
    if(argc!=3) return 2;
    const std::string root=argv[1],fixtures=argv[2];
    {
        std::ofstream out(root+"/synthetic_receiver.atx");
        record(out,"     1.4            M","ANTEX VERSION / SYST");
        record(out,"A","PCV TYPE / REFANT");record(out,"","END OF HEADER");
        for(bool satellite:{false,true}) {
            record(out,"","START OF ANTENNA");
            std::ostringstream type;type << std::left << std::setw(20) << (satellite?"BLOCK IIF":"TESTANT NONE")
                << std::setw(20) << (satellite?"G01":"RX001") << std::setw(10) << (satellite?"G063":"");
            record(out,type.str(),"TYPE / SERIAL NO");
            record(out,"     0.0","DAZI");record(out,"     0.0    90.0    10.0","ZEN1 / ZEN2 / DZEN");
            record(out,"     2","# OF FREQUENCIES");
            record(out,"  2020     1     1     0     0    0.0000000","VALID FROM");
            record(out,"  2030     1     1     0     0    0.0000000","VALID UNTIL");
            for(int f=0;f<2;++f) {
                const std::string code=f?"G02":"G01";record(out,"   "+code,"START OF FREQUENCY");
                record(out,satellite?"      0.00      0.00   2000.00":(f?"     -5.00     18.00    125.00":"     12.00     -8.00     90.00"),"NORTH / EAST / UP");
                out << "   NOAZI";
                for(int i=0;i<10;++i) out << std::right << std::fixed << std::setprecision(2) << std::setw(8)
                    << (satellite?0.0:(f?0.3*i*i-2.0:0.2*i*i-1.0));
                out << '\n';record(out,"   "+code,"END OF FREQUENCY");
            }
            record(out,"","END OF ANTENNA");
        }
    }
    std::array<Vector3d,2> enu;std::array<libgnss::ReceiverPcvGrid,2> grids;
    if(!nativeReceiverCalibration(root+"/synthetic_receiver.atx",enu,grids)) return 3;
    const Vector3d delta(0.025,-0.035,0.8);
    std::ofstream model(root+"/upstream_receiver_antenna.csv");model << std::setprecision(17)
        << "frequency,east,north,up,azimuth,elevation,pcv,range_model_m,observation_add_m\n";
    for(int f=0;f<2;++f) for(double az:{-M_PI,0.0,0.3,M_PI/2,2.9})
        for(double deg:{-5.0,0.0,7.3,10.0,35.25,89.0,90.0,95.0}) {
            const double el=deg*M_PI/180,v=pcv(grids[f],el);
            const Vector3d neu(enu[f].y(),enu[f].x(),enu[f].z());
            const double range=libgnss::clasReceiverAntennaCorrectionMeters(delta,neu,v,az,el);
            model << (f?"G02":"G01") << ',' << enu[f].x() << ',' << enu[f].y() << ',' << enu[f].z()
                << ',' << az << ',' << el << ',' << v << ',' << range << ',' << -range << '\n';
        }
    PreciseProducts products;
    if(!products.loadSP3File(fixtures+"/synthetic_ppp.sp3") || !products.loadClockFile(fixtures+"/synthetic_ppp.clk")) return 4;
    const auto rx=geodetic2ecef(35*M_PI/180,139*M_PI/180,45.0);
    std::ofstream trace(root+"/upstream_ppp_receiver_antenna.csv");trace << std::setprecision(17)
        << "epoch,sat,azimuth,elevation,add_l1_m,add_l2_m,if_add_m\n";
    const double f1=constants::GPS_L1_FREQ,f2=constants::GPS_L2_FREQ,a=f1*f1/(f1*f1-f2*f2),b=1-a;
    for(bool noisy:{false,true}) {
        std::ifstream input(fixtures+(noisy?"/synthetic_ppp_windup_noisy.obs":"/synthetic_ppp_windup.obs"));
        std::ofstream output(root+(noisy?"/synthetic_ppp_antenna_noisy.obs":"/synthetic_ppp_antenna.obs"));
        int index=-1;std::string line;
        while(std::getline(input,line)) {
            if(!line.empty() && line[0]=='>') ++index;
            if(index<0 || line.empty() || line[0]!='G') {output << line << '\n';continue;}
            const auto id=line.substr(0,3);const SatelliteId sat(GNSSSystem::GPS,std::stoi(id.substr(1)));
            const auto tx=transmit(products,sat,GNSSTime(2300,346200+30*index),rx,false);
            const auto los=ecef2enu(tx.p-rx,35*M_PI/180,139*M_PI/180);
            const double az=std::atan2(los.x(),los.y()),el=std::atan2(los.z(),std::hypot(los.x(),los.y()));
            std::array<double,2> add;
            for(int f=0;f<2;++f) {
                const Vector3d neu(enu[f].y(),enu[f].x(),enu[f].z());
                add[f]=-libgnss::clasReceiverAntennaCorrectionMeters(delta,neu,pcv(grids[f],el),az,el);
                for(int j=0;j<2;++j) {
                    const int start=3+16*(2*f+j);const double raw=std::stod(line.substr(start,14));
                    const double wavelength=constants::SPEED_OF_LIGHT/(f?f2:f1);
                    std::ostringstream value;value << std::fixed << std::setprecision(j?4:3) << std::setw(14)
                        << raw-add[f]/(j?wavelength:1.0);
                    line.replace(start,14,value.str());
                }
            }
            if(!noisy) trace << index << ',' << id << ',' << az << ',' << el << ',' << add[0] << ',' << add[1] << ',' << a*add[0]+b*add[1] << '\n';
            output << line << '\n';
        }
        if(index!=63 || !input.eof() || !output) return 5;
    }
    return model && trace?0:6;
}
