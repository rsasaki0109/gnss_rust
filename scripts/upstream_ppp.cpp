// Independent initial IF PPP model/observation reference, not PPPProcessor.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#define GNSS_RUST_PRECISE_TRANSMIT_MAIN precise_transmit_fixture_main
#include "upstream_precise_transmit.cpp"
#include <filesystem>
#include <libgnss++/models/troposphere.hpp>

int main(int argc,char** argv) {
    if (argc!=2) return 2;
    const std::string root=argv[1];
    const double rad=std::acos(-1.0)/180.0;
    const Vector3d receiver=geodetic2ecef(35*rad,139*rad,45.0),seed=receiver+Vector3d(5,-3,2);
    const auto satellites=visible(receiver);
    spp_products(root,satellites);
    std::filesystem::rename(root+"/synthetic_precise_spp.sp3",root+"/synthetic_ppp.sp3");
    std::filesystem::rename(root+"/synthetic_precise_spp.clk",root+"/synthetic_ppp.clk");
    PreciseProducts products;
    if(!products.loadSP3File(root+"/synthetic_ppp.sp3") || !products.loadClockFile(root+"/synthetic_ppp.clk")) return 3;
    std::ofstream mapping(root+"/upstream_niell.csv"),model(root+"/upstream_ppp_model.csv"),truth(root+"/synthetic_ppp_truth.csv"),bias(root+"/synthetic_ppp_bias.csv");
    mapping << "latitude_rad,height_m,elevation_rad,day,mapping\n" << std::setprecision(17);
    int index=0;
    for(double latitude : {-90.0,-60.0,-35.0,0.0,15.0,30.0,60.0,90.0}) {
        const double height=index%3==0?-100.0:(index%3==1?45.0:2000.0);++index;
        for(double elevation:{0.0,5*rad,90*rad}) for(int day:{1,39,200,366})
            mapping << latitude*rad << ',' << height << ',' << elevation << ',' << day << ',' << models::niellHydrostaticMapping(latitude*rad,height,elevation,day) << '\n';
    }
    model << "sat,week,tow,rx_x,rx_y,rx_z,clock_m,zenith_m,range_m,total_clock_s,elevation,mapping,predicted_m,h_x,h_y,h_z,h_clock,h_zenith\n" << std::setprecision(17);
    for(int i:{0,31,63}) for(const Vector3d& rx:{receiver,seed}) for(const auto& e:satellites) {
        const GNSSTime time(2300,346200.0+30*i);
        const auto tx=transmit(products,e.satellite,time,rx,false);
        double latitude,longitude,height;ecef2geodetic(rx,latitude,longitude,height);
        const auto difference=tx.p-rx;const auto los=difference.normalized();const auto enu=ecef2enu(difference,latitude,longitude);
        const double elevation=std::atan2(enu(2),std::hypot(enu(0),enu(1))),nmf=models::niellHydrostaticMapping(latitude,height,elevation,39);
        const double range=geodist(tx.p,rx),clock=(2e-6+i*5e-10)*constants::SPEED_OF_LIGHT;
        const double predicted=range+clock-tx.total*constants::SPEED_OF_LIGHT+nmf*2.4;
        model << e.satellite.toString() << ",2300," << time.tow;
        for(int j=0;j<3;++j) model << ',' << rx(j);
        model << ',' << clock << ",2.4," << range << ',' << tx.total << ',' << elevation << ',' << nmf << ',' << predicted
            << ',' << -los(0)-constants::OMEGA_E*tx.p(1)/constants::SPEED_OF_LIGHT
            << ',' << -los(1)+constants::OMEGA_E*tx.p(0)/constants::SPEED_OF_LIGHT
            << ',' << -los(2) << ",1," << nmf << '\n';
    }
    std::ofstream admission(root+"/upstream_ppp_admission.csv");
    admission << "epoch,sat,elevation_rad,admitted\n" << std::setprecision(17);
    for(int i=0;i<64;++i) for(const auto& e:satellites) {
        const auto tx=transmit(products,e.satellite,GNSSTime(2300,346200.0+30*i),receiver,false);
        const auto enu=ecef2enu(tx.p-receiver,35*rad,139*rad);
        const double elevation=std::atan2(enu(2),std::hypot(enu(0),enu(1)));
        admission << i << ',' << e.satellite.toString() << ',' << elevation << ',' << (elevation>=10*rad) << '\n';
    }
    truth << "# receiverXYZ,initialXYZ,initial_receiver_clock_s,zenith_m\n" << std::setprecision(17);
    for(const auto& p:{receiver,seed}) for(int j=0;j<3;++j) truth << p(j) << ',';
    truth << "0,2.4\n# satellite,N1,N2,IF_ambiguity_m\n";
    const double a=constants::IFLC_F1,b=constants::IFLC_F2,lam1=constants::GPS_L1_WAVELENGTH,lam2=constants::GPS_L2_WAVELENGTH;
    bias << "# satellite,tracking,source,from_week,from_tow,until_week,until_tow,code_add_m,phase_add_m,clock_reference,phase_reference\n" << std::setprecision(17);
    for(const auto& e:satellites) {
        const int prn=e.satellite.prn,n1=12+prn,n2=-7+2*prn;
        truth << e.satellite.toString() << ',' << n1 << ',' << n2 << ',' << a*lam1*n1+b*lam2*n2 << '\n';
        for(int j=0;j<2;++j) {
            const double code=j==0?3.2*(prn%3)-2.0:-1.1*(prn%4)+0.7;
            const double phase=j==0?0.015*(prn%3)-0.007:0.012*(prn%4)-0.020;
            bias << e.satellite.toString() << ',' << (j==0?"1C":"2W") << ",RinexClk,2300,346200,2300,348090," << -code << ',' << -phase << ",synthetic-precise-clock-v1,synthetic-phase-v1\n";
        }
    }
    for(bool noisy:{false,true}) {
        std::ofstream obs(root+(noisy?"/synthetic_ppp_noisy.obs":"/synthetic_ppp.obs"));
        header(obs,"     3.04           O                   G","RINEX VERSION / TYPE");
        std::ostringstream pos;pos << std::fixed << std::setprecision(4);for(int j=0;j<3;++j)pos << std::setw(14) << seed(j);
        header(obs,pos.str(),"APPROX POSITION XYZ");header(obs,"G    4 C1C L1C C2W L2W","SYS / # / OBS TYPES");header(obs,"","END OF HEADER");
        for(int i=0;i<64;++i) {
            const int seconds=600+30*i;
            obs << "> " << epoch(seconds) << "  0" << std::setw(3) << satellites.size() << '\n';
            for(const auto& e:satellites) {
                const int prn=e.satellite.prn;
                const auto tx=transmit(products,e.satellite,GNSSTime(2300,345600.0+seconds),receiver,false);
                const auto enu=ecef2enu(tx.p-receiver,35*rad,139*rad);
                const double elevation=std::atan2(enu(2),std::hypot(enu(0),enu(1)));
                const double common=geodist(tx.p,receiver)+(2e-6+i*5e-10)*constants::SPEED_OF_LIGHT-tx.total*constants::SPEED_OF_LIGHT
                    +models::niellHydrostaticMapping(35*rad,45.0,elevation,39)*2.4;
                const double ion1=4+0.01*prn+0.002*(seconds-600),ion2=ion1*constants::GPS_L1_FREQ*constants::GPS_L1_FREQ/(constants::GPS_L2_FREQ*constants::GPS_L2_FREQ);
                obs << e.satellite.toString();
                for(int j=0;j<2;++j) {
                    const double lambda=j==0?lam1:lam2,ion=j==0?ion1:ion2,n=j==0?12+prn:-7+2*prn;
                    const double code_bias=j==0?3.2*(prn%3)-2.0:-1.1*(prn%4)+0.7;
                    const double phase_bias=j==0?0.015*(prn%3)-0.007:0.012*(prn%4)-0.020;
                    const double code_noise=noisy?0.05*std::sin(0.37*i+0.53*prn+j):0;
                    const double phase_noise=noisy?0.001*std::sin(0.29*i+0.31*prn+j):0;
                    const double code=common+ion+code_bias+code_noise;
                    const double phase=(common-ion+phase_bias+phase_noise)/lambda+n;
                    obs << std::fixed << std::setprecision(5) << std::setw(14) << code << "  "
                        << std::setprecision(4) << std::setw(14) << phase << "07";
                }
                obs << '\n';
            }
        }
    }
    return mapping && model && admission && truth && bias ? 0 : 4;
}
