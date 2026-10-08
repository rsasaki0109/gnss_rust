// Actual pinned legacy bodyTideDisplacement and Sun/Moon helpers. NOT IERS.
// Complete synthetic native geometry/PCO/receiver/wind-up before one rounding.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#define GNSS_RUST_SATELLITE_ANTENNA_MAIN satellite_antenna_fixture_main
#include "upstream_satellite_antenna.cpp"
Vector3d nativeLegacyMoon(const GNSSTime&);
Vector3d nativeLegacyBodyTide(const Vector3d&,const Vector3d&,double);
Vector3d tide(const Vector3d& nominal,const GNSSTime& time) {
    return nativeLegacyBodyTide(nominal,nativeWindupSun(time),1.32712440018e20)
        +nativeLegacyBodyTide(nominal,nativeLegacyMoon(time),4.902801e12);
}
struct Model {
    Tx tx;
    Vector3d displacement,receiver,satellite;
    double range,elevation,mapping,prediction;
};
Model model(const PreciseProducts& products,const SatelliteId& sat,const GNSSTime& time,
            const Vector3d& nominal,const libgnss::SatelliteAntexEntry& calibration,double clock) {
    Model m;m.displacement=tide(nominal,time);m.receiver=nominal+m.displacement;
    m.tx=transmit(products,sat,time,m.receiver,false);
    const auto body=body_if(calibration,libgnss::SignalType::GPS_L1CA,libgnss::SignalType::GPS_L2C,constants::GPS_L1_FREQ,constants::GPS_L2_FREQ);
    m.satellite=apc(m.tx.p,nativeWindupSun(time),body);
    double latitude,longitude,height;ecef2geodetic(m.receiver,latitude,longitude,height);
    const auto enu=ecef2enu(m.satellite-m.receiver,latitude,longitude);
    m.elevation=std::atan2(enu.z(),std::hypot(enu.x(),enu.y()));
    m.mapping=libgnss::models::niellHydrostaticMapping(latitude,height,m.elevation,39);
    m.range=geodist(m.satellite,m.receiver);
    m.prediction=m.range+clock*constants::SPEED_OF_LIGHT-m.tx.total*constants::SPEED_OF_LIGHT+m.mapping*2.4;
    return m;
}
int main(int argc,char**argv) {
    if(argc!=3) return 2;
    const std::string root=argv[1],fixtures=argv[2];
    std::ofstream moon(root+"/upstream_legacy_moon.csv"),components(root+"/upstream_legacy_solid_tide.csv");
    moon << std::setprecision(17) << "week,tow,x,y,z\n";
    components << std::setprecision(17) << "week,tow,nom_x,nom_y,nom_z,sun_x,sun_y,sun_z,moon_x,moon_y,moon_z,sun_dx,sun_dy,sun_dz,moon_dx,moon_dy,moon_dz,dx,dy,dz,instant_x,instant_y,instant_z\n";
    for(int week:{0,1042,2300,2301}) for(double tow:{0.0,346200.25,604799.75}) {
        const GNSSTime time(week,tow);const auto s=nativeWindupSun(time),m=nativeLegacyMoon(time);
        moon << week << ',' << tow << ',' << m.x() << ',' << m.y() << ',' << m.z() << '\n';
        int index=0;
        for(double latitude:{-90.0,-60.0,-35.0,0.0,35.0,60.0,90.0}) for(double longitude:{-170.0,0.0,139.0}) {
            const double height=std::array<double,3>{-50.0,45.0,1500.0}[index++%3];
            const Vector3d rx=geodetic2ecef(latitude*M_PI/180,longitude*M_PI/180,height);
            const auto ds=nativeLegacyBodyTide(rx,s,1.32712440018e20),dm=nativeLegacyBodyTide(rx,m,4.902801e12);
            const Vector3d displacement=ds+dm,instant=rx+displacement;
            components << week << ',' << tow;
            for(const auto& v:{rx,s,m,ds,dm,displacement,instant}) for(int i=0;i<3;++i) components << ',' << v(i);
            components << '\n';
        }
    }
    PreciseProducts products;
    if(!products.loadSP3File(fixtures+"/synthetic_ppp_com.sp3") || !products.loadClockFile(fixtures+"/synthetic_ppp.clk")) return 3;
    std::vector<libgnss::SatelliteAntexEntry> entries;
    if(!nativeSatelliteCalibrations(fixtures+"/synthetic_satellite.atx",entries)) return 4;
    std::map<std::string,libgnss::SatelliteAntexEntry> calibrations;for(const auto& e:entries) calibrations.emplace(e.satellite.toString(),e);
    std::array<Vector3d,2> rx_enu;std::array<libgnss::ReceiverPcvGrid,2> grids;
    if(!nativeReceiverCalibration(fixtures+"/synthetic_receiver.atx",rx_enu,grids)) return 5;
    const Vector3d nominal=geodetic2ecef(35*M_PI/180,139*M_PI/180,45.0),seed=nominal+Vector3d(5,-3,2);
    std::ofstream models(root+"/upstream_ppp_legacy_tide_model.csv");models << std::setprecision(17)
        << "sat,week,tow,nom_x,nom_y,nom_z,instant_x,instant_y,instant_z,com_x,com_y,com_z,apc_x,apc_y,apc_z,range,elevation,mapping,prediction,hx,hy,hz,clock_s,tau_s\n";
    for(int i:{0,31,63}) for(const auto& rx:{nominal,seed}) for(const auto& [id,c]:calibrations) {
        const GNSSTime time(2300,346200+30*i);const auto m=model(products,c.satellite,time,rx,c,2e-6+i*5e-10);
        models << id << ",2300," << time.tow;
        for(const auto& v:{rx,m.receiver,m.tx.p,m.satellite}) for(int j=0;j<3;++j) models << ',' << v(j);
        const Vector3d los=(m.satellite-m.receiver).normalized();
        models << ',' << m.range << ',' << m.elevation << ',' << m.mapping << ',' << m.prediction
            << ',' << -los.x()-constants::OMEGA_E*m.satellite.y()/constants::SPEED_OF_LIGHT
            << ',' << -los.y()+constants::OMEGA_E*m.satellite.x()/constants::SPEED_OF_LIGHT << ',' << -los.z()
            << ',' << m.tx.total << ',' << m.tx.tau << '\n';
    }
    std::ofstream truth(root+"/synthetic_ppp_legacy_tide_truth.csv"),trace(root+"/upstream_ppp_legacy_tide.csv"),exact(root+"/upstream_ppp_legacy_tide_unrounded.csv");
    truth << std::setprecision(17) << "epoch,week,tow,nom_x,nom_y,nom_z,dx,dy,dz,instant_x,instant_y,instant_z\n";
    trace << std::setprecision(17) << "epoch,sat,windup,receiver_if_add,admitted\n";
    exact << std::setprecision(17) << "epoch,sat,code1_m,phase1_cycles,code2_m,phase2_cycles\n";
    const double f1=constants::GPS_L1_FREQ,f2=constants::GPS_L2_FREQ;
    const auto coefficients=libgnss::ppp_utils::getIonosphereFreeCoefficients(f1,f2);
    for(bool noisy:{false,true}) {
        std::ifstream input(fixtures+(noisy?"/synthetic_ppp_noisy.obs":"/synthetic_ppp.obs"));
        std::ofstream output(root+(noisy?"/synthetic_ppp_legacy_tide_noisy.obs":"/synthetic_ppp_legacy_tide.obs"));
        std::map<std::string,double> previous;int index=-1;std::string line;
        while(std::getline(input,line)) {
            if(!line.empty() && line[0]=='>') {
                ++index;if(!noisy) {
                    const auto displacement=tide(nominal,GNSSTime(2300,346200+30*index));const Vector3d instant=nominal+displacement;
                    truth << index << ",2300," << 346200+30*index;
                    for(const auto& v:{nominal,displacement,instant}) for(int j=0;j<3;++j) truth << ',' << v(j);
                    truth << '\n';
                }
            }
            if(index<0 || line.empty() || line[0]!='G') {output << line << '\n';continue;}
            const auto id=line.substr(0,3);const auto& calibration=calibrations.at(id);
            const GNSSTime time(2300,346200+30*index);const auto m=model(products,calibration.satellite,time,nominal,calibration,2e-6+index*5e-10);
            double latitude,longitude,height;ecef2geodetic(m.receiver,latitude,longitude,height);
            const auto enu=ecef2enu(m.satellite-m.receiver,latitude,longitude);const double az=std::atan2(enu.x(),enu.y());
            const double cycles=libgnss::ppp_utils::calculatePhaseWindup(m.receiver,m.satellite,nativeWindupSun(time),previous[id]);previous[id]=cycles;
            const int prn=calibration.satellite.prn;std::array<double,2> add;
            if(!noisy) exact << index << ',' << id;
            for(int f=0;f<2;++f) {
                const Vector3d neu(rx_enu[f].y(),rx_enu[f].x(),rx_enu[f].z());
                add[f]=-libgnss::clasReceiverAntennaCorrectionMeters(Vector3d(0.025,-0.035,0.8),neu,pcv(grids[f],m.elevation),az,m.elevation);
                const double lambda=constants::SPEED_OF_LIGHT/(f?f2:f1),ion1=4+0.01*prn+0.002*(30*index),ion=f?ion1*f1*f1/(f2*f2):ion1;
                const double code_bias=f?-1.1*(prn%4)+0.7:3.2*(prn%3)-2.0,phase_bias=f?0.012*(prn%4)-0.020:0.015*(prn%3)-0.007;
                const double code_noise=noisy?0.05*std::sin(0.37*index+0.53*prn+f):0.0,phase_noise=noisy?0.001*std::sin(0.29*index+0.31*prn+f):0.0;
                for(int j=0;j<2;++j) {
                    const double raw=j?(m.prediction-add[f]-ion+phase_bias+phase_noise)/lambda+(f?-7+2*prn:12+prn)+cycles
                        :m.prediction-add[f]+ion+code_bias+code_noise;
                    if(!noisy) exact << ',' << raw;
                    std::ostringstream value;value << std::fixed << std::setprecision(j?4:5) << std::setw(14) << raw;
                    line.replace(3+16*(2*f+j),14,value.str());
                }
            }
            if(!noisy) exact << '\n';
            if(!noisy) trace << index << ',' << id << ',' << cycles << ',' << coefficients.first*add[0]+coefficients.second*add[1] << ',' << (m.elevation>=10*M_PI/180) << '\n';
            output << line << '\n';
        }
        if(index!=63 || !input.eof() || !output) return 6;
    }
    return moon && components && models && truth && trace && exact?0:7;
}
