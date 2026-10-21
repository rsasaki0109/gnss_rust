// Actual pinned PPP PVA prediction and independently generated 5 Hz IF physics.
// Private access is confined to this adapter; no native source is patched.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <bits/stdc++.h>
#include <Eigen/Dense>
#define private public
#include <libgnss++/algorithms/ppp.hpp>
#undef private
#include <libgnss++/io/rinex.hpp>
#define GNSS_RUST_PRECISE_TRANSMIT_MAIN acceleration_unused_transmit_main
#include "upstream_precise_transmit.cpp"
#include <libgnss++/models/troposphere.hpp>

static PPPProcessor::PPPConfig pva_config() {
    PPPProcessor::PPPConfig c;c.kinematic_mode=true;c.low_dynamics_mode=false;
    c.use_dynamics_model=true;c.use_pva_dynamics=true;c.estimate_ionosphere=false;
    return c;
}
static PositionSolution header_seed(const Vector3d& xyz,GNSSTime time) {
    PositionSolution s;s.time=time;s.position_ecef=xyz;s.receiver_clock_bias=0;
    s.status=SolutionStatus::SPP;s.num_satellites=6;return s;
}
static void values(std::ostream& o,const Eigen::VectorXd& x,const std::vector<int>& ix) {
    for(int i:ix)o << ',' << x[i];
}
static void values(std::ostream& o,const Eigen::MatrixXd& p,const std::vector<int>& ix) {
    for(int i:ix)for(int j:ix)o << ',' << p(i,j);
}
static void prediction(const std::string& root) {
    std::ofstream out(root+"/upstream_ppp_acceleration_predict.csv");
    out << "# id,n,dt,qpos,qtrop,qamb,sigma_h,sigma_v,x,P,native_x,native_P; Rust XYZ/clock/ZTD/V/A/IFN layout\n" << std::setprecision(17);
    const double rad=std::acos(-1.0)/180;
    const std::vector<double> lat={0,35,-60,89.999,-89.999,45},lon={-179.9,139,0},gaps={0,.2,1,10,120},hs={0,.1,1,3},vs={0,.1,.3};
    for(int id=0;id<48;++id) {
        auto c=pva_config();c.process_noise_position=.003*(id%3);
        c.process_noise_troposphere=1e-8*(1+id%3);c.process_noise_ambiguity=1e-8*(1+id%4);
        c.process_noise_acceleration_horizontal=hs[id%hs.size()];c.process_noise_acceleration_vertical=vs[id%vs.size()];
        PPPProcessor processor(c);processor.precise_products_loaded_=true;
        ObservationData observation;observation.clear();observation.time=GNSSTime(2300,346200);
        const auto seed=header_seed(geodetic2ecef(lat[id%lat.size()]*rad,lon[id%lon.size()]*rad,id%2?45:-50),observation.time);
        NavigationData nav;if(!processor.initializeFilter(observation,nav,&seed))throw std::runtime_error("native PVA init failed");
        processor.filter_initialized_=true;auto& f=processor.filter_state_;const int n=14+id%4;
        f.total_states=n;f.state=Eigen::VectorXd::Zero(n);Eigen::MatrixXd matrix(n,n);
        for(int i=0;i<n;++i) {f.state[i]=.23*(i+1)+.01*id;for(int j=0;j<n;++j)matrix(i,j)=.1*std::sin(.7*(i+1)+1.3*(j+1)+id);}
        f.state.head<3>()=seed.position_ecef;f.covariance=matrix*matrix.transpose()+.2*Eigen::MatrixXd::Identity(n,n);
        std::vector<int> ix={0,1,2,9,11,3,4,5,6,7,8};for(int i=12;i<n;++i)ix.push_back(i);
        const double dt=gaps[id%gaps.size()];out << id << ',' << ix.size() << ',' << dt << ',' << c.process_noise_position << ',' << c.process_noise_troposphere << ',' << c.process_noise_ambiguity << ',' << c.process_noise_acceleration_horizontal << ',' << c.process_noise_acceleration_vertical;
        values(out,f.state,ix);values(out,f.covariance,ix);processor.predictState(dt,nullptr);values(out,f.state,ix);values(out,f.covariance,ix);out << '\n';
    }
}
static void initialization(const std::string& dataset,const std::string& root) {
    std::ofstream out(root+"/upstream_ppp_acceleration_initial.csv");out << "case,week,tow,x,y,z,clock_m,ZTD_m,pos_var,vel_var,acc_var,clock_var,ZTD_var,Vx,Vy,Vz,Ax,Ay,Az\n" << std::setprecision(17);
    for(const std::string city:{"nagoya","tokyo"})for(const std::string run:{"run1","run2","run3"}) {
        io::RINEXReader reader;if(!reader.open(dataset+"/"+city+"/"+run+"/rover.obs"))throw std::runtime_error("PPC rover unavailable");
        io::RINEXReader::RINEXHeader h;ObservationData obs;if(!reader.readHeader(h)||!h.has_approximate_position||!reader.readObservationEpoch(obs))throw std::runtime_error("PPC header unavailable");
        PPPProcessor p(pva_config());p.precise_products_loaded_=true;NavigationData nav;const auto seed=header_seed(h.approximate_position,obs.time);
        if(!p.initializeFilter(obs,nav,&seed))throw std::runtime_error("PVA header initialization failed");
        const auto& f=p.filter_state_;out << city << '-' << run << ',' << obs.time.week << ',' << obs.time.tow;
        for(int i:std::vector<int>{0,1,2,f.clock_index,f.trop_index})out << ',' << f.state[i];
        for(int i:std::vector<int>{0,f.vel_index,f.accel_index,f.clock_index,f.trop_index})out << ',' << f.covariance(i,i);
        for(int i=3;i<9;++i) {out << ',' << f.state[i];}
        out << '\n';
    }
}
static std::string fractional_epoch(double seconds) {
    std::ostringstream out;out << "2024 02 08 00 " << std::setfill('0') << std::setw(2) << int(seconds)/60 << ' ' << std::setfill(' ') << std::fixed << std::setprecision(8) << std::setw(11) << std::fmod(seconds,60);return out.str();
}
static void car_physics(const std::string& inputs,const std::string& root) {
    PreciseProducts products;if(!products.loadSP3File(inputs+"/synthetic_ppp.sp3")||!products.loadClockFile(inputs+"/synthetic_ppp.clk"))throw std::runtime_error("synthetic precise inputs unavailable");
    const double rad=std::acos(-1.0)/180;const Vector3d origin=geodetic2ecef(35*rad,139*rad,45),seed=origin+Vector3d(5,-3,2);const auto satellites=visible(origin);
    std::ofstream truth(root+"/synthetic_ppp_acceleration_truth.csv"),unrounded(root+"/upstream_ppp_acceleration_unrounded.csv");
    truth << "epoch,week,tow,x,y,z,vx,vy,vz,ax,ay,az,clock_s,ZTD_m\n" << std::setprecision(17);unrounded << "epoch,satellite,code1_m,phase1_cycles,code2_m,phase2_cycles\n" << std::setprecision(17);
    for(bool noisy:{false,true}) {
        std::ofstream out(root+(noisy?"/synthetic_ppp_acceleration_noisy.obs":"/synthetic_ppp_acceleration.obs"));header(out,"     3.04           O                   G","RINEX VERSION / TYPE");
        std::ostringstream pos;pos << std::fixed << std::setprecision(4);for(int j=0;j<3;++j)pos << std::setw(14) << seed[j];header(out,pos.str(),"APPROX POSITION XYZ");header(out,"G    4 C1C L1C C2W L2W","SYS / # / OBS TYPES");header(out,"","END OF HEADER");
        for(int i=0;i<64;++i) {
            const double t=.2*i,time=346200+t,clock=2e-6+i*5e-10;
            const Vector3d rx=origin+enu2ecef(Vector3d(20*t+.5*t*t,-8*t+4*std::sin(t/3),0),35*rad,139*rad);
            const Vector3d v=enu2ecef(Vector3d(20+t,-8+4.0/3*std::cos(t/3),0),35*rad,139*rad),a=enu2ecef(Vector3d(1,-4.0/9*std::sin(t/3),0),35*rad,139*rad);
            if(!noisy) {truth << i << ",2300," << time;for(const auto& x:{rx,v,a})for(int j=0;j<3;++j)truth << ',' << x[j];truth << ',' << clock << ",2.4\n";}
            out << "> " << fractional_epoch(600+t) << "  0" << std::setw(3) << satellites.size() << '\n';double lat,lon,height;ecef2geodetic(rx,lat,lon,height);
            for(const auto& e:satellites) {
                const int prn=e.satellite.prn;const auto tx=transmit(products,e.satellite,GNSSTime(2300,time),rx,false);const auto enu=ecef2enu(tx.p-rx,lat,lon);const double elev=std::atan2(enu[2],std::hypot(enu[0],enu[1]));
                const double common=geodist(tx.p,rx)+clock*constants::SPEED_OF_LIGHT-tx.total*constants::SPEED_OF_LIGHT+models::niellHydrostaticMapping(lat,height,elev,39)*2.4;
                const double ion1=4+.01*prn+.002*t,ion2=ion1*constants::GPS_L1_FREQ*constants::GPS_L1_FREQ/(constants::GPS_L2_FREQ*constants::GPS_L2_FREQ);
                out << e.satellite.toString();if(!noisy)unrounded << i << ',' << e.satellite.toString();
                for(int j=0;j<2;++j) {
                    const double lambda=j==0?constants::GPS_L1_WAVELENGTH:constants::GPS_L2_WAVELENGTH,ion=j==0?ion1:ion2;
                    const double code_bias=j==0?3.2*(prn%3)-2:-1.1*(prn%4)+.7,phase_bias=j==0?.015*(prn%3)-.007:.012*(prn%4)-.020;
                    const double code=common+ion+code_bias+(noisy?.05*std::sin(.37*i+.53*prn+j):0),phase=(common-ion+phase_bias+(noisy?.001*std::sin(.29*i+.31*prn+j):0))/lambda+(j==0?12+prn:-7+2*prn);
                    if(!noisy)unrounded << ',' << code << ',' << phase;
                    out << std::fixed << std::setprecision(5) << std::setw(14) << code << "  " << std::setprecision(4) << std::setw(14) << phase << "07";
                }
                if(!noisy) {unrounded << '\n';}
                out << '\n';
            }
        }
    }
}
int main(int argc,char**argv) {
    if(argc!=4)return 2;
    try {prediction(argv[3]);initialization(argv[2],argv[3]);car_physics(argv[1],argv[3]);return 0;}
    catch(const std::exception& e) {std::cerr << e.what() << '\n';return 1;}
}
