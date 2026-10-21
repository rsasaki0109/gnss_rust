// Actual native PPP predictState component plus independent moving IF physics.
// Private access is exposed only in this adapter; upstream source is unmodified.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <bits/stdc++.h>
#include <Eigen/Dense>
#define private public
#include <libgnss++/algorithms/ppp.hpp>
#undef private
#define GNSS_RUST_PRECISE_TRANSMIT_MAIN kinematic_unused_transmit_main
#include "upstream_precise_transmit.cpp"
#include <libgnss++/models/troposphere.hpp>

static void emit_vector(std::ostream& out,const Eigen::VectorXd& x,const std::vector<int>& indices) {
    for(int i:indices)out << ',' << x[i];
}
static void emit_matrix(std::ostream& out,const Eigen::MatrixXd& p,const std::vector<int>& indices) {
    for(int i:indices)for(int j:indices)out << ',' << p(i,j);
}
static void predict_cases(const std::string& root) {
    std::ofstream out(root+"/upstream_ppp_kinematic_predict.csv");
    out << "# id,n,dt,qpos,qtrop,qamb,qvelocity,x,P,native_predicted_x,native_predicted_P; row-major, Rust layout pos/clock/trop/velocity/ambiguities\n" << std::setprecision(17);
    for(int id=0;id<24;++id) {
        const int na=2+id%4,n=9+na;
        PPPProcessor::PPPConfig config;
        config.kinematic_mode=true;config.low_dynamics_mode=false;
        config.use_pva_dynamics=false;config.estimate_ionosphere=false;
        config.process_noise_position=.003*(id%3);
        config.process_noise_troposphere=1e-8*(1+id%3);
        config.process_noise_ambiguity=1e-8*(1+id%4);
        config.process_noise_velocity=1e-4*(1+id%3);
        PPPProcessor processor(config);
        processor.filter_initialized_=true;processor.precise_products_loaded_=true;
        auto& state=processor.filter_state_;state.total_states=n;
        state.state=Eigen::VectorXd::Zero(n);Eigen::MatrixXd matrix(n,n);
        for(int i=0;i<n;++i) {state.state[i]=.23*(i+1)+.01*id;
            for(int j=0;j<n;++j)matrix(i,j)=.1*std::sin(.7*(i+1)+1.3*(j+1)+id);}
        state.covariance=matrix*matrix.transpose()+.2*Eigen::MatrixXd::Identity(n,n);
        const std::vector<double> gaps={0.0,.2,1.0,30.0,120.0};
        const double dt=gaps[id%gaps.size()];
        std::vector<int> indices={0,1,2,6,8,3,4,5};
        for(int i=9;i<n;++i)indices.push_back(i);
        out << id << ',' << indices.size() << ',' << dt << ',' << config.process_noise_position << ','
            << config.process_noise_troposphere << ',' << config.process_noise_ambiguity << ',' << config.process_noise_velocity;
        emit_vector(out,state.state,indices);emit_matrix(out,state.covariance,indices);
        processor.predictState(dt,nullptr);
        emit_vector(out,state.state,indices);emit_matrix(out,state.covariance,indices);out << '\n';
    }
}
int main(int argc,char**argv) {
    if(argc!=3)return 2;
    const std::string inputs=argv[1],root=argv[2];
    predict_cases(root);
    PreciseProducts products;
    if(!products.loadSP3File(inputs+"/synthetic_ppp.sp3") || !products.loadClockFile(inputs+"/synthetic_ppp.clk"))return 3;
    const double rad=std::acos(-1.0)/180.0;
    const Vector3d origin=geodetic2ecef(35*rad,139*rad,45.0),seed=origin+Vector3d(5,-3,2);
    const auto satellites=visible(origin);
    std::ofstream truth(root+"/synthetic_ppp_kinematic_truth.csv");
    std::ofstream unrounded(root+"/upstream_ppp_kinematic_unrounded.csv");
    unrounded << "epoch,satellite,code1_m,phase1_cycles,code2_m,phase2_cycles\n" << std::setprecision(17);
    truth << "epoch,week,tow,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,clock_s,zenith_m\n" << std::setprecision(17);
    for(bool noisy:{false,true}) {
        std::ofstream obs(root+(noisy?"/synthetic_ppp_kinematic_noisy.obs":"/synthetic_ppp_kinematic.obs"));
        header(obs,"     3.04           O                   G","RINEX VERSION / TYPE");
        std::ostringstream pos;pos << std::fixed << std::setprecision(4);
        for(int j=0;j<3;++j)pos << std::setw(14) << seed[j];
        header(obs,pos.str(),"APPROX POSITION XYZ");header(obs,"G    4 C1C L1C C2W L2W","SYS / # / OBS TYPES");header(obs,"","END OF HEADER");
        for(int i=0;i<64;++i) {
            const int seconds=600+30*i;const double t=30.0*i;
            const Vector3d rx=origin+Vector3d(.12*t+2*std::sin(t/100.0),-.08*t+std::sin(t/130.0),.03*t);
            const Vector3d velocity(.12+.02*std::cos(t/100.0),-.08+std::cos(t/130.0)/130.0,.03);
            const double receiver_clock=2e-6+i*5e-10;
            double latitude,longitude,height;ecef2geodetic(rx,latitude,longitude,height);
            if(!noisy) {truth << i << ",2300," << 345600+seconds;
                for(int j=0;j<3;++j) {truth << ',' << rx[j];}
                for(int j=0;j<3;++j) {truth << ',' << velocity[j];}
                truth << ',' << receiver_clock << ",2.4\n";}
            obs << "> " << epoch(seconds) << "  0" << std::setw(3) << satellites.size() << '\n';
            for(const auto& e:satellites) {
                const int prn=e.satellite.prn;const auto tx=transmit(products,e.satellite,GNSSTime(2300,345600.0+seconds),rx,false);
                const auto enu=ecef2enu(tx.p-rx,latitude,longitude);
                const double elevation=std::atan2(enu[2],std::hypot(enu[0],enu[1]));
                const double common=geodist(tx.p,rx)+receiver_clock*constants::SPEED_OF_LIGHT-tx.total*constants::SPEED_OF_LIGHT
                    +models::niellHydrostaticMapping(latitude,height,elevation,39)*2.4;
                const double ion1=4+.01*prn+.002*t,ion2=ion1*constants::GPS_L1_FREQ*constants::GPS_L1_FREQ/(constants::GPS_L2_FREQ*constants::GPS_L2_FREQ);
                obs << e.satellite.toString();
                if(!noisy) {unrounded << i << ',' << e.satellite.toString();}
                for(int j=0;j<2;++j) {
                    const double lambda=j==0?constants::GPS_L1_WAVELENGTH:constants::GPS_L2_WAVELENGTH;
                    const double ion=j==0?ion1:ion2,ambiguity=j==0?12+prn:-7+2*prn;
                    const double code_bias=j==0?3.2*(prn%3)-2.0:-1.1*(prn%4)+.7;
                    const double phase_bias=j==0?.015*(prn%3)-.007:.012*(prn%4)-.020;
                    const double code_noise=noisy?.05*std::sin(.37*i+.53*prn+j):0.0;
                    const double phase_noise=noisy?.001*std::sin(.29*i+.31*prn+j):0.0;
                    const double code=common+ion+code_bias+code_noise;
                    const double phase=(common-ion+phase_bias+phase_noise)/lambda+ambiguity;
                    if(!noisy) {unrounded << ',' << code << ',' << phase;}
                    obs << std::fixed << std::setprecision(5) << std::setw(14) << code << "  "
                        << std::setprecision(4) << std::setw(14) << phase << "07";
                }
                if(!noisy) {unrounded << '\n';}
                obs << '\n';
            }
        }
    }
    return truth?0:4;
}
