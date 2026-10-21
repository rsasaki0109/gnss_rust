// Actual pinned native PPP initializer with independent RINEX header seeds.
// Component diagnostic only; does not execute the whole native PPP solver.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <bits/stdc++.h>
#include <Eigen/Dense>
#define private public
#include <libgnss++/algorithms/ppp.hpp>
#undef private
#include <libgnss++/io/rinex.hpp>
using namespace libgnss;
int main(int argc,char**argv) {
    if(argc!=2)return 2;
    std::cout << "case,week,tow,x_m,y_m,z_m,clock_m,zenith_m,pos_qxx_m2,pos_qyy_m2,pos_qzz_m2,velocity_qxx_m2_s2,clock_q_m2,zenith_q_m2\n" << std::setprecision(17);
    for(const std::string city:{"nagoya","tokyo"})for(const std::string run:{"run1","run2","run3"}) {
        io::RINEXReader reader;
        if(!reader.open(std::string(argv[1])+"/"+city+"/"+run+"/rover.obs"))return 3;
        io::RINEXReader::RINEXHeader header;
        ObservationData observation;
        if(!reader.readHeader(header)||!header.has_approximate_position||!reader.readObservationEpoch(observation))return 4;
        PositionSolution seed;seed.time=observation.time;seed.status=SolutionStatus::SPP;
        seed.position_ecef=header.approximate_position;seed.receiver_clock_bias=0.0;seed.num_satellites=6;
        PPPProcessor::PPPConfig config;config.kinematic_mode=true;config.low_dynamics_mode=false;
        config.use_pva_dynamics=false;config.estimate_ionosphere=false;config.estimate_troposphere=true;
        PPPProcessor processor(config);processor.precise_products_loaded_=true;
        NavigationData navigation;
        if(!processor.initializeFilter(observation,navigation,&seed))return 5;
        const auto& f=processor.filter_state_;
        std::cout << city << '-' << run << ',' << observation.time.week << ',' << observation.time.tow;
        for(int j=0;j<3;++j)std::cout << ',' << f.state[f.pos_index+j];
        std::cout << ',' << f.state[f.clock_index] << ',' << f.state[f.trop_index];
        for(int j=0;j<3;++j)std::cout << ',' << f.covariance(f.pos_index+j,f.pos_index+j);
        std::cout << ',' << f.covariance(f.vel_index,f.vel_index) << ',' << f.covariance(f.clock_index,f.clock_index)
                  << ',' << f.covariance(f.trop_index,f.trop_index) << '\n';
    }
    return 0;
}
