// Direct calls to pinned PPPProcessor::initializeAmbiguityState (precise IF).
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <bits/stdc++.h>
#include <Eigen/Dense>
#define private public
#include <libgnss++/algorithms/ppp.hpp>
#undef private
#include <libgnss++/core/coordinates.hpp>
using namespace libgnss;
int main(int argc,char**argv) {
    if(argc!=2)return 2;
    std::ofstream out(std::string(argv[1])+"/upstream_ppp_if_ambiguity_prior.csv");out<<std::setprecision(17);
    out<<"id,corrected_phase_if_m,predicted_range_m,initialized_ambiguity_m,variance_m2\n";
    for(int id=0;id<32;++id){
        PPPProcessor p;p.precise_products_loaded_=true;p.ppp_config_.use_ionosphere_free=true;p.ppp_config_.estimate_ionosphere=false;p.ppp_config_.estimate_troposphere=id%2==0;p.ppp_config_.initial_ambiguity_variance=3600+id;
        auto& f=p.filter_state_;f.pos_index=0;f.vel_index=3;f.accel_index=6;f.clock_index=9;f.glo_clock_index=10;f.trop_index=11;f.total_states=13;f.state=Eigen::VectorXd::Zero(13);f.covariance=Eigen::MatrixXd::Identity(13,13);
        const Vector3d receiver(-3947481.0,3431493.0,3637893.0);f.state.segment<3>(0)=receiver;f.state(9)=600+.3*id;f.state(11)=2.3+.01*id;
        PPPProcessor::IonosphereFreeObs o{};o.satellite=SatelliteId(GNSSSystem::GPS,1+id);o.valid=true;o.has_carrier_phase=true;o.receiver_position=id%3?receiver:Vector3d::Zero();o.satellite_position=receiver+Vector3d(1.7e7+1000*id,8e6-2000*id,1.2e7+3000*id);o.satellite_clock_bias=1e-5+.1e-6*id;o.trop_mapping=1.1+.1*id;o.modeled_trop_delay_m=2.5*o.trop_mapping;o.ambiguity_scale_m=.107;
        const double prediction=geodist(o.satellite_position,receiver)+f.state(9)-constants::SPEED_OF_LIGHT*o.satellite_clock_bias+(p.ppp_config_.estimate_troposphere?o.trop_mapping*f.state(11):o.modeled_trop_delay_m);
        o.carrier_phase_if=prediction+(id%5==0?0.0:(id-16)*.713);
        p.initializeAmbiguityState(o,12);
        out<<id<<','<<o.carrier_phase_if<<','<<prediction<<','<<f.state(12)<<','<<f.covariance(12,12)<<'\n';
    }
    return out?0:4;
}
