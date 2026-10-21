// Actual pinned native PVA restart, seed-disagreement and velocity-update calls.
// Private access is confined to this adapter; native source is unchanged.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <bits/stdc++.h>
#include <Eigen/Dense>
#define private public
#include <libgnss++/algorithms/ppp.hpp>
#undef private
#include "ppp_internal.hpp"
using namespace libgnss;

static void vector_values(std::ostream& out,const Eigen::VectorXd& x,const std::vector<int>& ix) {
    for(int i:ix)out << ',' << x[i];
}
static void matrix_values(std::ostream& out,const Eigen::MatrixXd& p,const std::vector<int>& ix) {
    for(int i:ix)for(int j:ix)out << ',' << p(i,j);
}
static Eigen::MatrixXd dense(int n,int id) {
    Eigen::MatrixXd b(n,n);for(int i=0;i<n;++i)for(int j=0;j<n;++j)b(i,j)=.1*std::sin(.31*i+.53*j+id);
    return b*b.transpose()+.2*Eigen::MatrixXd::Identity(n,n);
}
int main(int argc,char**argv) {
    if(argc!=2)return 2;
    const std::string root=argv[1];
    std::ofstream reset(root+"/upstream_ppp_motion_reset.csv"),velocity(root+"/upstream_ppp_motion_velocity.csv"),disagreement(root+"/upstream_ppp_motion_disagreement.csv");
    reset << "# id,n,position_var,gate,seed_valid,seed_has_velocity,seed_cov_present,seedXYZ,seedV,seedCov9,x,P,native_x,native_P\n" << std::setprecision(17);
    velocity << "# id,n,gate,measurementV,measurementCov9,x,P,native_applied,native_x,native_P\n" << std::setprecision(17);
    disagreement << "# id,predictedXYZ,seedXYZ,P9,C9,native_nis\n" << std::setprecision(17);
    const Vector3d origin(-3950000,3430000,3630000);
    for(int id=0;id<32;++id) {
        PPPProcessor::PPPConfig c;c.kinematic_mode=true;c.low_dynamics_mode=false;c.use_dynamics_model=true;c.use_pva_dynamics=true;c.estimate_ionosphere=false;c.pva_seed_disagreement_gate=id%7?16.27:0;
        PPPProcessor p(c);p.precise_products_loaded_=true;ObservationData obs;obs.clear();obs.time=GNSSTime(2300,346200);
        PositionSolution seed;seed.time=obs.time;seed.status=SolutionStatus::SPP;seed.num_satellites=6;seed.position_ecef=origin;seed.receiver_clock_bias=0;
        NavigationData nav;if(!p.initializeFilter(obs,nav,&seed))return 3;p.filter_initialized_=true;auto& f=p.filter_state_;const int n=14+id%3;
        f.total_states=n;f.state=Eigen::VectorXd::Zero(n);for(int i=0;i<n;++i)f.state[i]=.23*(i+1)+.1*id;f.state.head<3>()=origin;f.covariance=dense(n,id);
        const std::vector<double> position_vars={100,3600,std::nextafter(3600.0,std::numeric_limits<double>::infinity()),10000};
        const double pv=position_vars[id%position_vars.size()];for(int i=0;i<3;++i) {f.covariance(i,i)=pv;for(int j=0;j<3;++j)if(i!=j)f.covariance(i,j)=0;}
        const std::vector<double> deltas={0,5,50,500};seed.position_ecef+=Vector3d(deltas[(id/4)%4],0,0);seed.velocity_ecef=Vector3d(1+.1*id,-.5,.25);seed.has_velocity=id%5!=0;seed.position_covariance=Matrix3d::Identity()*(id%3?100:1);
        const bool cov_present=id%6!=0;if(!cov_present||id%11==0)seed.position_covariance(0,0)=std::numeric_limits<double>::quiet_NaN();const bool valid=id%9!=0;if(!valid)seed.status=SolutionStatus::NONE;
        std::vector<int> ix={0,1,2,9,11,3,4,5,6,7,8};for(int i=12;i<n;++i)ix.push_back(i);
        reset << id << ',' << ix.size() << ',' << c.initial_position_variance << ',' << c.pva_seed_disagreement_gate << ',' << valid << ',' << seed.has_velocity << ',' << cov_present;
        for(int i=0;i<3;++i) {reset << ',' << seed.position_ecef[i];}
        for(int i=0;i<3;++i) {reset << ',' << seed.velocity_ecef[i];}
        for(int i=0;i<3;++i)for(int j=0;j<3;++j) {reset << ',' << seed.position_covariance(i,j);}
        vector_values(reset,f.state,ix);matrix_values(reset,f.covariance,ix);p.predictState(0,&seed);vector_values(reset,f.state,ix);matrix_values(reset,f.covariance,ix);reset << '\n';
        // Actual helper uses a chosen velocity block; use the Rust head layout
        // directly here so every cross-covariance slot is compared.
        ppp_shared::PPPState state;state.vel_index=5;const int nv=13+id%3;state.state=Eigen::VectorXd::Zero(nv);for(int i=0;i<nv;++i)state.state[i]=.4*i+.07*id;state.covariance=dense(nv,id+40);
        Matrix3d R=dense(3,id+80);Vector3d measurement=state.state.segment<3>(5)+Vector3d(id%3?.2:30,-.3,.1);const double gate=id%4?16.27:0;
        std::vector<int> vi;for(int i=0;i<nv;++i)vi.push_back(i);
        velocity << id << ',' << nv << ',' << gate;for(int i=0;i<3;++i)velocity << ',' << measurement[i];for(int i=0;i<3;++i)for(int j=0;j<3;++j)velocity << ',' << R(i,j);vector_values(velocity,state.state,vi);matrix_values(velocity,state.covariance,vi);
        const bool applied=ppp_internal::applyVelocityMeasurement(state,measurement,R,gate);velocity << ',' << applied;vector_values(velocity,state.state,vi);matrix_values(velocity,state.covariance,vi);velocity << '\n';
        Vector3d predicted=origin+Vector3d(id,2*id,-id),other=origin;Matrix3d P=dense(3,id+120),C=dense(3,id+160);
        if(id==28) {predicted[0]=std::numeric_limits<double>::quiet_NaN();}
        if(id==29) {C(1,1)=std::numeric_limits<double>::infinity();}
        if(id==30) {P(0,0)=-10000;}
        if(id==31) {P.setZero();C.setZero();}
        disagreement << id;for(int i=0;i<3;++i)disagreement << ',' << predicted[i];for(int i=0;i<3;++i)disagreement << ',' << other[i];for(int i=0;i<3;++i)for(int j=0;j<3;++j)disagreement << ',' << P(i,j);for(int i=0;i<3;++i)for(int j=0;j<3;++j)disagreement << ',' << C(i,j);disagreement << ',' << ppp_internal::seedPositionDisagreement(predicted,P,other,C) << '\n';
    }
    return reset&&velocity&&disagreement?0:4;
}
