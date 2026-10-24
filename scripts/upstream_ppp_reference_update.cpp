// Native measurement assembly and residual screen; unchanged pinned library.
// Joseph update expressions copied from ppp_state.cpp::updateFilter.
// This is a component harness, not a call to the full updateFilter/processEpoch.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <bits/stdc++.h>
#include <Eigen/Dense>
#define private public
#include <libgnss++/algorithms/ppp.hpp>
#undef private
#include "ppp_internal.hpp"
#include <libgnss++/core/coordinates.hpp>
using namespace libgnss;
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    std::ofstream out(std::string(argv[1]) + "/upstream_ppp_reference_update.csv");
    out << std::setprecision(17);
    out << "id,n,m,kinematic,converged,success,screen_count,nis,mean_P_sat_phase_H_r_R_prefitmask_keepmask_screens_posterior_P\n";
    for (int id = 0; id < 32; ++id) {
        const int satellites = 4 + id % 6, n = 8 + satellites;
        PPPProcessor p;
        p.precise_products_loaded_ = true;
        p.ppp_config_.use_ionosphere_free = true;
        p.ppp_config_.estimate_ionosphere = false;
        p.ppp_config_.estimate_troposphere = true;
        p.ppp_config_.kinematic_mode = id % 4 != 0;
        p.converged_ = id % 3 == 0;
        auto& f = p.filter_state_;
        f.pos_index = 0; f.clock_index = 3; f.trop_index = 4;
        f.vel_index = 5; f.glo_clock_index = 3; f.accel_index = -1;
        f.total_states = n; f.state = Eigen::VectorXd::Zero(n);
        f.state.head<3>() << -3947481.0,3431493.0,3637893.0;
        f.state(3) = 600 + .3 * id; f.state(4) = 2.3;
        Eigen::MatrixXd b(n,n);
        for (int i=0;i<n;++i) for (int j=0;j<n;++j) b(i,j)=.1*std::sin(.3*i+.7*j+id);
        f.covariance = b*b.transpose();
        for (int i=0;i<n;++i) f.covariance(i,i) += i==3?400:i==4?.36:i<8?100:25;
        std::vector<PPPProcessor::IonosphereFreeObs> observations;
        for (int j=0;j<satellites;++j) {
            const SatelliteId sat(GNSSSystem::GPS,j+1);
            f.ambiguity_indices[sat] = 8+j; f.state(8+j) = (j-4)*.713;
            auto& a = p.ambiguity_states_[sat];
            a.lock_count = id % 5 == 0 ? 0 : 2;
            a.needs_reinitialization = false;
            PPPProcessor::IonosphereFreeObs o{};
            o.satellite=sat; o.valid=true; o.has_carrier_phase=true;
            o.receiver_position=f.state.head<3>();
            Vector3d direction(std::cos(j*1.7),std::sin(j*1.7),.2+.1*j);
            o.satellite_position=o.receiver_position+direction.normalized()*2.1e7;
            o.satellite_clock_bias=1e-5; o.trop_mapping=1.1+.3*j;
            o.elevation=.3+.05*j; o.variance_pr=2+.2*j; o.variance_cp=.001+.0001*j;
            const double predicted=geodist(o.satellite_position,o.receiver_position)+f.state(3)-constants::SPEED_OF_LIGHT*o.satellite_clock_bias+o.trop_mapping*f.state(4);
            double cr=.1*std::cos(id+j), pr=.0001*std::sin(id+j);
            // Native prefit code/phase exclusions; startup and converged floors.
            if (id%7==1 && j==satellites-1) cr=100000;
            if (id%7==2 && j==satellites-1) pr=250;
            // Screen one or multiple satellites, including exhaustion at four.
            if (id%7==3 && j==satellites-1) cr=100;
            if (id%7==4 && j>=satellites-2) cr=80+10*j;
            if (id==30) cr=100000;
            o.pseudorange_if=predicted+cr;
            o.carrier_phase_if=predicted+f.state(8+j)+pr;
            observations.push_back(o);
        }
        const NavigationData nav; const GNSSTime time(2300,100.0);
        const auto all=p.formMeasurementEquations(observations,nav,time,false,false);
        const auto gated=p.formMeasurementEquations(observations,nav,time,true,false);
        const int m=all.residuals.size();
        std::vector<int> prefit(m,0), keep(m,0), screens;
        for (int i=0;i<m;++i) for (int j=0;j<gated.residuals.size();++j) {
            if (all.row_satellites[i]==gated.row_satellites[j] && all.row_is_phase[i]==gated.row_is_phase[j]) prefit[i]=1;
        }
        keep=prefit;
        Eigen::VectorXd posterior=f.state;
        Eigen::MatrixXd covariance=f.covariance;
        bool success=false; double nis=0;
        for (;;) {
            std::vector<int> indices;
            for (int i=0;i<m;++i) if (keep[i]) indices.push_back(i);
            const int rows=indices.size();
            if (rows<4) break;
            Eigen::MatrixXd h(rows,n),rnoise=Eigen::MatrixXd::Zero(rows,rows);
            Eigen::VectorXd r(rows);
            for (int i=0;i<rows;++i) {
                const int k=indices[i];h.row(i)=all.design_matrix.row(k);
                r(i)=all.residuals(k);rnoise(i,i)=all.weight_matrix(k,k);
            }
            const Eigen::MatrixXd s=h*f.covariance*h.transpose()+rnoise;
            const Eigen::MatrixXd inv=s.ldlt().solve(Eigen::MatrixXd::Identity(rows,rows));
            const int worst=p.ppp_config_.kinematic_mode?ppp_internal::worstStandardizedResidualRow(r,inv,4.0):-1;
            if (worst>=0) {
                const auto sat=all.row_satellites[indices[worst]];
                screens.push_back(sat.prn);
                std::set<SatelliteId> remaining;
                for (int k:indices) if (!(all.row_satellites[k]==sat)) remaining.insert(all.row_satellites[k]);
                if (remaining.size()<4) break;
                for (int i=0;i<m;++i) if (all.row_satellites[i]==sat) keep[i]=0;
                continue;
            }
            const Eigen::MatrixXd gain=f.covariance*h.transpose()*inv;
            posterior=f.state+gain*r;
            const Eigen::MatrixXd identity=Eigen::MatrixXd::Identity(n,n);
            const Eigen::MatrixXd ikh=identity-gain*h;
            covariance=ikh*f.covariance*ikh.transpose()+gain*rnoise*gain.transpose();
            nis=r.dot(inv*r);success=true;break;
        }
        out<<id<<','<<n<<','<<m<<','<<p.ppp_config_.kinematic_mode<<','<<p.converged_<<','<<success<<','<<screens.size()<<','<<nis;
        for (int i=0;i<n;++i) out<<','<<f.state(i);
        for (int i=0;i<n;++i) for (int j=0;j<n;++j) out<<','<<f.covariance(i,j);
        for (int i=0;i<m;++i) out<<','<<static_cast<int>(all.row_satellites[i].prn);
        for (int i=0;i<m;++i) out<<','<<all.row_is_phase[i];
        for (int i=0;i<m;++i) for (int j=0;j<n;++j) out<<','<<all.design_matrix(i,j);
        for (int i=0;i<m;++i) out<<','<<all.residuals(i);
        for (int i=0;i<m;++i) out<<','<<all.weight_matrix(i,i);
        for (int value:prefit) out<<','<<value;
        for (int value:keep) out<<','<<value;
        for (int i=0;i<m;++i) out<<','<<(i<static_cast<int>(screens.size())?screens[i]:0);
        for (int i=0;i<n;++i) out<<','<<posterior(i);
        for (int i=0;i<n;++i) for (int j=0;j<n;++j) out<<','<<covariance(i,j);
        out<<'\n';
    }
    return out?0:4;
}
