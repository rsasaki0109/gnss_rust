// Finite-prior startup reference: native decoded observations/orbits/DD R/Kalman.
// This composes components; it is not the whole native RTKProcessor policy.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/io/rinex.hpp>
#include <libgnss++/algorithms/kalman.hpp>
#include <libgnss++/algorithms/rtk_measurement.hpp>
#include <libgnss++/core/coordinates.hpp>
#include <libgnss++/models/troposphere.hpp>
#include "navigation_internal.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace libgnss;
static double elevation(const Vector3d&p,const Vector3d&r) {
    double lat,lon,height;ecef2geodetic(r,lat,lon,height);
    const auto enu=ecef2enu(p-r,lat,lon);return std::atan2(enu(2),std::hypot(enu(0),enu(1)));
}
static Vector3d satellite_position(const Observation&o,GNSSTime t,const NavigationData&nav) {
    auto tx=t-o.pseudorange/constants::SPEED_OF_LIGHT;const auto*e=nav.getEphemeris(o.satellite,tx);
    if(!e)throw std::runtime_error("missing ephemeris");Vector3d p;double clock,drift;
    if(!navigation_internal::computeBroadcastState(*e,tx,p,clock,drift,true))throw std::runtime_error("state");
    tx=tx-clock;
    if(!navigation_internal::computeBroadcastState(*e,tx,p,clock,drift,true))throw std::runtime_error("state");
    return p;
}
struct Sample {SatelliteId id;Vector3d rp,bp;double code,phase;};

int main(int argc,char**argv) {
    if(argc!=2)return 2;const std::string root=argv[1];
    std::cout<<std::setprecision(17)<<"# case,epoch,posterior_xyz[3],position_covariance[9]; prior position=900m2, SD ambiguity=900cycles2\n";
    for(int scenario=0;scenario<3;++scenario) {
        const int epoch=scenario==2?1:0;
        io::RINEXReader rr,br,nr;rr.open(root+"/synthetic_rtk_mixed_rover.obs");br.open(root+"/synthetic_rtk_mixed_base.obs");nr.open(root+"/synthetic_rtk_mixed.nav");
        io::RINEXReader::RINEXHeader rh,bh;rr.readHeader(rh);br.readHeader(bh);NavigationData nav;nr.readNavigationData(nav);
        ObservationData rover,base;for(int i=0;i<=epoch;++i){if(!rr.readObservationEpoch(rover)||!br.readObservationEpoch(base))throw std::runtime_error("epoch");}
        const Vector3d seed=rh.approximate_position,base_position=bh.approximate_position;
        std::vector<Sample> data;std::map<GNSSSystem,int> counts;
        for(const auto&[key,ro]:rover.rinex_tracking_observations) {
            const bool gal=key.first.system==GNSSSystem::Galileo;
            if((!gal&&key.first.system!=GNSSSystem::GPS)||key.second!="1C"||(!gal&&scenario==2))continue;
            const auto found=base.rinex_tracking_observations.find({key.first,gal?"1X":"1C"});
            if(found==base.rinex_tracking_observations.end())throw std::runtime_error("pair");
            if(scenario==1 && counts[key.first.system]>=3)continue;
            ++counts[key.first.system];const auto&ba=found->second;
            data.push_back({key.first,satellite_position(ro,rover.time,nav),satellite_position(ba,base.time,nav),
                ro.pseudorange-ba.pseudorange,ro.carrier_phase-ba.carrier_phase});
        }
        std::map<GNSSSystem,int> references;
        for(int i=0;i<static_cast<int>(data.size());++i) {
            auto found=references.find(data[i].id.system);
            if(found==references.end()||elevation(data[i].rp,seed)>elevation(data[found->second].rp,seed))references[data[i].id.system]=i;
        }
        std::vector<std::pair<int,int>> pairs;std::vector<int> block_sizes;
        for(const auto&[sys,ref]:references) {
            const int start=static_cast<int>(pairs.size());
            for(int i=0;i<static_cast<int>(data.size());++i)if(i!=ref&&data[i].id.system==sys)pairs.emplace_back(ref,i);
            block_sizes.push_back(static_cast<int>(pairs.size())-start);
        }
        const auto blocks=block_sizes;block_sizes.insert(block_sizes.end(),blocks.begin(),blocks.end());
        const int n=3+static_cast<int>(data.size()),m=static_cast<int>(pairs.size());
        VectorXd prior=VectorXd::Zero(n);prior.head<3>()=seed-base_position;
        for(int i=0;i<static_cast<int>(data.size());++i)prior(3+i)=data[i].phase-data[i].code/constants::GPS_L1_WAVELENGTH;
        const MatrixXd p0=MatrixXd::Identity(n,n)*900;
        VectorXd estimate=prior;MatrixXd posterior;
        for(int iteration=0;iteration<8;++iteration) {
            const Vector3d position=base_position+estimate.head<3>();
            std::vector<Vector3d>los;std::vector<double>geom,el;
            for(const auto&d:data) {
                const double re=elevation(d.rp,position),be=elevation(d.bp,base_position);
                el.push_back(re);los.push_back((d.rp-position).normalized());
                geom.push_back(geodist(d.rp,position)-geodist(d.bp,base_position)
                    +models::tropDelaySaastamoinen(position,re)-models::tropDelaySaastamoinen(base_position,be));
            }
            MatrixXd h=MatrixXd::Zero(2*m,n);VectorXd v=VectorXd::Zero(2*m);
            std::vector<double>rv(2*m),sv(2*m);
            for(int row=0;row<m;++row) {
                const auto[ref,sat]=pairs[row];const double ddgeom=geom[ref]-geom[sat];
                h.block<1,3>(row,0)=(-los[ref]+los[sat]).transpose();h.block<1,3>(m+row,0)=h.block<1,3>(row,0);
                h(row,3+ref)=constants::GPS_L1_WAVELENGTH;h(row,3+sat)=-constants::GPS_L1_WAVELENGTH;
                v(row)=(data[ref].phase-data[sat].phase)*constants::GPS_L1_WAVELENGTH-ddgeom
                    -constants::GPS_L1_WAVELENGTH*(estimate(3+ref)-estimate(3+sat));
                v(m+row)=data[ref].code-data[sat].code-ddgeom;
                for(int kind=0;kind<2;++kind) {
                    const double sigma=kind==0?0.002:0.3;
                    auto variance=[&](int i){const double s=std::max(std::sin(el[i]),0.1);return 2*sigma*sigma*(1+1/(s*s));};
                    rv[kind*m+row]=variance(ref);sv[kind*m+row]=variance(sat);
                }
            }
            const MatrixXd r=rtk_measurement::buildDoubleDifferenceCovariance(block_sizes,rv,sv,2*m);
            const VectorXd innovation=v+h*(estimate-prior);
            VectorXd updated=prior;posterior=p0;
            if(kalmanFilter(updated,posterior,h,innovation,r,std::vector<bool>(n,true))!=0)throw std::runtime_error("Kalman");
            const double change=(updated.head<3>()-estimate.head<3>()).norm();estimate=updated;
            if(change<=1e-4)break;
            if(iteration==7)throw std::runtime_error("convergence");
        }
        const Vector3d result=base_position+estimate.head<3>();
        std::cout<<(scenario==0?"full":scenario==1?"three_each":"galileo_only")<<','<<epoch;
        for(int i=0;i<3;++i)std::cout<<','<<result(i);
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)std::cout<<','<<posterior(i,j);
        std::cout<<'\n';
    }
}
