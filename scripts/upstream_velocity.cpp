// Independent cases from the pinned, unmodified spp_velocity::solveVelocity.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/algorithms/spp_velocity.hpp>
#include <libgnss++/core/constants.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
using namespace libgnss;

static std::vector<spp_velocity::DopplerObservation> observations(
    const Vector3d& receiver,const Vector3d& velocity,double drift,int n,double noise,bool poor=false) {
    std::vector<spp_velocity::DopplerObservation> rows;
    for(int i=0;i<n;++i) {
        const double az=poor?0.01*i:0.93*i;
        const double el=poor?0.5+0.0001*i:0.2+0.1*(i%6);
        const Vector3d los(std::cos(el)*std::cos(az),std::cos(el)*std::sin(az),std::sin(el));
        spp_velocity::DopplerObservation d;
        d.satellite_position_ecef=receiver+los*(2.1e7+i*1e5);
        d.satellite_velocity_ecef=Vector3d(1800*std::sin(i+0.1),-2300*std::cos(i+0.2),1200*std::sin(i*0.5));
        d.satellite_clock_drift=(i-3)*1e-11;
        d.frequency_hz=constants::GPS_L1_FREQ;
        const double rate=(d.satellite_velocity_ecef-velocity).dot(los)
            +constants::SPEED_OF_LIGHT*(drift-d.satellite_clock_drift)+noise*std::sin(i*1.7+0.2);
        d.doppler_hz=-rate/constants::GPS_L1_WAVELENGTH;d.elevation_rad=el;
        rows.push_back(d);
    }
    return rows;
}
static void write(std::ostream& out,const std::string& name,
    const std::vector<spp_velocity::DopplerObservation>& rows,const Vector3d& receiver,double sigma=0.5,int minimum=4) {
    out<<name<<','<<rows.size();
    for(int j=0;j<3;++j) out<<','<<receiver(j);
    out<<','<<sigma<<','<<minimum;
    for(const auto& d:rows) {
        for(int j=0;j<3;++j) out<<','<<d.satellite_position_ecef(j);
        for(int j=0;j<3;++j) out<<','<<d.satellite_velocity_ecef(j);
        out<<','<<d.satellite_clock_drift<<','<<d.doppler_hz<<','<<d.frequency_hz<<','<<d.elevation_rad;
    }
    const auto r=spp_velocity::solveVelocity(rows,receiver,sigma,minimum);
    out<<','<<r.ok<<','<<r.num_satellites_used;
    for(int j=0;j<3;++j) out<<','<<r.velocity_ecef(j);
    out<<','<<r.receiver_clock_drift;
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) out<<','<<r.velocity_covariance(i,j);
    out<<','<<r.clock_drift_variance<<','<<r.residual_rms_mps<<'\n';
}
int main(int argc,char** argv) {
    if(argc!=2)return 2;
    std::ofstream out(argv[1]);out<<std::setprecision(17);
    out<<"# case,n,receiver[3],sigma,minimum,observations[n*(position3+velocity3+clock+doppler+frequency+elevation)],native_ok,used,velocity3,clock,covariance9,clock_variance,rms\n";
    const Vector3d r(-3817680.73,3562839.52,3650159.24),v(-3.2,4.5,1.2);
    write(out,"stationary",observations(r,Vector3d::Zero(),2e-10,8,0),r);
    write(out,"moving",observations(r,v,-1e-8,8,0),r);
    for(double noise:{0.1,0.5,1.0,2.0}) write(out,"noise_"+std::to_string(noise),observations(r,v,1e-9,9,noise),r);
    write(out,"four",observations(r,v,2e-9,4,0),r);
    write(out,"insufficient",observations(r,v,0,3,0),r);
    write(out,"poor_geometry",observations(r,v,0,4,0,true),r);
    auto bad=observations(r,v,0,8,0);bad[3].doppler_hz+=10000;
    write(out,"bad_fit",bad,r);
    auto singular=observations(r,v,0,6,0);
    for(auto& d:singular)d.satellite_position_ecef=singular[0].satellite_position_ecef;
    write(out,"singular",singular,r);
    auto zero=observations(r,Vector3d::Zero(),0,8,0);
    for(auto& d:zero) {d.doppler_hz=0;d.satellite_velocity_ecef.setZero();d.satellite_clock_drift=0;}
    write(out,"native_zero_sentinel",zero,r);
    return out?0:1;
}
