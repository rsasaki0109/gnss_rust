// Independent Eigen reference for the multi-clock extension of the pinned
// spp_velocity model. This is not an upstream multi-clock solver entry point.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#define main single_clock_generator_main
#include "upstream_velocity.cpp"
#undef main

static void multi_case(std::ostream& out,const std::string& name,int count,int clocks,
    const Vector3d& receiver,double noise,bool singular=false,bool bad=false) {
    const Vector3d velocity(-3.2,4.5,1.2);
    auto data=observations(receiver,velocity,0,count,noise);
    const double drifts[]={0.4,-2.0,4.7};
    const char prefixes[]={'G','E','J'};
    const int n=3+clocks;
    Eigen::MatrixXd h=Eigen::MatrixXd::Zero(count,n),weighted(count,n);
    Eigen::VectorXd y(count),wy(count),weights(count);
    for(int i=0;i<count;++i) {
        const int group=i%clocks;
        auto& d=data[i];
        d.doppler_hz-=drifts[group]/constants::GPS_L1_WAVELENGTH;
        if(singular)d.satellite_position_ecef=data[0].satellite_position_ecef;
        if(bad && i==3)d.doppler_hz+=10000;
        const Vector3d los=(d.satellite_position_ecef-receiver).normalized();
        h.block<1,3>(i,0)=-los.transpose();h(i,3+group)=1;
        y(i)=-d.doppler_hz*constants::GPS_L1_WAVELENGTH-d.satellite_velocity_ecef.dot(los)+constants::SPEED_OF_LIGHT*d.satellite_clock_drift;
        weights(i)=std::max(0.1,std::sin(d.elevation_rad))/0.5;
        weighted.row(i)=h.row(i)*weights(i);wy(i)=y(i)*weights(i);
    }
    out<<name<<','<<count<<','<<clocks;
    for(int j=0;j<3;++j)out<<','<<receiver(j);
    for(int i=0;i<count;++i) {
        const auto& d=data[i];
        out<<','<<prefixes[i%clocks]<<std::setfill('0')<<std::setw(2)<<(i/clocks+1)<<std::setfill(' ');
        for(int j=0;j<3;++j)out<<','<<d.satellite_position_ecef(j);
        for(int j=0;j<3;++j)out<<','<<d.satellite_velocity_ecef(j);
        out<<','<<d.satellite_clock_drift<<','<<d.doppler_hz<<','<<d.frequency_hz<<','<<d.elevation_rad;
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(weighted);
    if(count<n || qr.rank()<n) {out<<",0\n";return;}
    const Eigen::VectorXd x=qr.solve(wy),r=y-h*x;
    const double chi=(r.array()*weights.array()).square().sum()/std::max(1,count-n);
    const Eigen::MatrixXd q=(weighted.transpose()*weighted).ldlt().solve(Eigen::MatrixXd::Identity(n,n));
    if(!q.allFinite() || !x.allFinite() || chi>25 || q.diagonal().head<3>().maxCoeff()>400) {out<<",0\n";return;}
    const double scale=std::max(1.0,chi);
    out<<",1";for(int j=0;j<n;++j)out<<','<<x(j);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)out<<','<<q(i,j)*scale;
    out<<','<<chi<<','<<scale<<','<<std::sqrt(r.squaredNorm()/count)<<'\n';
}
int main(int argc,char**argv) {
    if(argc!=2)return 2;
    std::ofstream out(argv[1]);out<<std::setprecision(17);
    out<<"# case,n,clocks,receiver3,observations[n*(satellite+position3+velocity3+clock+doppler+frequency+elevation)],ok,[mean[3+clocks],covariance[(3+clocks)^2],chi,scale,rms if ok]\n";
    const std::vector<Vector3d> sites={{-3817680.73,3562839.52,3650159.24},{6378137,0,0}};
    int i=0;
    for(const auto& receiver:sites)for(double noise:{0.0,0.1,0.5,1.0})multi_case(out,"two_"+std::to_string(i++),8,2,receiver,noise);
    for(double noise:{0.0,0.1,0.5})multi_case(out,"three_"+std::to_string(i++),9,3,sites[0],noise);
    multi_case(out,"three_each",6,2,sites[0],0);
    multi_case(out,"minimum",5,2,sites[0],0);
    multi_case(out,"too_few",4,2,sites[0],0);
    multi_case(out,"singular",7,2,sites[0],0,true);
    multi_case(out,"bad_fit",8,2,sites[0],0,false,true);
    return out?0:1;
}
