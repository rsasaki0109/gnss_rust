// Actual pinned SPP variance, robust mask and RAIM/FDE selection utilities.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/algorithms/spp.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
using namespace libgnss;
int main(int argc,char**argv){
    if(argc!=2)return 2;
    std::ofstream o(argv[1]);if(!o)return 3;
    o<<std::setprecision(17)<<"# tagged native SPP quality utility cases; V variance / M inlier mask / F FDE selection\n";
    const double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
    for(int id=0;id<32;++id){
        spp_utils::MeasurementVarianceInputs v;v.elevation_rad=(id%7)*.2;v.snr_dbhz=id%6==0?0.:15.+5*(id%8);v.pseudorange_sigma_m=id%5==0?.01:3.;
        v.ionosphere_corrected=id%2;v.troposphere_corrected=id%3;v.ionosphere_delay_m=id-8.;v.troposphere_delay_m=id*.7;
        if(id==28)v.snr_dbhz=nan;
        if(id==29)v.elevation_rad=nan;
        if(id==30){v.snr_reference_dbhz=nan;v.ionosphere_sigma_scale=-1.;v.troposphere_sigma_scale=nan;}
        if(id==31){v.elevation_rad=std::acos(-1.)/2;v.snr_dbhz=33.;v.ionosphere_delay_m=100.;v.troposphere_delay_m=100.;}
        o<<"V,"<<id<<','<<v.elevation_rad<<','<<v.snr_dbhz<<','<<v.pseudorange_sigma_m<<','<<v.ionosphere_delay_m<<','<<v.troposphere_delay_m<<','<<v.ionosphere_corrected<<','<<v.troposphere_corrected<<','<<v.snr_reference_dbhz<<','<<v.ionosphere_sigma_scale<<','<<v.troposphere_sigma_scale<<','<<spp_utils::calculatePseudorangeVariance(v)<<'\n';
    }
    for(int id=0;id<24;++id){
        std::vector<double> r;double threshold=3.,floor=3.;
        for(int i=0;i<id%10;++i)r.push_back(.3*(i%3-1));
        if(id%3==0)r.push_back(20.);
        if(id==10)r={0,0,0,0,9.};
        if(id==11)r={0,0,0,0,std::nextafter(9.,inf)};
        if(id==12)r={nan,nan,0,1,2,3};
        if(id==13)r={nan,0,1,2,3,4};
        if(id==14)r={-100,-50,0,50,100};
        if(id==15){r={0,0,0,0,20};threshold=nan;}
        if(id==16)floor=nan;
        if(id==17){r={0,0,0,0,2};floor=-1;}
        if(id==18){r={0,0,0,1,2,3};threshold=.01;}
        const auto mask=spp_utils::calculateResidualInlierMask(r,threshold,floor);
        o<<"M,"<<id<<','<<r.size()<<','<<threshold<<','<<floor;for(auto v:r)o<<','<<v;for(auto v:mask)o<<','<<v;o<<'\n';
    }
    for(int id=0;id<24;++id){
        double baseline=10.,ratio=.25,minimum=1.;std::vector<double> r={9.,7.5,8.};
        if(id==1)r={9.,std::nextafter(7.5,inf),8.};
        if(id==2){baseline=2.;r={1.,1.,2.};}
        if(id==3){baseline=2.;r={std::nextafter(1.,inf),2.};}
        if(id==4)r={nan,inf,11.};
        if(id==5)r={0.,0.,3.};
        if(id==6)baseline=0.;
        if(id==7)baseline=nan;
        if(id==8){ratio=-1.;minimum=nan;r={9.999,10.};}
        if(id==9)r.clear();
        if(id>=10){baseline=1.3*id;ratio=(id%5)*.1;minimum=id%4;for(int i=0;i<4;++i)r.push_back(baseline-.7*(i+id%3));}
        o<<"F,"<<id<<','<<baseline<<','<<ratio<<','<<minimum<<','<<r.size();for(auto v:r)o<<','<<v;o<<','<<spp_utils::selectRaimFdeCandidate(baseline,r,ratio,minimum)<<'\n';
    }
    return o?0:1;
}
