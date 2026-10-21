// Native finite-variance prefilter and progressive variance drop paths.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/algorithms/rtk_ar_selection.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace libgnss;
int main(int argc,char** argv){
    if(argc!=2)return 2;
    std::ofstream out(argv[1]);
    if(!out)return 3;
    out<<std::setprecision(17)<<"# id,n,min_pairs,max_drop,multiplier,min_threshold,variances,selected_count,selected_indices,drop_count,(count,indices)*\n";
    for(int id=0;id<40;++id){
        std::vector<double> v;int minimum=4,steps=6;double multiplier=10,threshold=1e-4;
        switch(id){
            case 0:break;
            case 1:v={1,1,1,1,10};break;
            case 2:v={1,1,1,1,std::nextafter(10.,INFINITY)};break;
            case 3:v={1,2,3,29};break;
            case 4:v={0,0,0,1e-4,std::nextafter(1e-4,INFINITY)};break;
            case 5:v={0,0,0,1};break;
            case 6:v={0,0,0,0,0,0};threshold=0;break;
            case 7:v={4,4,3,2,1,1,1};break;
            case 8:v={1,2,3};break;
            case 9:v={1,1,1,1,100,10};steps=0;break;
            case 10:v={.3};minimum=1;break;
            case 11:v={1,1,1,1,1,1};minimum=1;steps=64;break;
            default:
                for(int i=0;i<4+id%13;++i)v.push_back(std::pow(10.,(i+id)%7-4)*(.5+.1*(i%4)));
                minimum=1+id%5;steps=id%9;multiplier=1+id%11;threshold=id%2?1e-4:0;
        }
        std::vector<rtk_ar_selection::PairDescriptor> pairs;
        for(double variance:v){rtk_ar_selection::PairDescriptor p;p.variance=variance;pairs.push_back(p);}
        auto selected=rtk_ar_selection::filterPairsByRelativeVariance(pairs,multiplier,threshold,minimum);
        auto dropped=rtk_ar_selection::buildProgressiveVarianceDropSubsets(pairs,minimum,steps);
        out<<id<<','<<v.size()<<','<<minimum<<','<<steps<<','<<multiplier<<','<<threshold;
        for(auto x:v)out<<','<<x;
        out<<','<<selected.size();for(auto x:selected)out<<','<<x;
        out<<','<<dropped.size();for(auto subset:dropped){out<<','<<subset.size();for(auto x:subset)out<<','<<x;}
        out<<'\n';
    }
    return out?0:1;
}
