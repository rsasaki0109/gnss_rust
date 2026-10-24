// Pinned native ZTD prior-variance helper; source remains unchanged.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <bits/stdc++.h>
#include <Eigen/Dense>
#include "ppp_internal.hpp"
int main(int argc,char**argv) {
    if(argc!=2)return 2;
    std::ofstream out(std::string(argv[1])+"/upstream_ppp_ztd_prior.csv");out<<std::setprecision(17);
    out<<"madoca_per_frequency,broadcast_model,configured_variance_m2,returned_variance_m2\n";
    for(bool madoca:{false,true})for(bool broadcast:{false,true})for(double configured:{.36,25.0,100.0})
        out<<madoca<<','<<broadcast<<','<<configured<<','<<libgnss::ppp_internal::initialTroposphereVariance(madoca,broadcast,configured)<<'\n';
    return out?0:4;
}
