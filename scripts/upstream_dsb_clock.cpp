// Independent Eigen gauge solve; not a native DSB-to-clock conversion API.
// Uses pinned native IF constants and additive measurement-correction sign.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <Eigen/Dense>
#include <libgnss++/core/constants.hpp>
#include <libgnss++/algorithms/ppp_correction_contract.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
    const double a=libgnss::constants::IFLC_F1,b=libgnss::constants::IFLC_F2;
    Eigen::Matrix4d h;
    h << 1,-1,0,0, -1,0,1,0, 0,1,0,-1, a,b,0,0;
    const auto decomposition=h.fullPivLu();
    if(decomposition.rank()!=4) return 2;
    const Eigen::Matrix4d inverse=decomposition.solve(Eigen::Matrix4d::Identity());
    const double sign=libgnss::algorithms::ppp_correction_contract::measurementCorrectionSign(false,false);
    std::cout << "case,d12_m,d1c_m,d2l_m,c1w_add_m,c2w_add_m,c1c_add_m,c2l_add_m,c1w_sigma_bound_m,c2w_sigma_bound_m,c1c_sigma_bound_m,c2l_sigma_bound_m\n" << std::setprecision(17);
    for(int i=0;i<24;++i) {
        const double d12=i==0?0.0:(i==1?1e-12:12.0*std::sin(0.31*i));
        const double d1c=i==0?0.0:3.0*std::cos(0.47*i),d2l=i==0?0.0:-1.7*std::sin(0.73*i);
        const Eigen::Vector4d rhs(d12,d1c,d2l,0.0),bias=decomposition.solve(rhs);
        if((h*bias-rhs).norm()>1e-11) return 3;
        std::cout << i << ',' << d12 << ',' << d1c << ',' << d2l;
        for(int j=0;j<4;++j)std::cout << ',' << sign*bias[j];
        for(int j=0;j<4;++j)std::cout << ',' << std::abs(inverse(j,0))*0.01+std::abs(inverse(j,1))*0.02+std::abs(inverse(j,2))*0.03;
        std::cout << '\n';
    }
}
