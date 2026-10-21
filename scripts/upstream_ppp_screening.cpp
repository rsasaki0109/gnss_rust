// Actual pinned native standardized-residual helper; source is unchanged.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <bits/stdc++.h>
#include <Eigen/Dense>
#include "ppp_internal.hpp"
using namespace libgnss;
int main(int argc,char**argv) {
    if(argc!=2)return 2;
    std::ofstream out(std::string(argv[1])+"/upstream_ppp_screening.csv");out<<std::setprecision(17);
    std::ofstream gaussian(std::string(argv[1])+"/upstream_ppp_screening_gaussian.csv");gaussian<<std::setprecision(17);
    out<<"id,m,threshold,worst,innovations_then_S_then_inverse_then_statistics\n";
    gaussian<<"id,m,threshold,worst,innovations_then_S_then_inverse_then_statistics\n";
    for(int id=0;id<64;++id) {
        int m=id<2?5:id<7?4:7+id%10;double gate=id%3?4:2.5;
        Eigen::MatrixXd b(m,m);for(int i=0;i<m;++i)for(int j=0;j<m;++j)b(i,j)=.2*std::sin(.3*i+.7*j+id);
        Eigen::MatrixXd S=b*b.transpose()+Eigen::MatrixXd::Identity(m,m)*(.2+.1*(id%5));
        Eigen::VectorXd r(m);for(int i=0;i<m;++i)r[i]=(id%4?.3:7)*std::cos(.5*i+id);
        if(id<2) {S=Eigen::MatrixXd::Ones(m,m)*100+Eigen::MatrixXd::Identity(m,m);r<<.3,-.4,.2,-.1,id==0?20:.5;gate=4;}
        if(id>=2&&id<7) {S.setIdentity();r.setZero();gate=4;r[1]=id==2?4:id==3?std::nextafter(4.0,5.0):id==4?std::nextafter(4.0,3.0):6;if(id>=5)r[3]=id==5?-6:-7;}
        if(id>=7&&id%4==0) { // large shared-state uncertainty is not small measurement noise
            S+=Eigen::MatrixXd::Ones(m,m)*1e4;
        }
        Eigen::MatrixXd inverse=S.inverse();
        const int worst=ppp_internal::worstStandardizedResidualRow(r,inverse,gate);
        const Eigen::VectorXd weighted=inverse*r;
        out<<id<<','<<m<<','<<gate<<','<<worst;
        for(int i=0;i<m;++i)out<<','<<r[i];
        for(int i=0;i<m;++i)for(int j=0;j<m;++j)out<<','<<S(i,j);
        for(int i=0;i<m;++i)for(int j=0;j<m;++j)out<<','<<inverse(i,j);
        for(int i=0;i<m;++i)out<<','<<weighted[i]/std::sqrt(inverse(i,i));
        out<<'\n';
        // Production PPP computes S^-1 with LDLT solve, unlike .inverse()
        // in the upstream helper unit test. Keep both paths explicitly.
        const Eigen::MatrixXd production_inverse=S.ldlt().solve(Eigen::MatrixXd::Identity(m,m));
        const int production_worst=ppp_internal::worstStandardizedResidualRow(r,production_inverse,gate);
        const Eigen::VectorXd production_weighted=production_inverse*r;
        gaussian<<id<<','<<m<<','<<gate<<','<<production_worst;
        for(int i=0;i<m;++i)gaussian<<','<<r[i];
        for(int i=0;i<m;++i)for(int j=0;j<m;++j)gaussian<<','<<S(i,j);
        for(int i=0;i<m;++i)for(int j=0;j<m;++j)gaussian<<','<<production_inverse(i,j);
        for(int i=0;i<m;++i)gaussian<<','<<production_weighted[i]/std::sqrt(production_inverse(i,i));
        gaussian<<'\n';
    }
    return out&&gaussian?0:4;
}
