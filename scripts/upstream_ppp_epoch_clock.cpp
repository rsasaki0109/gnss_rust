// Actual pinned native scalar-state covariance reset; not whole precise PPP.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <bits/stdc++.h>
#include <Eigen/Dense>
#define private public
#include <libgnss++/algorithms/ppp.hpp>
#undef private
using namespace libgnss;
int main(int argc,char**argv) {
    if(argc!=2)return 2;
    std::ofstream out(std::string(argv[1])+"/upstream_ppp_epoch_clock.csv");
    out << "# id,n,new_variance,P_before,P_after; row-major; clock slot 3\n" << std::setprecision(17);
    for(int id=0;id<16;++id) {
        const int n=5+id;Eigen::MatrixXd m(n,n);
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)m(i,j)=.1*std::sin(.8*i+1.2*j+id);
        PPPProcessor processor;auto& f=processor.filter_state_;f.total_states=n;f.state=Eigen::VectorXd::Zero(n);
        f.covariance=m*m.transpose()+.2*Eigen::MatrixXd::Identity(n,n);
        const double variance=id%2?1e8:100.0+id;
        out << id << ',' << n << ',' << variance;
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)out << ',' << f.covariance(i,j);
        processor.reinitializeScalarState(3,0.0,variance);
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)out << ',' << f.covariance(i,j);
        out << '\n';
    }
    return out?0:3;
}
