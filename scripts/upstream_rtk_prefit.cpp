// Native row classifications only; the Rust epoch gate takes a different action.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/algorithms/rtk_measurement.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
    std::cout << std::setprecision(17);
    std::cout << "# case,n_phase,threshold,phase[n],code[n],native_mask[2n],native_count\n";
    auto run = [](const char* name, std::vector<double> phase, std::vector<double> code, double threshold) {
        if(phase.size()!=code.size()) throw std::runtime_error("shape");
        const int n=static_cast<int>(phase.size()), m=2*n;
        Eigen::VectorXd v(m); Eigen::MatrixXd h(m,1);
        std::cout << name << ',' << n << ',' << threshold;
        for(int i=0;i<m;++i) {
            v(i)=i<n?phase[i]:code[i-n]; h(i,0)=i+1;
            std::cout << ',' << v(i);
        }
        const int count=libgnss::rtk_measurement::suppressOutlierRows(v,h,threshold);
        for(int i=0;i<m;++i) std::cout << ',' << (h(i,0)==0.0);
        std::cout << ',' << count << '\n';
    };
    run("small",{0.001,-0.002,0.003},{0.1,-0.2,0.3},3);
    run("boundary",{3,-3,std::nextafter(3.0,INFINITY)},
        {-std::nextafter(3.0,INFINITY),std::nextafter(3.0,0.0),0},3);
    run("phase_only",{10,0,-10},{0,0,0},3);
    run("code_only",{0,0,0},{10,0,-10},3);
    run("mixed",{0,4,0},{-4,0,4},3);
    run("all",{4,-4,100},{-100,4,-4},3);
    run("negative",{-3,-2,-1},{-4,-5,-6},3);
    run("other_threshold",{0.01,-0.1,0.3},{-0.3,0.2,-0.01},0.2);
    run("zero",{0},{0},3);
    run("many",{0,1,2,3,4,5,6,7},{7,6,5,4,3,2,1,0},3);
}
