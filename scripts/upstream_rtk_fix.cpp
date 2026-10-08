// Native fixed-head and CP/PR comparisons, independent of the Rust algorithms.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#define GNSS_RUST_RTK_FLOAT_HELPERS_ONLY
#include "upstream_rtk_float.cpp"
#include <libgnss++/algorithms/rtk_ar_evaluation.hpp>
#include <libgnss++/algorithms/rtk_cp_pr_gate.hpp>

int main(int argc,char** argv) {
    if(argc!=2) {std::cerr << "usage: upstream-rtk-fix OUTPUT_DIRECTORY\n";return 2;}
    const std::string root=argv[1];
    std::ofstream out(root+"/upstream_rtk_fix.csv"),gate(root+"/upstream_rtk_code_phase.csv");
    if(!out || !gate) return 1;
    out << "# name,h,k,x[h],Qaa[h*h],Qab[h*k],Qbb[k*k],b_float[k],b_integer[k],x_fixed[h],Qfixed[h*h],cost; row-major matrices\n" << std::setprecision(17);
    for(int sample=0;sample<5;++sample) {
        const int h=sample==4?2:3,k=sample==4?5:4;
        MatrixXd l=MatrixXd::Zero(h+k,h+k);
        for(int i=0;i<k;++i) for(int j=0;j<=i;++j)
            l(i,j)=i==j?0.2+0.13*i:0.07*std::sin((i+1)*(j+1));
        for(int i=0;i<h;++i) {
            for(int j=0;j<k;++j) l(k+i,j)=(sample==3?1.2:0.03)*std::cos((i+1)*1.3+(j+1)*0.7);
            for(int j=0;j<=i;++j) l(k+i,k+j)=i==j?(sample==2?1e-4:0.005)*(i+1):0.001*std::sin(i+j+1);
        }
        const MatrixXd joint=l*l.transpose();
        const MatrixXd qbb=joint.topLeftCorner(k,k),qaa=joint.bottomRightCorner(h,h),qab=joint.bottomLeftCorner(h,k);
        VectorXd x(h),bf(k),bi(k),xf;
        for(int i=0;i<h;++i) x(i)=sample==1?0.0:-1.3+2.7*i;
        for(int i=0;i<k;++i) {bi(i)=(i%2?-1:1)*(i*3+1);bf(i)=bi(i)+(sample==1?0.0:0.07*std::cos(i+sample+1));}
        if(!rtk_ar_evaluation::solveFixedHeadState(x,qab,qbb,bf,bi,xf)) return 1;
        // Independent Eigen Schur complement: upstream exports the fixed head
        // mean, but has no corresponding conditional-covariance API.
        const MatrixXd cf=qaa-qab*qbb.ldlt().solve(qab.transpose());
        const VectorXd delta=bf-bi;
        out << "conditional_" << sample << ',' << h << ',' << k;
        vector_csv(out,x);matrix_csv(out,qaa);matrix_csv(out,qab);matrix_csv(out,qbb);
        vector_csv(out,bf);vector_csv(out,bi);vector_csv(out,xf);matrix_csv(out,cf);
        out << ',' << delta.dot(qbb.ldlt().solve(delta)) << '\n';
    }
    gate << "# name,k,threshold,triples[3*k code/phase/ambiguity_m],consistent,bad_pairs,rms,max_abs\n" << std::setprecision(17);
    for(int sample=0;sample<3;++sample) {
        std::vector<rtk_cp_pr_gate::Observation> observations;
        for(int i=0;i<4;++i) {
            const double innovation=sample==0?(i%2?-0.125:0.125):sample==1?0.5:(i==2?0.625:0.125);
            observations.push_back({double(i),double(i)+2.0-innovation,2.0});
        }
        rtk_cp_pr_gate::Config config;config.innovation_threshold_m=0.5;config.min_pairs=4;
        const auto result=rtk_cp_pr_gate::evaluate(observations,config);if(!result.valid) return 1;
        gate << "code_phase_" << sample << ',' << observations.size() << ',' << config.innovation_threshold_m;
        for(const auto& o:observations) gate << ',' << o.dd_pseudorange_m << ',' << o.dd_carrier_m << ',' << o.fixed_ambiguity_m;
        gate << ',' << result.consistent << ',' << result.bad_pairs << ',' << result.rms_innovation_m << ',' << result.max_abs_innovation_m << '\n';
    }
    std::cerr << "generated 5 native conditional-head and 3 code/phase gate cases\n";
}
