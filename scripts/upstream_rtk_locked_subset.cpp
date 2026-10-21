// Native lock-qualified DD selection and marginal ambiguity conditioning.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/algorithms/rtk_selection.hpp>
#include <libgnss++/algorithms/rtk_ar_evaluation.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <stdexcept>
using namespace libgnss;
static void vector_csv(std::ostream& o,const Eigen::VectorXd& x){for(int i=0;i<x.size();++i)o<<','<<x(i);}
static void matrix_csv(std::ostream& o,const Eigen::MatrixXd& x){for(int i=0;i<x.rows();++i)for(int j=0;j<x.cols();++j)o<<','<<x(i,j);}
int main(int argc,char** argv){
    if(argc!=2)return 2;
    std::ofstream out(argv[1]);
    if(!out)return 3;
    out<<std::setprecision(17)<<"# id,satellites,states,pairs,locks,elevations,x,P,ref_target_indices,dd,Qbb,Qab,fixed_head,fixed_covariance; row-major\n";
    const int count=9,n=12;
    for(int id=0;id<24;++id){
        std::vector<rtk_selection::SatelliteSelectionData> sats(count);
        Eigen::VectorXd x(n);Eigen::MatrixXd b(n,n);
        for(int i=0;i<n;++i){x(i)=.23*(i+1)+.01*id;for(int j=0;j<n;++j)b(i,j)=.1*std::sin(.7*(i+1)+1.3*(j+1)+id);}
        Eigen::MatrixXd p=b*b.transpose()+.2*Eigen::MatrixXd::Identity(n,n);
        for(int i=0;i<count;++i){
            auto& s=sats[i];s.satellite={i<5?GNSSSystem::GPS:GNSSSystem::Galileo,static_cast<std::uint8_t>(i<5?i+1:i-4)};
            s.has_l1=s.n1_active=true;s.l1_wavelength=299792458./1575.42e6;s.elevation=.2+.06*i;s.lock_count_l1=5;
            if(id==1&&i==4)s.lock_count_l1=4;
            if(id==2&&i==8)s.lock_count_l1=4;
            if(id==3&&i>=6)s.lock_count_l1=4;
            if(id==4&&i>=5)s.lock_count_l1=4;
            if(id==5)s.lock_count_l1=0;
            if(id==6)s.lock_count_l1=i==0?5:0;
            if(id==7)s.elevation=.5;
            if(id>=8)s.lock_count_l1=(3*i+id)%8;
        }
        std::vector<rtk_selection::SelectionPair> selected,full;
        for(auto system:{GNSSSystem::GPS,GNSSSystem::Galileo}){
            auto sp=rtk_selection::buildDoubleDifferencePairsForSystem(sats,system,5,true);
            selected.insert(selected.end(),sp.begin(),sp.end());
            SatelliteId reference;
            if(!rtk_selection::selectSystemReferenceSatellite(sats,system,5,reference))
                if(!rtk_selection::selectSystemReferenceSatellite(sats,system,0,reference))throw std::runtime_error("no reference");
            auto fp=rtk_selection::buildDoubleDifferencePairsForSystem(sats,system,0,true,&reference);
            full.insert(full.end(),fp.begin(),fp.end());
        }
        auto index=[&](SatelliteId s){for(int i=0;i<count;++i)if(sats[i].satellite==s)return i;throw std::runtime_error("unknown satellite");};
        Eigen::MatrixXd d=Eigen::MatrixXd::Zero(full.size(),n);std::vector<int> subset;
        for(int i=0;i<(int)full.size();++i){d(i,3+index(full[i].ref_sat))=1;d(i,3+index(full[i].sat))=-1;}
        for(auto pair:selected){int found=-1;for(int i=0;i<(int)full.size();++i)if(full[i].sat==pair.sat&&full[i].ref_sat==pair.ref_sat)found=i;if(found<0)throw std::runtime_error("missing DD");subset.push_back(found);}
        auto m=rtk_ar_evaluation::extractSubset(d*x,d*p*d.transpose(),p.topRows(3)*d.transpose(),subset);
        out<<id<<','<<count<<','<<n<<','<<selected.size();
        for(auto s:sats)out<<','<<s.lock_count_l1;
        for(auto s:sats)out<<','<<s.elevation;
        vector_csv(out,x);matrix_csv(out,p);
        for(auto pair:selected)out<<','<<index(pair.ref_sat)<<','<<index(pair.sat);
        vector_csv(out,m.dd_float);matrix_csv(out,m.Qb);matrix_csv(out,m.Qab);
        if(!selected.empty()){
            Eigen::VectorXd integers=m.dd_float.unaryExpr([](double v){return std::round(v);}),fixed;
            if(!rtk_ar_evaluation::solveFixedHeadState(x.head(3),m.Qab,m.Qb,m.dd_float,integers,fixed))throw std::runtime_error("native condition failed");
            vector_csv(out,fixed);matrix_csv(out,p.topLeftCorner(3,3)-m.Qab*m.Qb.ldlt().solve(m.Qab.transpose()));
        }
        out<<'\n';
    }
    return out?0:1;
}
