// Execute pinned upstream measurement/ambiguity/Kalman components and generate
// independent base/rover observations with upstream broadcast/geometry models.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#define main navigation_fixture_main
#include "upstream_navigation.cpp"
#undef main
#include <libgnss++/algorithms/kalman.hpp>
#include <libgnss++/algorithms/rtk_measurement.hpp>

static void vector_csv(std::ostream& out,const VectorXd& v) {
    for (int i=0;i<v.size();++i) out << ',' << v(i);
}
static void matrix_csv(std::ostream& out,const MatrixXd& a) {
    for (int i=0;i<a.rows();++i) for (int j=0;j<a.cols();++j) out << ',' << a(i,j);
}
static void update_case(std::ostream& out,const std::string& name,VectorXd x,MatrixXd p,
                        const MatrixXd& h,const VectorXd& v,const MatrixXd& r) {
    const int n=x.size(),m=v.size();
    out << name << ',' << n << ',' << m;
    vector_csv(out,x);matrix_csv(out,p);matrix_csv(out,h);vector_csv(out,v);matrix_csv(out,r);
    VectorXd weighted,hph;
    if (kalmanFilter(x,p,h,v,r,std::vector<bool>(n,true),&weighted,&hph)!=0) throw std::runtime_error("native filter failed");
    vector_csv(out,x);matrix_csv(out,p);out << ',' << v.dot(weighted);
    vector_csv(out,VectorXd(hph+r.diagonal()));
    std::vector<rtk_measurement::AmbiguityDifference> pairs;
    for (int i=4;i<n;++i) pairs.push_back({3,i});
    const auto a=rtk_measurement::buildAmbiguityTransform(x,p,3,pairs);
    out << ',' << pairs.size();
    for (const auto& pair:pairs) out << ',' << pair.reference_state_index << ',' << pair.satellite_state_index;
    vector_csv(out,a.dd_float);matrix_csv(out,a.ambiguity_covariance);matrix_csv(out,a.head_ambiguity_covariance);out << '\n';
}
static void write_observation_header(std::ostream& out,const Vector3d& position) {
    header(out,"     3.04           O                   M","RINEX VERSION / TYPE");
    header(out,"gnss-rust RTK FLOAT reference generator","PGM / RUN BY / DATE");
    std::ostringstream text;text << std::fixed << std::setprecision(4);
    for (int i=0;i<3;++i) text << std::setw(14) << position(i);
    header(out,text.str(),"APPROX POSITION XYZ");
    header(out,"G    2 C1C L1C","SYS / # / OBS TYPES");
    header(out,"","END OF HEADER");
}
#ifndef GNSS_RUST_RTK_FLOAT_HELPERS_ONLY
int main(int argc,char** argv) {
    if (argc!=2) {std::cerr << "usage: upstream-rtk-float OUTPUT_DIRECTORY\n";return 2;}
    const std::string root=argv[1];
    std::ofstream updates(root+"/upstream_rtk_update.csv");
    updates << "# name,n,m,x[n],P[n*n],H[m*n],v[m],R[m*m],xp[n],Pp[n*n],NIS,innovation_variances[m],k,pairs[2*k],DD[k],Qdd[k*k],QheadDD[3*k]; row-major matrices\n" << std::setprecision(17);
    {
        const int n=6,m=5;VectorXd x(n),v(m);MatrixXd b(n,n),h(m,n);
        x << 0,2.5,-3.1,0,12.3,-4.7;
        for(int i=0;i<n;++i) for(int j=0;j<n;++j) b(i,j)=0.3*std::sin((i+1)*1.2+(j+1)*2.1);
        MatrixXd p=b*b.transpose()+MatrixXd::Identity(n,n)*0.5;
        for(int i=0;i<m;++i) {v(i)=0.07*std::cos(i+1);for(int j=0;j<n;++j) h(i,j)=std::sin((i+1)*0.7+(j+1)*1.3);}
        MatrixXd r=MatrixXd::Identity(m,m)*0.02+MatrixXd::Constant(m,m,0.01);
        update_case(updates,"dense_zero_states",x,p,h,v,r);
    }
    {
        const int n=8;VectorXd x=VectorXd::Zero(n);x << 0,-8,4,0,3.2,-4.1,5.6,7.3;
        MatrixXd p=MatrixXd::Identity(n,n)*900;
        std::vector<rtk_measurement::MeasurementBlock> blocks(2);
        for(int kind=0;kind<2;++kind) {
            blocks[kind].kind=kind==0?rtk_measurement::MeasurementKind::PHASE:rtk_measurement::MeasurementKind::CODE;
            for(int i=0;i<4;++i) {
                rtk_measurement::MeasurementRow row;row.baseline_coefficients << 0.2*(i+1),std::sin(i+1),std::cos(i+1);
                row.residual=kind==0?0.02*(i+1):0.3*(i-1);
                row.reference_variance=kind==0?0.00002:0.3;
                row.satellite_variance=(kind==0?0.00002:0.3)*(i+1);
                if(kind==0) row.state_coefficients={{3,constants::GPS_L1_WAVELENGTH},{4+i,-constants::GPS_L1_WAVELENGTH}};
                blocks[kind].rows.push_back(row);
            }
        }
        const auto system=rtk_measurement::assembleMeasurementSystem(blocks,n);
        update_case(updates,"dd_phase_code",x,p,system.design_matrix,system.residuals,system.covariance);
    }
    {
        const int n=6;VectorXd x=VectorXd::Constant(n,1),v(2);v << 1e-4,-2e-4;
        MatrixXd p=MatrixXd::Identity(n,n)*1e-3,h=MatrixXd::Zero(2,n),r=MatrixXd::Identity(2,2)*1e-8;
        h(0,0)=1;h(0,3)=0.19;h(0,4)=-0.19;h(1,1)=1;h(1,3)=0.19;h(1,5)=-0.19;
        update_case(updates,"precise_phase",x,p,h,v,r);
    }
    const double rad=std::acos(-1.0)/180;
    const Vector3d base=geodetic2ecef(35*rad,139*rad,45);
    const Vector3d rover=base+enu2ecef(Vector3d(12,-7,3),35*rad,139*rad);
    const Vector3d seed=rover+Vector3d(1,-0.5,0.3);
    std::ofstream nav(root+"/synthetic_rtk.nav"),truth(root+"/synthetic_rtk_truth.csv");
    std::ofstream base_obs(root+"/synthetic_rtk_base.obs"),rover_obs(root+"/synthetic_rtk_rover.obs");
    std::ofstream noisy_base(root+"/synthetic_rtk_base_noisy.obs"),noisy_rover(root+"/synthetic_rtk_rover_noisy.obs");
    if(!updates || !nav || !truth || !base_obs || !rover_obs || !noisy_base || !noisy_rover) return 1;
    header(nav,"     3.04           N                   M","RINEX VERSION / TYPE");header(nav,"","END OF HEADER");
    write_observation_header(base_obs,base);write_observation_header(noisy_base,base);
    write_observation_header(rover_obs,seed);write_observation_header(noisy_rover,seed);
    truth << std::setprecision(17) << "# base_xyz,rover_xyz,seed_xyz\n";
    for(const auto& position:{base,rover,seed}) for(int i=0;i<3;++i) truth << position(i) << (position==seed && i==2?'\n':',');
    std::vector<Ephemeris> visible;
    for(int prn=1;prn<=32;++prn) {
        auto e=make_eph(GNSSSystem::GPS,prn);e.omega0=0.31*prn;e.m0=0.77*prn;
        Vector3d p;double clock;state(e,e.toe,p,clock);
        const auto enu=ecef2enu(p-base,35*rad,139*rad);
        if(std::atan2(enu(2),std::hypot(enu(0),enu(1)))>15*rad) {
            visible.push_back(e);write_nav(nav,e);
            truth << e.satellite.toString() << ',' << -10000+prn+10*(prn%3) << '\n';
        }
    }
    for(int epoch=0;epoch<32;++epoch) {
        const std::time_t label=315964800LL+2300LL*604800+345600+epoch*30;
        char date[30];std::strftime(date,sizeof(date),"%Y %m %d %H %M",std::gmtime(&label));
        for(bool is_rover:{false,true}) for(bool noisy:{false,true}) {
            auto& out=is_rover?(noisy?noisy_rover:rover_obs):(noisy?noisy_base:base_obs);
            const Vector3d receiver=is_rover?rover:base;
            double lat,lon,height;ecef2geodetic(receiver,lat,lon,height);
            out << "> " << date << ' ' << std::fixed << std::setprecision(7) << std::setw(10)
                << static_cast<double>(epoch%2*30) << "  0" << std::setw(3) << visible.size() << '\n';
            for(const auto& e:visible) {
                const auto receive=GNSSTime(2300,345600+epoch*30.0);
                double code=26500000,phase=0,clock,initial_clock;Vector3d p;
                for(int iteration=0;iteration<12;++iteration) {
                    auto tx=receive-code/constants::SPEED_OF_LIGHT;state(e,tx,p,initial_clock);
                    tx=tx-initial_clock;state(e,tx,p,clock);
                    const auto enu=ecef2enu(p-receiver,lat,lon);
                    const double el=std::atan2(enu(2),std::hypot(enu(0),enu(1)));
                    const double range=geodist(p,receiver)+(is_rover?2e-6:-1e-6)*constants::SPEED_OF_LIGHT
                        -clock*constants::SPEED_OF_LIGHT+models::tropDelaySaastamoinen(receiver,el);
                    code=range+e.tgd*constants::SPEED_OF_LIGHT;
                    const int ambiguity=is_rover?10000+3*e.satellite.prn+10*(e.satellite.prn%3):20000+2*e.satellite.prn;
                    phase=range/constants::GPS_L1_WAVELENGTH+ambiguity;
                }
                if(noisy) {
                    code+=0.05*std::sin((epoch+1)*0.93+e.satellite.prn*1.37+(is_rover?0.2:2.1));
                    phase+=0.0005/constants::GPS_L1_WAVELENGTH*std::cos((epoch+1)*1.13+e.satellite.prn*0.77+(is_rover?0.8:2.7));
                }
                out << e.satellite.toString() << std::fixed << std::setprecision(5)
                    << std::setw(14) << code << "  " << std::setprecision(4) << std::setw(14) << phase << "  \n";
            }
        }
    }
    std::cerr << "generated 3 native RTK component cases and 32 static epochs, " << visible.size() << " GPS satellites\n";
}
#endif
