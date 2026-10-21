// Native-model moving GPS/Galileo observations with distinct receiver clocks.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#define GNSS_RUST_RTK_FLOAT_HELPERS_ONLY
#include "upstream_rtk_float.cpp"

static void mixed_header(std::ostream& out,const Vector3d& position,bool rover) {
    header(out,"     3.04           O                   M","RINEX VERSION / TYPE");
    header(out,"gnss-rust mixed-system reference generator","PGM / RUN BY / DATE");
    std::ostringstream text;text<<std::fixed<<std::setprecision(4);
    for(int i=0;i<3;++i)text<<std::setw(14)<<position(i);
    header(out,text.str(),"APPROX POSITION XYZ");
    header(out,"G    2 C1C L1C","SYS / # / OBS TYPES");
    header(out,rover?"E    2 C1C L1C":"E    2 C1X L1X","SYS / # / OBS TYPES");
    header(out,"","END OF HEADER");
}

int main(int argc,char**argv) {
    if(argc!=2)return 2;
    const std::string root=argv[1];const double rad=std::acos(-1.0)/180;
    const Vector3d base=geodetic2ecef(35*rad,139*rad,45);
    const Vector3d origin=base+enu2ecef(Vector3d(12,-7,3),35*rad,139*rad);
    const Vector3d seed=origin+Vector3d(1,-0.5,0.3);
    std::ofstream nav(root+"/synthetic_rtk_mixed.nav"),truth(root+"/synthetic_rtk_mixed_truth.csv"),
        rover_obs(root+"/synthetic_rtk_mixed_rover.obs"),base_obs(root+"/synthetic_rtk_mixed_base.obs");
    if(!nav||!truth||!rover_obs||!base_obs)return 1;
    header(nav,"     3.04           N                   M","RINEX VERSION / TYPE");
    header(nav,"","END OF HEADER");
    mixed_header(rover_obs,seed,true);mixed_header(base_obs,base,false);
    truth<<std::setprecision(17)<<"week,tow,x_m,y_m,z_m\n";
    std::vector<Ephemeris> visible;
    for(auto system:{GNSSSystem::GPS,GNSSSystem::Galileo}) {
        int count=0;
        for(int prn=1;prn<=32 && count<4;++prn) {
            auto e=make_eph(system,prn);
            const bool gal=system==GNSSSystem::Galileo;
            e.omega0=0.31*prn+(gal?0.7:0);e.m0=0.77*prn+(gal?1.1:0);
            Vector3d p;double clock;state(e,e.toe,p,clock);
            const auto enu=ecef2enu(p-base,35*rad,139*rad);
            if(std::atan2(enu(2),std::hypot(enu(0),enu(1)))>15*rad) {
                visible.push_back(e);write_nav(nav,e);++count;
            }
        }
        if(count!=4)throw std::runtime_error("insufficient synthetic sky coverage");
    }
    for(int epoch=0;epoch<32;++epoch) {
        const Vector3d rover=origin+enu2ecef(Vector3d(1.2*epoch,-0.8*epoch,0.1*epoch),35*rad,139*rad);
        truth<<2300<<','<<345600+epoch;
        for(int i=0;i<3;++i)truth<<','<<rover(i);truth<<'\n';
        const std::time_t label=315964800LL+2300LL*604800+345600+epoch;
        char date[30];std::strftime(date,sizeof(date),"%Y %m %d %H %M",std::gmtime(&label));
        for(bool is_rover:{false,true}) {
            auto&out=is_rover?rover_obs:base_obs;const Vector3d receiver=is_rover?rover:base;
            double lat,lon,height;ecef2geodetic(receiver,lat,lon,height);
            out<<"> "<<date<<' '<<std::fixed<<std::setprecision(7)<<std::setw(10)<<static_cast<double>(epoch)<<"  0"<<std::setw(3)<<visible.size()<<'\n';
            for(const auto&e:visible) {
                const bool gal=e.satellite.system==GNSSSystem::Galileo;
                const double receiver_clock=gal?(is_rover?-3e-6:4e-6):(is_rover?2e-6:-1e-6);
                const auto receive=GNSSTime(2300,345600+epoch);
                double code=26500000,phase=0,clock,initial;Vector3d p;
                for(int iteration=0;iteration<12;++iteration) {
                    auto tx=receive-code/constants::SPEED_OF_LIGHT;state(e,tx,p,initial);tx=tx-initial;state(e,tx,p,clock);
                    const auto enu=ecef2enu(p-receiver,lat,lon);
                    const double el=std::atan2(enu(2),std::hypot(enu(0),enu(1)));
                    const double range=geodist(p,receiver)+receiver_clock*constants::SPEED_OF_LIGHT-clock*constants::SPEED_OF_LIGHT
                        +models::tropDelaySaastamoinen(receiver,el);
                    code=range+e.tgd*constants::SPEED_OF_LIGHT;
                    phase=range/constants::GPS_L1_WAVELENGTH+(is_rover?10000+3*e.satellite.prn:20000+2*e.satellite.prn)
                        +(gal && !is_rover?0.25:0.0);
                }
                out<<e.satellite.toString()<<std::fixed<<std::setprecision(5)<<std::setw(14)<<code<<"  "
                    <<std::setprecision(4)<<std::setw(14)<<phase<<"  \n";
            }
        }
    }
    std::cerr<<"32 moving epochs; GPS/Galileo clocks differ by 10 microseconds in SD; Galileo base phase +0.25 cycle\n";
}
