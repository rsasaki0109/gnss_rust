// Actual pinned Dehant functions and libgnss++ SOFA/IERS wrappers.
// Explicit-body component fixtures; no complete native PPP processor.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <fstream>
#include <iomanip>
#include <cmath>
#include "libgnss++/iers/tides.hpp"
#include "libgnss++/iers/ephemeris.hpp"
#include "libgnss++/core/coordinates.hpp"
// Include the unmodified native implementation to expose its per-term helpers.
#include "dehanttideinel/dehanttide_all.cpp"
extern "C" {
#include "sofa.h"
}

void emit(std::ostream& out, int id, double t, double hour,
          const Eigen::Vector3d& station, const Eigen::Vector3d& sun,
          const Eigen::Vector3d& moon) {
    std::vector<Eigen::Vector3d> result;
    if(iers2010::dehanttideinel_impl(t,hour,sun,moon,{station},result)!=0 || result.size()!=1) throw std::runtime_error("native tide failed");
    const TideAux aux(station,sun,moon);const Step2Angles angles(t,hour);
    const auto diu=st1idiu(aux),sem=st1isem(aux),lat=st1l1(aux),diu2=step2diu(angles,aux),lon2=step2lon(angles,aux);
    const Eigen::Vector3d degree23=result[0]-diu-sem-lat-diu2-lon2,instant=station+result[0];
    out << id << ',' << t << ',' << hour;
    for(const auto& v:{station,sun,moon,degree23,diu,sem,lat,diu2,lon2,result[0],instant})
        for(int i=0;i<3;++i) out << ',' << v(i);
    out << '\n';
}

int main(int argc,char**argv) {
    if(argc!=2) return 2;
    std::ofstream kernel(std::string(argv[1])+"/upstream_iers2010_kernel.csv");
    kernel << std::setprecision(17)
        << "id,tt_centuries,ut_hour,nom_x,nom_y,nom_z,sun_x,sun_y,sun_z,moon_x,moon_y,moon_z,degree23_x,degree23_y,degree23_z,diu_x,diu_y,diu_z,sem_x,sem_y,sem_z,lat_x,lat_y,lat_z,diu2_x,diu2_y,diu2_z,lon2_x,lon2_y,lon2_z,dx,dy,dz,instant_x,instant_y,instant_z\n";
    const Eigen::Vector3d sun(137859926952.015,54228127881.4350,23509422341.6960),moon(-179996231.920342,-312468450.131567,-169288918.592160);
    emit(kernel,0,(54934.0+66.184/86400.0-51544.5)/36525.0,0.0,
         Eigen::Vector3d(4075578.385,931852.890,4801570.154),sun,moon);
    int id=1;
    for(double t:{-0.5,-0.2,0.0,0.0928,0.24,0.5}) for(double hour:{0.0,6.25,12.5,23.999999})
        for(double latitude:{-89.999,-60.0,-35.0,0.0,35.0,60.0,89.999}) for(double longitude:{-170.0,0.0,139.0}) {
            const double angle=(hour*15.0+t*37.0)*M_PI/180.0;
            Eigen::Matrix3d rotation;rotation << std::cos(angle),std::sin(angle),0,-std::sin(angle),std::cos(angle),0,0,0,1;
            const double height=std::array<double,3>{-50.0,45.0,1500.0}[id%3];
            emit(kernel,id++,t,hour,libgnss::geodetic2ecef(latitude*M_PI/180.0,longitude*M_PI/180.0,height),rotation*sun,rotation*moon);
        }
    std::ofstream epochs(std::string(argv[1])+"/upstream_iers2010_epochs.csv");
    epochs << std::setprecision(17)
        << "id,mjd_utc,tt_centuries,ut_hour,ut1_minus_utc_s,xp_arcsec,yp_arcsec,nom_x,nom_y,nom_z,sun_x,sun_y,sun_z,moon_x,moon_y,moon_z,dx,dy,dz,instant_x,instant_y,instant_z\n";
    id=0;
    for(double start:{54934.0,60348.0,60908.0}) for(int i=0;i<25;++i) for(bool nonzero_eop:{false,true}) {
        const double utc=start+i/24.0;
        libgnss::iers::EarthOrientationParams eop;
        if(nonzero_eop) {eop.ut1_minus_utc_seconds=-0.1234;eop.xp_arcsec=0.1;eop.yp_arcsec=-0.2;}
        double tai1,tai2,tt1,tt2;
        if(iauUtctai(2400000.5,utc,&tai1,&tai2)<0 || iauTaitt(tai1,tai2,&tt1,&tt2)!=0) return 3;
        const double t=((tt1-2451545.0)+tt2)/36525.0,hour=(utc-std::floor(utc))*24.0;
        const auto rotation=libgnss::iers::icrsToItrs(utc,eop);
        const Eigen::Vector3d sun=rotation*libgnss::iers::sunPositionIcrs(utc),moon=rotation*libgnss::iers::moonPositionIcrs(utc);
        const auto station=libgnss::geodetic2ecef((id%2?-35.0:35.0)*M_PI/180.0,139.0*M_PI/180.0,45.0);
        const auto displacement=libgnss::iers::solidEarthTideDisplacementAt(utc,station,eop);
        const Eigen::Vector3d instant=station+displacement;
        epochs << id++ << ',' << utc << ',' << t << ',' << hour << ',' << eop.ut1_minus_utc_seconds << ',' << eop.xp_arcsec << ',' << eop.yp_arcsec;
        for(const auto& v:{station,sun,moon,displacement,instant}) for(int j=0;j<3;++j) epochs << ',' << v(j);
        epochs << '\n';
    }
    return kernel && epochs ? 0:4;
}
