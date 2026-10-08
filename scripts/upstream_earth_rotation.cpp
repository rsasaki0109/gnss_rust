// Pinned native wrapper and unchanged SOFA time/ERA/polar/CIO references.
// Correct GPS->TAI->UTC oracle is retained separately from native one-pass UTC.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <cmath>
#include <fstream>
#include <iomanip>
#include <array>
#include "libgnss++/iers/earth_rotation.hpp"
extern "C" {
#include "sofa.h"
}
std::pair<int,double> split(double a,double b) {
    const double base=a-2400000.5;
    const double whole=std::floor(base),fraction=base-whole+b;
    const double shift=std::floor(fraction);
    return {static_cast<int>(whole+shift),(fraction-shift)*86400.0};
}
void emit(std::ostream& out,int id,int week,double tow,int mode) {
    const int day=44244+week*7+static_cast<int>(std::floor(tow/86400.0));
    double tai1=2400000.5+day,tai2=(std::fmod(tow,86400.0)+19.0)/86400.0;
    // Normalize to the same day+seconds partition used by Rust, preserving fractions.
    if(tai2>=1) {tai1+=1;tai2-=1;}
    double utc1,utc2,tt1,tt2;
    if(iauTaiutc(tai1,tai2,&utc1,&utc2)<0 || iauTaitt(tai1,tai2,&tt1,&tt2)!=0) throw std::runtime_error("SOFA time failed");
    int year,month,date,hms[4];
    if(iauD2dtf("UTC",9,utc1,utc2,&year,&month,&date,hms)<0) throw std::runtime_error("SOFA calendar failed");
    double mjd0,utcday,dat;
    iauCal2jd(year,month,date,&mjd0,&utcday);iauDat(year,month,date,0.0,&dat);
    double fraction;int ny,nm,nd;iauJd2cal(2400000.5,utcday+1,&ny,&nm,&nd,&fraction);
    double nextdat;iauDat(ny,nm,nd,0.0,&nextdat);
    const double daylength=86400+nextdat-dat;
    const double utcsec=hms[0]*3600+hms[1]*60+hms[2]+hms[3]*1e-9;
    libgnss::iers::EarthOrientationParams eop;
    eop.ut1_minus_utc_seconds=mode==0?0.0:(mode==1?-0.1234:0.3345);
    eop.xp_arcsec=mode==0?0.0:0.1*mode;eop.yp_arcsec=mode==0?0.0:-0.2*mode;
    double ut11,ut12;iauUtcut1(utc1,utc2,eop.ut1_minus_utc_seconds,&ut11,&ut12);
    const double era=iauEra00(ut11,ut12),sp=iauSp00(tt1,tt2);
    double ci[3][3],pom[3][3],c2t[3][3];iauC2i06a(tt1,tt2,ci);
    iauPom00(eop.xp_arcsec*M_PI/180.0/3600.0,eop.yp_arcsec*M_PI/180.0/3600.0,sp,pom);
    iauC2tcio(ci,era,pom,c2t);
    libgnss::GNSSTime gps;gps.week=week;gps.tow=tow;
    const double nativeutc=libgnss::iers::gnssTimeToMjdUtc(gps);
    const auto tt=split(tt1,tt2),ut1=split(ut11,ut12);
    out << id << ',' << week << ',' << tow << ',' << utcday << ',' << utcsec << ',' << daylength << ',' << dat
        << ',' << (utc1-2400000.5)+utc2 << ',' << nativeutc << ',' << tt.first << ',' << tt.second
        << ',' << ut1.first << ',' << ut1.second << ',' << eop.ut1_minus_utc_seconds << ',' << eop.xp_arcsec << ',' << eop.yp_arcsec << ',' << era << ',' << sp;
    for(const auto& matrix:{ci,pom,c2t}) for(int i=0;i<3;++i) for(int j=0;j<3;++j) out << ',' << matrix[i][j];
    out << ',' << ((tt1-2451545.0)+tt2)/36525.0 << ',' << (utcday-std::floor(utcday)+utcsec/daylength)*24.0 << '\n';
}
int main(int argc,char**argv) {
    if(argc!=2) return 2;
    std::ofstream out(std::string(argv[1])+"/upstream_earth_rotation.csv");
    out << std::setprecision(17) << "id,week,tow,utc_day,utc_seconds,day_length,tai_utc,sofa_utc_mjd,native_utc_mjd,tt_day,tt_seconds,ut1_day,ut1_seconds,dut1_s,xp_arcsec,yp_arcsec,era_rad,sp_rad";
    for(const char* prefix:{"ci","pom","c2t"}) for(int i=0;i<3;++i) for(int j=0;j<3;++j) out << ',' << prefix << i << j;
    out << ",tt_centuries,utc_hour\n";
    int id=0;
    for(int week:{0,1042,1527,1930,2300,2380}) for(double tow:{0.0,1.125,346200.25,604799.75}) for(int mode:{0,1,2}) emit(out,id++,week,tow,mode);
    const std::array<std::array<int,3>,18> transitions={{{1981,7,1},{1982,7,1},{1983,7,1},{1985,7,1},{1988,1,1},{1990,1,1},{1991,1,1},{1992,7,1},{1993,7,1},{1994,7,1},{1996,1,1},{1997,7,1},{1999,1,1},{2006,1,1},{2009,1,1},{2012,7,1},{2015,7,1},{2017,1,1}}};
    int index=0;
    for(const auto& date:transitions) {
        double mjd0,utcday;iauCal2jd(date[0],date[1],date[2],&mjd0,&utcday);
        const double effective_gps=(utcday-44244)*86400+(++index);
        for(double delta:{-86400.25,-43200.25,-20.0,-1.25,-0.75,-0.25,0.0,0.25,20.0}) for(int mode:{0,1,2}) {
            const double total=effective_gps+delta;
            const int week=static_cast<int>(std::floor(total/604800.0));
            emit(out,id++,week,total-week*604800.0,mode);
        }
    }
    return out?0:3;
}
