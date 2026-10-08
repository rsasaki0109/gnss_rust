// Generate navigation parity values and a deterministic RINEX-to-SPP fixture
// using unmodified libgnss++ broadcast/geometry/atmosphere functions.
#include "navigation_internal.hpp"
#include <libgnss++/core/coordinates.hpp>
#include <libgnss++/models/ionosphere.hpp>
#include <libgnss++/models/troposphere.hpp>
#include <array>
#include <cassert>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace libgnss;

Ephemeris make_eph(GNSSSystem system, int prn) {
    Ephemeris e;
    e.satellite = SatelliteId(system, prn);
    e.valid = true;
    e.toe = e.toc = GNSSTime(2300, 345600.0);
    e.toes = system == GNSSSystem::BeiDou ? 345586.0 : 345600.0;
    const bool geo = system == GNSSSystem::BeiDou && prn == 1;
    const double a = geo || system == GNSSSystem::QZSS ? 42164000.0 :
        system == GNSSSystem::Galileo ? 29600000.0 :
        system == GNSSSystem::BeiDou ? 27800000.0 : 26560000.0;
    e.sqrt_a = std::sqrt(a); e.e = 0.0012;
    e.i0 = geo ? 0.08 : 0.98;
    e.omega0 = 1.25; e.omega = 0.61; e.m0 = 0.42;
    e.delta_n = 1e-9; e.idot = -2e-10; e.omega_dot = -2.6e-9;
    e.cuc = 1e-6; e.cus = -2e-6; e.crc = 120.0; e.crs = -30.0;
    e.cic = 2e-7; e.cis = -3e-7;
    e.af0 = 2.3e-5; e.af1 = -1.2e-12; e.af2 = 2e-20;
    e.tgd = 2e-9; e.tgd_secondary = -3e-9;
    e.iode = e.iodc = 42; e.data_source_code = system == GNSSSystem::Galileo ? 513 : 0;
    return e;
}

void state(const Ephemeris& e, GNSSTime t, Vector3d& p, double& clock) {
    double drift;
    if (!navigation_internal::computeBroadcastState(e, t, p, clock, drift, true)) {
        throw std::runtime_error("upstream broadcast computation failed");
    }
}

void header(std::ostream& out, std::string content, const char* label) {
    if (content.size() > 60) throw std::runtime_error("header content overflow");
    out << std::left << std::setw(60) << content << label << '\n' << std::right;
}

void write_nav(std::ostream& out, const Ephemeris& e) {
    const bool bds = e.satellite.system == GNSSSystem::BeiDou;
    const std::time_t label = 315964800LL + 2300LL * 604800 + 345600 - (bds ? 14 : 0);
    char date[30]; std::strftime(date,sizeof(date),"%Y %m %d %H %M %S",std::gmtime(&label));
    out << e.satellite.toString() << ' ' << date << std::scientific << std::setprecision(12);
    for (double v : {e.af0,e.af1,e.af2}) out << std::setw(19) << v;
    out << '\n';
    const std::array<std::array<double,4>,7> rows = {{
        {42,e.crs,e.delta_n,e.m0}, {e.cuc,e.e,e.cus,e.sqrt_a},
        {e.toes,e.cic,e.omega0,e.cis}, {e.i0,e.crc,e.omega,e.omega_dot},
        {e.idot,static_cast<double>(e.data_source_code),bds ? 944.0 : 2300.0,0},
        {2.0,0,e.tgd,bds || e.satellite.system == GNSSSystem::Galileo ? e.tgd_secondary : 42.0},
        {e.toes,bds ? 42.0 : 4.0,0,0}
    }};
    for (const auto& row : rows) {
        out << "    "; for(double v : row) out << std::setw(19) << v; out << '\n';
    }
}

int main(int argc, char** argv) {
    if (argc != 2) { std::cerr << "usage: upstream-navigation OUTPUT_DIRECTORY\n"; return 2; }
    const std::string root = argv[1];
    std::ofstream parity(root+"/upstream_navigation.csv"), nav(root+"/upstream_broadcast.nav");
    std::ofstream synthetic_nav(root+"/synthetic_spp.nav"), obs(root+"/synthetic_spp.obs");
    std::ofstream atmosphere_obs(root+"/synthetic_spp_atmosphere.obs");
    std::ofstream atmosphere(root+"/upstream_atmosphere.csv");
    if (!parity || !nav || !synthetic_nav || !obs || !atmosphere || !atmosphere_obs) return 1;
    for (auto* out : {&nav,&synthetic_nav}) {
        header(*out,"     3.04           N                   M","RINEX VERSION / TYPE");
        header(*out,"gnss-rust upstream fixture generator","PGM / RUN BY / DATE");
        header(*out,"","END OF HEADER");
    }
    parity << "# upstream 72f3b7c2c4c088dc575286a6f0edf4e407bc499a\n"
           << "# satellite,offset_s,x,y,z,vx,vy,vz,clock_bias_s,clock_drift\n" << std::setprecision(17);
    for (const auto& identity : std::array<std::pair<GNSSSystem,int>,5>{{
        {GNSSSystem::GPS,1},{GNSSSystem::Galileo,27},{GNSSSystem::QZSS,1},
        {GNSSSystem::BeiDou,1},{GNSSSystem::BeiDou,19}}}) {
        const auto e = make_eph(identity.first,identity.second);
        write_nav(nav,e);
        for (double offset : {-300.0,0.0,120.0,7200.0}) {
            Vector3d p, forward, backward; double clock, cf, cb;
            const auto t = e.toe + offset;
            state(e,t,p,clock); state(e,t+0.5,forward,cf); state(e,t-0.5,backward,cb);
            const Vector3d velocity = forward-backward;
            parity << e.satellite.toString() << ',' << offset;
            for (int i=0;i<3;++i) parity << ',' << p(i);
            for (int i=0;i<3;++i) parity << ',' << velocity(i);
            parity << ',' << clock << ',' << cf-cb << '\n';
        }
    }
    const double rad = std::acos(-1.0)/180;
    const Vector3d receiver = geodetic2ecef(35*rad,139*rad,45);
    std::vector<Ephemeris> visible;
    for (int prn=1;prn<=32;++prn) {
        auto e = make_eph(GNSSSystem::GPS,prn);
        e.omega0 = 0.31*prn; e.m0 = 0.77*prn;
        Vector3d p; double clock; state(e,e.toe,p,clock);
        const auto enu = ecef2enu(p-receiver,35*rad,139*rad);
        const double elevation = std::atan2(enu(2),std::hypot(enu(0),enu(1)));
        if (elevation > 15*rad) { visible.push_back(e); write_nav(synthetic_nav,e); }
    }
    if (visible.size()<6) {std::cerr << "insufficient synthetic geometry\n";return 1;}
    for (auto* output : {&obs, &atmosphere_obs}) {
    header(*output,"     3.04           O                   M","RINEX VERSION / TYPE");
    header(*output,"gnss-rust upstream fixture generator","PGM / RUN BY / DATE");
    {
        std::ostringstream position;
        position << std::fixed << std::setprecision(4) << std::setw(14) << receiver(0)+100
            << std::setw(14) << receiver(1)-50 << std::setw(14) << receiver(2)+30;
        header(*output,position.str(),"APPROX POSITION XYZ");
    }
    header(*output,"G    4 C1C L1C C2W L2W","SYS / # / OBS TYPES");
    header(*output,"","END OF HEADER");
    }
    for (int epoch=0;epoch<8;++epoch) {
        const std::time_t label=315964800LL+2300LL*604800+345600+epoch*30;
        char date[30];std::strftime(date,sizeof(date),"%Y %m %d %H %M",std::gmtime(&label));
        for (bool with_atmosphere : {false, true}) {
        auto& output = with_atmosphere ? atmosphere_obs : obs;
        output << "> " << date << ' ' << std::fixed << std::setprecision(7) << std::setw(10)
            << static_cast<double>(epoch%2*30) << "  0" << std::setw(3) << visible.size() << '\n';
        for (const auto& e:visible) {
            const auto receive=GNSSTime(2300,345600.0+epoch*30);
            double pseudorange=26500000.0, clock, initial_clock; Vector3d p;
            double phase_range=0, phase_secondary=0, pseudorange_secondary=0;
            for (int iteration=0;iteration<12;++iteration) {
                auto tx=receive-pseudorange/constants::SPEED_OF_LIGHT;
                state(e,tx,p,initial_clock);tx=tx-initial_clock;state(e,tx,p,clock);
                const double base=geodist(p,receiver)+2e-6*constants::SPEED_OF_LIGHT-clock*constants::SPEED_OF_LIGHT;
                double trop=0.0, iono=0.0;
                if (with_atmosphere) {
                    const auto enu=ecef2enu(p-receiver,35*rad,139*rad);
                    const double elevation=std::atan2(enu(2),std::hypot(enu(0),enu(1)));
                    const double azimuth=std::atan2(enu(0),enu(1));
                    trop=models::tropDelaySaastamoinen(receiver,elevation);
                    iono=models::ionoDelayKlobuchar(35*rad,139*rad,azimuth,elevation,receive.tow,nullptr,nullptr);
                }
                const double gamma=std::pow(constants::GPS_L1_FREQ/constants::GPS_L2_FREQ,2);
                pseudorange=base+e.tgd*constants::SPEED_OF_LIGHT+trop+iono;
                phase_range=base+trop-iono;
                pseudorange_secondary=base+gamma*e.tgd*constants::SPEED_OF_LIGHT+trop+gamma*iono;
                phase_secondary=base+trop-gamma*iono;
            }
            output << e.satellite.toString() << std::fixed << std::setprecision(3)
                << std::setw(14) << pseudorange << "  "
                << std::setw(14) << phase_range/constants::GPS_L1_WAVELENGTH+10000+e.satellite.prn << "  "
                << std::setw(14) << pseudorange_secondary << "  "
                << std::setw(14) << phase_secondary/constants::GPS_L2_WAVELENGTH+20000+e.satellite.prn << "  " << '\n';
        }
        }
    }
    atmosphere << "# latitude,longitude,height,azimuth,elevation,tow,klobuchar_m,saastamoinen_m\n" << std::setprecision(17);
    for (const auto& sample:std::array<std::array<double,6>,4>{{
        {35,139,45,120,30,16000},{-33,151,500,30,70,45000},
        {0,-120,0,250,10,600000},{80,0,100,0,90,0}}}) {
        const double lat=sample[0]*rad,lon=sample[1]*rad,az=sample[3]*rad,el=sample[4]*rad;
        atmosphere << lat << ',' << lon << ',' << sample[2] << ',' << az << ',' << el << ',' << sample[5]
            << ',' << models::ionoDelayKlobuchar(lat,lon,az,el,sample[5],nullptr,nullptr)
            << ',' << models::tropDelaySaastamoinen(geodetic2ecef(lat,lon,sample[2]),el) << '\n';
    }
    std::cerr << "generated 20 broadcast cases, 4 atmosphere cases, 8 epochs with " << visible.size() << " GPS satellites\n";
    return 0;
}
