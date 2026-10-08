// Precise light-time component and code-SPP fixture generator.
// Uses unmodified native interpolation/geometry/atmosphere, plus the literal
// receive-time Taylor statements from upstream spp.cpp/ppp_corrections.cpp.
// This does not execute the full native SPPProcessor.
#define main precise_products_fixture_main
#include "upstream_precise.cpp"
#undef main
#define main navigation_fixture_main
#include "upstream_navigation.cpp"
#undef main

struct Tx {
    GNSSTime emission;
    double tau,published,drift,periodic,total;
    Vector3d p,v;
    int iterations;
};
static Tx transmit(const PreciseProducts& products,SatelliteId sat,GNSSTime receive,const Vector3d& receiver,bool emission_method) {
    Tx r;bool available=false;double bias;
    if(!products.interpolateOrbitClock(sat,receive,r.p,r.v,bias,r.drift,&available) || !available) throw std::runtime_error("receive state unavailable");
    r.tau=(r.p-receiver).norm()/constants::SPEED_OF_LIGHT;
    r.emission=receive-r.tau;r.iterations=1;
    if(emission_method) {
        bool converged=false;
        for(int i=1;i<=8;++i) {
            if(!products.interpolateOrbitClock(sat,r.emission,r.p,r.v,bias,r.drift,&available) || !available) throw std::runtime_error("emission unavailable");
            const double next=(r.p-receiver).norm()/constants::SPEED_OF_LIGHT;
            r.iterations=i;
            if(std::abs(next-r.tau)<=1e-11) {converged=true;break;}
            r.tau=next;r.emission=receive-r.tau;
        }
        if(!converged) throw std::runtime_error("emission did not converge");
        r.published=bias;
    } else {
        r.p-=r.v*r.tau;
        r.published=bias-r.drift*r.tau;
    }
    r.periodic=preciseClockRelativisticCorrection(r.p,r.v);
    r.total=r.published+r.periodic;
    return r;
}
static std::vector<Ephemeris> visible(const Vector3d& receiver) {
    std::vector<Ephemeris> result;
    const double rad=std::acos(-1.0)/180.0;
    for(int prn=1;prn<=32;++prn) {
        auto e=make_eph(GNSSSystem::GPS,prn);e.omega0=0.31*prn;e.m0=0.77*prn;
        Vector3d p;double clock;state(e,e.toe,p,clock);
        const auto enu=ecef2enu(p-receiver,35*rad,139*rad);
        if(std::atan2(enu(2),std::hypot(enu(0),enu(1)))>15*rad) result.push_back(e);
    }
    return result;
}
static void spp_products(const std::string& root,const std::vector<Ephemeris>& satellites) {
    std::ofstream orbit(root+"/synthetic_precise_spp.sp3"),clock(root+"/synthetic_precise_spp.clk");
    orbit << "#dP" << epoch(0) << ' ' << std::setw(7) << 12 << " ORBIT IGS20 FIT TEST\n";
    orbit << "## 2300 " << std::fixed << std::setprecision(8) << std::setw(15) << 345600.0 << ' ' << std::setw(14) << 300.0 << " 60348 0.0000000000000\n";
    orbit << "+  " << std::setw(3) << satellites.size() << "   ";
    for(const auto& e:satellites) orbit << e.satellite.toString();
    orbit << std::string(3*(17-satellites.size()),' ') << '\n';
    orbit << "++       ";for(int i=0;i<17;++i) orbit << std::setw(3) << (i<int(satellites.size())?5:0);orbit << '\n';
    orbit << "%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n";
    orbit << "%f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n";
    orbit << "%i    0    0    0    0      0      0      0      0         0\n";
    orbit << "/* Synthetic code model positions: ANTENNA PHASE CENTER\n";
    for(int i=0;i<12;++i) {
        orbit << "*  " << epoch(300*i) << '\n';
        for(const auto& e:satellites) {
            Vector3d p;double unused;state(e,GNSSTime(2300,345600.0+300*i),p,unused);
            orbit << 'P' << e.satellite.toString() << std::fixed << std::setprecision(6);
            for(int j=0;j<3;++j) orbit << std::setw(14) << p(j)/1000.0;
            orbit << std::setw(14) << 999999.999999 << '\n';
        }
    }
    orbit << "EOF\n";
    clock_header(clock,"     3.00           C","RINEX VERSION / TYPE");
    clock_header(clock,"GPS","TIME SYSTEM ID");clock_header(clock,"","END OF HEADER");
    for(int i=0;i<40;++i) for(const auto& e:satellites) {
        const double t=90.0*i;
        const double bias=5e-5+e.satellite.prn*1e-7+2e-11*t+3e-15*t*t;
        clock << "AS " << e.satellite.toString() << ' ' << epoch(int(t)) << " 2 " << std::scientific << std::setprecision(12) << bias << " 1.000000000000E-12\n";
    }
}
#ifndef GNSS_RUST_PRECISE_TRANSMIT_MAIN
#define GNSS_RUST_PRECISE_TRANSMIT_MAIN main
#endif
int GNSS_RUST_PRECISE_TRANSMIT_MAIN(int argc,char** argv) {
    if(argc!=2) {std::cerr << "usage: upstream-precise-transmit OUTPUT_DIRECTORY\n";return 2;}
    if(precise_products_fixture_main(argc,argv)!=0) return 1;
    const std::string root=argv[1];
    const double rad=std::acos(-1.0)/180.0;
    const Vector3d receiver=geodetic2ecef(35*rad,139*rad,45.0),seed=receiver+Vector3d(100,-50,30);
    std::ofstream out(root+"/upstream_precise_transmit.csv");
    out << "# grid,method,satellite,week,receive_tow,receiverXYZ,emit_tow,tau,iterations,positionXYZ,velocityXYZ,published_s,drift_s_s,relativity_s,total_s,geodist_m\n" << std::setprecision(17);
    int cases=0;
    for(int grid=0;grid<2;++grid) {
        PreciseProducts products;
        if(!products.loadSP3File(root+"/synthetic_precise.sp3")) return 1;
        if(grid==1) {
            for(auto& [sat,entries]:products.orbit_clock_data) if(sat.system==GNSSSystem::GPS)
                for(auto& e:entries) e.clock_valid=false;
            if(!products.loadClockFile(root+"/synthetic_precise.clk")) return 1;
        }
        for(bool exact:{false,true}) for(int s=0;s<2;++s) for(double t:{30.0,150.0,900.0,1350.125,3150.0}) for(const Vector3d& rx:{receiver,seed}) {
            const SatelliteId sat(s?GNSSSystem::Galileo:GNSSSystem::GPS,s?12:1);
            const auto r=transmit(products,sat,GNSSTime(2300,345600.0+t),rx,exact);
            out << grid << ',' << exact << ',' << sat.toString() << ",2300," << 345600.0+t;
            for(int j=0;j<3;++j) out << ',' << rx(j);
            out << ',' << r.emission.tow << ',' << r.tau << ',' << r.iterations;
            for(int j=0;j<3;++j) out << ',' << r.p(j);
            for(int j=0;j<3;++j) out << ',' << r.v(j);
            out << ',' << r.published << ',' << r.drift << ',' << r.periodic << ',' << r.total << ',' << geodist(r.p,rx) << '\n';
            ++cases;
        }
    }
    const auto satellites=visible(receiver);
    if(satellites.size()<6 || satellites.size()>17) return 1;
    spp_products(root,satellites);
    PreciseProducts spp;
    if(!spp.loadSP3File(root+"/synthetic_precise_spp.sp3") || !spp.loadClockFile(root+"/synthetic_precise_spp.clk")) return 1;
    std::ofstream truth(root+"/synthetic_precise_spp_truth.csv"),biases(root+"/synthetic_precise_spp_bias.csv");
    truth << "# receiverXYZ,initialXYZ,receiver_clock_s\n" << std::setprecision(17);
    for(const auto& p:{receiver,seed}) for(int j=0;j<3;++j) truth << p(j) << ',';
    truth << 2e-6 << '\n';
    biases << "# satellite,tracking,source,from_week,from_tow,until_week,until_tow,correction_m,clock_reference\n" << std::setprecision(17);
    for(const auto& e:satellites) biases << e.satellite.toString() << ",1C,RinexClk,2300,346200,2300,346410," << -(3.2*(e.satellite.prn%3)-2.0) << ",synthetic-precise-clock-v1\n";
    for(bool exact:{false,true}) {
        std::ofstream obs(root+(exact?"/synthetic_precise_spp_emission.obs":"/synthetic_precise_spp_taylor.obs"));
        header(obs,"     3.04           O                   G","RINEX VERSION / TYPE");
        std::ostringstream pos;pos << std::fixed << std::setprecision(4);for(int j=0;j<3;++j) pos << std::setw(14) << seed(j);
        header(obs,pos.str(),"APPROX POSITION XYZ");header(obs,"G    1 C1C","SYS / # / OBS TYPES");header(obs,"","END OF HEADER");
        for(int i=0;i<8;++i) {
            const int t=600+30*i;
            obs << "> " << epoch(t) << "  0" << std::setw(3) << satellites.size() << '\n';
            for(const auto& e:satellites) {
                const auto tx=transmit(spp,e.satellite,GNSSTime(2300,345600.0+t),receiver,exact);
                const auto enu=ecef2enu(tx.p-receiver,35*rad,139*rad);
                const double elevation=std::atan2(enu(2),std::hypot(enu(0),enu(1)));
                const double raw=geodist(tx.p,receiver)+2e-6*constants::SPEED_OF_LIGHT-tx.total*constants::SPEED_OF_LIGHT
                    +models::tropDelaySaastamoinen(receiver,elevation)+(3.2*(e.satellite.prn%3)-2.0);
                obs << e.satellite.toString() << std::fixed << std::setprecision(5) << std::setw(14) << raw << "  \n";
            }
        }
    }
    std::cerr << "generated " << cases << " precise transmit cases and 2 x 8 code-SPP epochs, " << satellites.size() << " satellites\n";
    return 0;
}
