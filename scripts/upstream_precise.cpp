// Independent native SP3/CLK/interpolation reference generator.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/core/navigation.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace libgnss;

static std::string epoch(int seconds) {
    std::ostringstream out;
    out << "2024 02 08 " << std::setfill('0') << std::setw(2) << seconds / 3600
        << ' ' << std::setw(2) << (seconds / 60) % 60 << ' ' << std::fixed
        << std::setprecision(8) << std::setfill(' ') << std::setw(11) << double(seconds % 60);
    return out.str();
}
static void clock_header(std::ostream& out, const std::string& value, const std::string& label) {
    out << std::left << std::setw(60) << value << label << '\n' << std::right;
}
static void write_sp3(const std::string& path, int count) {
    std::ofstream out(path);
    out << "#dP" << epoch(0) << ' ' << std::setw(7) << count << " ORBIT IGS20 FIT TEST\n";
    out << "## 2300 " << std::fixed << std::setprecision(8) << std::setw(15) << 345600.0
        << ' ' << std::setw(14) << 300.0 << " 60348 0.0000000000000\n";
    out << "+    2   G01E12" << std::string(45,' ') << '\n';
    out << "++         5  6" << std::string(45,' ') << '\n';
    out << "%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n";
    out << "%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n";
    out << "%f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n";
    out << "%f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n";
    out << "%i    0    0    0    0      0      0      0      0         0\n";
    out << "/* Synthetic moving tracks; ANTENNA PHASE CENTER\n";
    for(int i=0;i<count;++i) {
        const double t=300.0*i;
        out << "*  " << epoch(int(t)) << '\n';
        for(int s=0;s<2;++s) {
            // Curved polynomial track with a degree-nine term exercises the
            // same native Neville window and central difference near boundaries.
            const double u=t/3600.0;
            const double x=20200000.0+2500*t-0.03*t*t+20*std::pow(u,9)+s*100000.0;
            const double y=-14100000.0+1300*t+0.02*t*t-15*std::pow(u,7)-s*80000.0;
            const double z=21700000.0-900*t-0.01*t*t+9*std::pow(u,5)+s*50000.0;
            const double clk=(s?0.0:20.0+0.00002*t+0.00000001*t*t);
            out << 'P' << (s?"E12":"G01") << std::fixed << std::setprecision(6)
                << std::setw(14) << x/1000.0 << std::setw(14) << y/1000.0
                << std::setw(14) << z/1000.0 << std::setw(14) << clk
                << "  5  6  7   8       " << '\n';
        }
    }
    out << "EOF\n";
}
static void write_clk(const std::string& path) {
    std::ofstream out(path);
    clock_header(out,"     3.00           C", "RINEX VERSION / TYPE");
    clock_header(out,"GPS", "TIME SYSTEM ID");
    clock_header(out,"     1    AS", "# / TYPES OF DATA");
    clock_header(out,"Synthetic clock datum deliberately differs from SP3", "COMMENT");
    clock_header(out,"", "END OF HEADER");
    for(int i=0;i<40;++i) {
        const double t=i*90.0;
        const double clk=5e-5+3e-11*t-1e-15*t*t+3e-10*std::pow(t/3600.0,5);
        out << "AS G01 " << epoch(int(t)) << "  2 " << std::scientific << std::setprecision(12)
            << clk << " 1.000000000000E-12\n";
    }
}
int main(int argc,char** argv) {
    if(argc!=2) { std::cerr << "usage: upstream-precise OUTPUT_DIRECTORY\n";return 2; }
    const std::string root=argv[1];
    write_sp3(root+"/synthetic_precise.sp3",12);
    write_sp3(root+"/synthetic_precise_two.sp3",2);
    write_clk(root+"/synthetic_precise.clk");
    std::ofstream out(root+"/upstream_precise.csv");
    out << "# mode,satellite,week,tow,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,clock_available,bias_s,drift_s_s,relativity_s\n" << std::setprecision(17);
    int cases=0;
    for(int mode=0;mode<4;++mode) {
        PreciseProducts products;
        if(!products.loadSP3File(root+(mode==2?"/synthetic_precise_two.sp3":"/synthetic_precise.sp3"))) return 1;
        // Mode 1 compares the native interpolator on the same CLK-only grid
        // selected by Rust. Mode 3 records the upstream default mixed grid
        // separately, to expose its intentional source-selection difference.
        if(mode==1) for(auto& [sat, entries]:products.orbit_clock_data) {
            if(sat.system==GNSSSystem::GPS && sat.prn==1)
                for(auto& entry:entries) entry.clock_valid=false;
        }
        if((mode==1 || mode==3) && !products.loadClockFile(root+"/synthetic_precise.clk")) return 1;
        for(int s=0;s<2;++s) for(double offset:{0.0,0.125,37.5,150.0,300.0,450.0,1350.125,2850.0,3300.0}) {
            if(mode==2 && offset>300.0) continue;
            const SatelliteId satellite(s?GNSSSystem::Galileo:GNSSSystem::GPS,s?12:1);
            const GNSSTime time(2300,345600.0+offset);
            Vector3d p,v;double c=0.0,d=0.0;bool clock=false;
            if(!products.interpolateOrbitClock(satellite,time,p,v,c,d,&clock)) return 1;
            out << mode << ',' << (s?"E12":"G01") << ',' << 2300 << ',' << time.tow;
            for(int axis=0;axis<3;++axis) out << ',' << p(axis);
            for(int axis=0;axis<3;++axis) out << ',' << v(axis);
            out << ',' << clock << ',' << c << ',' << d << ',' << preciseClockRelativisticCorrection(p,v) << '\n';
            ++cases;
        }
    }
    std::cerr << "generated SP3/CLK tracks and " << cases << " native interpolation cases\n";
    return 0;
}
