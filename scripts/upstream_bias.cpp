// Native product-loader and literal PPP conversion/sign reference, not PPPProcessor.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/core/navigation.hpp>
#include <libgnss++/algorithms/ppp_correction_contract.hpp>
#include <libgnss++/core/constants.hpp>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include <algorithm>
using namespace libgnss;

// Literal conversion from ppp_corrections.cpp, anonymous namespace. The native
// product loader below is compiled unmodified. Phase cycles are a Rust extension.
double native_m(const DCBEntry& e) {
    std::string unit=e.unit;
    std::transform(unit.begin(),unit.end(),unit.begin(),[](unsigned char c){return std::toupper(c);});
    if(unit=="NS") return e.bias*constants::SPEED_OF_LIGHT*1e-9;
    if(unit=="M"||unit=="METER"||unit=="METERS") return e.bias;
    return std::numeric_limits<double>::quiet_NaN();
}
void header(std::ostream& out,int count) {
    out << "%=BIA 1.00 SYN 2024:039:00000 SYN 2024:039:00000 2024:040:00000 A "
        << std::setfill('0') << std::setw(8) << count << std::setfill(' ') << '\n'
        << "+BIAS/DESCRIPTION\n TIME_SYSTEM G\n BIAS_MODE ABSOLUTE\n-BIAS/DESCRIPTION\n"
        << "+BIAS/SOLUTION\n*BIAS SVN_ PRN STATION__ OBS1 OBS2 BEGIN__________ END____________ UNIT VALUE STD_DEV\n";
}
int main(int argc,char**argv) {
    if(argc!=2) return 2;
    const std::string root=argv[1];
    std::ofstream standard(root+"/synthetic_ppp_osb.bsx"),compatible(root+"/upstream_bias.bsx");
    header(standard,36);header(compatible,39);
    // The same satellite list as the independently generated nine-SV PPP data.
    for(int prn:{1,2,7,8,9,18,19,25,26}) {
        std::ostringstream id,svn;id << 'G' << std::setfill('0') << std::setw(2) << prn;
        svn << 'G' << std::setfill('0') << std::setw(3) << 80+prn;
        for(int j=0;j<2;++j) for(bool phase:{false,true}) {
            const std::string obs=std::string(phase?"L":"C")+(j==0?"1C":"2W");
            const double metres=phase?(j==0?0.015*(prn%3)-0.007:0.012*(prn%4)-0.020)
                :(j==0?3.2*(prn%3)-2.0:-1.1*(prn%4)+0.7);
            // Mixed ns/metre units; both are supported by native conversion.
            const bool ns=j==0;const std::string unit=ns?"ns":"m";
            const double value=ns?metres/constants::SPEED_OF_LIGHT*1e9:metres;
            for(bool native:{false,true}) {
                auto& out=native?compatible:standard;
                // Native whitespace loader requires OBS2, even for OSB. Its
                // repeated observable is omitted in the standard-style file.
                out << " OSB " << svn.str() << ' ' << id.str() << "           " << obs
                    << ' ' << (native?obs:"   ") << " 2024:039:00600 2024:039:02520 " << unit
                    << ' ' << std::setprecision(17) << value << ' ' << (phase?1e-5:1e-3) << '\n';
            }
        }
    }
    // Unselected native-readable DSB, alternate system, open interval and exponent.
    compatible << " DSB E215 E11 C1C C5Q 2024:039:00000 2024:040:00000 ns 1.234 0.1\n"
        << " OSB C201 C01 C2I C2I 0000:000:00000 0000:000:00000 meters -2.5 0.2\n"
        << " OSB J001 J01 C1C C1C 2024:039:00000 2024:040:00000 m 2.0E-3 0.01\n";
    standard << "-BIAS/SOLUTION\n%=ENDBIA\n";compatible << "-BIAS/SOLUTION\n%=ENDBIA\n";
    standard.close();compatible.close();
    DCBProducts products;
    if(!products.loadFile(root+"/upstream_bias.bsx") || products.entries.size()!=39) return 3;
    std::ofstream reference(root+"/upstream_bias.csv");
    reference << "kind,satellite,observation1,observation2,unit,value,sigma,metres,code_add_m\n" << std::setprecision(17);
    for(const auto&e:products.entries) {
        double value=0,sigma=0;
        if(!products.getBias(e.satellite,e.bias_type,e.observation_1,e.observation_2,value,&sigma)) return 4;
        const double m=native_m(e);
        reference << e.bias_type << ',' << e.satellite.toString() << ',' << e.observation_1 << ',' << e.observation_2
            << ',' << e.unit << ',' << value << ',' << sigma << ',' << m << ','
            << algorithms::ppp_correction_contract::measurementCorrectionSign(false,false)*m << '\n';
    }
    return reference ? 0 : 5;
}
