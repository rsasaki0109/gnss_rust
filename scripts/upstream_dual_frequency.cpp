// Independent dual-frequency reference using unmodified native PPP utilities.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/algorithms/ppp.hpp>
#include <libgnss++/core/constants.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace libgnss;
namespace c = libgnss::constants;
struct Pair { const char* satellite; const char* first; const char* second; int channel; double f1; double f2; };
int main(int argc, char** argv) {
    if (argc != 3) return 2;
    std::ofstream out(argv[1]);
    std::ofstream obs(argv[2]);
    if (!out || !obs) return 3;
    const std::array<Pair, 12> pairs{{
        {"G01","1C","2W",99,c::GPS_L1_FREQ,c::GPS_L2_FREQ},
        {"G02","1W","5Q",99,c::GPS_L1_FREQ,c::GPS_L5_FREQ},
        {"G03","2L","5I",99,c::GPS_L2_FREQ,c::GPS_L5_FREQ},
        {"J01","1C","2X",99,c::GPS_L1_FREQ,c::GPS_L2_FREQ},
        {"J02","1X","5Q",99,c::GPS_L1_FREQ,c::GPS_L5_FREQ},
        {"E12","1C","5Q",99,c::GAL_E1_FREQ,c::GAL_E5A_FREQ},
        {"E13","1X","7Q",99,c::GAL_E1_FREQ,c::GAL_E5B_FREQ},
        {"C06","2I","7I",99,c::BDS_B1I_FREQ,c::BDS_B2I_FREQ},
        {"C19","2I","6I",99,c::BDS_B1I_FREQ,c::BDS_B3I_FREQ},
        {"C20","1P","5P",99,c::BDS_B1C_FREQ,c::BDS_B2A_FREQ},
        {"R01","1C","2C",-7,c::GLO_L1_BASE_FREQ-7*c::GLO_L1_STEP_FREQ,c::GLO_L2_BASE_FREQ-7*c::GLO_L2_STEP_FREQ},
        {"R02","1P","2P",6,c::GLO_L1_BASE_FREQ+6*c::GLO_L1_STEP_FREQ,c::GLO_L2_BASE_FREQ+6*c::GLO_L2_STEP_FREQ}
    }};
    auto header = [&](const std::string& value, const std::string& label) {
        obs << std::left << std::setw(60) << value << label << '\n' << std::right;
    };
    header("     3.04           O                   M", "RINEX VERSION / TYPE");
    header("gnss-rust native dual-frequency synthetic fixture", "PGM / RUN BY / DATE");
    header("    18", "LEAP SECONDS");
    std::map<char,std::vector<std::string>> types;
    for (const auto& pair : pairs) {
        auto& codes = types[pair.satellite[0]];
        for (const auto* tracking : {pair.first,pair.second}) {
            for (char kind : {'C','L'}) {
                const std::string code = std::string(1,kind)+tracking;
                if (std::find(codes.begin(),codes.end(),code)==codes.end()) codes.push_back(code);
            }
        }
    }
    for (const auto& [system,codes] : types) {
        std::ostringstream line;
        line << system << "  " << std::setw(3) << codes.size() << ' ';
        for (const auto& code : codes) line << code << ' ';
        header(line.str(), "SYS / # / OBS TYPES");
    }
    header("  2 R01 -7 R02  6", "GLONASS SLOT / FRQ #");
    header("", "END OF HEADER");
    out << "sat,first,second,channel,case,f1,f2,p1,p2,l1_cycles,l2_cycles,a,b,p_if,l_if,gf,mw,code_variance,phase_variance,ambiguity_if_m\n";
    out << std::setprecision(17);
    for (const auto& pair : pairs) {
        const auto [a,b] = ppp_utils::getIonosphereFreeCoefficients(pair.f1,pair.f2);
        const double lam1=c::SPEED_OF_LIGHT/pair.f1, lam2=c::SPEED_OF_LIGHT/pair.f2;
        for (int i=0;i<4;++i) {
            const double range=23000000.125+i*17.5;
            const double ion=i*4.25;
            const double ion2=ion*pair.f1*pair.f1/(pair.f2*pair.f2);
            const double n1=12+i, n2=-7+i*2;
            const double p1=range+ion, p2=range+ion2;
            const double l1=(range-ion)/lam1+n1, l2=(range-ion2)/lam2+n2;
            const double lm1=l1*lam1,lm2=l2*lam2;
            out << pair.satellite << ',' << pair.first << ',' << pair.second << ',' << pair.channel << ',' << i << ','
                << pair.f1 << ',' << pair.f2 << ',' << p1 << ',' << p2 << ',' << l1 << ',' << l2 << ','
                << a << ',' << b << ',' << a*p1+b*p2 << ',' << a*lm1+b*lm2 << ',' << lm1-lm2 << ','
                << ppp_utils::calculateMelbourneWubbena(l1,l2,p1,p2,pair.f1,pair.f2) << ','
                << (a*a+b*b)*0.25 << ',' << (a*a+b*b)*0.0001 << ',' << a*lam1*n1+b*lam2*n2 << '\n';
        }
    }
    for (int i=0;i<4;++i) {
        obs << "> 2024 02 08 00 " << std::setfill('0') << std::setw(2) << (i*30)/60
            << ' ' << std::setfill(' ') << std::fixed << std::setprecision(7)
            << std::setw(10) << double((i*30)%60) << "  0 12\n";
        for (const auto& pair : pairs) {
            const double range=23000000.125+i*17.5, ion=i*4.25;
            const double ion2=ion*pair.f1*pair.f1/(pair.f2*pair.f2);
            const double p1=range+ion,p2=range+ion2;
            const double l1=(range-ion)*pair.f1/c::SPEED_OF_LIGHT+12+i;
            const double l2=(range-ion2)*pair.f2/c::SPEED_OF_LIGHT-7+i*2;
            const auto& codes=types.at(pair.satellite[0]);
            obs << pair.satellite;
            for (size_t j=0;j<codes.size();++j) {
                if (j && j%4==0) obs << "\n   ";
                const auto& code=codes[j];
                const bool first=code.substr(1)==pair.first,second=code.substr(1)==pair.second;
                if (!first && !second) obs << std::string(16,' ');
                else if (code[0]=='C') obs << std::setprecision(5) << std::setw(14) << (first?p1:p2) << "  ";
                else obs << std::setprecision(4) << std::setw(14) << (first?l1:l2) << "07";
            }
            obs << '\n';
        }
    }
    return out && obs ? 0 : 4;
}
