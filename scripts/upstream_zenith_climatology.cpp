// Actual pinned inline native compact climatology, without altered coefficients.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/models/troposphere.hpp>
#include <fstream>
#include <iomanip>
int main(int argc,char**argv) {
    if(argc!=2)return 2;
    std::ofstream out(std::string(argv[1])+"/upstream_zenith_climatology.csv");
    out << "latitude_deg,height_m,day,pressure_hpa,temperature_k,vapor_hpa,hydrostatic_m,wet_m,total_m\n" << std::setprecision(17);
    for(double latitude:{-90.,-75.,-60.,-45.,-30.,-15.,0.,15.,30.,45.,60.,75.,90.})
        for(double height:{-500.,0.,45.,10000.})for(int day:{1,28,91,182,216,366}) {
            const auto delay=libgnss::models::estimateZenithTroposphereClimatology(latitude*std::acos(-1.)/180.,height,day);
            out << latitude << ',' << height << ',' << day << ',' << delay.pressure_hpa << ',' << delay.temperature_k << ','
                << delay.water_vapor_pressure_hpa << ',' << delay.hydrostatic_delay_m << ',' << delay.wet_delay_m << ',' << delay.totalDelayMeters() << '\n';
        }
    return out?0:3;
}
