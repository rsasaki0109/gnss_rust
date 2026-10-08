// Generate reference data using the actual, unmodified libgnss++ header.
// See tests/fixtures/README.md for pinned revisions and reproduction commands.
#include <libgnss++/core/coordinates.hpp>
#include <array>
#include <iomanip>
#include <iostream>

int main() {
    const double rad = std::acos(-1.0) / 180.0;
    const std::array<std::array<double, 3>, 8> samples = {{
        {0.0, 0.0, 0.0}, {0.0, 90.0, 100.0},
        {35.0, 139.0, 45.0}, {-33.8, 151.2, 500.0},
        {89.999999, -170.0, -100.0}, {90.0, 0.0, 0.0},
        {-90.0, 180.0, 20.0}, {10.0, -170.0, 20200000.0}
    }};
    const Eigen::Vector3d difference(100.0, -200.0, 300.0);
    const Eigen::Vector3d satellite(20200000.0, 14000000.0, 21700000.0);
    std::cout << "# upstream 72f3b7c2c4c088dc575286a6f0edf4e407bc499a; Eigen 3.4.0\n";
    std::cout << "# lat,lon,h,x,y,z,inverse_lat,inverse_lon,inverse_h,east,north,up,dx,dy,dz,range\n";
    std::cout << std::setprecision(17);
    for (const auto& sample : samples) {
        const double lat = sample[0] * rad, lon = sample[1] * rad, h = sample[2];
        const auto ecef = libgnss::geodetic2ecef(lat, lon, h);
        double inverse_lat, inverse_lon, inverse_h;
        libgnss::ecef2geodetic(ecef, inverse_lat, inverse_lon, inverse_h);
        const auto enu = libgnss::ecef2enu(difference, lat, lon);
        const auto delta = libgnss::enu2ecef(difference, lat, lon);
        std::cout << lat << ',' << lon << ',' << h;
        for (int i = 0; i < 3; ++i) std::cout << ',' << ecef(i);
        std::cout << ',' << inverse_lat << ',' << inverse_lon << ',' << inverse_h;
        for (int i = 0; i < 3; ++i) std::cout << ',' << enu(i);
        for (int i = 0; i < 3; ++i) std::cout << ',' << delta(i);
        std::cout << ',' << libgnss::geodist(satellite, ecef) << '\n';
    }
}
