// Call pinned native common-clock DOP; Q is an independent Eigen normal inverse.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/core/observation.hpp>
#include <libgnss++/core/coordinates.hpp>
#include <libgnss++/algorithms/spp.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
    using namespace libgnss;
    const double radians = std::acos(-1.0) / 180;
    std::cout << "# pinned SPP default max_gdop=" << SPPProcessor::SPPConfig{}.max_gdop << '\n';
    std::cout << "# case,n,receiver[3],satellite_positions[3*n],native_gdop,pdop,tdop,Eigen_Q[16]\n";
    std::cout << std::setprecision(17);
    ObservationData observation;
    for (int frame : {0, 1}) for (int count : {4, 5, 7, 9}) {
        const double lat = (frame == 0 ? 35 : -33) * radians;
        const double lon = (frame == 0 ? 139 : -71) * radians;
        const Vector3d receiver = geodetic2ecef(lat, lon, 45 + frame * 200);
        for (bool weak : {false, true}) {
            std::vector<Vector3d> satellites;
            for (int index = 0; index < count; ++index) {
                const double azimuth = (13 + 360.0 * index / count) * radians;
                const double elevation = (weak ? 20 + 0.1 * index
                                                : 18 + std::fmod(31.0 * index, 65.0)) * radians;
                Vector3d enu(std::cos(elevation) * std::sin(azimuth),
                             std::cos(elevation) * std::cos(azimuth), std::sin(elevation));
                satellites.push_back(receiver + enu2ecef(enu * (21000000 + 10000 * index), lat, lon));
            }
            const auto dop = observation.calculateDOP(receiver, satellites);
            const MatrixXd h = observation.calculateGeometryMatrix(receiver, satellites);
            const MatrixXd covariance = (h.transpose() * h).inverse();
            std::cout << (weak ? "weak" : "ordinary") << frame << '_' << count << ',' << count;
            for (int i = 0; i < 3; ++i) std::cout << ',' << receiver(i);
            for (const auto& satellite : satellites)
                for (int i = 0; i < 3; ++i) std::cout << ',' << satellite(i);
            std::cout << ',' << dop.gdop << ',' << dop.pdop << ',' << dop.tdop;
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) std::cout << ',' << covariance(i, j);
            std::cout << '\n';
        }
    }
}
