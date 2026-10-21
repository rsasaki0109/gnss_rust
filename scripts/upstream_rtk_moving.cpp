// Independent moving receiver observations from pinned libgnss++ geometry.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#define GNSS_RUST_RTK_FLOAT_HELPERS_ONLY
#include "upstream_rtk_float.cpp"

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const std::string root = argv[1];
    const double rad = std::acos(-1.0) / 180;
    const Vector3d base = geodetic2ecef(35 * rad, 139 * rad, 45);
    const Vector3d origin = base + enu2ecef(Vector3d(12, -7, 3), 35 * rad, 139 * rad);
    const Vector3d seed = origin + Vector3d(1, -0.5, 0.3);
    std::ofstream truth(root + "/synthetic_rtk_moving_truth.csv");
    std::ofstream base_obs(root + "/synthetic_rtk_moving_base.obs");
    std::ofstream rover_obs(root + "/synthetic_rtk_moving_rover.obs");
    std::ofstream noisy_base(root + "/synthetic_rtk_moving_base_noisy.obs");
    std::ofstream noisy_rover(root + "/synthetic_rtk_moving_rover_noisy.obs");
    if (!truth || !base_obs || !rover_obs || !noisy_base || !noisy_rover) return 1;
    write_observation_header(base_obs, base);
    write_observation_header(noisy_base, base);
    write_observation_header(rover_obs, seed);
    write_observation_header(noisy_rover, seed);
    truth << "week,tow,x_m,y_m,z_m\n" << std::setprecision(17);
    std::vector<Ephemeris> visible;
    for (int prn = 1; prn <= 32; ++prn) {
        auto e = make_eph(GNSSSystem::GPS, prn);
        e.omega0 = 0.31 * prn; e.m0 = 0.77 * prn;
        Vector3d p; double clock; state(e, e.toe, p, clock);
        const auto enu = ecef2enu(p - base, 35 * rad, 139 * rad);
        if (std::atan2(enu(2), std::hypot(enu(0), enu(1))) > 15 * rad) visible.push_back(e);
    }
    for (int epoch = 0; epoch < 64; ++epoch) {
        // A known ENU trajectory; no position from a Rust solver is used.
        const Vector3d rover = origin + enu2ecef(
            Vector3d(1.2 * epoch, -0.8 * epoch, 0.1 * epoch), 35 * rad, 139 * rad);
        truth << 2300 << ',' << 345600 + epoch;
        for (int i = 0; i < 3; ++i) truth << ',' << rover(i);
        truth << '\n';
        const std::time_t label = 315964800LL + 2300LL * 604800 + 345600 + epoch;
        char date[30]; std::strftime(date, sizeof(date), "%Y %m %d %H %M", std::gmtime(&label));
        for (bool is_rover : {false, true}) for (bool noisy : {false, true}) {
            auto& out = is_rover ? (noisy ? noisy_rover : rover_obs) : (noisy ? noisy_base : base_obs);
            const Vector3d receiver = is_rover ? rover : base;
            double lat, lon, height; ecef2geodetic(receiver, lat, lon, height);
            out << "> " << date << ' ' << std::fixed << std::setprecision(7)
                << std::setw(10) << static_cast<double>(epoch % 60)
                << "  0" << std::setw(3) << visible.size() << '\n';
            for (const auto& e : visible) {
                const auto receive = GNSSTime(2300, 345600 + epoch);
                double code = 26500000, phase = 0, clock, initial_clock; Vector3d p;
                for (int iteration = 0; iteration < 12; ++iteration) {
                    auto tx = receive - code / constants::SPEED_OF_LIGHT;
                    state(e, tx, p, initial_clock); tx = tx - initial_clock; state(e, tx, p, clock);
                    const auto enu = ecef2enu(p - receiver, lat, lon);
                    const double el = std::atan2(enu(2), std::hypot(enu(0), enu(1)));
                    const double range = geodist(p, receiver)
                        + (is_rover ? 2e-6 : -1e-6) * constants::SPEED_OF_LIGHT
                        - clock * constants::SPEED_OF_LIGHT
                        + models::tropDelaySaastamoinen(receiver, el);
                    code = range + e.tgd * constants::SPEED_OF_LIGHT;
                    const int ambiguity = is_rover ? 10000 + 3 * e.satellite.prn + 10 * (e.satellite.prn % 3)
                                                  : 20000 + 2 * e.satellite.prn;
                    phase = range / constants::GPS_L1_WAVELENGTH + ambiguity;
                }
                if (noisy) {
                    code += 0.05 * std::sin((epoch + 1) * 0.93 + e.satellite.prn * 1.37 + (is_rover ? 0.2 : 2.1));
                    phase += 0.0005 / constants::GPS_L1_WAVELENGTH
                        * std::cos((epoch + 1) * 1.13 + e.satellite.prn * 0.77 + (is_rover ? 0.8 : 2.7));
                }
                out << e.satellite.toString() << std::fixed << std::setprecision(5)
                    << std::setw(14) << code << "  " << std::setprecision(4)
                    << std::setw(14) << phase << "  \n";
            }
        }
    }
    std::cerr << "generated 64 moving epochs and " << visible.size()
              << " GPS satellites; native geometry components, not whole RTKProcessor\n";
}
