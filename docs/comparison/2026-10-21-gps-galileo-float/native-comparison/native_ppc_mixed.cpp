// Pinned whole RTKProcessor adapter; no upstream solver source is modified.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/algorithms/rtk.hpp>
#include <libgnss++/io/rinex.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace libgnss;


static bool eligible(const std::pair<SatelliteId,std::string>& key) {
    return (key.first.system==GNSSSystem::GPS && key.second=="1C") ||
        (key.first.system==GNSSSystem::Galileo && (key.second=="1C" || key.second=="1X"));
}
static ObservationData selected_l1(const ObservationData& raw) {
    ObservationData selected(raw.time);
    selected.receiver_position=raw.receiver_position;selected.receiver_clock_bias=raw.receiver_clock_bias;
    std::map<SatelliteId,std::pair<std::string,Observation>> choices;
    auto rank=[](const std::string& code,const Observation& o) {
        bool complete=o.has_pseudorange && std::isfinite(o.pseudorange) && o.pseudorange>0 &&
            o.has_carrier_phase && std::isfinite(o.carrier_phase);
        return std::pair<bool,bool>(!complete,code!="1C");
    };
    for(const auto&[key,o]:raw.rinex_tracking_observations) {
        if(!eligible(key))continue;
        const auto old=choices.find(key.first);
        if(old==choices.end() || rank(key.second,o)<rank(old->second.first,old->second.second))
            choices[key.first]={key.second,o};
    }
    for(const auto&[sat,choice]:choices) {
        (void)sat;selected.observations.push_back(choice.second);
        selected.addRinexTrackingObservation(choice.first,choice.second);
    }
    selected.setRinexFrequencySlot(GNSSSystem::GPS,0,"1C");
    return selected;
}
static bool continuity(const ObservationData& raw) {
    for(const auto&[key,value]:raw.rinex_tracking_observations)
        if(eligible(key) && (value.lli&3))return true;
    return false;
}
static void ledger(std::ostream& out,const char* side,const ObservationData& epoch,bool reset) {
    for(const auto&[key,o]:epoch.rinex_tracking_observations) {
        if(!o.has_pseudorange && !o.has_carrier_phase && !o.has_doppler && o.snr==0.0 && o.lli==0)continue;
        out<<side<<','<<epoch.time.week<<','<<epoch.time.tow<<','<<o.satellite.toString()<<','<<key.second<<','
           <<o.has_pseudorange<<','<<o.pseudorange<<','<<o.has_carrier_phase<<','<<o.carrier_phase<<','
           <<o.has_doppler<<','<<o.doppler<<','<<o.snr<<','<<static_cast<int>(o.lli)<<','<<reset<<'\n';
    }
}
static const char* status_name(SolutionStatus status) {
    switch (status) {
        case SolutionStatus::NONE: return "NONE";
        case SolutionStatus::SPP: return "SPP";
        case SolutionStatus::DGPS: return "DGPS";
        case SolutionStatus::FLOAT: return "FLOAT";
        case SolutionStatus::FIXED: return "FIX";
        case SolutionStatus::PPP_FLOAT: return "PPP_FLOAT";
        case SolutionStatus::PPP_FIXED: return "PPP_FIX";
        case SolutionStatus::PROPAGATED: return "PROPAGATED";
    }
    throw std::runtime_error("unknown native status");
}
int main(int argc, char** argv) {
    try {
        if (argc != 8) {
            std::cerr << "usage: native_ppc_rtk CASE_DIR BASE_X BASE_Y BASE_Z OUTPUT_CSV INPUT_LEDGER_CSV previous|spp\n";
            return 2;
        }
        for (int index : {5, 6})
            if (std::filesystem::exists(argv[index])) throw std::runtime_error("output already exists");
        const std::string root = argv[1];
        const Vector3d base_position(std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]));
        const std::string seed_policy=argv[7]; if(seed_policy!="previous" && seed_policy!="spp") throw std::runtime_error("seed policy"); const std::size_t maximum=0;
        io::RINEXReader rover_reader, base_reader, nav_reader;
        if (!rover_reader.open(root+"/rover.obs") || !base_reader.open(root+"/base.obs") ||
            !nav_reader.open(root+"/base.nav")) throw std::runtime_error("cannot open input");
        io::RINEXReader::RINEXHeader rover_header, base_header;
        if (!rover_reader.readHeader(rover_header) || !base_reader.readHeader(base_header))
            throw std::runtime_error("cannot read observation header");
        if (!rover_header.has_approximate_position) throw std::runtime_error("missing header seed");
        NavigationData navigation;
        if (!nav_reader.readNavigationData(navigation)) throw std::runtime_error("cannot read navigation");
        std::size_t gps_navigation_records=0;
        for (const auto& [satellite,records] : navigation.ephemeris_data)
            if (satellite.system==GNSSSystem::GPS) gps_navigation_records+=records.size();
        std::ofstream out(argv[5]), inputs(argv[6]);
        if (!out || !inputs) throw std::runtime_error("cannot create outputs");
        out << std::scientific << std::setprecision(17);
        inputs << std::scientific << std::setprecision(17);
        inputs << "side,week,tow,satellite,tracking_code,has_code,code_m,has_phase,phase_cycles,has_doppler,doppler_hz,snr_db_hz,lli,derived_reset\n";
        out << "kind,week,tow,side,status,x_m,y_m,z_m,satellites,is_valid,native_status,ratio,position_qxx_m2,position_qyy_m2,position_qzz_m2,code_slips,doppler_slips,lli_slips,derived_reset,solution_week,solution_tow,seed_x_m,seed_y_m,seed_z_m\n";
        // Library diagnostic stdout must not be mixed with adapter output.
        std::cout.rdbuf(std::cerr.rdbuf());
        RTKProcessor::RTKConfig config;
        config.position_mode = RTKProcessor::RTKConfig::PositionMode::KINEMATIC;
        config.elevation_mask = 10.0 * std::acos(-1.0) / 180;
        config.max_baseline_length = 10000.0;
        config.prefer_rover_position_seed = seed_policy=="previous";
        RTKProcessor processor(config);
        ProcessorConfig general;
        general.mode = PositioningMode::RTK_FIXED;
        general.elevation_mask = 10.0;
        processor.initialize(general);
        processor.setBasePosition(base_position);
        Vector3d previous_seed = rover_header.approximate_position;
        ObservationData rover, base;
        bool have_rover = rover_reader.readObservationEpoch(rover);
        bool have_base = base_reader.readObservationEpoch(base);
        bool reset_pending = false;
        std::size_t pairs=0, rover_only=0, base_only=0, notices=0, valid=0, fixed=0;
        GNSSTime last_rover, last_base;
        bool seen_rover=false, seen_base=false;
        auto advance = [](io::RINEXReader& reader, ObservationData& epoch,
                          bool& seen, GNSSTime& previous) {
            const bool available = reader.readObservationEpoch(epoch);
            if (available) {
                if (seen && epoch.time-previous <= 0.0) throw std::runtime_error("non-increasing epoch");
                previous=epoch.time; seen=true;
            }
            return available;
        };
        if (have_rover) {last_rover=rover.time;seen_rover=true;}
        if (have_base) {last_base=base.time;seen_base=true;}
        while (have_rover || have_base) {
            if (have_rover && have_base && std::abs(rover.time-base.time)<=1e-7) {
                auto r=selected_l1(rover), b=selected_l1(base);
                ledger(inputs,"rover",r,reset_pending);
                ledger(inputs,"base",b,reset_pending);
                r.receiver_position=previous_seed;
                const Vector3d supplied_seed=previous_seed;
                b.receiver_position=base_position;
                if (reset_pending) {
                    for (auto* epoch : {&r,&b}) {
                        for (auto& observation : epoch->observations) {observation.lli|=1;observation.loss_of_lock=true;}
                        for (auto& [key,observation] : epoch->rinex_tracking_observations) {
                            (void)key;observation.lli|=1;observation.loss_of_lock=true;
                        }
                    }
                }
                const auto solution=processor.processRTKEpoch(r,b,navigation);
                const bool usable=solution.isValid() && solution.position_ecef.allFinite();
                const auto& telemetry=processor.getLastDebugTelemetry();
                out << "pair," << rover.time.week << ',' << rover.time.tow << ",both,"
                    << (usable?status_name(solution.status):"INVALID") << ',';
                if (usable) out << solution.position_ecef(0) << ',' << solution.position_ecef(1) << ',' << solution.position_ecef(2);
                else out << ",,";
                out << ',' << solution.num_satellites << ',' << usable << ',' << status_name(solution.status)
                    << ',' << solution.ratio << ',';
                for (int i=0;i<3;++i) {
                    if (std::isfinite(solution.position_covariance(i,i))) out << solution.position_covariance(i,i);
                    out << ',';
                }
                out << telemetry.code_slip_l1_count << ',' << telemetry.doppler_slip_l1_count << ','
                    << telemetry.lli_slip_l1_count << ',' << reset_pending
                    << ',' << solution.time.week << ',' << solution.time.tow
                    << ',' << supplied_seed(0) << ',' << supplied_seed(1) << ',' << supplied_seed(2) << '\n';
                if (usable) ++valid;
                if (usable && solution.isFixed()) ++fixed;
                if (usable && (solution.status==SolutionStatus::FLOAT || solution.status==SolutionStatus::FIXED))
                    previous_seed=solution.position_ecef;
                reset_pending=false; ++pairs;
                have_rover=advance(rover_reader,rover,seen_rover,last_rover);
                have_base=advance(base_reader,base,seen_base,last_base);
                if (maximum && pairs>=maximum) break;
            } else {
                const bool is_rover=have_rover && (!have_base || rover.time-base.time<0);
                const auto& epoch=is_rover?rover:base;
                const bool broken=continuity(epoch);
                if (broken) {reset_pending=true;++notices;}
                out << "unmatched," << epoch.time.week << ',' << epoch.time.tow << ','
                    << (is_rover?"rover":"base") << ",UNMATCHED";
                for (int column=5;column<24;++column) out << ',';
                out << '\n';
                if (is_rover) {++rover_only;have_rover=advance(rover_reader,rover,seen_rover,last_rover);}
                else {++base_only;have_base=advance(base_reader,base,seen_base,last_base);}
            }
        }
        std::cerr << "native_whole_RTKProcessor=true GPS_1C_Galileo_E1_only=true pairs=" << pairs
                  << " valid=" << valid << " FIX=" << fixed << " rover_only=" << rover_only
                  << " base_only=" << base_only << " unmatched_continuity_notices=" << notices
                  << " gps_navigation_records=" << gps_navigation_records
                  << " input_drained=" << (!have_rover&&!have_base)
                  << " seed_policy=" << seed_policy << " AR_default_continuous=true performance_accepted=false\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "adapter error: " << error.what() << '\n';
        return 1;
    }
}
