// Call the pinned, unmodified RTKProcessor::updateBias through its existing
// unit-test friend. Isolates L1 datum alignment from slip detectors and Q noise.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/algorithms/rtk.hpp>
#include <libgnss++/core/constants.hpp>
#include <fstream>
#include <iomanip>
#include <optional>
#include <string>
#include <vector>

namespace libgnss {
struct DatumCaseInput {
    SatelliteId satellite;
    double raw;
    std::optional<double> retained;
    bool slip = false;
};
struct RTKProcessorTestAccess {
    static void run(std::ostream& out,const std::string& name,
                    const std::vector<DatumCaseInput>& inputs) {
        RTKProcessor::RTKConfig config;
        config.enable_cycle_slip_detection=false;
        config.enable_code_slip_detection=false;
        config.enable_doppler_slip_detection=false;
        config.process_noise_ambiguity=0.0;
        RTKProcessor processor(config);
        auto& state=processor.filter_state_;
        state.state=VectorXd::Zero(RTKProcessor::LEGACY_NX);
        state.covariance=MatrixXd::Zero(RTKProcessor::LEGACY_NX,RTKProcessor::LEGACY_NX);
        std::map<SatelliteId,RTKProcessor::SatelliteData> data;
        for(const auto& input:inputs) {
            auto& d=data[input.satellite];
            d.satellite=input.satellite;
            d.has_l1=true;
            d.l1_wavelength=constants::GPS_L1_WAVELENGTH;
            d.rover_l1_phase=input.raw;
            d.l1_lli=input.slip?1:0;
            if(input.retained) processor.getOrCreateN1Index(input.satellite,*input.retained);
        }
        // Nontrivial cross-covariance tests that datum translation is not a
        // measurement update. New/slipped arcs start with upstream 900 cycles².
        for(const auto& a:inputs) if(a.retained) {
            const int i=state.n1_indices.at(a.satellite);
            state.covariance(i,i)=4.0;
            for(const auto& b:inputs) if(b.retained && a.satellite!=b.satellite)
                state.covariance(i,state.n1_indices.at(b.satellite))=0.125;
        }
        processor.updateBias(data,1.0);
        for(const auto& input:inputs) {
            const int i=state.n1_indices.at(input.satellite);
            out<<name<<','<<input.satellite.toString()<<','<<input.raw<<',';
            if(input.retained) out<<*input.retained;
            out<<','<<input.slip<<','<<state.state(i)<<','<<state.covariance(i,i);
            double cross=0;
            for(const auto& other:inputs) if(other.satellite!=input.satellite)
                cross+=state.covariance(i,state.n1_indices.at(other.satellite));
            out<<','<<cross<<'\n';
        }
    }
};
}

int main(int argc,char** argv) {
    if(argc!=2) return 2;
    using namespace libgnss;
    const auto g=[](int p){return SatelliteId{GNSSSystem::GPS,static_cast<uint8_t>(p)};};
    const auto e=[](int p){return SatelliteId{GNSSSystem::Galileo,static_cast<uint8_t>(p)};};
    std::ofstream out(argv[1]);
    out<<"case,satellite,raw_cycles,retained_cycles,slip,aligned_cycles,variance_cycles2,sum_cross_cycles2\n"<<std::setprecision(17);
    using Peer=RTKProcessorTestAccess;
    Peer::run(out,"new",{{g(1),10,{}},{g(2),-4,{}},{e(1),3.25,{}}});
    Peer::run(out,"retained",{{g(1),20,10},{g(2),17,-4},{g(3),-1,2}});
    Peer::run(out,"join",{{g(1),20,10},{g(2),17,-4},{g(3),-1,{}}});
    Peer::run(out,"one_retained",{{g(1),-20,10},{g(2),17,{}},{g(3),-1,{}}});
    Peer::run(out,"mixed",{{g(1),20,10},{g(2),17,-4},{g(3),-1,{}},{e(1),-25,3},{e(2),-29,-5},{e(3),15,{}}});
    Peer::run(out,"slip",{{g(1),20,10},{g(2),17,-4,true},{g(3),-1,2}});
    Peer::run(out,"all_slipped",{{g(1),20,10,true},{g(2),17,-4,true},{e(1),-25,3,true}});
    Peer::run(out,"zero_sum",{{g(1),15,10},{g(2),-9,-4},{g(3),0,{}}});
    Peer::run(out,"large_common",{{g(1),1000010.25,10},{g(2),999996.25,-4},{g(3),1000033,{}}});
    // Explicitly records the upstream zero sentinel behavior. Rust uses Option
    // and tests this intentional difference, rather than dropping valid zero.
    Peer::run(out,"native_zero_sentinel",{{g(1),10,0},{g(2),5,2},{g(3),8,{}}});
    return out?0:1;
}
