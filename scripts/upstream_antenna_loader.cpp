// Execute the pinned anonymous-namespace ANTEX receiver loader.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include "ppp.cpp"
bool nativeSatelliteCalibrations(const std::string& path,
                                std::vector<libgnss::SatelliteAntexEntry>& entries) {
    return libgnss::loadSatelliteAntexOffsets(path,entries);
}
bool nativeReceiverCalibration(const std::string& path,
                               std::array<libgnss::Vector3d,2>& offsets_enu,
                               std::array<libgnss::ReceiverPcvGrid,2>& grids) {
    using namespace libgnss;
    std::map<std::string,std::map<SignalType,Vector3d>> offsets;
    std::map<std::string,std::map<SignalType,ReceiverPcvGrid>> pcv;
    if(!loadReceiverAntexOffsets(path,offsets,pcv)) return false;
    int i=0;
    for(auto signal:{SignalType::GPS_L1CA,SignalType::GPS_L2C}) {
        offsets_enu[i]=offsets.at("TESTANT NONE").at(signal);
        grids[i]=pcv.at("TESTANT NONE").at(signal);++i;
    }
    return true;
}
