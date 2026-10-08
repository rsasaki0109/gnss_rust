// Expose the actual pinned anonymous-namespace solar helper in this TU.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include "ppp_corrections.cpp"
libgnss::Vector3d nativeWindupSun(const libgnss::GNSSTime& time) {
    return libgnss::approximateSunPositionEcef(time);
}
libgnss::Vector3d nativeLegacyMoon(const libgnss::GNSSTime& time) {
    return libgnss::approximateMoonPositionEcef(time);
}
libgnss::Vector3d nativeLegacyBodyTide(const libgnss::Vector3d& station,
                                     const libgnss::Vector3d& body,double gm) {
    return libgnss::bodyTideDisplacement(station,body,gm);
}
