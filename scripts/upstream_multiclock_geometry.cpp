// Independent Eigen geometry with one receiver-clock range per constellation.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/core/coordinates.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
    using namespace libgnss;
    const double rad=std::acos(-1.0)/180;
    std::cout << std::setprecision(17) << "# case,systems,per_system,receiver[3],satellite_xyz[3*m],Eigen_Q[(3+systems)^2]\n";
    for(int frame=0;frame<2;++frame) for(int groups:{1,2,3}) for(int per:{4,5,7,9}) {
        const double lat=(frame?-33:35)*rad,lon=(frame?-71:139)*rad;
        const Vector3d receiver=geodetic2ecef(lat,lon,45+frame*200);
        const int m=groups*per,n=3+groups;
        Eigen::MatrixXd h=Eigen::MatrixXd::Zero(m,n);
        std::cout << "frame" << frame << '_' << groups << '_' << per << ',' << groups << ',' << per;
        for(int j=0;j<3;++j) std::cout << ',' << receiver(j);
        for(int i=0;i<m;++i) {
            const double az=(13+137*i)*rad,el=(18+std::fmod(31.0*i,65.0))*rad;
            Vector3d direction(std::cos(el)*std::sin(az),std::cos(el)*std::cos(az),std::sin(el));
            const Vector3d p=receiver+enu2ecef(direction*(21000000+10000*i),lat,lon);
            for(int j=0;j<3;++j)std::cout << ',' << p(j);
            h.block<1,3>(i,0)=-(p-receiver).normalized().transpose();
            h(i,3+i/per)=1;
        }
        const Eigen::MatrixXd q=(h.transpose()*h).inverse();
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)std::cout << ',' << q(i,j);
        std::cout << '\n';
    }
}
