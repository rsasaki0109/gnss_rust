// Native measurement/Kalman components under every reference of a fixed SD set.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <libgnss++/algorithms/kalman.hpp>
#include <libgnss++/algorithms/rtk_measurement.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>
using namespace libgnss;
using Eigen::Vector3d;
static void vector_csv(std::ostream& out, const VectorXd& value) {
    for (int i = 0; i < value.size(); ++i) out << ',' << value(i);
}
static void matrix_csv(std::ostream& out, const MatrixXd& value) {
    for (int i = 0; i < value.rows(); ++i)
        for (int j = 0; j < value.cols(); ++j) out << ',' << value(i, j);
}
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    std::ofstream out(argv[1]);
    if (!out) return 3;
    out << std::setprecision(17)
        << "# case,satellites,states,rows,reference_prn,x,P,H,v,R,xp,Pp,NIS,vRinvv,phase_RMS,code_RMS; matrices row-major\n";
    const double wavelength = 299792458.0 / 1575.42e6;
    for (int id = 0; id < 4; ++id) {
        const int count = 4 + id, n = 3 + count, k = count - 1;
        VectorXd x(n);
        MatrixXd b(n, n);
        for (int i = 0; i < n; ++i) {
            x(i) = 0.03 * (i - 2);
            for (int j = 0; j < n; ++j)
                b(i, j) = 0.2 * std::sin((i + 1) * .7 + (j + 1) * 1.3);
        }
        const MatrixXd p = b * b.transpose() + MatrixXd::Identity(n, n) * .5;
        std::vector<Vector3d> geometry(count);
        VectorXd sd_phase(count), sd_code(count), shared_phase(count), shared_code(count);
        std::vector<double> phase_variance(count), code_variance(count);
        for (int i = 0; i < count; ++i) {
            geometry[i] = Vector3d(std::sin(i + .3), std::cos(2 * i + .4), .5 + .1 * i).normalized();
            sd_phase(i) = .002 * std::sin(i + id) + (i == 0 ? .008 : 0);
            sd_code(i) = 2 * std::sin((i + 1) * (id + 1));
            shared_phase(i) = .0004 * (i + 1);
            shared_code(i) = .03 * std::cos(i + 1);
            phase_variance[i] = std::pow(.01 * (1 + .1 * i), 2);
            code_variance[i] = std::pow(.3 * (1 + .2 * i), 2);
        }
        for (int reference = 0; reference < count; ++reference) {
            std::vector<rtk_measurement::MeasurementBlock> blocks(2);
            VectorXd shared(2 * k);
            for (int kind = 0; kind < 2; ++kind) {
                blocks[kind].kind = kind == 0
                    ? rtk_measurement::MeasurementKind::PHASE
                    : rtk_measurement::MeasurementKind::CODE;
                int target = 0;
                for (int i = 0; i < count; ++i) {
                    if (i == reference) continue;
                    rtk_measurement::MeasurementRow row;
                    row.baseline_coefficients = geometry[reference] - geometry[i];
                    row.residual = kind == 0
                        ? sd_phase(reference) - sd_phase(i)
                        : sd_code(reference) - sd_code(i);
                    const auto& variance = kind == 0 ? phase_variance : code_variance;
                    row.reference_variance = variance[reference];
                    row.satellite_variance = variance[i];
                    if (kind == 0) row.state_coefficients = {
                        {3 + reference, wavelength}, {3 + i, -wavelength}};
                    blocks[kind].rows.push_back(row);
                    shared(kind * k + target++) = kind == 0
                        ? shared_phase(reference) - shared_phase(i)
                        : shared_code(reference) - shared_code(i);
                }
            }
            const auto system = rtk_measurement::assembleMeasurementSystem(blocks, n);
            // One actual shared SD noise source couples satellites and phase/code.
            const MatrixXd r = system.covariance + shared * shared.transpose();
            VectorXd updated = x, weighted, hph;
            MatrixXd covariance = p;
            if (kalmanFilter(updated, covariance, system.design_matrix,
                    system.residuals, r, std::vector<bool>(n, true), &weighted, &hph) != 0)
                throw std::runtime_error("native Kalman update failed");
            out << id << ',' << count << ',' << n << ',' << 2 * k << ',' << reference + 1;
            vector_csv(out, x); matrix_csv(out, p);
            matrix_csv(out, system.design_matrix); vector_csv(out, system.residuals);
            matrix_csv(out, r); vector_csv(out, updated); matrix_csv(out, covariance);
            out << ',' << system.residuals.dot(weighted)
                << ',' << system.residuals.dot(r.ldlt().solve(system.residuals))
                << ',' << std::sqrt(system.residuals.head(k).squaredNorm() / k)
                << ',' << std::sqrt(system.residuals.tail(k).squaredNorm() / k) << '\n';
        }
    }
    return out ? 0 : 1;
}
