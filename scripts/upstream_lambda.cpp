// Independent reference: execute the unmodified pinned libgnss++ LAMBDA.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE and LICENSE-RTKLIB.
#include <libgnss++/algorithms/lambda.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using Eigen::MatrixXd;
using Eigen::VectorXd;

static void emit(const std::string& name, const VectorXd& a, const MatrixXd& q, int m) {
    libgnss::LambdaCandidateDiagnostics out;
    if (!libgnss::lambdaSearchTopK(a, q, m, out)) throw std::runtime_error(name);
    double ratio = 0;
    if (m >= 2 && out.squared_residuals(0) > 0) ratio = out.squared_residuals(1) / out.squared_residuals(0);
    const int n = a.size();
    for (int candidate = 0; candidate < m; ++candidate) {
        std::cout << name << ',' << n << ',' << m << ',' << candidate;
        for (int i = 0; i < n; ++i) std::cout << ',' << a(i);
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) std::cout << ',' << q(i,j);
        for (int i = 0; i < n; ++i) std::cout << ',' << out.candidates(i,candidate);
        std::cout << ',' << out.squared_residuals(candidate);
        for (int i = 0; i < n; ++i) std::cout << ',' << out.conditional_variances(i);
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) std::cout << ',' << out.decorrelation_transform(i,j);
        for (int i = 0; i < n; ++i) std::cout << ',' << out.decorrelated_float(i);
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) std::cout << ',' << out.decorrelated_covariance(i,j);
        std::cout << ',' << ratio << '\n';
    }
}

int main() {
    std::cout << std::setprecision(17);
    std::cout << "# name,n,m,candidate,a[n],Q[n*n row-major],fixed[n],residual,D[n],Z[n*n row-major],z[n],Qz[n*n row-major],ratio\n";
    VectorXd a(1); MatrixXd q(1,1);
    a << -1.5; q << 0.04; emit("negative_half", a, q, 4);
    a << 7.0; emit("exact_integer", a, q, 2);
    a << 0.49; emit("single_candidate", a, q, 1);
    a.resize(2); q.resize(2,2);
    a << 1.4, 2.4; q << 1., 0.99, 0.99, 1.; emit("correlated", a, q, 5);
    a << -12.32, 43.74; q << 0.02, -0.018, -0.018, 0.02; emit("negative_correlation", a, q, 3);
    a << 1000000.12, -2000000.34; q << 0.05, 0.049999, 0.049999, 0.05; emit("near_singular", a, q, 2);
    a.resize(3); q.resize(3,3);
    a << 5.45, 3.1, 2.97;
    q << 6.29, 5.978, 0.544, 5.978, 6.292, 2.34, 0.544, 2.34, 6.288;
    emit("three_dimensions", a, q, 4);
    a << -0.5, 0.5, 1.5; q = MatrixXd::Identity(3,3) * 0.25; emit("ties", a, q, 4);
    for (int n : {4, 8, 12}) {
        MatrixXd b(n,n);
        a.resize(n);
        for (int i = 0; i < n; ++i) {
            a(i) = 10 * i - 13.0 + 0.37 * std::sin(1.3 * (i + 1));
            for (int j = 0; j < n; ++j) b(i,j) = 0.2 * std::sin((i+1)*2.3 + (j+1)*1.7) + (i == j ? 0.05 : 0.0);
        }
        q = b * b.transpose() + MatrixXd::Identity(n,n) * 0.001;
        emit("dense_" + std::to_string(n), a, q, 3);
    }
}
