#pragma once
#include <random>
#include <limits>
#include <cmath>

namespace kcored {

// Geometric sampler matching Go implementation of Google DP Laplace mechanism
class GeomSampler {
public:
    explicit GeomSampler(double lambda, uint64_t seed = std::random_device{}());
    void reseed(uint64_t seed);
    long long twoSidedGeometric();

private:
    long long geometric();

    double lambda_;
    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
};

} // namespace kcored
