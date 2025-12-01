#include "noise.h"

namespace kcored {

namespace {
constexpr long long kMax = std::numeric_limits<long long>::max();
}

GeomSampler::GeomSampler(double lambda, uint64_t seed)
    : lambda_(lambda > 0 ? lambda : 1e-12), rng_(seed) {}

void GeomSampler::reseed(uint64_t seed) {
    rng_.seed(seed);
}

long long GeomSampler::geometric() {
    if (uniform_(rng_) > -1.0 * std::expm1(-1.0 * lambda_ * static_cast<double>(kMax))) {
        return kMax;
    }
    long long left = 0;
    long long right = kMax;
    while (left + 1 < right) {
        double numerator = std::log(0.5) + std::log1p(std::exp(lambda_ * static_cast<double>(left - right)));
        long long mid = left - static_cast<long long>(std::floor(numerator / lambda_));
        if (mid <= left) {
            mid = left + 1;
        } else if (mid >= right) {
            mid = right - 1;
        }
        double q = std::expm1(lambda_ * static_cast<double>(left - mid)) /
                   std::expm1(lambda_ * static_cast<double>(left - right));
        if (uniform_(rng_) <= q) {
            right = mid;
        } else {
            left = mid;
        }
    }
    return right;
}

long long GeomSampler::twoSidedGeometric() {
    long long sample = 0;
    long long sign = -1;
    while (sample == 0 && sign == -1) {
        sample = geometric() - 1;
        sign = (uniform_(rng_) < 0.5) ? -1 : 1;
    }
    return sample * sign;
}

} // namespace kcored
