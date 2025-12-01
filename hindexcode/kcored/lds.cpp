#include "lds.h"
#include <stdexcept>

namespace kcored {

LDS::LDS(int n, double levelsPerGroup) {
    reset(n, levelsPerGroup);
}

void LDS::reset(int n, double levelsPerGroup) {
    levels_.assign(n, 0u);
    levelsPerGroup_ = levelsPerGroup > 0 ? levelsPerGroup : 1.0;
}

uint32_t LDS::getLevel(int idx) const {
    if (idx < 0 || idx >= static_cast<int>(levels_.size())) {
        throw std::out_of_range("LDS::getLevel index out of range");
    }
    return levels_[idx];
}

void LDS::levelIncrease(int idx) {
    if (idx < 0 || idx >= static_cast<int>(levels_.size())) {
        throw std::out_of_range("LDS::levelIncrease index out of range");
    }
    ++levels_[idx];
}

int LDS::groupForLevel(int level) const {
    if (levelsPerGroup_ <= 0) return 0;
    return static_cast<int>(std::floor(static_cast<double>(level) / levelsPerGroup_));
}

} // namespace kcored
