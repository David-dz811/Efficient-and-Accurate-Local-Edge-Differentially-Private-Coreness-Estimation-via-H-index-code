#pragma once
#include <vector>
#include <cstdint>
#include <cmath>

namespace kcored {

class LDS {
public:
    LDS() = default;
    LDS(int n, double levelsPerGroup);

    void reset(int n, double levelsPerGroup);

    uint32_t getLevel(int idx) const;
    void levelIncrease(int idx);
    int groupForLevel(int level) const;
    int size() const { return static_cast<int>(levels_.size()); }
    double levelsPerGroup() const { return levelsPerGroup_; }

private:
    std::vector<uint32_t> levels_;
    double levelsPerGroup_ = 1.0;
};

} // namespace kcored
