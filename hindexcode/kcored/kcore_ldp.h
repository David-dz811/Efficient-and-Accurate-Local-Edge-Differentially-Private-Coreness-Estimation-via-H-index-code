#pragma once
#include "lds.h"
#include "noise.h"
#include <string>
#include <vector>
#include <cstdint>

namespace kcored {

struct RunConfig {
    int n = 0;
    double psi = 0.5;
    double epsilon = 1.0;
    double split = 0.8;
    bool bias = true;
    int biasFactor = 8;
    bool noise = true;
    int numWorkers = 1;
    std::string partitionDir; // path containing files 0.txt ...
    double lambdaCoordinator = 0.5;
    // Optional result writing
    bool writePerRun = false;
    bool writeAverage = false;
    std::string resultDir;
    std::string resultPrefix;
    std::string resultTag;
    int runOffset = 0;
    std::vector<std::string> perRunFilenames;
    std::string meanFilename;
};

struct RunStatistics {
    double preprocessingSeconds = 0.0;
    double algorithmSeconds = 0.0;
};

struct SingleRunResult {
    std::vector<double> coreNumbers;
    RunStatistics stats;
};

SingleRunResult RunSingle(const RunConfig& config, uint64_t seed = 0);
std::vector<double> RunMultiple(const RunConfig& config, int runs,
                                std::vector<RunStatistics>* stats = nullptr,
                                uint64_t seed = 0);

} // namespace kcored
