#include "kcore_ldp.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <filesystem>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace kcored {
namespace {

using Clock = std::chrono::steady_clock;

inline double log_a_to_base_b(int a, double b) {
    if (a <= 0 || b <= 0.0 || b == 1.0) return 0.0;
    return std::log(static_cast<double>(a)) / std::log(b);
}

struct Vertex {
    int globalId = 0;
    bool active = true;  // corresponds to permanent_zero in Go code
    int roundThreshold = 1;
    std::vector<int> neighbours;
};

class WorkerContext {
public:
    bool load(const std::string& path,
              int offset,
              int workLoad,
              double lambda1,
              double levelsPerGroup,
              bool bias,
              int biasFactor,
              bool noise,
              std::function<uint64_t()> nextSeed,
              int n);

    std::vector<int> runRound(int round,
                              double lambda2,
                              double psi,
                              int groupIndex,
                              const LDS& lds,
                              bool noise,
                              bool bias,
                              std::function<uint64_t()> nextSeed) {
        std::vector<int> promoted;
        promoted.reserve(vertices_.size());
        for (auto& v : vertices_) {
            if (v.roundThreshold == round) {
                v.active = false;
            }
            uint32_t currentLevel = 0;
            if (v.globalId >= lds.size()) {
                continue;
            }
            currentLevel = lds.getLevel(v.globalId);
            if (static_cast<int>(currentLevel) != round || !v.active) {
                continue;
            }
            int neighborCount = 0;
            for (int ngh : v.neighbours) {
                if (ngh < 0 || ngh >= lds.size()) continue;
                if (static_cast<int>(lds.getLevel(ngh)) == round) {
                    ++neighborCount;
                }
            }
            long long noised = neighborCount;
            if (noise) {
                int threshold = std::max(1, v.roundThreshold);
                double scale = lambda2 / (2.0 * static_cast<double>(threshold));
                GeomSampler sampler(scale, nextSeed());
                noised += sampler.twoSidedGeometric();
                if (bias) {
                    double denom = std::exp(2.0 * scale) - 1.0;
                    if (std::abs(denom) > 1e-12) {
                        double extra = 3.0 * (2.0 * std::exp(scale)) / std::pow(denom, 3.0);
                        noised += static_cast<long long>(extra);
                    }
                }
            }
            double rhs = std::pow(1.0 + psi, static_cast<double>(std::max(0, groupIndex)));
            if (static_cast<double>(noised) > rhs) {
                promoted.push_back(v.globalId);
            } else {
                v.active = false;
            }
        }
        return promoted;
    }

    int maxRoundThreshold() const { return maxRoundThreshold_; }

private:
    std::vector<Vertex> vertices_;
    int maxRoundThreshold_ = 0;
};

bool WorkerContext::load(const std::string& path,
                         int offset,
                         int workLoad,
                         double lambda1,
                         double levelsPerGroup,
                         bool bias,
                         int biasFactor,
                         bool noise,
                         std::function<uint64_t()> nextSeed,
                         int n) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
        throw std::runtime_error("Failed to open partition file: " + path);
    }
    std::vector<std::vector<int>> adjacency(workLoad);
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        int u, v;
        if (!(iss >> u >> v)) continue;
        int local = u - offset;
        if (local < 0 || local >= workLoad) continue;
        adjacency[local].push_back(v);
    }
    fin.close();

    GeomSampler sampler(lambda1 / 2.0, nextSeed());
    vertices_.clear();
    vertices_.reserve(workLoad);
    maxRoundThreshold_ = 0;

    for (int i = 0; i < workLoad; ++i) {
        int globalId = offset + i;
        if (globalId >= n) {
            continue;
        }
        Vertex vertex;
        vertex.globalId = globalId;
        vertex.neighbours = std::move(adjacency[i]);
        long long degree = static_cast<long long>(vertex.neighbours.size());
        long long noisedDegree = degree;
        if (noise) {
            sampler.reseed(nextSeed());
            noisedDegree += sampler.twoSidedGeometric();
            if (bias) {
                double numerator = 2.0 * std::exp(lambda1);
                double denominator = std::exp(2.0 * lambda1) - 1.0;
                if (std::abs(denominator) > 1e-12) {
                    double cap = static_cast<double>(biasFactor) * (numerator / denominator);
                    double boundedCap = std::min(cap, static_cast<double>(noisedDegree));
                    noisedDegree -= static_cast<long long>(std::floor(boundedCap));
                }
            }
            noisedDegree = std::max<long long>(0, noisedDegree);
            noisedDegree += 1; // ensure positive
        }
        double thresholdVal = 1.0;
        if (noisedDegree > 0) {
            thresholdVal = std::ceil(log_a_to_base_b(static_cast<int>(noisedDegree), 2.0)) * levelsPerGroup;
        }
        vertex.roundThreshold = static_cast<int>(thresholdVal) + 1;
        vertex.active = true;
        vertices_.push_back(std::move(vertex));
        maxRoundThreshold_ = std::max(maxRoundThreshold_, vertices_.back().roundThreshold);
    }
    return true;
}

std::string JoinPath(const std::string& base, const std::string& file) {
    if (base.empty()) return file;
    if (base.back() == '/' || base.back() == '\\') {
        return base + file;
    }
    return base + "/" + file;
}

std::vector<double> EstimateCoreNumbers(const LDS& lds,
                                        int n,
                                        double levelsPerGroup,
                                        double lambda,
                                        double phi) {
    std::vector<double> core(n, 0.0);
    double twoPlusLambda = 2.0 + lambda;
    double onePlusPhi = 1.0 + phi;
    for (int i = 0; i < n; ++i) {
        double nodeLevel = static_cast<double>(lds.getLevel(i));
        double fracNumerator = nodeLevel + 1.0;
        double power = std::max(std::floor(fracNumerator / levelsPerGroup) - 1.0, 0.0);
        core[i] = twoPlusLambda * std::pow(onePlusPhi, power);
    }
    return core;
}

} // namespace

void WriteRunOutput(const RunConfig& config,
                    const SingleRunResult& result,
                    const std::vector<double>& values,
                    int runIndex,
                    const std::string* explicitName = nullptr,
                    const char* defaultSuffix = "") {
    if (!config.writePerRun && !config.writeAverage) return;
    if (config.resultDir.empty()) return;

    std::filesystem::path dir(config.resultDir);
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);

    std::string filename;
    if (explicitName && !explicitName->empty()) {
        filename = *explicitName;
    } else {
        std::string prefix = config.resultPrefix.empty() ? "kcored" : config.resultPrefix;
        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss.precision(6);
        oss << prefix;
        if (!config.resultTag.empty()) {
            oss << "_" << config.resultTag;
        }
        if (defaultSuffix && *defaultSuffix) {
            oss << defaultSuffix;
        } else {
            oss << "_run" << (config.runOffset + runIndex);
        }
        oss << ".txt";
        filename = oss.str();
    }

    auto filepath = dir / filename;
    std::ofstream fout(filepath);
    if (!fout.is_open()) return;

    fout.setf(std::ios::fixed);
    fout.precision(8);
    fout << "Preprocessing Time: " << result.stats.preprocessingSeconds << "\n";
    fout << "Algorithm Time: " << result.stats.algorithmSeconds << "\n";
    fout.precision(4);
    for (size_t i = 0; i < values.size(); ++i) {
        fout << i << ": " << values[i] << "\n";
    }
}

SingleRunResult RunSingle(const RunConfig& config, uint64_t seed) {
    if (config.n <= 0) {
        throw std::invalid_argument("RunSingle requires positive n");
    }
    if (config.numWorkers <= 0) {
        throw std::invalid_argument("RunSingle requires positive numWorkers");
    }

    std::mt19937_64 seedRng(seed ? seed : std::random_device{}());
    auto nextSeed = [&]() -> uint64_t { return seedRng(); };

    auto start = Clock::now();
    double levelsPerGroup = std::ceil(log_a_to_base_b(std::max(2, config.n), 1.0 + config.psi)) / 4.0;
    double logTerm = log_a_to_base_b(std::max(2, config.n), 1.0 + config.psi);
    int roundsParam = static_cast<int>(std::ceil(4.0 * std::pow(logTerm, 1.2)));
    int numberOfRounds = std::max(0, roundsParam - 2);

    double lambda1 = config.epsilon * config.split;
    double lambda2 = config.epsilon * (1.0 - config.split);

    int chunk = config.n / config.numWorkers;
    int extra = config.n % config.numWorkers;

    std::vector<WorkerContext> workers;
    workers.reserve(config.numWorkers);
    int maxPublicRoundThreshold = 0;
    for (int i = 0; i < config.numWorkers; ++i) {
        int workLoad = (i == config.numWorkers - 1) ? (chunk + extra) : chunk;
        int offset = i * chunk;
        WorkerContext worker;
        std::string filename = JoinPath(config.partitionDir, std::to_string(i) + ".txt");
        worker.load(filename, offset, workLoad, lambda1, levelsPerGroup,
                    config.bias, config.biasFactor, config.noise, nextSeed, config.n);
        maxPublicRoundThreshold = std::max(maxPublicRoundThreshold, worker.maxRoundThreshold());
        workers.push_back(std::move(worker));
    }

    int totalRounds = std::min(numberOfRounds, maxPublicRoundThreshold);
    if (totalRounds < 0) totalRounds = 0;

    LDS lds(config.n, levelsPerGroup);
    auto preprocessingDone = Clock::now();

    for (int round = 0; round < totalRounds; ++round) {
        int groupIndex = lds.groupForLevel(round);
        std::vector<int> promoted;
        for (auto& worker : workers) {
            auto local = worker.runRound(round, lambda2, config.psi, groupIndex,
                                         lds, config.noise, config.bias, nextSeed);
            promoted.insert(promoted.end(), local.begin(), local.end());
        }
        for (int id : promoted) {
            if (id >= 0 && id < lds.size()) {
                lds.levelIncrease(id);
            }
        }
    }

    auto algorithmDone = Clock::now();

    SingleRunResult result;
    result.coreNumbers = EstimateCoreNumbers(lds, config.n, levelsPerGroup, config.lambdaCoordinator, config.psi);
    result.stats.preprocessingSeconds = std::chrono::duration<double>(preprocessingDone - start).count();
    result.stats.algorithmSeconds = std::chrono::duration<double>(algorithmDone - preprocessingDone).count();
    if (config.writePerRun) {
        const std::string* name = nullptr;
        if (!config.perRunFilenames.empty()) {
            name = &config.perRunFilenames.front();
        }
        WriteRunOutput(config, result, result.coreNumbers, /*runIndex=*/0, name);
    }
    return result;
}

std::vector<double> RunMultiple(const RunConfig& config,
                                int runs,
                                std::vector<RunStatistics>* stats,
                                uint64_t seed) {
    if (runs <= 0) {
        auto single = RunSingle(config, seed);
        if (stats) stats->push_back(single.stats);
        return single.coreNumbers;
    }
    std::vector<long double> accum(config.n, 0.0L);
    std::vector<double> last;
    std::vector<RunStatistics> localStats;
    localStats.reserve(runs);
    std::mt19937_64 seedRng(seed ? seed : std::random_device{}());
    for (int run = 0; run < runs; ++run) {
        uint64_t runSeed = seedRng();
        RunConfig runConfig = config;
        runConfig.writePerRun = false;
        runConfig.writeAverage = false;
        auto single = RunSingle(runConfig, runSeed);
        last = single.coreNumbers;
        if (stats) stats->push_back(single.stats);
        localStats.push_back(single.stats);
        if (config.writePerRun) {
            RunConfig outConfig = config;
            outConfig.writeAverage = false;
            const std::string* name = nullptr;
            if (run < config.perRunFilenames.size()) {
                name = &config.perRunFilenames[run];
            }
            WriteRunOutput(outConfig, single, single.coreNumbers, run, name);
        }
        for (int i = 0; i < config.n; ++i) {
            accum[i] += static_cast<long double>(single.coreNumbers[i]);
        }
    }
    std::vector<double> mean(config.n, 0.0);
    for (int i = 0; i < config.n; ++i) {
        mean[i] = static_cast<double>(accum[i] / static_cast<long double>(runs));
    }
    if (config.writeAverage) {
        RunConfig outConfig = config;
        outConfig.writePerRun = false;
        std::vector<double> meanCopy = mean;
        SingleRunResult avgResult;
        avgResult.coreNumbers = meanCopy;
        if (!localStats.empty()) {
            double sumPre = 0.0, sumAlgo = 0.0;
            for (const auto& st : localStats) {
                sumPre += st.preprocessingSeconds;
                sumAlgo += st.algorithmSeconds;
            }
            avgResult.stats.preprocessingSeconds = sumPre / localStats.size();
            avgResult.stats.algorithmSeconds = sumAlgo / localStats.size();
        }
        const std::string* meanName = config.meanFilename.empty() ? nullptr : &config.meanFilename;
        WriteRunOutput(outConfig, avgResult, meanCopy, /*runIndex=*/0, meanName, "_mean");
    }
    return mean;
}

} // namespace kcored
