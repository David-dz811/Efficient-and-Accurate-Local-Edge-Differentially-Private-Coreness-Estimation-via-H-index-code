#include "counting.h"
#include "kcored/kcore_ldp.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>
#include <thread>
using namespace std;

// helper: simple parallel_for wrapper
template<typename Func>
void parallel_for(int start, int end, Func func, int num_threads = std::thread::hardware_concurrency()) {
    if (num_threads <= 1) {
        for (int i = start; i < end; ++i) {
            func(i);
        }
        return;
    }
    
    std::vector<std::future<void>> futures;
    int chunk_size = (end - start) / num_threads;
    
    for (int t = 0; t < num_threads; ++t) {
        int thread_start = start + t * chunk_size;
        int thread_end = (t == num_threads - 1) ? end : start + (t + 1) * chunk_size;
        
        futures.push_back(std::async(std::launch::async, [=]() {
            for (int i = thread_start; i < thread_end; ++i) {
                func(i);
            }
        }));
    }
    
    for (auto& future : futures) {
        future.wait();
    }
}

static int sample_symmetric_geometric(double epsilon, std::mt19937 &rng)
{
    std::geometric_distribution<int> geo(1.0 - std::exp(-epsilon));
    int g = geo(rng);
    return (rng() & 1) ? g : -g;
}

// compute H-index of a degree multiset

int compute_h_index(const std::vector<int> &degs)
{
    if (degs.empty())
        return 0;

    // copy and sort descending
    std::vector<int> sorted = degs;
    std::sort(sorted.begin(), sorted.end(), std::greater<int>());

    int h = 0;
    for (int d : sorted)
    {
        if (d >= h + 1)
            ++h;
        else
            break;
    }
    return h;
}

// Algorithm 3 (kcoreH): iterative H-index with per-round noise split
std::vector<int> estimate_coreness_LDP_hindex_iteration(Graph &g, double epsilon, int seed, int max_iter)
{
    const int n = g.num_nodes();
    std::mt19937 rng(seed);

    // use KCORE_H_INIT_RATIO (default 0) to split epsilon into init vs iter noise
    auto parse_ratio = [](const char* s, double defv){
        if (!s) return defv;
        try { return std::stod(s); } catch(...) { return defv; }
    };
    double r_init = parse_ratio(std::getenv("KCORE_H_INIT_RATIO"), 0.0);
    if (r_init < 0.0) r_init = 0.0; if (r_init > 1.0) r_init = 1.0;
    double eps_init = r_init * epsilon;
    double eps_iter = (max_iter > 0) ? (epsilon - eps_init) / std::max(1, max_iter) : 0.0;

    // initialization: H-index over neighbor degrees plus optional noise
    std::vector<int> coreness(n);
    for (int u = 0; u < n; ++u)
    {
        std::vector<int> degs;
        degs.reserve(g.neighbor[u].size());
        for (int v : g.neighbor[u])
            degs.push_back((int)g.neighbor[v].size());
        int h0 = compute_h_index(degs);
        if (eps_init > 0.0) h0 += sample_symmetric_geometric(eps_init, rng);
        coreness[u] = std::max(0, h0);
    }

    // pre-sample per-round noise so total budget stays within epsilon
    std::vector<int> noise_values(n * std::max(0, max_iter));
    for (int i = 0; i < (int)noise_values.size(); ++i)
        noise_values[i] = sample_symmetric_geometric(eps_iter, rng);
    
    // update by recomputing neighbor H-index each round plus noise
    for (int it = 0; it < max_iter; ++it)
    {
        std::vector<int> next(n);

        parallel_for(0, n, [&](int u) {
            std::vector<int> neigh_cor;
            neigh_cor.reserve(g.neighbor[u].size());
            for (int v : g.neighbor[u])
                neigh_cor.push_back(coreness[v]);
            int noise_idx = it * n + u;
            int new_c = compute_h_index(neigh_cor) + noise_values[noise_idx];
            next[u] = std::max(0, new_c);
        });

        coreness.swap(next);
    }
    return coreness;
}

// H-index iteration where low-degree nodes stay fixed and high-degree nodes update
static std::vector<int> estimate_coreness_LDP_hindex_high_only(
    Graph &g,
    double eps_per_iter,
    int seed,
    int max_iter,
    const std::vector<int>& init_coreness,
    const std::vector<char>& is_low_fixed)
{
    const int n = g.num_nodes();
    std::mt19937 rng(seed);
    std::vector<int> coreness = init_coreness;

    for (int it = 0; it < max_iter; ++it)
    {
        std::vector<int> next = coreness;
        for (int u = 0; u < n; ++u)
        {
            if (is_low_fixed[u]) continue; // keep low-degree nodes fixed
            std::vector<int> neigh_cor;
            neigh_cor.reserve(g.neighbor[u].size());
            for (int v : g.neighbor[u])
                neigh_cor.push_back(coreness[v]);
            int new_c = compute_h_index(neigh_cor);
            if (eps_per_iter > 0.0)
                new_c += sample_symmetric_geometric(eps_per_iter, rng);
            next[u] = std::max(0, new_c);
        }
        coreness.swap(next);
    }
    return coreness;
}

// Algorithm 6 variant: low-degree nodes stop early
std::vector<int> estimate_coreness_LDP_hindex_iter_threshold(
    Graph &g,
    double epsilon,
    int seed,
    int max_iter,
    const std::vector<std::pair<int,int>>& rules_input)
{
    const int n = g.num_nodes();
    std::mt19937 rng(seed);

    std::vector<int> degrees(n);
    for (int u = 0; u < n; ++u)
        degrees[u] = static_cast<int>(g.neighbor[u].size());

    std::vector<int> coreness(n);
    for (int u = 0; u < n; ++u)
    {
        std::vector<int> degs;
        degs.reserve(g.neighbor[u].size());
        for (int v : g.neighbor[u])
            degs.push_back(static_cast<int>(g.neighbor[v].size()));

        coreness[u] = static_cast<int>(g.neighbor[u].size());
        coreness[u] += sample_symmetric_geometric(epsilon, rng);
    }

    std::vector<int> noise_values(n * max_iter);
    for (int i = 0; i < n * max_iter; ++i)
        noise_values[i] = sample_symmetric_geometric(epsilon, rng);

    std::vector<std::pair<int,int>> rules = rules_input;
    if (rules.empty()) {
        rules = {{8,1},{15,2},{30,3},{80,4}};
    }
    std::sort(rules.begin(), rules.end(), [](const auto& a, const auto& b){ return a.first < b.first; });

    std::vector<int> freeze_after(n, max_iter);
    for (int u = 0; u < n; ++u)
    {
        int limit = max_iter;
        int deg = degrees[u];
        for (const auto& rule : rules) {
            if (deg < rule.first) {
                limit = std::clamp(rule.second, 0, max_iter);
                break;
            }
        }
        freeze_after[u] = limit;
    }

    for (int it = 0; it < max_iter; ++it)
    {
        std::vector<int> next(n);

        parallel_for(0, n, [&](int u) {
            if (it >= freeze_after[u]) {
                next[u] = coreness[u];
                return;
            }

            std::vector<int> neigh_cor;
            neigh_cor.reserve(g.neighbor[u].size());
            for (int v : g.neighbor[u])
                neigh_cor.push_back(coreness[v]);

            int noise_idx = it * n + u;
            int new_c = compute_h_index(neigh_cor) + noise_values[noise_idx];
            next[u] = new_c;
        });

        coreness.swap(next);
    }

    return coreness;
}



std::vector<int> compute_coreness(Graph &g)
{
    int n = g.num_nodes();
    // 1) initial degree of every node
    std::vector<int> degree(n);
    for (int v = 0; v < n; ++v)
    {
        degree[v] = (int)g.neighbor[v].size();
    }
    // 2) find maximum degree
    int max_deg = *std::max_element(degree.begin(), degree.end());
    // 3) bucketize degrees in bin[d]
    std::vector<int> bin(max_deg + 1, 0);
    for (int d : degree)
        ++bin[d];
    // 4) convert bin to prefix sums storing bucket start indices
    int start = 0;
    for (int d = 0; d <= max_deg; ++d)
    {
        int cnt = bin[d];
        bin[d] = start;
        start += cnt;
    }
    // 5) build vert/pos arrays describing the degeneracy ordering
    std::vector<int> vert(n), pos(n);
    for (int v = 0; v < n; ++v)
    {
        int d = degree[v];
        vert[bin[d]] = v;
        pos[v] = bin[d];
        bin[d]++;
    }
    // 6) reset bin[d] to bucket starts
    for (int d = max_deg; d > 0; --d)
    {
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;
    // 7) peeling loop that updates neighbors once a node is removed
    std::vector<int> coreness(n);
    for (int i = 0; i < n; ++i)
    {
        int v = vert[i];
        coreness[v] = degree[v];
        // decrease neighbor degrees when they are still higher than current bucket
        for (int u : g.neighbor[v])
        {
            if (degree[u] > degree[v])
            {
                int du = degree[u];
                int pu = pos[u];
                int pw = bin[du];
                int w = vert[pw];
                if (u != w)
                {
                    // swap positions in vert/pos to keep buckets compact
                    vert[pu] = w;
                    pos[w] = pu;
                    vert[pw] = u;
                    pos[u] = pw;
                }
                bin[du]++;
                degree[u]--;
            }
        }
    }
    return coreness;
}

struct KcoredPartitionInfo
{
    std::string path;
    int num_workers = 1;
    bool is_temporary = false;
};

static bool parse_partition_index(const std::string& filename, int& idx_out)
{
    auto dot = filename.find('.');
    if (dot == std::string::npos)
        return false;
    if (filename.substr(dot) != ".txt")
        return false;
    std::string prefix = filename.substr(0, dot);
    if (prefix.empty())
        return false;
    try
    {
        size_t used = 0;
        int value = std::stoi(prefix, &used);
        if (used != prefix.size() || value < 0)
            return false;
        idx_out = value;
        return true;
    }
    catch (...)
    {
        return false;
    }
}

static std::vector<int> discover_partition_indices(const std::string& dir)
{
    namespace fs = std::filesystem;
    std::vector<int> indices;
    std::error_code ec;
    if (!fs::exists(dir, ec) || !fs::is_directory(dir, ec))
        return indices;

    for (const auto& entry : fs::directory_iterator(dir, ec))
    {
        if (!entry.is_regular_file())
            continue;
        int idx = -1;
        if (parse_partition_index(entry.path().filename().string(), idx))
            indices.push_back(idx);
    }

    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    return indices;
}

static int contiguous_worker_count(const std::vector<int>& indices)
{
    if (indices.empty())
        return 0;
    int max_idx = indices.back();
    for (int i = 0; i <= max_idx; ++i)
    {
        if (!std::binary_search(indices.begin(), indices.end(), i))
            return 0;
    }
    return max_idx + 1;
}

static std::string materialize_partitions_from_graph(const Graph& g, int num_workers)
{
    namespace fs = std::filesystem;
    const int n = g.num_nodes();
    int workers = std::max(1, std::min(num_workers, std::max(1, n)));
    auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    fs::path base = fs::temp_directory_path() / ("kcored_partitions_" + std::to_string(stamp));
    fs::create_directories(base);

    int chunk = n / workers;
    int extra = n % workers;
    for (int i = 0; i < workers; ++i)
    {
        int workLoad = (i == workers - 1) ? (chunk + extra) : chunk;
        int offset = i * chunk;
        std::ofstream fout(base / (std::to_string(i) + ".txt"));
        for (int j = 0; j < workLoad; ++j)
        {
            int u = offset + j;
            if (u >= n)
                break;
            for (int v : g.neighbor[u])
                fout << u << " " << v << "\n";
        }
    }
    return base.string();
}

static KcoredPartitionInfo prepare_kcored_partitions(const Graph& g, const std::string& partition_dir_env)
{
    KcoredPartitionInfo info;
    info.path = partition_dir_env;
    info.num_workers = 1;
    info.is_temporary = false;

    if (!partition_dir_env.empty())
    {
        auto indices = discover_partition_indices(partition_dir_env);
        int workers = contiguous_worker_count(indices);
        if (workers > 0)
        {
            info.num_workers = workers;
            return info;
        }
    }

    info.path = materialize_partitions_from_graph(g, 1);
    info.num_workers = 1;
    info.is_temporary = true;
    return info;
}

static std::vector<int> round_coreness_to_int(const std::vector<double>& vals)
{
    std::vector<int> out(vals.size());
    for (size_t i = 0; i < vals.size(); ++i)
    {
        int v = static_cast<int>(std::lround(vals[i]));
        out[i] = std::max(0, v);
    }
    return out;
}

struct BaselineOutcome
{
    std::vector<int> coreness;
    std::string partition_dir;
    int num_workers = 1;
    bool used_temp_partitions = false;
    kcored::RunStatistics stats{};
    bool success = false;
};

static BaselineOutcome run_kcored_baseline(const Graph& g,
                                           double epsilon,
                                           double split,
                                           double psi,
                                           double lambda_coord,
                                           int bias_factor,
                                           const std::string& partition_dir_env)
{
    BaselineOutcome out;
    auto parts = prepare_kcored_partitions(g, partition_dir_env);
    out.partition_dir = parts.path;
    out.num_workers = parts.num_workers;
    out.used_temp_partitions = parts.is_temporary;

    kcored::RunConfig cfg;
    cfg.n = g.num_nodes();
    cfg.psi = psi;
    cfg.epsilon = epsilon;
    cfg.split = split;
    cfg.bias = true;
    cfg.biasFactor = bias_factor;
    cfg.noise = true;
    cfg.numWorkers = parts.num_workers;
    cfg.partitionDir = parts.path;
    cfg.lambdaCoordinator = lambda_coord;

    try
    {
        auto res = kcored::RunSingle(cfg, 42);
        out.stats = res.stats;
        out.coreness = round_coreness_to_int(res.coreNumbers);
        out.success = true;
    }
    catch (const std::exception& e)
    {
        std::cerr << "[baseline] kcored RunSingle failed: " << e.what() << std::endl;
    }

    if (parts.is_temporary)
    {
        std::error_code ec;
        std::filesystem::remove_all(parts.path, ec);
    }

    return out;
}

struct GlobalMetrics
{
    double mae = 0.0;
    double rmse = 0.0;
    double acc_k = 0.0;
    double mean_factor = std::numeric_limits<double>::quiet_NaN();
    double p80_factor = std::numeric_limits<double>::quiet_NaN();
    double p95_factor = std::numeric_limits<double>::quiet_NaN();
};

static bool compute_global_metrics(const std::vector<int> &true_vals,
                                   const std::vector<int> &pred_vals,
                                   int k,
                                   GlobalMetrics &out)
{
    if (true_vals.size() != pred_vals.size() || true_vals.empty())
    {
        return false;
    }

    const size_t n = true_vals.size();
    double mae = 0.0, mse = 0.0;
    size_t agree_k = 0;
    std::vector<double> factors;
    factors.reserve(n);

    for (size_t i = 0; i < n; ++i)
    {
        double err = std::abs(true_vals[i] - pred_vals[i]);
        mae += err;
        mse += err * err;
        if ((true_vals[i] >= k) == (pred_vals[i] >= k))
            ++agree_k;

        int true_c = std::max(true_vals[i], 1);
        int pred_c = std::max(pred_vals[i], 1);
        double max_val = static_cast<double>(std::max(pred_c, true_c));
        double min_val = static_cast<double>(std::min(pred_c, true_c));
        factors.push_back(max_val / min_val);
    }

    mae /= n;
    mse /= n;
    double rmse = std::sqrt(mse);
    double acc_k = static_cast<double>(agree_k) / n;

    std::sort(factors.begin(), factors.end());
    double sum = 0.0;
    size_t cnt = 0;
    for (double v : factors)
    {
        if (!std::isinf(v))
        {
            sum += v;
            ++cnt;
        }
    }
    double mean_f = cnt ? sum / cnt : std::numeric_limits<double>::quiet_NaN();

    auto percentile = [&](double q) {
        if (factors.empty()) return std::numeric_limits<double>::quiet_NaN();
        size_t idx = static_cast<size_t>(std::ceil(q * factors.size()));
        if (idx) --idx;
        return factors[idx];
    };

    out.mae = mae;
    out.rmse = rmse;
    out.acc_k = acc_k;
    out.mean_factor = mean_f;
    out.p80_factor = percentile(0.80);
    out.p95_factor = percentile(0.95);
    return true;
}


// Algorithm 2 variant: fixed mapping from pre-uploaded noisy arrays
std::vector<int> estimate_coreness_LDP_fixed_mapping_from_preuploads(
    Graph &g,
    const std::vector<int>& noisy_degree,
    const std::vector<int>& noisy_hindex,
    int max_iter,
    double momentum)
{
    const int n = g.num_nodes();
    if ((int)noisy_degree.size() != n || (int)noisy_hindex.size() != n) {
        std::cerr << "[ERROR] preuploads size mismatch in fixed_mapping_from_preuploads\n";
        return std::vector<int>();
    }

    // start from noisy H-index as the current estimate
    std::vector<double> coreness(n);
    for (int u = 0; u < n; ++u)
        coreness[u] = std::max(0, noisy_hindex[u]);

    int Dmax = 0;
    for (int d : noisy_degree) Dmax = std::max(Dmax, d);

    std::vector<double> next(n);
    std::vector<std::vector<double>> degree_bins(Dmax + 1);

    for (int it = 0; it < max_iter; ++it)
    {
        for (auto &bin : degree_bins) bin.clear();

        // (a) bucket by noisy degree
        for (int u = 0; u < n; ++u)
        {
            int d = noisy_degree[u];
            if (d < 0) d = 0; if (d > Dmax) d = Dmax;
            degree_bins[d].push_back(coreness[u]);
        }

        // (b) mapping via bucket averages
        std::vector<double> mapping(Dmax + 1, 0.0);
        for (int d = 0; d <= Dmax; ++d)
        {
            if (!degree_bins[d].empty())
            {
                double sum = 0.0;
                for (double v : degree_bins[d]) sum += v;
                mapping[d] = sum / degree_bins[d].size();
            }
        }

        // (c) apply mapping with momentum
        bool changed = false;
        for (int u = 0; u < n; ++u)
        {
            int d = noisy_degree[u];
            if (d < 0) d = 0; if (d > Dmax) d = Dmax;
            double new_est = mapping[d];
            next[u] = momentum * new_est + (1.0 - momentum) * coreness[u];
            if (std::fabs(next[u] - coreness[u]) > 1e-6) changed = true;
        }

        coreness.swap(next);
        if (!changed) break;
    }

    std::vector<int> result(n);
    for (int u = 0; u < n; ++u)
        result[u] = std::max(0, (int)std::round(coreness[u]));
    return result;
}


// Algorithm 2: iterative mapping (kcoreM)
std::vector<int> estimate_coreness_LDP_iterative_mapping(
    Graph &g,
    double eps1,
    double eps2,
    int seed,
    int max_iter,
    double momentum)
{
    const int n = g.num_nodes();
    std::mt19937 rng(seed);

    // round 1: collect noisy degree
    std::vector<int> real_deg(n), noisy_deg(n);
    for (int u = 0; u < n; ++u)
    {
        real_deg[u] = (int)g.neighbor[u].size();
        int noise = sample_symmetric_geometric(eps1, rng);
        noisy_deg[u] = std::max(0, real_deg[u] + noise);
    }

    // round 2: neighbors compute noisy H-index and upload
    std::vector<std::vector<int>> neigh_noisy_deg(n);
    for (int u = 0; u < n; ++u)
    {
        for (int v : g.neighbor[u])
            neigh_noisy_deg[u].push_back(noisy_deg[v]);
    }

    std::vector<double> coreness(n);
    for (int u = 0; u < n; ++u)
    {
        int h = compute_h_index(neigh_noisy_deg[u]);
        int noise = sample_symmetric_geometric(eps2, rng);
        coreness[u] = std::max(0, h + noise);
    }

    // step 3: dynamic mapping iteration based on rounded coreness keys
    const double smooth_alpha = 0.0;
    std::vector<double> prev_mapping;

    std::vector<int> keys(n);
    std::vector<double> next(n);
    
    for (int it = 0; it < max_iter; ++it)
    {
        // use rounded coreness as the key for buckets
        for (int u = 0; u < n; ++u)
            keys[u] = std::max(0, (int)std::round(coreness[u]));

        int Kmax = *std::max_element(keys.begin(), keys.end());
        std::vector<std::vector<double>> bins(Kmax + 1);

        // (a) bucket by key
        for (int u = 0; u < n; ++u)
            bins[keys[u]].push_back(coreness[u]);

        // (b) fill non-empty buckets with their mean
        std::vector<double> mapping(Kmax + 1, std::numeric_limits<double>::quiet_NaN());
        for (int k = 0; k <= Kmax; ++k)
        {
            if (!bins[k].empty())
            {
                // mean mapping; other robust stats also work
                double sum = 0.0;
                for (double v : bins[k])
                    sum += v;
                mapping[k] = sum / bins[k].size();
            }
        }

        // (b-1) fill empty buckets via nearest or interpolation
        auto is_nan = [](double x)
        { return std::isnan(x); };

        std::vector<int> Lnear(Kmax + 1, -1), Rnear(Kmax + 1, -1);
        int last = -1;
        for (int k = 0; k <= Kmax; ++k)
        {
            if (!is_nan(mapping[k]))
                last = k;
            Lnear[k] = last;
        }
        last = -1;
        for (int k = Kmax; k >= 0; --k)
        {
            if (!is_nan(mapping[k]))
                last = k;
            Rnear[k] = last;
        }

        for (int k = 0; k <= Kmax; ++k)
        {
            if (!is_nan(mapping[k]))
                continue;
            int L = Lnear[k], R = Rnear[k];
            if (L != -1 && R != -1 && L != R)
            {
                double w = double(k - L) / double(R - L);
                mapping[k] = (1.0 - w) * mapping[L] + w * mapping[R];
            }
            else if (L != -1)
            {
                mapping[k] = mapping[L];
            }
            else if (R != -1)
            {
                mapping[k] = mapping[R];
            }
            else
            {
                mapping[k] = 0.0;
            }
        }

        // (b-2) enforce mapping[k] <= k
        for (int k = 0; k <= Kmax; ++k)
            if (mapping[k] > k)
                mapping[k] = (double)k;

        // (b-3) optional cross-iteration smoothing
        if (smooth_alpha > 0.0 && !prev_mapping.empty())
        {
            int M = std::max((int)mapping.size(), (int)prev_mapping.size());
            mapping.resize(M, mapping.back());
            prev_mapping.resize(M, prev_mapping.back());
            for (int k = 0; k < M; ++k)
                mapping[k] = smooth_alpha * mapping[k] + (1.0 - smooth_alpha) * prev_mapping[k];
        }
        prev_mapping = mapping;

        // (c) apply mapping with momentum
        std::vector<double> next = coreness;
        bool changed = false;
        for (int u = 0; u < n; ++u)
        {
            int key = keys[u];
            if (key < 0)
                key = 0;
            if (key > (int)mapping.size() - 1)
                key = (int)mapping.size() - 1; // clamp
            double new_est = mapping[key];
            next[u] = momentum * new_est + (1.0 - momentum) * coreness[u];

            // use float threshold to decide convergence
            if (std::fabs(next[u] - coreness[u]) > 1e-6)
                changed = true;
        }

        coreness.swap(next);
        if (!changed)
            break;
    }

    // final rounding
    std::vector<int> result(n);
    for (int u = 0; u < n; ++u)
        result[u] = std::max(0, (int)std::round(coreness[u]));
    return result;
}





// degree-segment aggregation helpers
struct DegreeMetrics {
    double mae = 0.0;
    double rmse = 0.0;
    double mean_factor = std::numeric_limits<double>::quiet_NaN();
    double p80_factor = std::numeric_limits<double>::quiet_NaN();
    double p95_factor = std::numeric_limits<double>::quiet_NaN();
    size_t count = 0;
};

static bool compute_degree_metrics(const Graph &g,
                                   const std::vector<int> &true_vals,
                                   const std::vector<int> &pred_vals,
                                   int k,
                                   bool ge,
                                   DegreeMetrics &out)
{
    if (true_vals.size() != pred_vals.size()) return false;
    const size_t n = true_vals.size();
    double mae = 0.0, mse = 0.0;
    std::vector<double> factors;
    size_t cnt = 0;
    for (size_t i = 0; i < n; ++i) {
        if (ge) {
            if (g.degree[i] < k) continue;
        } else {
            if (g.degree[i] >= k) continue;
        }
        ++cnt;
        double err = std::abs(true_vals[i] - pred_vals[i]);
        mae += err;
        mse += err * err;
        int true_c = std::max(true_vals[i], 1);
        int pred_c = std::max(pred_vals[i], 1);
        double mx = static_cast<double>(std::max(true_c, pred_c));
        double mn = static_cast<double>(std::min(true_c, pred_c));
        factors.push_back(mx / mn);
    }
    if (cnt == 0) return false;
    mae /= cnt;
    mse /= cnt;
    double rmse = std::sqrt(mse);
    double sum = 0.0; size_t valid = 0;
    for (double v : factors) if (!std::isinf(v)) { sum += v; ++valid; }
    double mean_f = valid ? sum / valid : std::numeric_limits<double>::quiet_NaN();
    std::sort(factors.begin(), factors.end());
    auto percentile = [&](double q) {
        if (factors.empty()) return std::numeric_limits<double>::quiet_NaN();
        size_t idx = static_cast<size_t>(std::ceil(q * factors.size()));
        if (idx > 0) --idx;
        return factors[idx];
    };
    double p80 = percentile(0.80);
    size_t idx95 = static_cast<size_t>(std::ceil(0.95 * factors.size()));
    if (idx95 > 0) --idx95;
    double p95 = factors[idx95];
    out.mae = mae;
    out.rmse = rmse;
    out.mean_factor = mean_f;
    out.p80_factor = p80;
    out.p95_factor = p95;
    out.count = cnt;
    return true;
}

struct AlgoView { const char* name; const std::vector<int>* pred; };

static void print_degree_segment_comparison(const Graph &g,
                                            const std::vector<int> &real,
                                            const std::vector<AlgoView> &algos)
{
    std::cout << "\n=== 度数分段对比（算法2/3/5/8） ===\n";
    struct Segment { const char* label; bool ge; int k; };
    std::vector<Segment> segs = {
        {"degree < 30", false, 30},
        {"degree < 50", false, 50},
        {"degree < 80", false, 80},
        {"degree < 100", false, 100},
        {"degree < 120", false, 120},
        {"degree < 150", false, 150},
        {"degree < 200", false, 200},
        {"degree < 250", false, 250},
        {"degree < 300", false, 300},
        {"degree < 350", false, 350},
        {"degree < 400", false, 400},
        {"degree ≥ 30", true, 30},
        {"degree ≥ 50", true, 50},
        {"degree ≥ 80", true, 80},
        {"degree ≥ 100", true, 100},
        {"degree ≥ 120", true, 120},
        {"degree ≥ 150", true, 150},
        {"degree ≥ 200", true, 200},
        {"degree ≥ 250", true, 250},
        {"degree ≥ 300", true, 300},
        {"degree ≥ 350", true, 350},
        {"degree ≥ 400", true, 400},
    };

    for (const auto &seg : segs) {
        std::cout << "\n-- " << seg.label << " --" << std::endl;
        std::cout << std::left
                  << std::setw(14) << "算法"
                  << std::setw(12) << "MAE"
                  << std::setw(12) << "RMSE"
                  << std::setw(14) << "MeanFactor"
                  << std::setw(12) << "P80"
                  << std::setw(12) << "P95"
                  << std::setw(8)  << "N" << std::endl;
        std::cout << std::string(72, '-') << std::endl;

        // compute N for this degree range (independent of algorithms)
        size_t countN = 0;
        for (size_t i = 0; i < real.size(); ++i) {
            if (seg.ge ? (g.degree[i] >= seg.k) : (g.degree[i] < seg.k)) ++countN;
        }

        for (const auto &av : algos) {
            if (!av.pred || av.pred->size() != real.size()) continue;
            DegreeMetrics m; bool ok = compute_degree_metrics(g, real, *av.pred, seg.k, seg.ge, m);
            std::cout << std::left << std::setw(14) << av.name;
            if (ok) {
                std::cout << std::setw(12) << std::fixed << std::setprecision(4) << m.mae
                          << std::setw(12) << std::fixed << std::setprecision(4) << m.rmse
                          << std::setw(14) << std::fixed << std::setprecision(4) << m.mean_factor
                          << std::setw(12) << std::fixed << std::setprecision(4) << m.p80_factor
                          << std::setw(12) << std::fixed << std::setprecision(4) << m.p95_factor
                          << std::setw(8)  << m.count << std::endl;
            } else {
                std::cout << std::setw(12) << "N/A"
                          << std::setw(12) << "N/A"
                          << std::setw(14) << "N/A"
                          << std::setw(12) << "N/A"
                          << std::setw(12) << "N/A"
                          << std::setw(8)  << countN << std::endl;
            }
        }
    }
}

struct AlgoSummary {
    const char* name;
    const std::vector<int>* pred;
    double runtime_seconds = std::numeric_limits<double>::quiet_NaN();
};

static void print_global_summary_comparison(const std::vector<int> &real,
                                            const std::vector<AlgoSummary> &algos,
                                            int acc_k)
{
    if (algos.empty()) return;
    std::cout << "\n=== 全图指标对比（算法2/3/5/8） ===\n\n";
    std::cout << std::left
              << std::setw(14) << "算法"
              << std::setw(10) << "MAE"
              << std::setw(10) << "RMSE"
              << std::setw(12) << ("Acc@" + std::to_string(acc_k))
              << std::setw(14) << "MeanFactor"
              << std::setw(12) << "P80"
              << std::setw(12) << "P95"
              << std::setw(10) << "Time(s)" << std::endl;
    std::cout << std::string(94, '-') << std::endl;

    for (const auto &algo : algos)
    {
        if (!algo.pred || algo.pred->size() != real.size())
        {
            std::cout << std::left << std::setw(14) << algo.name
                      << std::setw(10) << "N/A"
                      << std::setw(10) << "N/A"
                      << std::setw(12) << "N/A"
                      << std::setw(14) << "N/A"
                      << std::setw(12) << "N/A"
                      << std::setw(12) << "N/A"
                      << std::setw(10) << "N/A" << std::endl;
            continue;
        }
        GlobalMetrics gm;
        if (!compute_global_metrics(real, *algo.pred, acc_k, gm))
        {
            std::cout << std::left << std::setw(14) << algo.name
                      << std::setw(10) << "N/A"
                      << std::setw(10) << "N/A"
                      << std::setw(12) << "N/A"
                      << std::setw(14) << "N/A"
                      << std::setw(12) << "N/A"
                      << std::setw(12) << "N/A"
                      << std::setw(10) << "N/A" << std::endl;
            continue;
        }
        std::cout << std::left << std::setw(14) << algo.name
                  << std::setw(10) << std::fixed << std::setprecision(4) << gm.mae
                  << std::setw(10) << std::fixed << std::setprecision(4) << gm.rmse
                  << std::setw(12) << std::fixed << std::setprecision(2) << gm.acc_k * 100.0
                  << std::setw(14) << std::fixed << std::setprecision(4) << gm.mean_factor
                  << std::setw(12) << std::fixed << std::setprecision(4) << gm.p80_factor
                  << std::setw(12) << std::fixed << std::setprecision(4) << gm.p95_factor
                  << std::setw(10) << std::fixed << std::setprecision(4) << algo.runtime_seconds
                  << std::endl;
    }
}


// Algorithm 7: kcoreHB-th (mix of algo6 and algo2 with degree threshold)
std::vector<int> estimate_coreness_LDP_algo26_hybrid(
    Graph &g,
    double epsilon,
    int seed,
    int max_iter,
    double momentum,
    int degree_threshold)
{
    const int n = g.num_nodes();
    std::mt19937 rng(seed);

    // share epsilon across the three internal stages
    double eps1 = epsilon / 3.0; // noisy degree
    double eps2 = epsilon / 3.0; // H-index based path (algo6)
    double eps_iter2 = epsilon - eps1 - eps2; // remaining to algo2 inner noise, used inside

    // Step 1: upload noisy degree once
    std::vector<int> degrees(n), noisy_degree(n);
    for (int u = 0; u < n; ++u) {
        degrees[u] = (int)g.neighbor[u].size();
        int noise = sample_symmetric_geometric(eps1, rng);
        noisy_degree[u] = std::max(0, degrees[u] + noise);
    }

    // Step 2a: run the algo6 path for low-degree nodes
    std::vector<std::pair<int,int>> rules = {{8,1},{15,2},{30,3},{80,4}};
    auto coreness_6 = estimate_coreness_LDP_hindex_iter_threshold(g, eps2, seed + 17, max_iter, rules);

    // Step 2b: run the algo2 path with eps_iter2 split in half
    auto coreness_2 = estimate_coreness_LDP_iterative_mapping(
        g, eps_iter2 / 2.0, eps_iter2 / 2.0, seed + 23, max_iter, momentum);

    // Step 3: pick which coreness to return per node
    std::vector<int> final_result(n);
    for (int u = 0; u < n; ++u) {
        if (noisy_degree[u] > degree_threshold) final_result[u] = coreness_2[u];
        else                                    final_result[u] = coreness_6[u];
    }
    return final_result;
}




// Algorithm 8: kcoreHB (mix of algo2 and algo3 with high-node iteration)
std::vector<int> estimate_coreness_LDP_algo23_hybrid(
    Graph &g,
    double epsilon,
    int seed,
    int max_iter,
    double momentum,
    int degree_threshold)
{
    // first two rounds upload global noisy degree / noisy H-index, then freeze lows and iterate highs

    const int n = g.num_nodes();
    std::mt19937 rng(seed);

    // epsilon split
    double eps1 = epsilon / 3.0;
    double eps2 = epsilon / 3.0;
    double eps3 = epsilon - eps1 - eps2;

    // Round 1: noisy degree upload
    std::vector<int> degrees(n), noisy_degree(n);
    for (int u = 0; u < n; ++u) {
        degrees[u] = (int)g.neighbor[u].size();
        int z = sample_symmetric_geometric(eps1, rng);
        noisy_degree[u] = std::max(0, degrees[u] + z);
    }

    // Round 2: noisy H-index based on noisy degrees
    std::vector<int> noisy_hindex(n);
    for (int u = 0; u < n; ++u) {
        std::vector<int> neigh_noisy;
        neigh_noisy.reserve(g.neighbor[u].size());
        for (int v : g.neighbor[u]) neigh_noisy.push_back(noisy_degree[v]);
        int h = compute_h_index(neigh_noisy);
        int z = sample_symmetric_geometric(eps2, rng);
        noisy_hindex[u] = std::max(0, h + z);
    }

    // Path A (kcoreM): full-graph mapping from uploaded arrays
    auto outM = estimate_coreness_LDP_fixed_mapping_from_preuploads(
        g, noisy_degree, noisy_hindex, max_iter, momentum);

    // High-degree path (kcoreH): keep lows fixed, iterate only highs
    std::vector<char> is_low(n, 0);
    for (int u = 0; u < n; ++u) is_low[u] = (noisy_degree[u] <= degree_threshold) ? 1 : 0;

    std::vector<int> init = outM;
    for (int u = 0; u < n; ++u) {
        if (!is_low[u]) init[u] = noisy_hindex[u];
    }

    double epsH = (max_iter > 0) ? eps3 / std::max(1, max_iter) : 0.0;
    auto high_updated = estimate_coreness_LDP_hindex_high_only(
        g, epsH, seed + 202, max_iter, init, is_low);

    std::vector<int> final_result(n);
    for (int u = 0; u < n; ++u)
        final_result[u] = is_low[u] ? outM[u] : high_updated[u];
    return final_result;
}



// run the full comparison suite (2/3/5/8)
void run_kcore_suite(double eps, const std::string& dataset, int rounds, double v_ratio, int algo_switch_)
{
    // set globals
    Eps = eps;
    num_rounds = rounds;
    vertex_ratio = v_ratio;
    algo_switch = algo_switch_;

    RR_time = 0; server_side_time = 0; naive_server_side = 0;

    bool is_bipartite = false;
    Graph g(dataset, is_bipartite);

    std::cout << "vertex_ratio = " << vertex_ratio << std::endl;
    std::cout << "n = " << g.num_nodes() << std::endl;

    std::unordered_map<std::string, long int> dataset_map = {
        {"../graphs/email.edges", 727044},
        {"../graphs/youtube.edges", 3056386},
        {"../graphs/wiki.edges", 9203519},
        {"../graphs/gow.edges", 2273138},
        {"../graphs/dblp.edges", 2224385},
        {"../graphs/gplus.edges", 1073677742},
        {"../graphs/skitter.edges", 28769868},
        {"../graphs/imdb.edges", 3856982376},
        {"../graphs/lj.edges", 177820130},
        {"../graphs/orkut.edges", 627584181}
    };

    tri_cnt.resize(g.num_nodes());
    fill(tri_cnt.begin(), tri_cnt.end(), 0);

    auto it = dataset_map.find(dataset);
    if (it != dataset_map.end()) {
        real_val = it->second;
    } else {
        sup_vector.resize(g.num_nodes());
        std::cout << "len = " << sup_vector.size() << std::endl;
        real_val = 0;
    }
    printf("coreness = %lld\n", real_val);

    const char* partition_dir_env = std::getenv("KCORE_PARTITION_DIR");
    std::string partition_dir = partition_dir_env ? partition_dir_env : "";

    double kcorem_momentum = 1.0;
    if (const char* env = std::getenv("KCOREM_MOMENTUM")) {
        try {
            kcorem_momentum = std::stod(env);
        } catch (...) {
            kcorem_momentum = 1.0;
        }
    }

    // exact coreness
    std::vector<int> real_coreness = compute_coreness(g);

    std::vector<int> H_est2 = real_coreness;
    std::vector<int> H_est3 = real_coreness;
    std::vector<int> H_est5 = real_coreness;
    std::vector<int> H_est8 = real_coreness;
    double time2 = 0.0, time3 = 0.0, time5 = 0.0, time8 = 0.0;

    std::cout << "\n=== 2. kcoreM ===" << std::endl;
    {
        double eps1 = eps / 2.0;
        double eps2 = eps - eps1;
        std::cout << "参数: ε1=" << eps1 << ", ε2=" << eps2
                  << ", seed=42, max_iter=100, momentum=" << kcorem_momentum << std::endl;
        auto t0 = std::chrono::high_resolution_clock::now();
        H_est2 = estimate_coreness_LDP_iterative_mapping(g, eps1, eps2, 42, 100, kcorem_momentum);
        auto t1 = std::chrono::high_resolution_clock::now();
        time2 = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1e6;
    }

    std::cout << "\n=== 3. kcoreH ===" << std::endl;
    {
        double eps_iter = eps;
        std::cout << "参数: ε=" << eps_iter << ", seed=42, max_iter=5" << std::endl;
        auto t0 = std::chrono::high_resolution_clock::now();
        H_est3 = estimate_coreness_LDP_hindex_iteration(g, eps_iter, 42, 5);
        auto t1 = std::chrono::high_resolution_clock::now();
        time3 = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1e6;
    }

    std::cout << "\n=== 5*. baseline (distributed) ===" << std::endl;
    {
        double eps_local = eps;
        double split = 0.8;
        double psi = 0.5;
        double eta = 3.625;
        int bias_factor = 8;
        auto t0 = std::chrono::high_resolution_clock::now();
        auto baseline = run_kcored_baseline(g, eps_local, split, psi, eta, bias_factor, partition_dir);
        auto t1 = std::chrono::high_resolution_clock::now();
        time5 = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1e6;
        if (baseline.success)
        {
            H_est5.swap(baseline.coreness);
            // Prefer kcored timing if available
            double kcored_time = baseline.stats.preprocessingSeconds + baseline.stats.algorithmSeconds;
            if (kcored_time > 0.0)
                time5 = kcored_time;
        }
        std::cout << "参数: ε=" << eps_local << ", split=" << split
                  << ", ψ=" << psi << ", η=" << eta
                  << ", bias_factor=" << bias_factor
                  << ", workers=" << baseline.num_workers
                  << ", partition_dir=" << (baseline.partition_dir.empty() ? "<auto>" : baseline.partition_dir);
        if (baseline.used_temp_partitions)
            std::cout << " (auto-generated from input graph)";
        std::cout << std::endl;
        if (baseline.success)
        {
            std::cout << "kcored耗时: 预处理 " << baseline.stats.preprocessingSeconds
                      << " 秒, 算法 " << baseline.stats.algorithmSeconds << " 秒" << std::endl;
        }
        else
        {
            std::cout << "[WARN] kcored baseline 运行失败，结果为空" << std::endl;
            H_est5.assign(real_coreness.size(), 0);
        }
    }

    std::cout << "\n=== 8. kcoreHB ===" << std::endl;
    {
        double eps_per = eps / 3.0;
        std::cout << "参数: ε_total=" << eps << " (ε1=ε2=ε3=" << eps_per << ")"
                  << ", seed=42, max_iter=4, momentum=1.0, degree_threshold=30" << std::endl;
        auto t0 = std::chrono::high_resolution_clock::now();
        H_est8 = estimate_coreness_LDP_algo23_hybrid(g, eps, 42, 4, 1.0, 30);
        auto t1 = std::chrono::high_resolution_clock::now();
        time8 = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1e6;
    }

    {
        std::vector<AlgoView> algos;
        algos.push_back({"kcoreM", &H_est2});
        algos.push_back({"kcoreH", &H_est3});
        algos.push_back({"baseline", &H_est5});
        algos.push_back({"kcoreHB", &H_est8});
        print_degree_segment_comparison(g, real_coreness, algos);
    }

    {
        std::vector<AlgoSummary> all_algos;
        all_algos.push_back({"kcoreM", &H_est2, time2});
        all_algos.push_back({"kcoreH", &H_est3, time3});
        all_algos.push_back({"baseline", &H_est5, time5});
        all_algos.push_back({"kcoreHB", &H_est8, time8});
        print_global_summary_comparison(real_coreness, all_algos, 4);
    }

    std::cout << "\n=== 运行时间总结（算法2/3/5/8） ===" << std::endl;
    std::cout << "2. kcoreM: " << time2 << " 秒" << std::endl;
    std::cout << "3. kcoreH: " << time3 << " 秒" << std::endl;
    std::cout << "5*. baseline(distributed): " << time5 << " 秒" << std::endl;
    std::cout << "8. kcoreHB: " << time8 << " 秒" << std::endl;
}
