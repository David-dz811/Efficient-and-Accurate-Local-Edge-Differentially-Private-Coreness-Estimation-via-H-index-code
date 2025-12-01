#pragma once

#include "graph.h"
#include <string>

// Entry point that runs algorithms 2 (kcoreM), 3 (kcoreH),
// 5 (baseline) and 8 (kcoreHB) on the provided dataset.
void run_kcore_suite(double eps,
                     const std::string& dataset,
                     int rounds,
                     double vertex_ratio,
                     int algo_switch);
