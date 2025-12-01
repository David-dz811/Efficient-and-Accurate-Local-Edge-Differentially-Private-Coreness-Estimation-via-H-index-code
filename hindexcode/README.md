# Minimal k-core LDP Suite (Algorithms 2/3/5/8)

This repository now contains only the code that is required to reproduce
the four coreness estimation algorithms used in the paper:

* **Algorithm 2 – kcoreM** (iterative mapping)
* **Algorithm 3 – kcoreH** (H-index iteration)
* **Algorithm 5 – baseline / distributed**
* **Algorithm 8 – kcoreHB** (hybrid of 2 and 3)

All legacy triangle-counting and unused helper code was removed so that
`counting.cpp`, `graph.cpp`, `utility.cpp`, and `main.cpp` are the only
source files needed to build and run the experiments.

## Build

```bash
make
```

The build produces a single binary `counting`. Run `make clean` to remove
any artefacts.

## Input format

Datasets use the standard edge-list format:

```
n m
u v
u v
...
```

`n` is the number of vertices and `m` the number of undirected edges. Each
subsequent line contains one edge (`0`-based vertex ids).

## Running algorithms 2/3/5/8

```
./counting <epsilon> <dataset> [num_rounds] [vertex_ratio] [algo_switch]
```

* `epsilon` – privacy budget used throughout the suite.
* `dataset` – path to the edge-list file.
* `num_rounds` – optional (default `40`), used by the iterative pieces.
* `vertex_ratio` – optional (default `0.1`), allows vertex sampling.
* `algo_switch` – accepted for compatibility; the current driver always
  runs algorithms 2/3/5/8 sequentially and prints their metrics.

Examples:

```
./counting 1.0 graphs/email.edges 10 1.0 0
```

```
KCORE_PARTITION_DIR=graphs/<name>/<name>_partitioned_80 \
./counting <epsilon> graphs/<name>/<name>.txt <rounds> <vertex_ratio> <algo_switch>
```

Environment variables:

* `KCORE_PARTITION_DIR` – optional path that gets reported alongside
  algorithm 5 runs (useful when you keep partitioned datasets elsewhere).
* `KCOREM_MOMENTUM` – overrides the default momentum (1.0) used by
  Algorithms 2 and 8.
* `KCORE_H_INIT_RATIO` – splits the privacy budget between initialization
  noise and iterative noise for Algorithm 3.

On completion the program prints per-algorithm summaries and per-degree
segment metrics so you can compare the estimators directly.
