

This repository contains the code that is required to reproduce
the four coreness estimation algorithms used in the paper:

* **Algorithm 2 – kcoreM** 
* **Algorithm 3 – kcoreH** 
* **Algorithm 5 – baseline / distributed**
* **Algorithm 8 – kcoreHB** 

All legacy triangle-counting and unused helper code was removed so that
`counting.cpp`, `graph.cpp`, `utility.cpp`, and `main.cpp` are the only
source files needed to build and run the experiments.

## Build

```bash
make clean && make
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


* `epsilon` – privacy budget used throughout the suite.
* `dataset` – path to the edge-list file.
* `num_rounds` – optional (default `40`), used by the iterative pieces.
* `vertex_ratio` – optional (default `0.1`), allows vertex sampling.
* `algo_switch` – accepted for compatibility; the current driver always
  runs algorithms 2/3/5/8 sequentially and prints their metrics.

Example:

```
KCORE_PARTITION_DIR=graphs/<name>/<name>_partitioned_80 \
./counting <epsilon> graphs/<name>/<name>.txt <rounds> <vertex_ratio> <algo_switch>




KCORE_PARTITION_DIR=graphs/facebook/facebook_partitioned_80 \
./counting 1.0 graphs/facebook/facebook.txt 1 1 1
```


Environment variables:

* `KCORE_PARTITION_DIR` – optional path that gets reported alongside
  algorithm 5 runs (useful when you keep partitioned datasets elsewhere).
* `KCORE_H_INIT_RATIO` – splits the privacy budget between initialization
  noise and iterative noise for Algorithm 3.
