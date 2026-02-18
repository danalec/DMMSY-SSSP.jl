# DMMSY-SSSP.jl

Experimental Julia implementation of the single-source shortest path (SSSP) algorithm from:

> Ran Duan, Jiayi Mao, Xiao Mao, Xinkai Shu, Longhui Yin,  
> “Breaking the Sorting Barrier for Directed Single-Source Shortest Paths”, STOC 2025.

The goal of this repository is to provide a readable and hackable version of the Duan–Mao–Mao–Shu–Yin algorithm for directed graphs with non‑negative edge weights, and to make it easy to compare against a standard Dijkstra baseline.

## Features

- `CSRGraph` type: cache-friendly compressed-sparse-row representation for directed graphs.
- `ssp_duan(g, s)`: DMMSY-style SSSP from a single source `s`.
- `reconstruct_path(pred, s, t)`: reconstruct an `s → t` path from the predecessor array.
- `dijkstra_ref(g, s)`: vanilla heap-based Dijkstra implementation for verification.
- `verify_correctness()`: small test suite comparing DMMSY against Dijkstra on several graphs.
- `benchmark_comparison(n, m, trials)`: micro-benchmark harness for random graphs.

The implementation is **research / prototype quality**: it mirrors the algorithmic ideas (BMSSP recursion, pivots, and blocked partial priority queues) but does not strictly match every invariant and constant-factor optimization in the paper.

## Installation

Clone or copy this file into a Julia project:

```bash
git clone https://github.com/danalec/DMMSY-SSSP.jl
```

## License

This project is dual-licensed under both the MIT License and the Apache License 2.0. You may choose either license at your option.

- [MIT License](LICENSE-MIT)
- [Apache License 2.0](LICENSE-APACHE)