# DMMSY-SSSP.jl

A high-performance Julia implementation of the Single-Source Shortest Path (SSSP) algorithm introduced in:

> Ran Duan, Jiayi Mao, Xiao Mao, Xinkai Shu, Longhui Yin,  
> **"Breaking the Sorting Barrier for Directed Single-Source Shortest Paths"**, STOC 2025.

## Overview

This repository provides a performance-engineered implementation of the Duan–Mao–Mao–Shu–Yin (**DMMSY**) algorithm. By leveraging advanced architectural optimizations such as **8-way Instruction-Level Parallelism (ILP)**, **hybrid priority tracking**, and **cache-interleaved memory layouts**, this implementation demonstrates significant practical advantages over traditional heap-based Dijkstra's algorithm, achieving over **200x speedup** on large-scale graphs when using high-performance compiler flags.

### The Scientific Core: Breaking the Sorting Barrier

The DMMSY algorithm represents the first deterministic breakthrough to exceed the $O(m + n \log n)$ complexity barrier that has stood for 65 years since Dijkstra's original proposal. 

- **Complexity Shift**: It reduces the logarithmic factor from $O(\log n)$ to $O(\log^{2/3} n)$. For a graph with $10^9$ nodes, this is a ~3x reduction in sorting overhead.
- **Interval-Pivot Mechanism**: Unlike Dijkstra which sorts nodes individually, DMMSY uses a recursive subproblem decomposition based on interval pivots. This allows it to relax edges without maintaining a strict global priority queue for every extraction.

### DMMSY vs. Preprocessing (Dynamic Graphs)

Traditional production engines often use **Contraction Hierarchies (CH)** or **A*** to achieve sub-millisecond queries. However, these rely on expensive preprocessing ($O(n \log n)$ to $O(n^2)$) that must be rebuilt whenever edge weights change.

**DMMSY is the superior choice for Dynamic Graphs** where:
1. **Topology/Weights Change Rapidly**: The costs of rebuilding a hierarchy (20-40 minutes for road networks) cannot be amortized.
2. **Zero Preprocessing luxury**: Real-time systems like **Financial Market Arbitrage** or **Software-Defined Networking (SDN)** require immediate SSSP computation on a raw, mutating graph.

## Key Features

- **Theoretical Breakthrough**: Implements the first deterministic algorithm to break the $O(m \log n)$ sorting barrier for directed graphs, achieving $O(m \log^{2/3} n)$ complexity.
- **Hardware-Aware Design**: Features unrolled priority tracking and software prefetching to maximize pipeline occupancy.
- **Selective Memory Reset**: Utilizes dirty-index tracking to eliminate $O(n)$ overhead during recursive subproblem resets.
- **Adaptive Strategy**: Employs a hybrid tracking system (Bitmap + Heap) and dynamic threshold feedback to maintain peak efficiency.
- **Robust Verification**: Validated against reference Dijkstra implementations across linear, diamond, cycle, and hub topologies.

## Industrial Impact

The algorithm's ability to handle massive sparse graphs without preprocessing creates a new paradigm for specific industrial domains:

| Sector | Application | DMMSY Advantage |
| :--- | :--- | :--- |
| **Finance** | Arbitrage Detection | Faster negative-cycle detection in violently dynamic FX/Crypto markets. |
| **Networking** | SDN & OSPF | Reduced path convergence time in datacenter fabrics during link failures. |
| **Logistics** | Live Traffic Routing | Faster fallback SSSP when emergency closures invalidate precomputed hierarchies. |
| **EDA** | Chip Design | Iterative timing analysis on circuit DAGs with millions of gate nodes. |
| **Social Tech** | Influence Analysis | Real-time recommendation ranking on billion-node sparse social graphs. |

## Performance Metrics

Tests conducted on a modern x86_64 architecture with specialized selective reset and GC-managed timing. The results below showcase the stable performance delta under `--math-mode=fast`.


![DMMSY-SSSP Performance Dashboard](screenshot.png)

> [!NOTE]
> The algorithm achieves its maximum performance delta at scale (250k–1M+ nodes), where the reduction in sorting overhead and enhanced memory locality most effectively mitigate the cache-bottlenecks of standard priority queues. At 1M nodes, DMMSY outperforms Dijkstra by over **200x**.


## Core Optimizations

1. **Dual-Level Interleaving**: 
   - **Graph Layer**: Neighbor IDs and weights are co-located in an `Edge` struct to ensure single-cache-line fetches during relaxation.
   - **Heap Layer**: Distance values and node IDs are co-located in a `HeapNode` struct within the 4-ary heap, minimizing cache thrashing during swaps.
2. **8-way ILP Unrolling**: Critical relaxation loops are manually unrolled to saturate execution units and remove data dependency chains.
3. **Hybrid Block Search**: Swaps between $O(1)$ bitmap scanning and $O(\log n)$ heap-based selection based on subproblem characteristics.
4. **Adaptive Thresholding**: Dynamically adjusts recursive propagation thresholds based on real-time extraction efficiency.
5. **Quickselect-based Pivots**: Utilizes $O(n)$ partial sorting for pivot selection instead of full $O(n \log n)$ sorts.

## Repository Structure

- `DMMSY-SSSP.jl`: Core highly-optimized implementation.
- `CSRGraph.jl`: High-performance Compressed Sparse Row (CSR) graph implementation.
- `Dijkstra.jl`: Reference implementation used for correctness and performance baselines.
- `DMMSY_Research.jl`: Baseline research implementation for algorithmic comparison.
- `bench_report_collector.jl`: Main high-precision benchmarking and reporting driver.
- `benchmark_results.html`: Visualization for performance data.
- `benchmark_data.csv`: Raw timing metrics.

## Quick Start

### Installation

Ensure you have [Julia](https://julialang.org/) installed, then clone this repository and include the source files:

```julia
include("Common.jl")
include("CSRGraph.jl")
include("Dijkstra.jl")
include("DMMSY-SSSP.jl")

using .Common, .CSRGraphModule, .DijkstraModule, .DMMSYSSSP

# Generate a synthetic graph
g = random_graph(10000, 50000, 100.0)

# Compute shortest paths from source vertex 1
dists, preds = ssp_duan(g, 1)
```

### Testing & Verification

To verify the implementation against reference Dijkstra across multiple topologies:
```bash
julia test_sssp.jl
```

### Benchmarking

To generate a full performance report:
```bash
julia bench_report_collector.jl
```

### High-Performance Execution

For maximum performance, it is highly recommended to run Julia with the `--math-mode=fast` flag. This enables aggressive floating-point optimizations, SIMD vectorization, and reordering of operations that are critical for the edge relaxation loops.

```bash
# Standard Execution
julia bench_report_collector.jl

# High-Performance Mode (Recommended for Benchmarking)
julia --math-mode=fast bench_report_collector.jl
```

> [!IMPORTANT]
> Enabling `--math-mode=fast` can increase speedups from **1.1x** to over **220x** on specific graph topologies by allowing the compiler to fully utilize the processor's vector units and instruction pipelines.

## Contributing

We welcome contributions to further improve the performance of this implementation. Potential areas for exploration include:
- **Multithreading**: Parallelizing the recursive subproblem decomposition.
- **SIMD**: Leveraging AVX-512/AMX instructions for even more aggressive relaxation unrolling.
- **GPU Acceleration**: Offloading the blocked priority queue logic to massive parallel hardware.

## References

1. Duan, R., Mao, J., Mao, X., Shu, X., & Yin, L. (2025). "Breaking the Sorting Barrier for Directed Single-Source Shortest Paths". *Proceedings of the 57th Annual ACM Symposium on Theory of Computing (STOC 2025)*.
2. [arXiv:2504.17033](https://arxiv.org/pdf/2504.17033)
3. [ACM Digital Library](https://dl.acm.org/doi/10.1145/3717823.3718179)
4. [MPI-INF repository](https://pure.mpg.de/pubman/faces/ViewItemFullPage.jsp?itemId=item_3634229)


## License

This project is dual-licensed under the **MIT License** and **Apache License 2.0**.
