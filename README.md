# PBTM-CORDIC: Prefix Binary Tree Mapping and Table Generation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18780760.svg)](https://doi.org/10.5281/zenodo.18780760)

This repository provides the MATLAB implementation for the **Prefix Binary Tree Mapping (PBTM)** algorithm, as described in the paper *"CORDIC-based Computation of Arcsine/Arccosine Integrated with Prefix Binary Tree Mapping"* (submitted to IEEE TCAS-II).

## Project Overview

The provided scripts automate the generation of non-uniform solution space boundaries and the construction of the prefix search tree. This process is essential for bypassing the first  iterations of the CORDIC algorithm while maintaining full precision.

### Core Components

1. **Rounding Boundary Search**: Implements the Hamming-weight optimization for decision boundaries.
2. **Tree Construction**: A recursive DFS algorithm that builds the Prefix Binary Tree.
3. **Hardware Parameter Extraction**: Generates  (prefix widths),  (offsets), and ROM initialization files (.COE).

---

## Algorithm Workflow & Code Mapping

The generation procedure follows a strict logical flow corresponding to the sections in the manuscript:

### Step 1: Solution Space & Boundary Optimization

**Script**: `search_rounding_boundary.m`

* **Mathematical Basis**: Equations (7), (8).
* **Function**:
* Calculates the discrete candidate vectors  and their convergence ranges .
* Identifies overlaps between adjacent solution spaces.
* **Boundary Selection**: The function `find_rounding_boundary` implements Eq. (7). It finds the bit position  (the most significant differing bit) and constructs  to maximize trailing zeros, minimizing hardware comparison cost.


* **Output**: `invc_param_kn.mat` (contains the optimized boundary table `intlv_table`).

### Step 2: Prefix Binary Tree Construction

**Script**: `create_PBT.m`

* **Mathematical Basis**: Section III-B (Prefix Binary Tree Mapping).
* **Function**:
* **Recursive DFS**: The `search_tree` function performs a Depth-First Search to segment the input domain.
* **Node Constraint**: It enforces the rule that each leaf node covers at most two solution spaces, ensuring single-comparator latency.
* **Parameter Extraction**: Calculates the variable prefix width and the level-dependent offset used in Eq. (9).


* **Output**:
* `pt.txt`: Tree structure for visualization.
* Verilog Code Snippets: Priority logic and subtraction units for FPGA implementation.
* `.coe` files: Memory initialization for the Solution Space Table.


---

## How to Run

1. **Environment**: MATLAB R2022b or later (requires Fixed-Point Designer).
2. **Execution Order**:
* Run `search_rounding_boundary.m` first. This generates the necessary boundary parameters.
* Run `create_PBT.m`. This will load the parameters, build the tree, and output the hardware configuration.


3. **Visualization**: The scripts will plot the solution space ranges and the resulting tree structure (`treeplot`).