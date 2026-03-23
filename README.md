# CLuMP — Cluster-based Message-Passing Optimization

**Cluster-based Message-Passing (CLuMP) Optimization in Rugged Energy Landscapes**

Paolo Rissone¹, Stefan Boettcher², Alfonso Amendola³, Simone Sala³, and Federico Ricci-Tersenghi¹·⁴·⁵

> ¹ Department of Physics, Sapienza University of Rome, Rome, Italy  
> ² Department of Physics, Emory University, Atlanta, GA, USA  
> ³ Eni S.p.A., San Donato Milanese, Italy  
> ⁴ CNR-Nanotec, Rome unit, Rome, Italy  
> ⁵ INFN-Rome 1 section, Rome, Italy

---

## Overview

This repository contains the **source code** for CLuMP and the benchmark algorithms presented in the paper. The associated **data** (pre-generated graphs and results) are archived separately on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19048238.svg)](https://doi.org/10.5281/zenodo.19048238).

CLuMP is a novel message-passing-based optimization algorithm designed for spin-glass-like problems defined on sparse graphs (random regular graphs and 2D/3D lattices). It identifies frustrated clusters via belief propagation and proposes collective moves that reduce energy, outperforming standard heuristics such as Simulated Annealing and Extremal Optimization in rugged energy landscapes.

---

## Repository Structure

```
/
├── README.md
├── LICENSE
└── code/
    ├── Makefile                  # Build all binaries with a single command
    ├── GraphGenerator.c          # Generate random regular graphs (RRG) and 2D/3D lattices
    ├── CLuMP.c                   # CLuMP optimization algorithm
    ├── RCLuMP.c                  # Population-based variant of CLuMP (R-CLuMP)
    ├── SimulatedAnnealing.c      # Simulated Annealing baseline
    ├── PopulationAnnealing.c     # Population Annealing baseline
    ├── ExtremalOpt.c             # Extremal Optimization (EO / τEO) baseline
    ├── ResultsAnalysis.py        # Post-processing and figure generation
    └── [header files]            # Shared C headers
```

> **Data** (graphs, raw results, pre-computed outputs for paper figures) are hosted on Zenodo — see the link above.

---

## Requirements

### C code
- **GCC** (GNU Compiler Collection) and **make**. Pre-installed on most Linux systems.
  ```bash
  # Ubuntu/Debian
  sudo apt install gcc make

  # macOS (via Homebrew)
  brew install gcc make
  ```

### Python (analysis)
- Python 3.x
- Standard library modules (no install needed): `math`, `sys`, `os`, `re`, `itertools`, `io`, `random`, `operator`
- Third-party packages:
  ```bash
  pip install numpy pandas matplotlib seaborn scipy
  ```

---

## Compilation

From the `code/` directory, build all binaries at once:

```bash
cd code/
make
```

To build a single binary:

```bash
make CLuMP
make GraphGenerator
# etc.
```

To remove all compiled binaries:

```bash
make clean
```

> The Makefile detects your OS automatically and applies the correct flags (`-std=c99` on macOS, `-O2 -Wall` on both platforms). No manual editing required.

---

## Usage

### Step 1 — Generate a graph

```bash
./GraphGenerator [TYPE] [J] [N] [C]
```

| Parameter | Options | Description |
|-----------|---------|-------------|
| `TYPE` | `2`, `3`, `RRG` | Lattice dimension or Random Regular Graph |
| `J` | `Gaussian`, `PM1` | Coupling distribution |
| `N` | integer | Number of nodes. **Note:** for d=2,3 lattices `N` = nodes per side |
| `C` | integer | Connectivity (degree) **Note:** for d=2,3 `C` must be 4,6, respectively |

**Example:**
```bash
./GraphGenerator RRG Gaussian 1000 3
```

---

### Step 2 — Run an optimization algorithm

#### CLuMP
```bash
./CLuMP [MAXFRUS] [beta] [log2_iters] [RUN] [graph.conf]
```
| Parameter | Description |
|-----------|-------------|
| `MAXFRUS` | Maximum frustration threshold for cluster selection |
| `beta` | Inverse effective temperature (move acceptance) |
| `log2_iters` | Total iterations as a power of 2 |
| `RUN` | Number of independent runs |
| `graph.conf` | Input graph file |

#### R-CLuMP
```bash
./RCLuMP [MAXFRUS] [R] [log2_iters] [RUN] [graph.conf]
```
*(Replace `beta` with population size `R`.)*

#### Simulated Annealing
```bash
./SimulatedAnnealing [betaf] [log2_iters] [RUN] [graph.conf]
```

#### Population Annealing
```bash
./PopulationAnnealing [betaf] [R] [log2_iters] [RUN] [graph.conf]
```

#### Extremal Optimization
```bash
./ExtremalOpt [log2_iters] [RUN] [graph.conf]
```

#### τ - Extremal Optimization

The code used to generate the τ-Extremal Optimization (τEO) results was not developed as part of this project. These results were kindly provided by Stefan Boettcher.
For a detailed description of the τEO method, please refer to the following publication:

Extremal optimization: heuristics via coevolutionary avalanches, S. Boettcher, Computing in Science & Engineering, 2(6), 75–82 (2002), DOI: https://ieeexplore.ieee.org/document/881710


---

### Step 3 — Analyse results

```bash
python3 ResultsAnalysis.py [OUTPUT] [tolerance] [data_file(s)]
```

| Parameter | Options / Notes |
|-----------|-----------------|
| `OUTPUT` | `INFO` (stats only) · `OVERVIEW` (stats + energy evolution) · `FULL` (all + cluster move analysis) |
| `tolerance` | `[0., 1.]` — energy tolerance for convergence. Default `0` (exact match) |
| `data_file(s)` | One or more result files — **all must be in the same folder** |

> **Note for CLuMP/R-CLuMP:** both the `.bp` and `EVO.bp` files must be present in the same folder.

---

## Reproducing Paper Figures

The exact algorithm parameters used to produce each figure are recorded in:
```
data/results/[graph]/[graph]_ALGOspecs.dat
```
*(Download the data archive from Zenodo first.)*

**Quick-start workflow:**
```bash
# 1. Generate graph
./GraphGenerator RRG Gaussian 1000 3

# 2. Run CLuMP
./CLuMP <MAXFRUS> <beta> <log2_iters> <RUN> graph.conf

# 3. Analyse
python3 ResultsAnalysis.py FULL 0 results/[graph]/*/*.bp results/[graph]/*/*.bpr results/[graph]/*/*.eo results/[graph]/*/*.taueo results/[graph]/*/*.sa
```

---

## License

**Code:** MIT License — see [`LICENSE`](LICENSE).  
**Data:** CC BY 4.0 — see the [Zenodo record](https://doi.org/10.5281/zenodo.19048238).

---

## Citation

If you use this code, please cite:

```bibtex
@article{Rissone2025CLuMP,
  author  = {Rissone, Paolo and Boettcher, Stefan and Amendola, Alfonso and Sala, Simone and Ricci-Tersenghi, Federico},
  title   = {Cluster-based Message-Passing ({CLuMP}) Optimization in Rugged Energy Landscapes},
  journal = {XXX},
  year    = {2025},
  doi     = {10.XXXX/XXXXXX}
}
```
*(TBD)*

---

## Contact

Dr. Paolo Rissone — [rissone.paolo@gmail.com](mailto:rissone.paolo@gmail.com)
