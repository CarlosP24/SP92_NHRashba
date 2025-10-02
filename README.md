# Non-Hermitian Skin Effect and Electronic Nonlocal Transport

[![arXiv](https://img.shields.io/badge/arXiv-2510.00921-b31b1b.svg)](https://arxiv.org/abs/2510.00921)
[![DOI](https://zenodo.org/badge/966776654.svg)](https://doi.org/10.5281/zenodo.17243672)
[![Julia v1.11+](https://img.shields.io/badge/Julia-v1.11+-blue.svg)](https://julialang.org/)
[![Quantica badge](https://raw.githubusercontent.com/pablosanjose/Quantica.jl/master/docs/src/assets/badge.svg)](https://github.com/pablosanjose/Quantica.jl)
[![License](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](LICENSE.md)

This repository contains the code used to perform numerical calculations and generate all figures for the manuscript "Non-Hermitian Skin Effect and Electronic Nonlocal Transport".

## Abstract

Open quantum systems governed by non-Hermitian effective Hamiltonians exhibit unique phenomena, such as the non-Hermitian skin effect, where eigenstates localize at system boundaries. We investigate this effect in a Rashba nanowire coupled to a ferromagnetic lead and demonstrate that it can be detected via non-local transport spectroscopy: while local conductance remains symmetric, the non-local conductance becomes non-reciprocal. We account for this behavior using both conventional transport arguments and the framework of non-Hermitian physics. Furthermore, we explain that exceptional points shift in parameter space when transitioning from periodic to open boundary conditions, a phenomenon observed in other non-Hermitian systems but so far not explained. Our results establish transport spectroscopy as a tool to probe non-Hermitian effects in open electronic systems.

## About

### Prerequisites

- Julia v1.11 or higher
- Recommended: 8+ CPU cores for parallel calculations

### Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/CarlosP24/SP92_NHRashba.git
   cd SP92_NHRashba
   ```

2. Install dependencies:

   ```bash
   julia --project=. -e "using Pkg; Pkg.instantiate()"
   ```

### Running Calculations

#### Local Execution

For quick local calculations with parallel processing:

```bash
julia bin/launch_local.jl base_wire_nh
```

#### HPC Cluster Execution

For large-scale calculations on computing clusters, use a script like `launch_esbirro.sh`:

```bash
# Submit to SLURM queue
bash bin/launch_esbirro.sh base_wire_nh
```

#### Available Systems

- `base_wire`: Basic Rashba nanowire
- `base_wire_nh`: Non-Hermitian Rashba nanowire
- `base_wire_long`: Extended length version
- `base_wire_nh_long`: Extended non-Hermitian version

#### Calculation Types

- Add `_cond` suffix for conductance calculations only
- Add `_spec` suffix for spectrum calculations only
- Use base name for both spectrum and conductance

### Generating Figures

1. Activate the plots environment:

   ```bash
   cd plots
   julia --project=. plots.jl
   ```

2. Generate individual figures:

   ```julia
   include("fig_bands.jl")      # Band structure and wavefunctions
   include("fig_conductance.jl") # Transport properties
   include("fig_spectrum.jl")   # Energy spectra
   ```

All figures are saved as PDF files in `plots/figures/`.

## Repository Structure

```text
├── src/                    # Main source code
│   ├── main.jl            # Main execution script
│   ├── builders/          # System construction modules
│   │   ├── Wire.jl        # Rashba nanowire builder
│   │   └── leads.jl       # Lead attachment utilities
│   ├── calculations/      # Numerical calculation modules
│   │   ├── Spectrum.jl    # Energy spectrum calculations
│   │   └── Conductance.jl # Transport calculations
│   ├── models/            # Physical system definitions
│   │   ├── systems.jl     # System parameter structures
│   │   └── wires.jl       # Predefined wire configurations
│   └── parallelizers/     # Parallel computing utilities
├── plots/                 # Figure generation
│   ├── plotters.jl        # Plotting utilities
│   ├── fig_*.jl          # Individual figure scripts
│   └── figures/          # Generated PDF figures
├── bin/                   # Execution scripts
│   ├── launch_local.jl    # Local parallel execution
│   └── launch_*.sh       # HPC cluster scripts
├── data/                  # Computational results (JLD2 format)
└── config/               # Configuration files
```

## Dependencies

### Core Packages

- **[Quantica.jl](https://github.com/pablosanjose/Quantica.jl)** (v1.2.0): Quantum tight-binding calculations
- **[JLD2.jl](https://github.com/JuliaIO/JLD2.jl)** (v0.5.13): Data serialization
- **[Arpack.jl](https://github.com/JuliaLinearAlgebra/Arpack.jl)** (v0.5.4): Eigenvalue computations
- **[CairoMakie.jl](https://github.com/JuliaPlots/Makie.jl)** (v0.12.18): Scientific plotting

### Computational Environment

- **[SlurmClusterManager.jl](https://github.com/JuliaParallel/SlurmClusterManager.jl)** (v1.1.0): HPC integration
- **[ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl)** (v1.10.4): Progress tracking
- **[Parameters.jl](https://github.com/mauro3/Parameters.jl)** (v0.12.3): Parameter management

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{Paya:25,
  title={Non-Hermitian Skin Effect and Electronic Nonlocal Transport}, 
  author={Carlos Payá and Oliver Solow and Elsa Prada and Ramón Aguado and Karsten Flensberg},
  year={2025},
  eprint={2510.00921},
  archivePrefix={arXiv},
  primaryClass={cond-mat.mes-hall},
  url={https://arxiv.org/abs/2510.00921}
}
```

You can also cite the Zenodo repository for the code:

[![DOI](https://zenodo.org/badge/966776654.svg)](https://doi.org/10.5281/zenodo.17243672)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details.

## Support

For questions about the code or manuscript, please open a new issue with a detailed description or contact the first author(s) listed in the manuscript.

## Acknowledgments

This work was supported by Grants No. PID2021-125343NB-I00, PRE2022-101362 and CEX2024-001445-S, funded by MICIU/AEI/10.13039/501100011033, ERDF A way of making Europe and ESF+; the European Research Council (Grant Agreement No. 856526) and by the DFG Collaborative Research Center 183, Project No. 277101999. 
