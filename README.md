# Oceananigans LES Tutorial

Tutorial materials for the **Large-Eddy Simulation (LES) with Oceananigans** session, part of the seminar:

> **Novas abordagens computacionais para mecânica dos fluidos: Julia e Chapel**
> Monday, April 13 2026 — ICMC/USP

## About the seminar

This free seminar explores how modern programming tools are transforming complex fluid simulations, bridging the gap between ease of use and the extreme performance required in the field.

## This repository

This repo contains the tutorial scripts used in the afternoon hands-on session. Each subdirectory is a self-contained example that can be run interactively in a Julia REPL or executed as a script.

| Directory | Description |
|---|---|
| `dry_atmosphere_les/` | Doubly-periodic 3D atmosphere LES driven by a pressure-gradient body force and quadratic surface drag, with a passive tracer |
| `hill_flow/` | 2D nonhydrostatic flow past a Gaussian hill using the immersed boundary method |
| `kelvin_helmholtz/` | 2D Kelvin-Helmholtz instability in a stratified shear layer; implicit LES via WENO dissipation |
| `free_convection/` | 3D atmospheric free convection heated from below with dynamic Smagorinsky SGS closure |
| `rain_over_ocean/` | 2D precipitating shallow cumulus convection (RICO case) using Breeze's anelastic model with one-moment cloud microphysics |

Each directory contains a simulation script (e.g. `hill_flow.jl`) and a matching plot script (e.g. `plot_hill_flow.jl`). Run the simulation first to produce the output file, then run the plot script to generate the animation.

## Prerequisites

- [Julia](https://julialang.org/downloads/) (1.12 or later recommended)
- All dependencies are listed in `Project.toml`. Install them in one step:

```julia
using Pkg
Pkg.instantiate()
```

This will install [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl), [Breeze.jl](https://github.com/NumericalEarth/Breeze.jl), [GLMakie.jl](https://github.com/MakieOrg/Makie.jl), [Oceanostics.jl](https://github.com/tomchor/Oceanostics.jl), and other dependencies automatically.

## Running the examples

Each script can be run from the Julia REPL (recommended for interactive exploration) or from the terminal:

```bash
julia kelvin_helmholtz/kelvin_helmholtz.jl   # run simulation
julia kelvin_helmholtz/plot_kelvin_helmholtz.jl  # produce animation
```

Because the scripts avoid module-level constants, you can freely modify parameters and re-`include` them in the same REPL session without restarting Julia.
