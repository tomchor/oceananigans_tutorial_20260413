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
| `channel_flow_les/` | Doubly-periodic LES channel flow driven by a pressure-gradient body force and damped by quadratic bottom drag, with a passive tracer |
| `seamount_flow/` | 2-D nonhydrostatic flow past a Gaussian seamount using the immersed boundary method |

## Prerequisites

- [Julia](https://julialang.org/downloads/) (1.10 or later recommended)
- [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl)
- [NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl)
- [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) (for plotting)

Install dependencies from the Julia REPL:

```julia
using Pkg
Pkg.add(["Oceananigans", "NCDatasets", "CairoMakie"])
```

## Running the examples

Each script can be run from the Julia REPL (recommended for interactive exploration) or from the terminal:

```bash
julia channel_flow_les/channel_flow_les.jl
julia seamount_flow/seamount_flow.jl
```

Because the scripts avoid module-level constants, you can freely modify parameters and re-`include` them in the same REPL session without restarting Julia.
