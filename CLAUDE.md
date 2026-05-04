# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a MATLAB-based Computational Fluid Dynamics (CFD) project that solves the **Shallow Water Equations (SWE)** with topography using the **MUSCL-Hancock finite volume method** and **HLLC approximate Riemann solver**. It was authored by Narong Batsuwan.

## Running the Code

Each case script is self-contained and run directly in MATLAB:

```matlab
% Run a specific case
run('case_1_1DHump_and_Exact.m')
run('case_2_circle_dam_break.m')
run('case_3_water_flowing_through_the_river.m')
run('case_4_dam_dry.m')
```

Each script uses `clear` and `close all` at the top, so they can be run independently. There is no build step, no test runner, and no package manager.

## Numerical Algorithm Architecture

All case scripts implement the same two-stage MUSCL-Hancock pipeline. Understanding this pipeline is essential before modifying any case:

### State Variables (per cell)
- `h` — water depth
- `eta` — water surface elevation (`eta = h + z`)
- `hu`, `hv` — discharge (momentum) in x and y directions
- `z` — bed topography elevation
- `zbx`, `zby` — single topography value at each cell interface (well-balancing)

### Grid Layout
Cells are indexed with **2 ghost layers on each side**. For `n` interior cells, arrays are sized `n+4`. Interior cells occupy indices `3:n+2`; flux loops run over `3:n+3` (interfaces).

### Two-Stage Update (per time step)
1. **Slope limiting** via `minmod()` for `h`, `eta`, `hu`/`hv` on both left and right sides of each interface.
2. **Data reconstruction** — piecewise-linear extrapolation to cell faces.
3. **Well-balancing fix** — single interface topography `zbx(i) = max(zl, zr)`, then recompute `hl = max(0, etal - zbx(i))` to prevent spurious flow at rest.
4. **HLLC flux computation** — wave speeds `Sl`, `Sr`, `Sm`; flux selected from left/middle/right region.
5. **Solution update** — flux divergence + bed slope source term, using factor `0.5*dt/dx` (each stage advances a half step).
6. Repeat stages 1–5 a second time using the updated values (completing the full Hancock two-stage predictor-corrector).

### Key Constants and Thresholds
- `grav = 9.806` m/s² (global variable)
- Dry-bed threshold: `h < 1.0e-06` — velocity is set to zero to avoid division by near-zero depth.
- Negative depth correction: `h` floored at `1.0e-06` after each update to maintain positivity.
- CFL factor: `dt = 0.5 * dx / max(|u| + sqrt(g*h))`

### Source Term (bed slope)
The gravitational source on momentum is discretized with central-differencing at the interface:
```
sox = ((hlx(i+1) + hrx(i)) / 2) * ((zbx(i+1) - zbx(i)) / dx)
```

## Shared Utility

`minmod.m` — slope limiter used in all cases:
- Returns `a` if `|a| ≤ |b|` and same sign
- Returns `b` if `|b| < |a|` and same sign
- Returns `0` if opposite signs (limiter activates)

## Cases Summary

| File | Dimensions | Domain | Cells | Steps | Scenario |
|------|-----------|--------|-------|-------|----------|
| `case_1_1DHump_and_Exact.m` | 1D | 1000 m | 41 | 2522 | Flow over cosine hump, compared to exact solution |
| `case_2_circle_dam_break.m` | 2D | 200×200 m | 85×85 | 200 | Circular dam break, radial propagation |
| `case_3_water_flowing_through_the_river.m` | 2D | 6×6 m | 60×60 | 600 | River flow with natural topography |
| `case_4_dam_dry.m` | 2D | 100×100 m | 200×200 | 350 | Dam break onto dry ground through a gate |

## 2D Extension

Cases 2–4 extend the 1D algorithm to 2D by **operator splitting** (sequential x-sweep then y-sweep within each time step). The x-direction sweep uses `zbx`/`hlx`/`hrx` arrays; the y-direction uses `zby`/`hly`/`hry`. Both sweeps share the same time step `dt` computed from the 2D CFL condition:
```matlab
dt = 0.5 / (max(max(abs(u)+sqrt(grav*h)))/dx + max(max(abs(v)+sqrt(grav*h)))/dy)
```
