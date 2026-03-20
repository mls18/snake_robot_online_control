# Data-Driven Adaptive Control for Planar Snake Robots

This repository accompanies the journal paper:

> **"A Data-Driven Adaptive Control Framework for Planar Snake Robots under Unknown Ground Friction"**

The code implements an online adaptive pipeline that:
1. Designs a robust periodic controller from offline data 
2. Detects ground friction changes online via tangential velocity monitoring
3. Estimates updated friction coefficients
4. Collects new closed-loop data and redesigns the controller automatically

Two model variants are validated: a **simplified** planar snake robot model and a **complex** model based on Liljeback et al. (2013).

---

## Requirements

| Toolbox / Software | Purpose |
|---|---|
| MATLAB R2020b or later | Main simulation environment |
| [CVX](http://cvxr.com/cvx/) | Convex optimisation (SDP solver) |
| [MOSEK](https://www.mosek.com/) | Recommended backend solver for CVX |
| MATLAB Optimization Toolbox | Used internally by `ode45` |


---

---

## Quick Start

### 1. Simplified Model (recommended entry point)
```matlab
% In MATLAB, from the repository root:
run('main_online_simplified.m')
```

This script runs the full adaptive pipeline on the simplified snake model:
- **Steps 1–3 000:** Initial ground (cₙ = 3, cₜ = 1), closed-loop tracking
- **Step 3 000:** Ground change to (cₙ = 16, cₜ = 4)
- **Steps ~3 000–5 400:** Detection → friction estimation → data collection → SDP redesign
- **Steps 5 400–20 000:** Resumed tracking under new controller

Three figures are produced automatically (joint tracking errors, error norm vs. bound, joint trajectories vs. reference).

### 2. Complex Model
```matlab
run('main_adaptive_complex.m')
```

Identical pipeline but uses `createComplexSnakeModelFast` for simulation — a numerically efficient implementation of the Liljeback et al. (2013) model that avoids symbolic matrix inversion.

---

## Key Parameters

| Parameter | Variable | Default | Meaning |
|---|---|---|---|
| Number of links | `N` | 5 | Snake robot link count |
| Undulation frequency | `omega` | 120°/s | Joint reference frequency |
| Sample time | `Ts` | 0.005 s | Discrete-time step |
| SDP upper bound | `rho_sdp` | 8000 | Upper bound on P[k] eigenvalues |
| SDP lower bound | `eta_sdp` | 1 | Lower bound on P[k] eigenvalues |
| Experiments | `L` | n+m = 9 | Number of parallel data-collection experiments |
| Exploration std | `sigma_explore` | 6.0 | Standard deviation of random exploration signal |
| PD gain (position) | `k_i` | 20 | Tracking PD gain — position |
| PD gain (velocity) | `k_d` | 5 | Tracking PD gain — velocity |

---

---

## Citing This Work

If you use this code, please cite the accompanying paper (citation details to be added upon publication).

The models are based on:

> P. Liljeback, K. Y. Pettersen, Ø. Stavdahl, and J. T. Gravdahl,
> *Snake Robots: Modelling, Mechatronics, and Control.*
> Springer, 2013.

---

## Troubleshooting

**SDP infeasible during online redesign**
- Try increasing `rho_sdp` (e.g. 10 000–20 000)
- Change `sigma_explore` slightly to improve Z[k] conditioning
- Change `R_bar_new` 

**CVX not found / solver error**
- Run `cvx_setup` in MATLAB after installing CVX
- Ensure MOSEK is installed and licensed; run `cvx_solver mosek` before the main scripts
