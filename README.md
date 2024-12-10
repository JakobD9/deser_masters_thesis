# Code Masters' Thesis Jakob Deser

## 3b_chaotic.jl
Implementation of AMOC 3-box model with optional hosing, T-forcing, and Lorenz forcing. Compares two different time series.

- Choose parameterset (see *parameters.jl*) to choose what model scenario to use.
- Setting Coupled_l enables correlated 2D noise, by standard 10xHadGEM-MM. Otherwise, only chaos in S_N.

## 3b_chaotic_rate.jl
AMOC 3-box model with hosing, T-forcing, and Lorenz forcing. Plots different time series with evenly spaced forcing rates.

## 5b_rescaled.jl
Implementation of AMOC 5-box model with optional hosing, T-forcing. Compares two time series. By standard adjusted to recreate one T-forcing hysteresis cycle.

## basin_of_attraction.jl
Plots equilibria and respective basins of attraction for a fixed hosing/T-forcing value on a specified grid in the transformed S_N, S_T plane for 3- and 5-box model.

## lorenz_precompute.jl
Precomputes N trajectories of Lorenz 63 system and saves data in CSV files. **Attention:** Dependent on setting this requires a lot of storage. Mainly used for *kramers_escape.jl*. Predefine saveat (fixed time steps) and final time dependent on needs.

## sigma_approximation.jl
Approximates σ² using the Green-Kubo formula, either solving one double integral or averaging over single integrals of drawn initial conditions (`sampled = true`).

## sigma_estimation.jl
Estimates σ² using the weak invariance principle.

## amoc_model_dissipativ.jl
Creates a plot of the maximum value of AMOC 3-box model on circles of radius R around the origin.

## kramers_escape.jl
Integrates `num_ic` trajectories for a long time and plots the logarithmic average time against the inverse strength of perturbation.  
Supports double-well potential (`double_well = true`) and AMOC 3-box model (else). Parameter `noise_limit` determines whether to use stochastic forcing or chaotic forcing. For chaotic forcing, it is possible to use precomputed Lorenz 63 trajectories (see *lorenz_precompute.jl*).  

**Notice:** The results are only valid if all trajectories tipped at least once; otherwise, the average might be too low.

## parameters.jl
Reads the model parameters for different fits to GCM from `parameters.txt`. Contains location of stable equilibria and defines structure for forcing the AMOC model and noise coefficients.

## utils.jl
Defines all model functions and other useful functions.
