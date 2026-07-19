"""
    SharkDentalDynamics

Clean, config-driven reimplementation of the shark tooth-drop metapopulation
simulation and empirical-comparison pipeline behind Kim, Yeakel et al. (2022,
Proc. R. Soc. B), "Decoding the dynamics of dental distributions." Reproduces the
behavior of the original `sharksims_*.jl`/`sharkcompare_*.jl` scripts and `src/*.jl`
in this repository, with three differences:

  1. Scenario parameters (temperature, distance, grid resolution, body size) are
     explicit `SiteConfig` fields with named presets (`modern`, `eocene_highlatitude`,
     `eocene_lowlatitude`) instead of copy-pasted across near-duplicate scripts.
  2. Parallelism is `Threads.@threads` instead of `@distributed`/`SharedArray`.
  3. The reps-averaging step in `error_surface` fixes a bug found in the legacy
     `sharkcompare_modern.jl` (4 of 7 shape-descriptor error terms were averaged
     across the sigtau grid axis instead of across reps; `sharkcompare_eocene.jl`/
     `sharkcompare_eocene_lowlatitude.jl` were unaffected -- see `error_surface`'s
     docstring for details, and note `sharkcompare_modern.jl` has since been fixed
     too). Everything else is a faithful, validated port; see
     `check_legacy_equivalence.jl`.

Public API: `SiteConfig`, `modern`, `eocene_highlatitude`, `eocene_lowlatitude`,
`build_library!`, `compare_site`, `simulate_metapopulation`, `precompute_growth`.
"""
module SharkDentalDynamics

using Random
using LinearAlgebra: dot
using Statistics: mean
using JLD2
using KernelDensity: kde
import Plots

include("config.jl")
include("growth.jl")
include("simulate.jl")
include("shape_descriptors.jl")
include("compare.jl")
include("library.jl")
include("plotting.jl")

export SiteConfig, modern, eocene_highlatitude, eocene_lowlatitude
export precompute_growth, simulate_metapopulation, tooth_length, SimulationTimeout
export shape_descriptors, detect_modes, local_maxima, local_minima
export compare_to_empirical, error_surface, best_fit, compare_site
export build_library!, save_settings, load_settings, save_simdata, load_simdata,
       settings_path, simdata_path, legacy_simdata_path
export plot_scenario, best_fit_densities, sigtau_mass

end # module
