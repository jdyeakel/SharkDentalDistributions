#!/usr/bin/env julia
# Seeded, numerical old-vs-new equivalence check.
#
# Run with (from the package root):  julia --project=. scripts/check_legacy_equivalence.jl
#
# This does NOT modify or depend on any file outside this folder being changed --
# it only *reads* the legacy `../src/*.jl` files, `include`-ing them into an isolated
# `Legacy` module so their function names don't clash with SharkDentalDynamics.
#
# Two things are checked:
#   1. The stochastic kernel: legacy `popgen_migrate_g` vs new
#      `simulate_metapopulation`, for all 3 presets at several representative
#      (sigtau, tau) grid points (including the grid corners), using the *literal
#      same* default RNG object (re-seeded identically immediately before each call)
#      -- this sidesteps any risk of two independently-constructed RNG types not
#      being bit-identical for the same integer seed.
#
#      This is run with `gen` overridden to a tiny fraction (see `CHECK_GEN`) rather
#      than each preset's real value. Equivalence between the two implementations is
#      a *code-path* property, not a statistical one: if the ported logic diverges
#      from the original in even a single random draw, every subsequent draw
#      diverges too (the dynamics are chaotic), so a handful of Gillespie events is
#      just as conclusive as a full-length run -- and vastly cheaper. This matters
#      because some grid corners (e.g. large sigtau + large tau in
#      eocene_highlatitude, where dispersal is maximally flexible in both directions)
#      let the metapopulation grow substantially within one full `gen` cycle, and
#      since Gillespie event size dt ~ 1/N, a growing population needs far more
#      discrete events to simulate the same time span -- one such corner took nearly
#      an hour and 10GB of RAM at the real `gen` before being killed. A separate,
#      clearly-labeled realistic-timing smoke test below runs a couple of *moderate*
#      (non-corner) grid points at each preset's real `gen`, to sanity-check actual
#      production timing without hitting that pathology.
#   2. The deterministic shape-descriptor/comparison pipeline: legacy
#      `toothdist_analysis`/`toothdist_emp_analysis`/`empirical_sim_comparison` vs new
#      `shape_descriptors`/`compare_to_empirical`, on fixed synthetic inputs (no RNG
#      involved at all).
#
# NOT checked here: the legacy `sharkcompare_*.jl` driver scripts' reps-averaging
# step, which has a confirmed bug (see `error_surface`'s docstring in
# src/compare.jl) that this package deliberately does not reproduce.

using Random

module Legacy
    using Random, LinearAlgebra, Statistics, KernelDensity
    # scripts/ -> SharkDentalDynamics/ -> 2018_sharks/ (where the legacy src/ lives)
    const _ROOT = normpath(joinpath(@__DIR__, "..", ".."))

    # `src/ts.jl`'s `ts(epsilon1,epsilon2,a,eta)` references `M` (asymptotic mass)
    # without taking it as a parameter -- it only worked in the original driver
    # scripts because `M` was already a top-level global by the time `ts`/
    # `popgen_migrate_g` were called. Replicate that here with an explicit setter.
    M = 0.0
    set_M!(v) = (global M = v)

    include(joinpath(_ROOT, "src", "ts.jl"))
    include(joinpath(_ROOT, "src", "popgen_migrate_g.jl"))
    include(joinpath(_ROOT, "src", "findlocalmaxima.jl"))
    include(joinpath(_ROOT, "src", "findlocalminima.jl"))
    include(joinpath(_ROOT, "src", "modality_analysis.jl"))
    include(joinpath(_ROOT, "src", "toothdist_analysis.jl"))
    include(joinpath(_ROOT, "src", "toothdist_emp_analysis.jl"))
    include(joinpath(_ROOT, "src", "empirical_sim_comparison.jl"))
end

using SharkDentalDynamics

const RESULTS = Tuple{String, Bool}[]

function check(name::AbstractString, cond::Bool; elapsed = nothing)
    push!(RESULTS, (name, cond))
    suffix = elapsed === nothing ? "" : "  ($(round(elapsed; digits = 2))s)"
    println(cond ? "  PASS  " : "  FAIL  ", name, suffix)
    flush(stdout)
    return cond
end

# ---------------------------------------------------------------------------
# 1. Stochastic kernel equivalence (cheap: tiny `gen` override)
# ---------------------------------------------------------------------------

const CHECK_GEN = 1e-4  # ~0.01% of a real growth cycle -> a handful of events per run

function kernel_matches(config::SiteConfig, precomp, sigtau, tau, seed::Integer;
        gen_override::Union{Nothing, Real} = nothing)
    distance_m = config.distance_km * 1000.0
    Legacy.set_M!(precomp.M)

    # `SiteConfig.gen::Int` (matching the original's `gen=1`) can't take a tiny
    # fractional override, so truncate `tmax` directly instead -- it's just
    # `mean(...) * gen` internally, so scaling it by `gen_override/config.gen` is
    # exactly equivalent to having built the config with that `gen` in the first
    # place. The legacy call still gets an explicit `gen` argument (a plain
    # function parameter with no type constraint) so both sides see the same tmax.
    legacy_gen = gen_override === nothing ? config.gen : gen_override
    new_precomp = gen_override === nothing ? precomp :
        merge(precomp, (; tmax = precomp.tmax * (gen_override / config.gen)))

    Random.seed!(seed)
    legacy_mass1, legacy_mass2, legacy_epsilonvec, legacy_clock, legacy_popstate,
        legacy_toothdrop, legacy_state = Legacy.popgen_migrate_g(
        precomp.m0, precomp.M, precomp.tempvec1, precomp.tempvec2,
        config.n0, legacy_gen, distance_m, config.velocity, config.D, sigtau, tau,
    )

    Random.seed!(seed)
    new_result = simulate_metapopulation(config, new_precomp, sigtau, tau; rng = Random.default_rng())

    return legacy_mass1 == new_result.mass1 &&
           legacy_mass2 == new_result.mass2 &&
           legacy_epsilonvec == new_result.epsilonvec &&
           legacy_clock == new_result.clock &&
           legacy_popstate == new_result.popstate &&
           legacy_toothdrop == new_result.toothdrop &&
           legacy_state == new_result.state
end

println("Checking stochastic kernel equivalence (gen=$CHECK_GEN override, all grid corners)...")
flush(stdout)
for (label, config) in ("modern" => modern(), "eocene_highlatitude" => eocene_highlatitude(),
                         "eocene_lowlatitude" => eocene_lowlatitude())
    precomp = precompute_growth(config)
    nsig, ntau = length(config.sigtauvec), length(config.tauvec)
    gridpoints = [(1, 1), (1, ntau), (nsig, 1), (nsig, ntau), (cld(nsig, 2), cld(ntau, 2))]
    for (si, ti) in gridpoints
        sigtau, tau = config.sigtauvec[si], config.tauvec[ti]
        for seed in (1, 2024)
            t0 = time()
            ok = kernel_matches(config, precomp, sigtau, tau, seed; gen_override = CHECK_GEN)
            check("$label  sigtau_pos=$si tau_pos=$ti  seed=$seed", ok; elapsed = time() - t0)
        end
    end
end

# ---------------------------------------------------------------------------
# 2. Realistic-timing smoke test (real `gen`, moderate/interior grid points only)
# ---------------------------------------------------------------------------

println("\nRealistic-timing smoke test (real gen, moderate grid points, legacy vs new)...")
flush(stdout)
for (label, config) in ("modern" => modern(), "eocene_highlatitude" => eocene_highlatitude(),
                         "eocene_lowlatitude" => eocene_lowlatitude())
    precomp = precompute_growth(config)
    nsig, ntau = length(config.sigtauvec), length(config.tauvec)
    # an interior point away from the sigtau-large/tau-large corner
    si, ti = max(1, cld(nsig, 3)), max(1, cld(ntau, 3))
    sigtau, tau = config.sigtauvec[si], config.tauvec[ti]
    t0 = time()
    ok = kernel_matches(config, precomp, sigtau, tau, 1)
    check("$label  sigtau_pos=$si tau_pos=$ti  seed=1  (real gen=$(config.gen))", ok; elapsed = time() - t0)
end

# ---------------------------------------------------------------------------
# 3. Deterministic shape-descriptor / comparison pipeline equivalence
# ---------------------------------------------------------------------------

println("\nChecking deterministic analysis pipeline (compare_to_empirical vs legacy empirical_sim_comparison)...")
flush(stdout)

function deterministic_pipeline_matches()
    rng = Random.Xoshiro(2024)
    toothlength = collect(0.0:0.5:40.0)
    n = length(toothlength)
    juv_density = exp.(-((toothlength .- 12.0) .^ 2) ./ 20) .+ 0.05 .* rand(rng, n)
    adult_density = exp.(-((toothlength .- 25.0) .^ 2) ./ 40) .+ 0.05 .* rand(rng, n)
    toothdrop = hcat(juv_density, adult_density)
    measures = 8.0 .+ 15.0 .* rand(rng, 150)

    legacy_vals = Legacy.empirical_sim_comparison(toothdrop, toothlength, measures)
    new_res = compare_to_empirical(toothdrop, toothlength, measures)
    new_vals = (
        new_res.mean_j, new_res.mean_a, new_res.mode_j, new_res.mode_a,
        new_res.modedist_j, new_res.modedist_a, new_res.sd_j, new_res.sd_a,
        new_res.median_j, new_res.median_a, new_res.q25_j, new_res.q25_a,
        new_res.q75_j, new_res.q75_a,
    )
    return all(isapprox.(legacy_vals, new_vals; atol = 1e-9))
end

check("compare_to_empirical vs legacy empirical_sim_comparison (synthetic bimodal input)",
      deterministic_pipeline_matches())

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

n_pass = count(r -> r[2], RESULTS)
n_total = length(RESULTS)
println("\n$(n_pass)/$(n_total) checks passed.")
if n_pass == n_total
    println("ALL CHECKS PASSED.")
else
    println("FAILURES:")
    for (name, ok) in RESULTS
        ok || println("  - ", name)
    end
    exit(1)
end
