#!/usr/bin/env julia
# Worked example: define a brand-new scenario from scratch, build its simulation
# library FRESH (via build_library!, not any pre-existing/legacy data), and compare
# it against a new empirical site -- here, Blackheath, a high-latitude locality
# added to empirical/SandTiger_all_2026.csv (column 6).
#
# This is meant as a template: edit the "FREE PARAMETERS" block below for your own
# new site, then run (from the package root):
#   julia --project=. -t auto scripts/example_new_site.jl

using SharkDentalDynamics

const PKG_ROOT = normpath(joinpath(@__DIR__))

# ============================== FREE PARAMETERS ===============================
# Edit these for a new site. See README.md's "Presets reference" table and
# src/config.jl's docstrings for what each field means and the values behind the
# existing presets, for comparison/reference.

# Body size (cm): birth length, asymptotic length. `L` sets the maximum tooth
# length the model can produce (via tooth_length) -- make sure it's large enough
# to cover your empirical data's largest specimen.
l0 = 50.0
L = 477.0

# Initial population size; number of growth cycles to simulate.
n0 = 1000
gen = 1

# Seasonal temperature range at each site (°C). Blackheath is a high-latitude
# locality, so these mirror `eocene_highlatitude()` as a starting point -- replace
# with site-specific estimates if you have them.
mintemp_j, maxtemp_j = 12.0, 24.0   # nursery/juvenile site
mintemp_a, maxtemp_a = 9.0, 17.0    # adult site

# Distance between the two sites (km).
distance_km = 400.0

# Dispersal-window grid to sweep, and stochastic replicates per grid cell.
# NOTE: kept small and away from large-sigtau/large-tau (both near their max
# simultaneously is the known expensive corner -- see README's note on pathological
# grid corners, and SimulationTimeout's docstring) so this example finishes quickly
# with no cells skipped. `build_library!` has a 120s-per-cell timeout by default
# regardless (skips + warns rather than hanging), so widening this grid towards
# eocene_highlatitude()'s full 50x50x25 is safe to try -- just expect a warning for
# any cell that lands on that corner, and a correspondingly longer total runtime.
sigtauvec = collect(0.5:2.0:12.5)
tauvec = collect(1.0:5.0:26.0)
reps = 3

# Where to save the simulation library, and which empirical column to compare
# against.
data_dir = joinpath(PKG_ROOT, "data", "example_blackheath")
csv_path = joinpath(PKG_ROOT, "empirical", "SandTiger_all_2026.csv")
csv_column = 6   # Blackheath
site_label = "Blackheath"
figure_path = joinpath(PKG_ROOT, "figures", "fig_example_blackheath_julia.pdf")
# ================================================================================

config = SiteConfig(;
    l0, L, n0, gen,
    mintemp_j, maxtemp_j, mintemp_a, maxtemp_a,
    distance_km, sigtauvec, tauvec, reps,
)

"""
    read_csv_column(path, col)

Minimal CSV column reader (no CSV.jl/DataFrames.jl dependency needed) -- matches
`scripts/reproduce_figures.jl`.
"""
function read_csv_column(path::AbstractString, col::Int)
    lines = readlines(path)
    values = Float64[]
    for line in lines[2:end]
        cell = split(line, ',')[col]
        isempty(strip(cell)) || push!(values, parse(Float64, cell))
    end
    return values
end

measures = read_csv_column(csv_path, csv_column)
println("Loaded $(length(measures)) measurements for \"$site_label\"."); flush(stdout)

its = length(sigtauvec) * length(tauvec) * reps
println("Building simulation library at $data_dir ($(length(sigtauvec))x$(length(tauvec))x$reps = $its runs)..."); flush(stdout)
@time build_library!(config, data_dir; seed = 1)

println("Comparing against \"$site_label\"..."); flush(stdout)
result = compare_site(config, data_dir, measures)
println("  juvenile-site hypothesis: error = $(round(result.best_j.error, digits = 2)) at sigtau=$(result.best_j.sigtau), tau=$(result.best_j.tau)")
println("  adult-site hypothesis:    error = $(round(result.best_a.error, digits = 2)) at sigtau=$(result.best_a.sigtau), tau=$(result.best_a.tau)")
verdict = result.best_j.error < result.best_a.error ? "juvenile/nursery site" : "adult site"
println("  => better supported as a $verdict"); flush(stdout)

mkpath(dirname(figure_path))
plot_scenario(config, data_dir, [(site_label, measures)]; filename = figure_path)
println("\nSaved figure -> $figure_path")
