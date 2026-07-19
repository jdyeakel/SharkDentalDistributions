#!/usr/bin/env julia
# Reproduces the three legacy comparison figures (modern, eocene high-latitude,
# eocene low-latitude) natively in Julia via `Plots.jl`, reading the already-built
# legacy-format simulation libraries in the parent repo's `data/` directory --
# no rebuilding needed. Doesn't touch or depend on the R/RCall plotting code.
#
# Run with (from the package root):  julia --project=. -t auto scripts/reproduce_figures.jl

using SharkDentalDynamics

const PKG_ROOT = normpath(joinpath(@__DIR__))          # SharkDentalDynamics/
const LEGACY_ROOT = normpath(joinpath(@__DIR__, "..","..","2018_sharks"))
# normpath(joinpath(@__DIR__, "..", ".."))  # 2018_sharks/ (for the legacy data/ libraries)

"""
    read_csv_column(path, col)

Minimal CSV column reader (no CSV.jl/DataFrames.jl dependency needed): returns the
non-empty numeric entries of column `col` (1-indexed) of the comma-separated file at
`path`, skipping the header row.
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

# Local copy (see empirical/) -- this script never reaches outside the package.
csv_path = joinpath(PKG_ROOT, "empirical", "SandTiger_all.csv")

println("Reproducing modern (Delaware Bay) comparison figure...")
modern_sites = [("Delaware Bay", read_csv_column(csv_path, 5))]
plot_scenario(modern(), joinpath(LEGACY_ROOT, "data", "sharks_modern2"), modern_sites;
    path_fn = legacy_simdata_path, filename = joinpath(PKG_ROOT, "fig_modern_julia.pdf"))
println("  -> fig_modern_julia.pdf")

println("Reproducing eocene high-latitude (Banks Island / Seymour Island) comparison figure...")
highlat_sites = [("Banks Island", read_csv_column(csv_path, 3)),
                  ("Seymour Island", read_csv_column(csv_path, 4))]
plot_scenario(eocene_highlatitude(), joinpath(LEGACY_ROOT, "data", "sharks_eocene2"), highlat_sites;
    path_fn = legacy_simdata_path, filename = joinpath(PKG_ROOT, "fig_eocene_highlatitude_julia.pdf"))
println("  -> fig_eocene_highlatitude_julia.pdf")

println("Reproducing eocene low-latitude (Red Hot Truck Stop / Whiskey Bridge) comparison figure...")
lowlat_sites = [("Red Hot Truck Stop", read_csv_column(csv_path, 1)),
                 ("Whiskey Bridge", read_csv_column(csv_path, 2))]
plot_scenario(eocene_lowlatitude(), joinpath(LEGACY_ROOT, "data", "sharks_eocene_lowlatitude"), lowlat_sites;
    path_fn = legacy_simdata_path, filename = joinpath(PKG_ROOT, "fig_eocene_lowlatitude_julia.pdf"))
println("  -> fig_eocene_lowlatitude_julia.pdf")

println("\nDone.")
