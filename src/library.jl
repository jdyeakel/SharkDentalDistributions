"""
    settings_path(data_dir)
    simdata_path(data_dir, rep, sigtau_pos, tau_pos)

Plain `joinpath`-based file naming, replacing `src/smartpath.jl`'s `homedir()`-hardcoded
string logic with an explicit `data_dir`. The per-cell naming scheme
(`simdata_<rep>_<sigtau_pos>_<tau_pos>.jld2`) matches the legacy scripts, so resumability
(skip a cell if its file already exists) works the same way.
"""
settings_path(data_dir::AbstractString) = joinpath(data_dir, "simsettings.jld2")
simdata_path(data_dir::AbstractString, rep::Integer, sigtau_pos::Integer, tau_pos::Integer) =
    joinpath(data_dir, "simdata_$(rep)_$(sigtau_pos)_$(tau_pos).jld2")

"""
    legacy_simdata_path(data_dir, rep, sigtau_pos, tau_pos)

Path function for reading a simulation library saved by the *legacy* `sharksims_*.jl`
scripts (`.jld` extension, no version suffix) -- e.g. the already-built
`data/sharks_modern2/`, `data/sharks_eocene2/`, `data/sharks_eocene_lowlatitude/`
directories in the parent repo. The internal JLD2 keys (`mass1`, `mass2`, `toothdrop`,
`toothlength1`, `toothlength2`) are identical to the new format -- only the filename
differs -- so `load_simdata(dir, r, s, t; path_fn = legacy_simdata_path)` reads a
legacy-built library directly, with no need to rebuild it.
"""
legacy_simdata_path(data_dir::AbstractString, rep::Integer, sigtau_pos::Integer, tau_pos::Integer) =
    joinpath(data_dir, "simdata_$(rep)_$(sigtau_pos)_$(tau_pos).jld")

function save_settings(config::SiteConfig, data_dir::AbstractString)
    mkpath(data_dir)
    JLD2.jldsave(settings_path(data_dir); config)
    return nothing
end

load_settings(data_dir::AbstractString) = JLD2.load(settings_path(data_dir), "config")

function save_simdata(data_dir::AbstractString, rep::Integer, sigtau_pos::Integer, tau_pos::Integer, result)
    path = simdata_path(data_dir, rep, sigtau_pos, tau_pos)
    JLD2.jldsave(path;
        mass1 = result.mass1, mass2 = result.mass2, toothdrop = result.toothdrop,
        toothlength1 = tooth_length(result.mass1), toothlength2 = tooth_length(result.mass2),
    )
    return path
end

function load_simdata(data_dir::AbstractString, rep::Integer, sigtau_pos::Integer, tau_pos::Integer;
        path_fn::Function = simdata_path)
    path = path_fn(data_dir, rep, sigtau_pos, tau_pos)
    d = JLD2.load(path)
    return (;
        mass1 = d["mass1"], mass2 = d["mass2"], toothdrop = d["toothdrop"],
        toothlength1 = d["toothlength1"], toothlength2 = d["toothlength2"],
    )
end

"""
    build_library!(config::SiteConfig, data_dir::AbstractString; seed=1, overwrite=false,
                   max_seconds_per_run=120.0)

Builds the full `(sigtau, tau, rep)` simulation library for `config`, saving one JLD2
file per grid cell under `data_dir` (same naming/resumability scheme as the legacy
`sharksims_*.jl` scripts -- skips a cell whose file already exists unless
`overwrite=true`) plus a `simsettings.jld2`. Replaces `@distributed`/`SharedArray`
with `Threads.@threads` over the flattened `(rep, sigtau_pos, tau_pos)` grid; growth
tables are computed once via `precompute_growth` (they don't depend on sigtau/tau) and
shared read-only across threads/cells, rather than recomputed for every cell.

Each grid cell gets its own `Xoshiro` RNG, seeded deterministically from
`(seed, rep, sigtau_pos, tau_pos)`, so a rerun with the same `seed` reproduces the
library exactly regardless of thread scheduling order -- a new property; the legacy
sweep was never seeded (relied on whatever the OS-seeded global RNG gave each worker
process at the time).

`max_seconds_per_run` (default 120s -- normal cells typically take 15-35s even at
real `gen`, so this leaves a wide margin) guards against the pathologically expensive
`(sigtau, tau)` corners described in `SimulationTimeout`'s docstring: a cell that
exceeds it is skipped with a `@warn` rather than blocking the whole sweep. Skipped
cells stay missing from `data_dir` and will be retried (and may time out again) on
the next `build_library!` call for the same `data_dir` -- rerun with a larger
`max_seconds_per_run`, or accept the gap (`error_surface`/`compare_site` will simply
have fewer reps at that grid cell) if you've confirmed it's a genuine pathological
corner rather than just a slow machine. Pass `max_seconds_per_run = nothing` to
disable the guard entirely.
"""
function build_library!(config::SiteConfig, data_dir::AbstractString; seed::Integer = 1, overwrite::Bool = false,
        max_seconds_per_run::Union{Nothing, Real} = 120.0)
    mkpath(data_dir)
    save_settings(config, data_dir)

    precomp = precompute_growth(config)
    nsig = length(config.sigtauvec)
    ntau = length(config.tauvec)
    reps = config.reps
    its = reps * nsig * ntau
    combos = [(r, s, t) for r in 1:reps for s in 1:nsig for t in 1:ntau]

    Threads.@threads for idx in 1:its
        r, s, t = combos[idx]
        path = simdata_path(data_dir, r, s, t)
        if overwrite || !isfile(path)
            rng = Random.Xoshiro(hash((seed, r, s, t)))
            sigtau = config.sigtauvec[s]
            tau = config.tauvec[t]
            try
                result = simulate_metapopulation(config, precomp, sigtau, tau; rng, max_seconds = max_seconds_per_run)
                save_simdata(data_dir, r, s, t, result)
            catch e
                e isa SimulationTimeout || rethrow()
                @warn "skipping grid cell after timeout" rep=r sigtau_pos=s tau_pos=t sigtau tau elapsed=e.elapsed n_events=e.n_events
            end
        end
    end
    return data_dir
end
