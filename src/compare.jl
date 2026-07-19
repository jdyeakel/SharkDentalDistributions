"""
    relative_error(pred, obs)

Symmetric relative-error measure used throughout: `abs(sum((pred .- obs) ./ obs))`.
Works elementwise for both scalars and vectors (e.g. comparing sorted mode-location
pairs). Ported verbatim from the `errorfunc` closure in `empirical_sim_comparison.jl`.
"""
relative_error(pred, obs) = abs(sum((pred .- obs) ./ obs))

"""
    compare_to_empirical(toothdrop, toothlength, measures)

For one simulated `(toothdrop, toothlength)` pair (juvenile-site density in column 1,
adult-site density in column 2) and an empirical sample `measures`, computes the 7
shape-descriptor relative errors (mean, mode value, mode distance, SD, median, 25th
and 75th percentile) for both the juvenile- and adult-site hypotheses -- i.e. equation
2.1 of Kim, Yeakel et al. 2022, evaluated at one (sigtau, tau) grid point.

Ported from `src/empirical_sim_comparison.jl`, using the unified `shape_descriptors`
in place of the separate `toothdist_analysis`/`toothdist_emp_analysis`.
"""
function compare_to_empirical(toothdrop::AbstractMatrix, toothlength::AbstractVector, measures::AbstractVector)
    juv = shape_descriptors(toothdrop[:, 1], toothlength)
    adult = shape_descriptors(toothdrop[:, 2], toothlength)

    U = kde(measures)
    emp = shape_descriptors(U.density, collect(U.x))

    mean_j = relative_error(juv.mean, emp.mean)
    mean_a = relative_error(adult.mean, emp.mean)
    sd_j = relative_error(sqrt(juv.var), sqrt(emp.var))
    sd_a = relative_error(sqrt(adult.var), sqrt(emp.var))
    median_j = relative_error(juv.median, emp.median)
    median_a = relative_error(adult.median, emp.median)
    q25_j = relative_error(juv.quartile25, emp.quartile25)
    q25_a = relative_error(adult.quartile25, emp.quartile25)
    q75_j = relative_error(juv.quartile75, emp.quartile75)
    q75_a = relative_error(adult.quartile75, emp.quartile75)

    if emp.peak_bin == 0
        mode_j = relative_error(maximum(juv.modes), maximum(emp.modes))
        mode_a = relative_error(maximum(adult.modes), maximum(emp.modes))
        modedist_j = juv.peak_bin == 0 ? 0.0 : 1.0
        modedist_a = adult.peak_bin == 0 ? 0.0 : 1.0
    else
        mode_j = relative_error(juv.modes, emp.modes)
        mode_a = relative_error(adult.modes, emp.modes)
        modedist_j = relative_error(juv.peak_dist, emp.peak_dist)
        modedist_a = relative_error(adult.peak_dist, emp.peak_dist)
    end

    return (;
        mean_j, mean_a, mode_j, mode_a, modedist_j, modedist_a,
        sd_j, sd_a, median_j, median_a, q25_j, q25_a, q75_j, q75_a,
    )
end

const ERROR_TERM_FIELDS = (
    :mean_j, :mean_a, :mode_j, :mode_a, :modedist_j, :modedist_a,
    :sd_j, :sd_a, :median_j, :median_a, :q25_j, :q25_a, :q75_j, :q75_a,
)

"""
    error_surface(config, data_dir, measures)

Loads every saved `(rep, sigtau_pos, tau_pos)` simulation result under `data_dir` for
`config`, runs `compare_to_empirical` against `measures` at each grid cell, and
averages the 14 error terms **across reps** (a plain `dims=1`-style mean over the
`reps` axis for every term). Parallelized with `Threads.@threads` over the flattened
`(rep, sigtau_pos, tau_pos)` grid.

**Historical note (now fixed in the legacy script too):** `sharkcompare_modern.jl`
had a bug where `mean_*`/`mode_*`/`modedist_*` were correctly averaged across reps
(`dims=1`, since its arrays are shaped `(reps, nsig, ntau)`), but `sd_*`/`median_*`/
`q25_*`/`q75_*` were averaged with `dims=2` -- i.e. across the *sigtau* grid axis,
not reps -- then only replicate #1's resulting `(1, ntau)`-shaped slice was kept,
broadcasting identically across every sigtau row when summed into the final error
surface. That meant 4 of the 7 shape-descriptor error terms contributed no
sigtau-axis information to the modern/Delaware Bay `qmatrixj`/`qmatrixa` surfaces
(and therefore to its best-fit search and heatmap) -- only mean/mode/mode-distance
actually varied across sigtau there. This was a copy-paste `dims=1`/`dims=2` slip
specific to `sharkcompare_modern.jl`'s 3D array layout; `sharkcompare_eocene.jl` and
`sharkcompare_eocene_lowlatitude.jl` use a 4D `(num, reps, nsig, ntau)` layout where
`dims=2` is the *correct* reps axis for all 14 terms, so the four Eocene site
comparisons were never affected. `sharkcompare_modern.jl` has since been corrected
to `dims=1` throughout, matching this implementation; see
`check_legacy_equivalence.jl` for how this was scoped during validation.
"""
function error_surface(config::SiteConfig, data_dir::AbstractString, measures::AbstractVector;
        path_fn::Function = simdata_path)
    nsig = length(config.sigtauvec)
    ntau = length(config.tauvec)
    reps = config.reps
    its = reps * nsig * ntau

    raw = Dict(f => Array{Float64}(undef, reps, nsig, ntau) for f in ERROR_TERM_FIELDS)
    combos = [(r, s, t) for r in 1:reps for s in 1:nsig for t in 1:ntau]

    Threads.@threads for idx in 1:its
        r, s, t = combos[idx]
        result = load_simdata(data_dir, r, s, t; path_fn)
        toothlength = result.toothlength1[1, :]
        res = compare_to_empirical(result.toothdrop, toothlength, measures)
        for f in ERROR_TERM_FIELDS
            raw[f][r, s, t] = getfield(res, f)
        end
    end

    avg(f) = dropdims(mean(raw[f], dims = 1), dims = 1)
    terms = NamedTuple{ERROR_TERM_FIELDS}(avg(f) for f in ERROR_TERM_FIELDS)

    qmatrix_j = terms.mean_j .+ terms.mode_j .+ terms.modedist_j .+ terms.sd_j .+ terms.median_j .+ terms.q25_j .+ terms.q75_j
    qmatrix_a = terms.mean_a .+ terms.mode_a .+ terms.modedist_a .+ terms.sd_a .+ terms.median_a .+ terms.q25_a .+ terms.q75_a

    return (; qmatrix_j, qmatrix_a, terms)
end

"""
    best_fit(qmatrix, config)

Grid point minimizing `qmatrix` (an (nsig, ntau) error surface), returning
`(; error, sigtau, tau)`.
"""
function best_fit(qmatrix::AbstractMatrix, config::SiteConfig)
    val, idx = findmin(qmatrix)
    return (; error = val, sigtau = config.sigtauvec[idx[1]], tau = config.tauvec[idx[2]])
end

"""
    compare_site(config, data_dir, measures; path_fn=simdata_path)

Full site-comparison pipeline (the numerical core of `sharkcompare_*.jl`, excluding
figure generation): builds the error surfaces via `error_surface`, then reports the
best-fit `(sigtau, tau)` and minimum error under both the juvenile-site and
adult-site hypotheses. A lower `best_j.error` vs `best_a.error` favors interpreting
the empirical sample as a nursery/juvenile site, and vice versa (see the paper's
eq. 2.1 and section 2c).

Pass `path_fn = legacy_simdata_path` to read an already-built legacy-format library
(e.g. `../data/sharks_modern2/`) instead of one built by `build_library!`.
"""
function compare_site(config::SiteConfig, data_dir::AbstractString, measures::AbstractVector;
        path_fn::Function = simdata_path)
    surf = error_surface(config, data_dir, measures; path_fn)
    best_j = best_fit(surf.qmatrix_j, config)
    best_a = best_fit(surf.qmatrix_a, config)
    return (; surf.qmatrix_j, surf.qmatrix_a, surf.terms, best_j, best_a)
end
