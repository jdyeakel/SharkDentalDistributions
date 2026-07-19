"""
    sigtau_mass(config)

Converts `config.sigtauvec` (juvenile dispersal-window steepness, in mass-bin units)
to kg, for x-axis labeling -- matches the `sigtauvecmass` computation in the legacy
`sharkcompare_*.jl` R-plotting sections (`((sigtauvec .* M) ./ 50) ./ 1000`, where 50
is `EPSILON_STEP` and `M` is the asymptotic mass in grams).
"""
function sigtau_mass(config::SiteConfig)
    M = (0.00013 * config.L^2.4) * 1000
    return ((config.sigtauvec .* M) ./ EPSILON_STEP) ./ 1000
end

"""
    best_fit_densities(config, data_dir, best, site_type; path_fn=simdata_path)

Loads every replicate's simulated tooth-length density at the `(sigtau, tau)` grid
cell identified by `best` (a NamedTuple from `best_fit`, e.g. `compare_site(...).best_j`),
for `site_type` (`:juv` -> toothdrop column 1, `:adult` -> column 2), each scaled so
its maximum is 1 (matching the legacy `plotcompare.jl`'s `scaledsimdensity`). Returns
`(toothlength, densities::Vector{Vector{Float64}})`, one density vector per replicate.
"""
function best_fit_densities(config::SiteConfig, data_dir::AbstractString, best, site_type::Symbol;
        path_fn::Function = simdata_path)
    s = findfirst(==(best.sigtau), config.sigtauvec)
    t = findfirst(==(best.tau), config.tauvec)
    col = site_type === :juv ? 1 : 2
    toothlength = Float64[]
    densities = Vector{Vector{Float64}}(undef, config.reps)
    for r in 1:config.reps
        d = load_simdata(data_dir, r, s, t; path_fn)
        toothlength = d.toothlength1[1, :]
        density = d.toothdrop[:, col] ./ sum(d.toothdrop[:, col])
        densities[r] = density ./ maximum(density)
    end
    return toothlength, densities
end

"""
    plot_scenario(config, data_dir, sites; path_fn=simdata_path, filename=nothing)

Julia-native replacement for the R/RCall plotting sections of `sharkcompare_*.jl`.
`sites` is a vector of `(label, measures)` pairs (one per empirical site being tested
against this `config`'s simulation library). Produces one row of 4 panels per site --
juvenile-hypothesis error-surface heatmap (best fit marked), adult-hypothesis
heatmap, juvenile density comparison (empirical vs. best-fit simulated, all reps),
adult density comparison -- laid out exactly like the legacy figures
(`fig_empirical_comp_*`), with a shared color scale across every panel/site so
error magnitudes are comparable at a glance. Returns the `Plots.jl` figure; also
saves it to `filename` if given.
"""
function plot_scenario(config::SiteConfig, data_dir::AbstractString, sites;
        path_fn::Function = simdata_path, filename::Union{Nothing, AbstractString} = nothing)
    n = length(sites)
    results = [compare_site(config, data_dir, measures; path_fn) for (_, measures) in sites]

    zmin = minimum(min(minimum(r.qmatrix_j), minimum(r.qmatrix_a)) for r in results)
    zmax = maximum(max(maximum(r.qmatrix_j), maximum(r.qmatrix_a)) for r in results)
    sigmass = sigtau_mass(config)

    # `best_j.sigtau`/`best_a.sigtau` are in raw `config.sigtauvec` units, but the
    # heatmap x-axis is in mass units (`sigmass`) -- convert by looking up the
    # matching grid index, rather than plotting the raw value on the mass axis.
    sigmass_at(sigtau) = sigmass[findfirst(==(sigtau), config.sigtauvec)]

    panels = Plots.Plot[]
    for (i, (label, measures)) in enumerate(sites)
        res = results[i]
        U = kde(measures)
        emp_density = U.density ./ maximum(U.density)

        hj = Plots.heatmap(sigmass, config.tauvec, res.qmatrix_j'; clims = (zmin, zmax),
            color = :Spectral, colorbar = false,
            xlabel = "", ylabel = i == 1 ? "Adult migration window" : "",
            title = i == 1 ? "Juvenile-site error" : "")
        Plots.scatter!(hj, [sigmass_at(res.best_j.sigtau)], [res.best_j.tau]; color = :white,
            markerstrokecolor = :black, markersize = 6, label = false)

        ha = Plots.heatmap(sigmass, config.tauvec, res.qmatrix_a'; clims = (zmin, zmax),
            color = :Spectral, colorbar = false,
            title = i == 1 ? "Adult-site error" : "")
        Plots.scatter!(ha, [sigmass_at(res.best_a.sigtau)], [res.best_a.tau]; color = :white,
            markerstrokecolor = :black, markersize = 6, label = false)

        # Mark whichever hypothesis (juvenile vs. adult) has the lower error for this
        # site with an asterisk -- computed from the actual values each time (not
        # hardcoded by panel position), matching the automatic-asterisk fix already
        # applied to the legacy R scripts.
        star_j = res.best_j.error < res.best_a.error ? "*" : ""
        star_a = res.best_a.error < res.best_j.error ? "*" : ""

        toothlength_j, densities_j = best_fit_densities(config, data_dir, res.best_j, :juv; path_fn)
        pj = Plots.plot(U.x, emp_density; color = :black, linewidth = 2, label = false,
            xlabel = i == n ? "Tooth length (mm)" : "", ylabel = "Scaled density",
            title = i == 1 ? "Juvenile-site fit" : "")
        for d in densities_j
            Plots.plot!(pj, toothlength_j, d; color = :gray, alpha = 0.4, label = false)
        end
        Plots.annotate!(pj, maximum(U.x) * 0.98, 0.92,
            Plots.text("$label\nerror = $(round(res.best_j.error, digits = 2))$star_j", 8, :right))

        toothlength_a, densities_a = best_fit_densities(config, data_dir, res.best_a, :adult; path_fn)
        pa = Plots.plot(U.x, emp_density; color = :black, linewidth = 2, label = false,
            xlabel = i == n ? "Tooth length (mm)" : "",
            title = i == 1 ? "Adult-site fit" : "")
        for d in densities_a
            Plots.plot!(pa, toothlength_a, d; color = :gray, alpha = 0.4, label = false)
        end
        Plots.annotate!(pa, maximum(U.x) * 0.98, 0.92,
            Plots.text("$label\nerror = $(round(res.best_a.error, digits = 2))$star_a", 8, :right))

        append!(panels, [hj, ha, pj, pa])
    end

    fig = Plots.plot(panels...; layout = (n, 4), size = (1200, 280 * n), legend = false)
    if filename !== nothing
        Plots.savefig(fig, filename)
    end
    return fig
end
