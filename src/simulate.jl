"""
    SimulationTimeout(sigtau, tau, elapsed, n_events)

Thrown by `simulate_metapopulation` when `max_seconds` is exceeded. Some
`(sigtau, tau)` corners let the metapopulation grow substantially within one `gen`
cycle -- since the Gillespie event size `dt ~ 1/N`, a growing population needs far
more discrete events to simulate the same time span, and such a corner can run for
an hour or more (one did during development: eocene_highlatitude at large
sigtau + large tau, ~10GB RAM before being killed manually). `build_library!` catches
this and skips the cell (with a warning) rather than blocking the whole sweep.
"""
struct SimulationTimeout <: Exception
    sigtau::Float64
    tau::Float64
    elapsed::Float64
    n_events::Int
end

function Base.showerror(io::IO, e::SimulationTimeout)
    print(io, "SimulationTimeout: sigtau=$(e.sigtau), tau=$(e.tau) exceeded $(round(e.elapsed, digits=1))s ",
        "after $(e.n_events) events -- likely a pathologically expensive grid corner")
end

"""
    simulate_metapopulation(config::SiteConfig, precomp, sigtau, tau; rng=Random.default_rng(), max_seconds=nothing)

Stochastic (Gillespie) simulation of the two-site shark metapopulation for one
particular `(sigtau, tau)` dispersal-window pair. This is a faithful port of the
original `popgen_migrate_g` in `src/popgen_migrate_g.jl` -- same algorithm, same order
of random draws (`rand(rng)` for the event-selection draw, then conditionally
`rand(rng, ...)` for the birth-state draw) -- restructured to:

  - take an explicit `rng` (needed for reproducible `Threads`-based parallelism --
    each library-build task gets its own RNG rather than sharing global state), and
  - take precomputed, config-only growth tables (`precompute_growth`) instead of
    recomputing them on every call (sigtau/tau never affected those tables).

Two purely-deterministic simplifications versus the original (verified not to touch
the random-draw sequence, and confirmed by `check_legacy_equivalence.jl`):
  - `findfirst(isodd, ...)` replaces `findall(isodd, ...)[1]` (same result, no
    intermediate array of every matching index).
  - the `Apos == 4` (migrate) branch drops the original's `growthovermigration`
    bookkeeping, which was hardcoded to `0` in both its branches (dead code -- the
    `if growthovermigration > 1 ...` condition could never be true), leaving only the
    single migrate-without-additional-growth update it always actually performed.

Returns a NamedTuple `(; mass1, mass2, epsilonvec, clock, popstate, toothdrop, state)`,
matching the original's return values.

`max_seconds` (disabled by default, `nothing`) throws `SimulationTimeout` if wall-clock
time exceeds it, checked every 5,000 events (cheap relative to the per-event cost, so
this doesn't measurably affect normal-run timing). Useful for guarding against the
pathologically expensive grid corners described in `SimulationTimeout`'s docstring;
`build_library!` sets this by default so one bad cell can't block an entire sweep.
"""
function simulate_metapopulation(config::SiteConfig, precomp, sigtau, tau;
        rng::Random.AbstractRNG = Random.default_rng(),
        max_seconds::Union{Nothing, Real} = nothing)
    n0 = config.n0
    distance = config.distance_km * 1000.0
    velocity = config.velocity
    D = config.D

    ltime = precomp.ltime
    juvpos = precomp.juvpos
    tempvec1 = precomp.tempvec1
    tempvec2 = precomp.tempvec2
    tint1 = precomp.tint1
    tint2 = precomp.tint2
    r_sizetemp1 = precomp.r_sizetemp1
    temptime1 = precomp.temptime1
    temptime2 = precomp.temptime2
    tmax = precomp.tmax
    toothlossrate = precomp.toothlossrate
    states_index = precomp.states_index
    cistaterate = precomp.cistaterate

    peakday = 180

    state = zeros(Int64, ltime, 2)
    state[1, 1] = n0
    state[juvpos, 2] = n0

    toothdrop = zeros(Float64, ltime, 2)

    pop1 = Array{Int64}(undef, 0)
    pop2 = Array{Int64}(undef, 0)
    clock = Array{Float64}(undef, 0)

    tcum = 0.0
    day = 1
    n_events = 0
    start_time = max_seconds === nothing ? 0.0 : time()
    while tcum < tmax
        n_events += 1
        if max_seconds !== nothing && n_events % 5000 == 0
            elapsed = time() - start_time
            elapsed > max_seconds && throw(SimulationTimeout(sigtau, tau, elapsed, n_events))
        end

        daytemp1 = temptime1[day]
        daytemp2 = temptime2[day]

        k1 = findmin((tempvec1 .- daytemp1) .^ 2)[2]
        k2 = findmin((tempvec2 .- daytemp2) .^ 2)[2]

        if sum(state) == 0
            break
        end

        N1 = state[:, 1]
        N2 = state[:, 2]

        # Reproduction rate per size class/site (nonlinear in temperature; no
        # pre-juvenile reproduction; none at the adult site).
        r1 = r_sizetemp1[k1, :]
        r1[1:juvpos] .= 0
        r2 = repeat([0], ltime)

        # Growth rate per size class/site.
        g1 = 1 ./ tint1[k1, :]
        g2 = 1 ./ tint2[k2, :]

        # Mortality rate per size class/site (elevated in the terminal class).
        d2 = repeat([MORTALITY_RATE], ltime)
        d2[ltime] = MORTALITY_RATE * 10
        d1 = copy(d2)

        # Migration rate: sigmoidal in mass for juveniles leaving the nursery,
        # a seasonal (Gaussian) pulse in time for adults returning to it.
        m1 = D * (velocity / distance) ./ (1 .+ exp.(-(1 / sigtau) .* (collect(1:ltime) .- juvpos)))
        m2 = repeat([D * (velocity / distance) * exp(-((day - peakday)^2) / (2 * tau^2))], ltime)

        Rate = dot(N1, r1 .+ g1 .+ d1 .+ m1) + dot(N2, r2 .+ g2 .+ d2 .+ m2)
        dt = 1 / Rate

        NAline_staterates = [
            state[:, 1] .* r1; state[:, 1] .* g1; state[:, 1] .* d1; state[:, 1] .* m1;
            state[:, 2] .* r2; state[:, 2] .* g2; state[:, 2] .* d2; state[:, 2] .* m2;
        ]
        NAline = cumsum(NAline_staterates / sum(NAline_staterates))

        NArand = rand(rng)
        NApos = findfirst(isodd, NArand .<= NAline)

        Nloc = cistaterate[NApos][2]
        Nstate = states_index[cistaterate[NApos][1]]
        Apos = findfirst(isodd, cistaterate[NApos][1] .<= (ltime .* [1, 2, 3, 4]))
        Altloc = setdiff([1, 2], Nloc)[1]

        if Apos == 1
            # Reproduce: a newborn appears in a random early mass class.
            rbirthstate = rand(rng, collect(1:Int64(EPSILON_STEP / 10)))
            state[rbirthstate, Nloc] += 1
        elseif Apos == 2
            # Grow: advance one mass class (clamped at the terminal class).
            state[Nstate, Nloc] -= 1
            if Nstate + 1 < ltime
                state[Nstate + 1, Nloc] += 1
            else
                state[ltime, Nloc] += 1
            end
        elseif Apos == 3
            # Die.
            state[Nstate, Nloc] -= 1
        elseif Apos == 4
            # Migrate to the other site (same mass class).
            state[Nstate, Nloc] -= 1
            state[Nstate, Altloc] += 1
        end

        tcum += dt
        push!(clock, tcum)

        day = round(Int64, DAYS_IN_YEAR * ((tcum / SECONDS_IN_YEAR) - floor(tcum / SECONDS_IN_YEAR)))
        day = max(1, day)

        pop = vec(sum(state, dims = 1))
        push!(pop1, pop[1])
        push!(pop2, pop[2])

        toothdrop .+= state .* (toothlossrate * dt)

        if any(state .< 0)
            @warn "negative population state encountered; aborting run" sigtau tau
            break
        end
    end

    popstate = [pop1 pop2]

    return (;
        mass1 = precomp.mass1, mass2 = precomp.mass2, epsilonvec = precomp.epsilonvec,
        clock, popstate, toothdrop, state,
    )
end

"""
    tooth_length(mass)

Anterior tooth crown height (mm) from body mass (g), Shimada (2004) allometry.
"""
tooth_length(mass) = 2.13337 .+ (0.187204 .* mass .^ (0.416667))
