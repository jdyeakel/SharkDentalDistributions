"""
    growth_time(epsilon1, epsilon2, a, eta)

Time (seconds) to grow from proportion `epsilon1` to `epsilon2` of asymptotic mass `M`,
under the West et al. (2001) ontogenetic growth model. Ported verbatim from the
original `ts.jl` (renamed for clarity; `M` is captured from the enclosing scope in the
original -- here it is threaded through explicitly via `growth_time(...,M,...)`).
"""
function growth_time(epsilon1, epsilon2, a, eta, M)
    return log((1 - epsilon1^(1 - eta)) / (1 - epsilon2^(1 - eta))) * (M^(1 - eta) / (a * (1 - eta)))
end

# Fixed biological constants used across all scenarios (not scenario parameters).
const METABOLIC_C = 18.47          # FISH normalization constant
const METABOLIC_EM = 5774.0        # energy density to synthesize a unit of mass (J/g)
const METABOLIC_ETA = 3 / 4        # metabolic scaling exponent
# NB: written as `X * 10^(-N.0)` rather than the equivalent-looking `Xe-N` literal
# deliberately -- Julia's runtime `^` for a Float64 exponent doesn't always round to
# the same bit pattern as the directly-parsed scientific-notation literal (e.g.
# `5.70776*10^(-9.0) != 5.70776e-9`), and the original source computes these at
# runtime this way, so matching the exact expression is needed for bit-for-bit
# equivalence with the legacy simulation (see check_legacy_equivalence.jl).
const MORTALITY_RATE = 5.70776 * 10^(-9.0)  # per-individual mortality rate (Schindler et al. 2002)
const RMAX_WARM = 1.47451 * 10^(-7.0)       # max per-capita reproductive rate (Cortes & Parsons 1996)
const RMAX_COLD = RMAX_WARM        # no cold/warm difference assumed (kept distinct for clarity)
const EPSILON_STEP = 50            # number of mass-class bins
const EPSILON_MAX = 0.99
const DAYS_IN_YEAR = 365
const SECONDS_IN_DAY = 24 * 60 * 60
const SECONDS_IN_YEAR = DAYS_IN_YEAR * SECONDS_IN_DAY
const UPPER_TEETH_LOST_PER_DAY = 1 / 40
const LOWER_TEETH_LOST_PER_DAY = 1 / 40

"""
    temperature_vector(tmin_c, tmax_c; n=100)

100-point (by default) interpolation of a site's seasonal temperature range, in Kelvin.
Matches the original inline logic in `popgen_migrate_g.jl` exactly, including the
degenerate case `tmin_c == tmax_c`.
"""
function temperature_vector(tmin_c, tmax_c; n::Int = 100)
    tmin = tmin_c + 273.15
    tmax = tmax_c + 273.15
    return tmin == tmax ? fill(tmin, n) : collect(range(tmin, tmax, length = n))
end

"""
    precompute_growth(config::SiteConfig)

Computes everything that depends only on `config` (birth/asymptotic mass, per-site
temperature vectors, growth-time and reproduction-rate tables, initial-state
bookkeeping) -- i.e. everything that does *not* depend on the (sigtau, tau) dispersal
parameters swept in a library build. The original scripts recomputed this fresh for
every one of the `its` grid cells; here it's computed once per config and reused,
which is a free, behavior-preserving optimization (sigtau/tau never affected these
tables).
"""
function precompute_growth(config::SiteConfig)
    m0 = (0.00013 * config.l0^2.4) * 1000
    M = (0.00013 * config.L^2.4) * 1000

    tempvec1 = temperature_vector(config.mintemp_j, config.maxtemp_j)
    tempvec2 = temperature_vector(config.mintemp_a, config.maxtemp_a)

    epsilonvec = collect(m0 / M:(EPSILON_MAX - m0 / M) / EPSILON_STEP:EPSILON_MAX)
    lspan = length(epsilonvec)
    ltime = lspan - 1
    ltemp = length(tempvec1)

    tint1 = Array{Float64}(undef, ltemp, ltime)
    tint2 = Array{Float64}(undef, ltemp, ltime)
    mass1 = Array{Float64}(undef, ltemp, ltime)
    mass2 = Array{Float64}(undef, ltemp, ltime)
    for k in 1:ltemp
        temp1 = tempvec1[k]
        temp2 = tempvec2[k]
        B01 = exp(METABOLIC_C) / exp(0.63 / (8.61733 * 10^(-5.0) * temp1))
        B02 = exp(METABOLIC_C) / exp(0.63 / (8.61733 * 10^(-5.0) * temp2))
        a1 = B01 / METABOLIC_EM
        a2 = B02 / METABOLIC_EM
        for i in 1:ltime
            epsilon1 = epsilonvec[i]
            epsilon2 = epsilonvec[i + 1]
            tint1[k, i] = growth_time(epsilon1, epsilon2, a1, METABOLIC_ETA, M)
            tint2[k, i] = growth_time(epsilon1, epsilon2, a2, METABOLIC_ETA, M)
            mass1[k, i] = (epsilon1 * M + epsilon2 * M) * 0.5
            mass2[k, i] = (epsilon1 * M + epsilon2 * M) * 0.5
        end
    end

    tvec1 = cumsum(tint1, dims = 2)
    tvec2 = cumsum(tint2, dims = 2)

    juvpos = Int64(floor(ltime / 4))

    mintemp1, maxtemp1 = extrema(tempvec1)
    mintemp2, maxtemp2 = extrema(tempvec2)
    r_sizetemp1 = Array{Float64}(undef, ltemp, ltime)
    r_sizetemp2 = Array{Float64}(undef, ltemp, ltime)
    for i in 1:ltemp
        rtemp1 = RMAX_COLD + (RMAX_WARM - RMAX_COLD) * ((tempvec1[i] - mintemp1) / (maxtemp1 - mintemp1))
        rtemp2 = RMAX_COLD + (RMAX_WARM - RMAX_COLD) * ((tempvec2[i] - mintemp2) / (maxtemp2 - mintemp2))
        for j in 1:ltime
            r_sizetemp1[i, j] = rtemp1
            r_sizetemp2[i, j] = rtemp2
        end
    end

    tmax = mean([maximum(tvec1), maximum(tvec2)]) * config.gen

    temptime1 = mean(tempvec1) .+ (maximum(tempvec1) .- mean(tempvec1)) .* sin.((pi / (DAYS_IN_YEAR / 2)) .* collect(0:1:DAYS_IN_YEAR))
    temptime2 = mean(tempvec2) .+ (maximum(tempvec2) .- mean(tempvec2)) .* sin.((pi / (DAYS_IN_YEAR / 2)) .* collect(0:1:DAYS_IN_YEAR))

    toothlossrate = (LOWER_TEETH_LOST_PER_DAY + UPPER_TEETH_LOST_PER_DAY) / 24 / 60 / 60

    # Bookkeeping arrays that depend only on `ltime` (constant across the whole
    # sigtau/tau/rep sweep for a given config) -- precomputed once rather than
    # rebuilt inside the stochastic loop.
    states_index = repeat(collect(1:ltime), outer = 4)
    cistaterate = CartesianIndices((4 * ltime, 2))

    return (;
        m0, M, tempvec1, tempvec2, epsilonvec, ltime, ltemp,
        tint1, tint2, mass1, mass2, r_sizetemp1, r_sizetemp2,
        juvpos, tmax, temptime1, temptime2, toothlossrate,
        states_index, cistaterate,
    )
end
