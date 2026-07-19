"""
Parameters for a two-site (juvenile/nursery + adult) shark metapopulation scenario.

`sigtauvec`/`tauvec` are the (juvenile dispersal window, adult dispersal window) grids
swept when building a simulation library; `reps` is the number of stochastic replicates
per grid cell. Temperatures are in °C, `distance_km` in km, `l0`/`L` (birth/asymptotic
body length) in cm.
"""
Base.@kwdef struct SiteConfig
    l0::Float64
    L::Float64
    n0::Int
    gen::Int
    mintemp_j::Float64
    maxtemp_j::Float64
    mintemp_a::Float64
    maxtemp_a::Float64
    distance_km::Float64
    velocity::Float64 = 1.0
    D::Float64 = 1.0
    sigtauvec::Vector{Float64}
    tauvec::Vector{Float64}
    reps::Int
end

# SiteConfig isn't `isbits` (it holds Vector fields), so Julia doesn't auto-derive a
# field-wise `==` for it and falls back to identity comparison. Define it explicitly
# so e.g. `load_settings(dir) == config` works as expected after a JLD2 round-trip.
function Base.:(==)(a::SiteConfig, b::SiteConfig)
    return all(getfield(a, f) == getfield(b, f) for f in fieldnames(SiteConfig))
end

function _preset(defaults::NamedTuple; kwargs...)
    overrides = NamedTuple(kwargs)
    return SiteConfig(; merge(defaults, overrides)...)
end

"""
    modern(; kwargs...)

Extant Delaware Bay sand tiger (*Carcharias taurus*) scenario. Historical parameters
from `sharksims_modern.jl`. Any field can be overridden via keyword, e.g.
`modern(reps=10)`.
"""
modern(; kwargs...) = _preset((
    l0 = 90.0, L = 300.0, n0 = 1000, gen = 1,
    mintemp_j = 17.0, maxtemp_j = 25.0, mintemp_a = 13.0, maxtemp_a = 23.0,
    distance_km = 700.0, velocity = 1.0, D = 1.0,
    sigtauvec = collect(0.5:1.0:20.0), tauvec = collect(1.0:2.0:40.0), reps = 5,
); kwargs...)

"""
    eocene_highlatitude(; kwargs...)

Eocene high-latitude scenario (Banks Island / Seymour Island). Historical parameters
from `sharksims_eocene.jl`, including the corrected `L=477` asymptotic length (see
`L=295` vs `L=477` history: the earlier `sharks_eocene` library used an asymptotic
length too small to reach the empirical max tooth length; `sharks_eocene2` fixed this).
"""
eocene_highlatitude(; kwargs...) = _preset((
    l0 = 50.0, L = 477.0, n0 = 1000, gen = 1,
    mintemp_j = 12.0, maxtemp_j = 24.0, mintemp_a = 9.0, maxtemp_a = 17.0,
    distance_km = 400.0, velocity = 1.0, D = 1.0,
    sigtauvec = collect(0.5:0.5:25.0), tauvec = collect(1.0:1.0:50.0), reps = 25,
); kwargs...)

"""
    eocene_lowlatitude(; kwargs...)

Eocene low-latitude scenario (Red Hot Truck Stop / Whiskey Bridge). Historical
parameters from `sharksims_eocene_lowlatitude.jl`. Note juvenile and adult site
temperature ranges are identical here (23-30°C) -- low-latitude sites have little
thermal contrast between nursery and adult habitat, unlike `eocene_highlatitude`
or `modern`.
"""
eocene_lowlatitude(; kwargs...) = _preset((
    l0 = 50.0, L = 477.0, n0 = 1000, gen = 1,
    mintemp_j = 23.0, maxtemp_j = 30.0, mintemp_a = 23.0, maxtemp_a = 30.0,
    distance_km = 400.0, velocity = 1.0, D = 1.0,
    sigtauvec = collect(0.5:1.0:20.0), tauvec = collect(1.0:2.0:40.0), reps = 5,
); kwargs...)
