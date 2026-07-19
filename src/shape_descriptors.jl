"""
    local_maxima(signal::AbstractVector)

Indices of local maxima in `signal` (interior points higher than both neighbors, plus
either endpoint if it's higher than its single neighbor). Ported verbatim from
`src/findlocalmaxima.jl`.
"""
function local_maxima(signal::AbstractVector)
    inds = Int[]
    if length(signal) > 1
        if signal[1] > signal[2]
            push!(inds, 1)
        end
        for i in 2:(length(signal) - 1)
            if signal[i - 1] < signal[i] > signal[i + 1]
                push!(inds, i)
            end
        end
        if signal[end] > signal[end - 1]
            push!(inds, length(signal))
        end
    end
    return inds
end

"""
    local_minima(signal::AbstractVector)

Indices of local minima in `signal`. Ported verbatim from `src/findlocalminima.jl`,
**including a pre-existing bug**: the last-point check is `signal[end] > signal[end-1]`
(the same condition used for a local *maximum*), not `<`, so a trailing local minimum
is never actually detected via that branch. Preserved as-is for exact behavioral
equivalence with the legacy pipeline (see `check_legacy_equivalence.jl`); worth fixing
in a later cleanup pass, but changing it here would change results relative to every
existing published figure.
"""
function local_minima(signal::AbstractVector)
    inds = Int[]
    if length(signal) > 1
        if signal[1] < signal[2]
            push!(inds, 1)
        end
        for i in 2:(length(signal) - 1)
            if signal[i - 1] > signal[i] < signal[i + 1]
                push!(inds, i)
            end
        end
        if signal[end] > signal[end - 1]  # NB: inherited bug, see docstring
            push!(inds, length(signal))
        end
    end
    return inds
end

"""
    detect_modes(pdist::AbstractVector, support::AbstractVector)

Detects whether `pdist` (a probability mass/density vector over `support`) has a
significant second mode, using a 10%-of-peak-height threshold on the trough between
the two tallest peaks. Ported verbatim from `src/modality_analysis.jl` (only the
active code path; the large commented-out duplicate block in the original is
dropped). Returns `(secondpeak, secondtooth, troughtooth, troughmin, maxtooth)`,
where `secondpeak == 0.0` means no significant second mode was found.
"""
function detect_modes(pdist::AbstractVector, support::AbstractVector)
    lmax = local_maxima(pdist)
    lmin = local_minima(pdist)
    peakprobs = pdist[lmax]
    troughprobs = pdist[lmin]
    sortlmax = lmax[sortperm(peakprobs)]
    sortlmin = lmin[sortperm(troughprobs)]
    sortedpeaks = pdist[sortlmax]
    maxpeak = last(sortedpeaks)
    sortedtooth = support[lmax[sortperm(peakprobs)]]
    maxtooth = last(sortedtooth)

    secondpeak = 0.0
    secondtooth = 0.0
    troughtooth = 0.0
    troughmin = 0.0

    if length(sortlmin) > 1
        troughpos1 = sortlmin[findall(x -> (x < sortlmax[end] && x > sortlmax[end - 1]), sortlmin)]
        troughpos2 = sortlmin[findall(x -> (x > sortlmax[end] && x < sortlmax[end - 1]), sortlmin)]
        troughpos = [troughpos1; troughpos2]

        troughposbackup = Int[]
        if length(sortlmax) > 2
            troughpos3 = sortlmin[findall(x -> (x < sortlmax[end] && x > sortlmax[end - 2]), sortlmin)]
            troughpos4 = sortlmin[findall(x -> (x > sortlmax[end] && x < sortlmax[end - 2]), sortlmin)]
            troughposbackup = [troughpos3; troughpos4]
        end

        trougha1 = pdist[troughpos][findmin(pdist[troughpos])[2]]
        trougha2 = Float64[]
        if length(troughposbackup) > 0
            trougha2 = pdist[troughposbackup][findmin(pdist[troughposbackup])[2]]
        end

        peakprop = [-1.0, -1.0]
        peakpos = 0
        if length(sortedpeaks) > 1
            peakprop[1] = (sortedpeaks[end - 1] - trougha1) / (maxpeak - trougha1)
            if length(troughposbackup) > 0
                peakprop[2] = (sortedpeaks[end - 2] - trougha2) / (maxpeak - trougha2)
            end
            if peakprop[1] > peakprop[2]
                peakpos = 1
            end
            if peakprop[2] > peakprop[1]
                peakpos = 2
            end
            if peakprop[peakpos] > 0.10
                secondpeak = sortedpeaks[end - peakpos]
                secondtooth = sortedtooth[end - peakpos]
            end
        end

        troughposmin = length(troughposbackup) > 0 ? [troughpos[1], troughposbackup[1]][peakpos] : troughpos[1]
        troughmin = [trougha1, trougha2][peakpos]
        troughtooth = support[troughposmin]
    end

    return secondpeak, secondtooth, troughtooth, troughmin, maxtooth
end

"""
    shape_descriptors(density::AbstractVector, support::AbstractVector)

Unified replacement for `src/toothdist_analysis.jl` (which computed this for both
juvenile and adult simulated densities inline) and `src/toothdist_emp_analysis.jl`
(the near-identical empirical version) -- both were ~90% duplicated code for mean,
variance, median, quartiles, and modality of a density-over-support pair. Call this
once per density (simulated juvenile, simulated adult, empirical) instead.

Returns a NamedTuple `(; mean, var, modes, peak_bin, peak_dist, median, quartile25,
quartile75)`. `modes` is `sort([secondtooth, maxtooth])` unconditionally (matching the
original: when there's no significant second mode, `secondtooth` is `0.0`, so `modes`
is `[0.0, primary_mode_location]`); `peak_bin` is `1` iff a significant second mode was
found, and `peak_dist` is the distance between the two modes in that case (else `0`).
"""
function shape_descriptors(density::AbstractVector, support::AbstractVector)
    p = density ./ sum(density)
    m = dot(support, p)
    v = dot(support .^ 2, p) - m^2

    cump = cumsum(p)
    med = support[findfirst(x -> x > 0.5, cump)]
    q25 = support[findfirst(x -> x > 0.25, cump)]
    q75 = support[findfirst(x -> x > 0.75, cump)]

    secondpeak, secondtooth, _, _, maxtooth = detect_modes(p, support)
    peak_bin = secondpeak == 0.0 ? 0 : 1
    peak_dist = secondpeak == 0.0 ? 0.0 : abs(secondtooth - maxtooth)
    modes = sort([secondtooth, maxtooth])

    return (; mean = m, var = v, modes, peak_bin, peak_dist, median = med, quartile25 = q25, quartile75 = q75)
end
