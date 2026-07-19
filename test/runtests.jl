using Test
using Random
using SharkDentalDynamics

@testset "presets" begin
    m = modern()
    @test m.l0 == 90.0 && m.L == 300.0 && m.n0 == 1000 && m.distance_km == 700.0
    @test length(m.sigtauvec) == 20 && length(m.tauvec) == 20 && m.reps == 5

    eh = eocene_highlatitude()
    @test eh.l0 == 50.0 && eh.L == 477.0 && eh.distance_km == 400.0
    @test length(eh.sigtauvec) == 50 && length(eh.tauvec) == 50 && eh.reps == 25

    el = eocene_lowlatitude()
    @test el.l0 == 50.0 && el.L == 477.0 && el.distance_km == 400.0
    @test el.mintemp_j == el.mintemp_a && el.maxtemp_j == el.maxtemp_a
    @test length(el.sigtauvec) == 20 && length(el.tauvec) == 20 && el.reps == 5

    # keyword overrides work
    tiny = modern(reps = 2, sigtauvec = [1.0, 2.0], tauvec = [5.0, 10.0])
    @test tiny.reps == 2 && tiny.sigtauvec == [1.0, 2.0]
    @test tiny.l0 == 90.0  # untouched fields keep the preset default
end

@testset "precompute_growth" begin
    for cfg in (modern(), eocene_highlatitude(), eocene_lowlatitude())
        precomp = precompute_growth(cfg)
        @test precomp.ltime == 50  # EPSILON_STEP
        @test size(precomp.tint1) == (100, precomp.ltime)
        @test size(precomp.r_sizetemp1) == (100, precomp.ltime)
        @test precomp.juvpos == Int64(floor(precomp.ltime / 4))
        @test precomp.m0 < precomp.M
    end
end

@testset "local_maxima / local_minima" begin
    signal = [1.0, 3.0, 2.0, 5.0, 1.0, 4.0, 0.5]
    @test local_maxima(signal) == [2, 4, 6]
    # index 1 is flagged via the (correct) first-point rule (1.0 < 3.0); index 7
    # is a real downward trend at the end but is NOT flagged, because of the
    # inherited endpoint bug (see docstring): signal[7]=0.5 > signal[6]=4.0 is
    # false, so the buggy ">" check never fires here either way.
    @test local_minima(signal) == [1, 3, 5]
end

@testset "shape_descriptors" begin
    support = collect(0.0:1.0:10.0)
    density = exp.(-((support .- 5.0) .^ 2) ./ 2)  # unimodal Gaussian-ish bump at 5
    d = shape_descriptors(density, support)
    @test isapprox(d.mean, 5.0; atol = 0.1)
    @test d.peak_bin == 0
    @test d.peak_dist == 0.0
    @test d.quartile25 < d.median < d.quartile75

    # bimodal density should be detected as such
    bimodal = exp.(-((support .- 2.0) .^ 2) ./ 0.5) .+ 0.8 .* exp.(-((support .- 8.0) .^ 2) ./ 0.5)
    dm = shape_descriptors(bimodal, support)
    @test dm.peak_bin == 1
    @test dm.peak_dist > 0
end

@testset "simulate_metapopulation determinism and shape" begin
    cfg = modern()
    precomp = precompute_growth(cfg)
    sigtau, tau = cfg.sigtauvec[1], cfg.tauvec[1]

    r1 = simulate_metapopulation(cfg, precomp, sigtau, tau; rng = Random.Xoshiro(42))
    r2 = simulate_metapopulation(cfg, precomp, sigtau, tau; rng = Random.Xoshiro(42))
    @test r1.toothdrop == r2.toothdrop
    @test r1.state == r2.state
    @test r1.clock == r2.clock

    @test size(r1.toothdrop) == (precomp.ltime, 2)
    @test all(r1.toothdrop .>= 0)
    @test sum(r1.state) >= 0

    r3 = simulate_metapopulation(cfg, precomp, sigtau, tau; rng = Random.Xoshiro(43))
    @test r1.toothdrop != r3.toothdrop  # different seed -> (almost certainly) different trajectory
end

@testset "compare_to_empirical" begin
    cfg = modern()
    precomp = precompute_growth(cfg)
    result = simulate_metapopulation(cfg, precomp, cfg.sigtauvec[1], cfg.tauvec[1]; rng = Random.Xoshiro(1))
    toothlength = tooth_length(result.mass1)[1, :]
    measures = 5.0 .+ 20.0 .* rand(Random.Xoshiro(7), 100)

    res = compare_to_empirical(result.toothdrop, toothlength, measures)
    @test all(isfinite, values(res))
    @test all(v -> v >= 0, values(res))
end

@testset "build_library! + compare_site end-to-end (tiny grid)" begin
    mktempdir() do tmp
        cfg = modern(sigtauvec = [2.0, 8.0], tauvec = [5.0, 20.0], reps = 2)
        build_library!(cfg, tmp; seed = 99)

        # every expected file exists
        for r in 1:cfg.reps, s in 1:length(cfg.sigtauvec), t in 1:length(cfg.tauvec)
            @test isfile(simdata_path(tmp, r, s, t))
        end
        @test isfile(settings_path(tmp))
        @test load_settings(tmp) == cfg

        measures = 5.0 .+ 20.0 .* rand(Random.Xoshiro(11), 100)
        res = compare_site(cfg, tmp, measures)

        @test size(res.qmatrix_j) == (2, 2)
        @test size(res.qmatrix_a) == (2, 2)
        @test all(isfinite, res.qmatrix_j)
        @test res.best_j.sigtau in cfg.sigtauvec
        @test res.best_j.tau in cfg.tauvec
        @test res.best_a.error >= 0

        # rebuild without overwrite should be a no-op (files untouched, same content)
        mtimes_before = [mtime(simdata_path(tmp, 1, 1, 1))]
        build_library!(cfg, tmp; seed = 99)
        @test mtime(simdata_path(tmp, 1, 1, 1)) == mtimes_before[1]
    end
end
