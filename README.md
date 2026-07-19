# SharkDentalDynamics

A clean, config-driven Julia package for simulating a two-site shark metapopulation
(juvenile/nursery site + adult site) and comparing its simulated tooth-drop size
distributions against empirical shark tooth measurements — the model behind Kim,
Yeakel et al. (2022), *"Decoding the dynamics of dental distributions,"* Proc. R.
Soc. B 289: 20220808.

It is a from-scratch reimplementation of the original `sharksims_*.jl` /
`sharkcompare_*.jl` scripts and `src/*.jl` helpers one level up in this repository —
not a fork of them. It is validated to reproduce their numerical output exactly (see
[Validation](#validation)) and fixes one real bug found along the way (see
[Known issues found and fixed](#known-issues-found-and-fixed)).

## Why this exists

The original pipeline works, but scenario parameters (temperature ranges, migration
distance, grid resolution, body size) were copy-pasted and hand-edited across
near-duplicate driver scripts (`sharksims_modern.jl`, `sharksims_eocene.jl`,
`sharksims_eocene_lowlatitude.jl`, ...), which is how a couple of real inconsistencies
crept in undetected. This package instead expresses every scenario as one
`SiteConfig` value with named presets, and swaps `@distributed`/`SharedArray`
parallelism for `Threads.@threads`.

## Repository layout

```
SharkDentalDynamics/
├── Project.toml                        # package manifest (deps: JLD2, KernelDensity, Plots, stdlibs)
├── Manifest.toml                       # locked dependency versions (commit this for reproducibility)
├── src/
│   ├── SharkDentalDynamics.jl          # module entry point; see its docstring for the full public API
│   ├── config.jl                       # SiteConfig struct + modern()/eocene_highlatitude()/eocene_lowlatitude() presets
│   ├── growth.jl                       # ontogenetic growth model (was ts.jl) + precompute_growth (config-only setup)
│   ├── simulate.jl                     # the Gillespie metapopulation kernel (was popgen_migrate_g.jl) + tooth_length() + SimulationTimeout
│   ├── shape_descriptors.jl            # mean/median/quartiles/modality of a density (was toothdist_*_analysis.jl, modality_analysis.jl)
│   ├── compare.jl                      # empirical-vs-simulated error terms + error surface + best-fit search (was empirical_sim_comparison.jl + the sharkcompare_*.jl reduction/best-fit logic)
│   ├── library.jl                      # build_library!() + JLD2 read/write (was the sharksims_*.jl driver loops + smartpath.jl)
│   └── plotting.jl                     # native Julia (Plots.jl) figures (replaces the RCall/R sections of sharkcompare_*.jl)
├── test/
│   └── runtests.jl                     # self-contained unit tests (no dependency on the legacy code)
├── scripts/
│   ├── check_legacy_equivalence.jl     # seeded numerical comparison against the legacy ../../src/*.jl (see Validation)
│   ├── reproduce_figures.jl            # regenerates the three published-comparison figures from the existing legacy data/ libraries
│   └── example_new_site.jl             # worked template: free parameters -> fresh build_library! -> compare -> plot (see Running a new empirical-site comparison)
├── empirical/
│   ├── SandTiger_all.csv               # local copy -- the five originally-published sites (used by reproduce_figures.jl)
│   └── SandTiger_all_2026.csv          # local copy -- adds Blackheath (used by example_new_site.jl)
├── fig_modern_julia.pdf                # ) output of scripts/reproduce_figures.jl, committed for convenience
├── fig_eocene_highlatitude_julia.pdf   # )
├── fig_eocene_lowlatitude_julia.pdf    # )
└── data/                               # simulation libraries land here by default (gitignored; see Data below)
```

Scripts are meant to be run from the package root (so `--project=.` resolves correctly),
e.g. `julia --project=. -t auto scripts/reproduce_figures.jl`.

**Deliberately out of scope for this package:** the tooth-seasonality variant
(`popgen_migrate_g_toothseasonality.jl`), and performance optimizations to the
Gillespie inner loop that were identified but not applied (a few redundant/dead
computations — see git history / prior discussion). Both are natural follow-ups but
weren't needed to reproduce the published results.

## Setup

Requires Julia ≥ 1.9. From this directory:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Everything in this package is `Threads`-parallel (library builds, error-surface
computation), so launch Julia with multiple threads to actually benefit:

```sh
julia --project=. -t auto
```

(`-t auto` picks up all available cores; use `-t N` for a specific count.)

## Quick start

```julia
using SharkDentalDynamics

# 1. A scenario is a SiteConfig. Use a preset...
config = modern()

# ...or override specific fields (anything not overridden keeps the preset's value):
config = modern(reps = 10, distance_km = 650.0)

# 2. Build the simulation library (one JLD2 file per (rep, sigtau, tau) grid cell).
build_library!(config, "data/my_run")

# 3. Load your empirical tooth-length measurements into a plain Vector{Float64}.
measures = Float64[...]

# 4. Compare.
result = compare_site(config, "data/my_run", measures)
result.best_j   # (; error, sigtau, tau) -- best fit under the juvenile/nursery hypothesis
result.best_a   # (; error, sigtau, tau) -- best fit under the adult hypothesis
# lower error => that hypothesis is the better-supported explanation for `measures`

# 5. Plot (heatmaps + density-comparison panels, one row per site).
plot_scenario(config, "data/my_run", [("My Site", measures)]; filename = "my_figure.pdf")
```

A lower `best_j.error` than `best_a.error` favors interpreting the empirical sample
as a nursery/juvenile site, and vice versa — see equation 2.1 and section 2c of the
paper for the underlying logic.

## Reproducing the published figures

The three published comparisons (Delaware Bay/modern, Banks Island + Seymour
Island/eocene high-latitude, Red Hot Truck Stop + Whiskey Bridge/eocene low-latitude)
can be regenerated **without rebuilding any simulation library** — they read the
already-built legacy-format libraries one level up in this repo
(`../data/sharks_modern2/`, `../data/sharks_eocene2/`,
`../data/sharks_eocene_lowlatitude/`) directly, via `legacy_simdata_path` (the
internal JLD2 keys are identical between formats; only the filename pattern
differs — `simdata_<rep>_<sigtau_pos>_<tau_pos>.jld` vs this package's own
`...jld2`).

```sh
julia --project=. -t auto scripts/reproduce_figures.jl
```

This reads real empirical measurements from the local `empirical/SandTiger_all.csv`
(a tiny manual column parser is used, so there's no CSV.jl/DataFrames.jl dependency —
see [Data](#data) for why this is a local copy rather than a `../` reference) and
writes `fig_modern_julia.pdf`, `fig_eocene_highlatitude_julia.pdf`, and
`fig_eocene_lowlatitude_julia.pdf` into the package root. Takes well under a minute on
a multi-core machine (`eocene_highlatitude` alone touches ~125,000 saved simulation
files across its 2 sites, so wall time scales with core count).

If you don't have the legacy `../data/sharks_*` folders this reads from (e.g. on a
fresh clone outside the `2018_sharks` monorepo), you'll need to either copy them over
or build fresh libraries with this package's own `build_library!` (dropping the
`path_fn = legacy_simdata_path` argument in `scripts/reproduce_figures.jl` once the
new-format library exists at whatever `data_dir` you point it at) — see
`scripts/example_new_site.jl` for exactly that pattern.

## Running a new empirical-site comparison

This is the main thing you'll want to do for a new fossil or modern locality. Three
steps: **pick/build a config**, **build (or reuse) a simulation library**, **compare**.

`scripts/example_new_site.jl` is a ready-to-edit, fully worked template covering all
three steps end-to-end — it defines a scenario from scratch (a "FREE PARAMETERS"
block at the top), builds its simulation library fresh with `build_library!`, and
compares it against Blackheath (a real high-latitude locality in
`empirical/SandTiger_all_2026.csv`, column 6). Copy it for your own new site and edit
the parameters block. Run with:

```sh
julia --project=. -t auto scripts/example_new_site.jl
```

The steps below explain what that template is doing and why.

### 1. Choose or build a `SiteConfig`

If your new site is a plausible fit for an existing scenario's biology (e.g. another
Eocene low-latitude locality with a similar climate to Red Hot Truck Stop/Whiskey
Bridge), start from that preset and override only what differs:

```julia
config = eocene_lowlatitude(distance_km = 350.0)
```

If it needs genuinely new parameters, construct a `SiteConfig` directly (see
`src/config.jl` for field docs and the exact values behind each preset):

```julia
config = SiteConfig(
    l0 = 55.0, L = 450.0,           # birth / asymptotic body length (cm)
    n0 = 1000, gen = 1,             # initial population size; number of growth cycles to simulate
    mintemp_j = 15.0, maxtemp_j = 26.0,   # nursery/juvenile site seasonal temperature range (°C)
    mintemp_a = 10.0, maxtemp_a = 18.0,   # adult site seasonal temperature range (°C)
    distance_km = 500.0,            # distance between sites (km)
    sigtauvec = collect(0.5:0.5:25.0),    # juvenile dispersal-window grid to sweep (mass-bin steepness)
    tauvec = collect(1.0:1.0:50.0),       # adult dispersal-window grid to sweep (days)
    reps = 10,                      # stochastic replicates per grid cell
)
```

A few things worth knowing when choosing `sigtauvec`/`tauvec`/`reps`:

- **Grid resolution and `reps` trade off directly against runtime.** `eocene_highlatitude`'s
  50×50×25 grid (62,500 simulation runs) is what took ~1 minute to *compare* against
  with 36 threads once already built, but building it from scratch is the expensive
  part — each individual simulation run can itself take anywhere from under a second
  to tens of minutes depending on where in `(sigtau, tau)` space it lands (see the note
  on pathological corners below). Start with a coarse grid (e.g. `modern()`/`eocene_lowlatitude()`'s
  20×20×5) to sanity-check a new site before committing to a fine one.
- **Some `(sigtau, tau)` corners are much more expensive than others.** Large
  `sigtau` + large `tau` (very flexible dispersal in both directions) can let the
  simulated population grow substantially within one `gen` cycle; since the
  underlying Gillespie algorithm's step size shrinks as population grows
  (`dt ~ 1/N`), a growing population needs far more discrete events to simulate the
  same time span — one such corner in `eocene_highlatitude` took nearly an hour and
  10GB of RAM before being killed during development (see
  `scripts/check_legacy_equivalence.jl`'s comments for the full story).
  `build_library!` now guards against this automatically: any single grid cell that
  runs longer than `max_seconds_per_run` (default 120s -- normal cells typically take
  15-35s even at real `gen`) is skipped with a `@warn` (`SimulationTimeout`) rather
  than blocking the whole sweep. A skipped cell just means one fewer rep at that
  particular `(sigtau, tau)` point; if you see the warning a lot, either your grid is
  probing a genuinely expensive region (consider narrowing it) or you want to raise
  `max_seconds_per_run` for a production run where you're willing to wait.
- `l0`/`L` set the tooth-length ceiling (via the `tooth_length` allometry) —
  make sure `L` is large enough that the simulation's maximum possible tooth length
  covers your empirical data's maximum, or the model structurally can't reproduce
  your largest specimens. (This was a real bug in the original codebase — see below.)

### 2. Build the simulation library

```julia
build_library!(config, "data/my_new_site"; seed = 1)
```

This is resumable: rerunning with the same `data_dir` skips any grid cell whose
file already exists (pass `overwrite = true` to force a rebuild). Each cell gets its
own seeded `Xoshiro` RNG derived from `(seed, rep, sigtau_pos, tau_pos)`, so — unlike
the original, never-seeded sweeps — rerunning with the same `seed` reproduces the
library exactly.

### 3. Compare against your empirical measurements and plot

```julia
measures = Float64[...]  # your tooth-length measurements, one per specimen

result = compare_site(config, "data/my_new_site", measures)
plot_scenario(config, "data/my_new_site", [("My New Site", measures)];
    filename = "figures/my_new_site.pdf")
```

`plot_scenario` accepts multiple `(label, measures)` pairs at once if you want
several empirical sites compared against the *same* simulation library in one
figure (as in the published Eocene figures, which each show two sites sharing one
temperature-regime library) — just make sure every site you group together is
scientifically comparable (i.e. built with the same `config`, at the same climate
regime).

## Presets reference

Exact historical parameters, verified against the scripts that produced the
published results:

| param | `modern()` | `eocene_highlatitude()` | `eocene_lowlatitude()` |
|---|---|---|---|
| `l0` (cm) | 90.0 | 50.0 | 50.0 |
| `L` (cm) | 300.0 | 477.0 | 477.0 |
| `n0` | 1000 | 1000 | 1000 |
| `gen` | 1 | 1 | 1 |
| `mintemp_j`/`maxtemp_j` (°C) | 17 / 25 | 12 / 24 | 23 / 30 |
| `mintemp_a`/`maxtemp_a` (°C) | 13 / 23 | 9 / 17 | 23 / 30 |
| `distance_km` | 700 | 400 | 400 |
| `sigtauvec` | `0.5:1:20` (20) | `0.5:0.5:25` (50) | `0.5:1:20` (20) |
| `tauvec` | `1:2:40` (20) | `1:1:50` (50) | `1:2:40` (20) |
| `reps` | 5 | 25 | 5 |
| Real empirical sites | Delaware Bay | Banks Island, Seymour Island | Red Hot Truck Stop, Whiskey Bridge |

Note `eocene_highlatitude` and `eocene_lowlatitude` use the *same* `distance_km` —
only temperature differs between them in the original scripts, despite the name.

## Validation

Two independent layers, both currently passing:

1. **`test/runtests.jl`** — ordinary unit tests with no dependency on the legacy
   code (presets, growth precomputation, shape descriptors, simulation determinism,
   a full tiny-grid `build_library!`/`compare_site` end-to-end run).

   ```sh
   julia --project=. -e 'using Pkg; Pkg.test()'
   ```

2. **`scripts/check_legacy_equivalence.jl`** — seeds the *literal same* default RNG
   before calling both the legacy `popgen_migrate_g` (`include`d fresh from
   `../../src/`, in an isolated module so nothing here touches those files) and this
   package's `simulate_metapopulation`, and asserts every output field is bit-for-bit
   identical, across every grid corner of all three presets, plus a couple of
   full-`gen`-scale realistic-timing checks. Also checks the deterministic
   `compare_to_empirical` pipeline against the legacy `empirical_sim_comparison` on
   synthetic input.

   ```sh
   julia --project=. -t auto scripts/check_legacy_equivalence.jl
   ```

   Read the comments at the top of that file before increasing the grid-corner
   coverage or the realistic-timing section's `gen` — see the pathological-corner
   note above.

## Known issues found and fixed

Building this package as an independent reimplementation surfaced two real,
previously-unnoticed issues in the *original* codebase (not just porting bugs —
these affect the legacy scripts too, and one has since been fixed there as well):

1. **`sharkcompare_modern.jl` reps-averaging bug (now fixed in both places).**
   4 of 7 shape-descriptor error terms (SD, median, Q25, Q75) were averaged across
   the wrong array axis (`dims=2`, the sigtau grid axis) instead of across
   replicates (`dims=1`), then only replicate #1's degenerate slice was kept and
   broadcast across every sigtau row. This meant those 4 terms contributed **no
   sigtau-axis information** to the Delaware Bay error surface/best-fit search —
   only mean/mode/mode-distance did. Using the real saved `../data/sharks_modern2/`
   library, the corrected values are ε_j = 1.35, ε_a = 0.30 (best-fit adult tau = 21
   days), vs. the originally published ε_j = 1.40, ε_a = 0.74 (best-fit adult tau =
   3 days) — the qualitative conclusion (adult site favored) holds either way and is
   actually *stronger* once corrected, but the specific adult dispersal-window
   estimate shifts substantially. **`sharkcompare_eocene.jl`/
   `sharkcompare_eocene_lowlatitude.jl` were never affected** — they use a 4D
   `(num, reps, sigtau, tau)` array layout where `dims=2` is the *correct* reps axis,
   so all four Eocene site comparisons are unaffected and reproduce their published
   values exactly (verified — see `scripts/reproduce_figures.jl`'s output). See
   `error_surface`'s docstring in `src/compare.jl` for the full technical detail.
   `sharkcompare_modern.jl` (one level up in this repo) has since been corrected to
   `dims=1` to match.

2. **`L=295` vs `L=477` (already fixed upstream long before this package).** An
   earlier version of the Eocene simulation library (`data/sharks_eocene/`, git
   history only, superseded ~2021) used an asymptotic body length too small to reach
   the empirical maximum tooth length, meaning the model couldn't structurally
   reproduce the largest observed fossil teeth. `data/sharks_eocene2/` (the current
   library) fixed this; this package's `eocene_highlatitude`/`eocene_lowlatitude`
   presets use the corrected `L=477`.

A few smaller, purely cosmetic/structural things inherited *by design* rather than
fixed (see the relevant docstrings for why): `local_minima`'s endpoint check has a
copy-paste bug inherited from the legacy `findlocalminima.jl` (checks the wrong
inequality, so a trailing local minimum is never actually detected via that branch)
— preserved for exact behavioral equivalence with the published figures rather than
silently changed.

## Data

Two very different kinds of data live near this package, and they're treated
differently on purpose:

- **`data/`** — generated simulation libraries (many JLD2 files per scenario). Not
  committed; excluded via this package's own `.gitignore`. Regenerable from a
  `SiteConfig` via `build_library!` — never hand-edit or commit its contents.
- **`empirical/`** — real tooth-measurement CSVs, committed. These are small,
  human-collected data, not generated output, so they belong in version control;
  scripts read them as **local copies** rather than reaching outside the package
  (`../SandTiger_all.csv` etc.) specifically so this package keeps working if it's
  ever cloned or pushed on its own, separate from the `2018_sharks` monorepo.

`Manifest.toml` *is* committed, so cloning + `Pkg.instantiate()` reproduces the exact
dependency versions this was built and validated against.

Note: `scripts/example_new_site.jl`'s Blackheath comparison uses a deliberately small,
narrow grid to run quickly as a demo (see the script's comments) — its best-fit points
land at the edges of the tested grid, a sign the true optimum lies outside the tested
range. Treat its printed result as a demonstration that the pipeline works end-to-end,
not as a real finding about Blackheath; widen the grid (see the pathological-corners
note above) for an actual analysis.

## Citation

If you use this code, please cite:

> Kim SL, Yeakel JD, Balk MA, Eberle JJ, Zeichner S, Fieman D, Kriwet J. 2022.
> Decoding the dynamics of dental distributions: insights from shark demography and
> dispersal. *Proc. R. Soc. B* 289: 20220808. https://doi.org/10.1098/rspb.2022.0808
