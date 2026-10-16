## Tooling plan

This is a sound Julia-first stack. I would organize it into a **minimal reproducible core** plus optional demonstration components, so the team can deliver even if one advanced piece proves difficult.

### Repository

https://github.com/gaelforget/CMEMS_DTOR_Hackathon_2026

A useful organization could be:

```text
cmems-dtor-2026/
├── README.md
├── Project.toml
├── Manifest.toml
├── CITATION.cff
├── LICENSE
├── src/
│   ├── download_cmems.jl
│   ├── load_glorys.jl
│   ├── sample_profiles.jl
│   ├── simulate_drifters.jl
│   ├── polygon_timeseries.jl
│   └── plotting.jl
├── notebooks/
│   ├── 01_explore_cmems_data.jl
│   ├── 02_argo_glorys_comparison.jl
│   ├── 03_virtual_drifter_experiment.jl
│   └── 04_regional_dashboard.jl
├── scripts/
│   ├── setup.jl
│   ├── fetch_data.jl
│   └── run_demo.jl
├── config/
│   ├── regions.geojson
│   └── experiments.toml
├── data/
│   └── README.md
├── results/
│   └── README.md
└── docs/
    ├── architecture.md
    ├── data-provenance.md
    └── team-workflow.md
```

The repository should contain **code, small configuration files, metadata, and reproducible download instructions**, but not the large Copernicus data products themselves.

## Environment and reproducibility

| Component | Decision |
|---|---|
| Julia version | **Julia 1.12 preferred**, with Julia 1.10 as a fallback if a dependency or cluster environment requires it |
| Packages | Use released versions only; avoid development branches during the hackathon |
| Dependency lock | Commit both `Project.toml` and `Manifest.toml` |
| Cluster use | Develop locally on a small test subset; move bulk downloads, interpolation, ensembles, and long particle integrations to Engaging |
| Secrets | Keep Copernicus Marine credentials out of Git; use environment variables or a local untracked credentials file |
| Data provenance | Record product IDs, download dates, variables, spatial/temporal extent, and preprocessing choices in each experiment configuration |

I would ask everyone to instantiate exactly the project environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()
```

And record the Julia version in `README.md`, for example:

```text
Julia: 1.12.x preferred
Fallback: Julia 1.10.x
```

## Data access and compute

### MyOcean Pro / Copernicus Marine

Use it initially to settle three things before building analysis code:

1. **Which exact GLORYS12 product(s)** provide the needed fields and temporal resolution.
2. **Which variables** the MVP genuinely needs:
   - `uo`, `vo` for advection;
   - temperature and salinity for Argo comparisons;
   - sea-surface or depth-resolved temperature for regional anomaly / marine-heatwave work;
   - optional mixed-layer or sea-level fields if they sharpen a regional story.
3. **The preferred automated download method** compatible with a reproducible Julia or shell workflow.

A very useful early deliverable is a small `fetch_data.jl` or shell wrapper that takes a dataset/product identifier, dates, variables, and bounding box from `config/experiments.toml`.

### Engaging cluster

The cluster is particularly appropriate for:

- staging sizeable GLORYS12 subsets;
- testing multiple regions or periods;
- running an ensemble of simulated drifters / release locations;
- collocating GLORYS12 with many Argo profiles;
- producing final results without requiring every team member to download the same large files.

Avoid making the demo dependent on live cluster access. A compact precomputed result bundle or a small regional NetCDF subset should allow the notebook/demo to run on a laptop.

## Julia package stack

### Core data and workflow packages

| Package | Main role |
|---|---|
| `NCDatasets.jl` | Read and write NetCDF data; inspect GLORYS12 fields and dimensions |
| `Climatology.jl` | Work with gridded ocean/climate data products and support field-oriented workflows |
| `MeshArrays.jl` | Grid-aware spatial operations, particularly useful if it fits the selected grid and interpolation needs |
| `GeoJSON.jl` | Store and load target regions, deployment domains, and analysis polygons |
| `OceanRobots.jl` | Access, organize, and analyze ocean-observing-platform observations |
| `Drifters.jl` | Advect virtual particles/drifters using GLORYS12 velocity fields |
| `STAC.jl` | SpatioTemporal Asset Catalog (STAC) julia client |


### Visual and interactive layer

| Package | Use |
|---|---|
| `Makie.jl` | Maps, trajectories, profile comparisons, animations, and publication/demo-quality figures |
| `Pluto.jl` | Preferred candidate for an interactive and reproducible hackathon demonstrator |
| `IJulia.jl` | Useful if contributors prefer Jupyter or need compatibility with an existing notebook workflow |

I would select **one primary presentation environment**, rather than maintaining the same demo in both Pluto and Jupyter. For this project, Pluto seems especially well matched: a user can change a region, polygon, launch date, or release strategy and immediately see the consequence.

### Event / extremes extension

| Package | Potential role |
|---|---|
| `ExtremeTracker.jl` | Define, identify, or track ocean extremes—particularly relevant if the marine-heatwave collaborators join and the final use case includes event-aware observing strategies |

This should initially remain an **optional module**, rather than a dependency of the first advection and Argo–GLORYS workflows. That protects the core MVP while making the project naturally expandable into a marine-heatwave demonstrator.

## Claude Code preparation

Using Claude Code could be worthwhile for accelerating routine implementation, but I would set it up **before** the hackathon and establish a few guardrails:

- use it to scaffold code, tests, documentation, configuration parsing, and plotting boilerplate;
- have it inspect package APIs and produce small isolated examples before it touches the main workflow;
- require generated code to run in the pinned project environment;
- never expose Copernicus credentials, private cluster information, or non-public data;
- keep commits small and reviewable;
- test interpolation, time handling, longitude conventions, depth coordinates, and units manually—these are likely to be the scientifically fragile areas.

## Suggested staged implementation

### Stage 0 — before the hackathon

- Create repository, pinned Julia environment, and basic team conventions.
- Confirm all required packages install under Julia 1.12; document a Julia 1.10 fallback only if needed.
- Identify a GLORYS12 product and produce one small regional/time subset.
- Verify access to MyOcean Pro and Engaging.
- Make a first map: GLORYS12 field + polygon + real observation locations.

### Stage 1 — core data products

- Load GLORYS12 variables from NetCDF.
- Extract a time series averaged or otherwise diagnosed within a GeoJSON polygon.
- Collocate GLORYS12 fields at Argo locations, times, and depths.
- Quantify basic differences: bias, RMSE, correlation, distributional summaries, and profile-level diagnostics.

### Stage 2 — digital-twin demonstration

- Seed virtual particles/drifters at a real or hypothetical observing-system location.
- Advect them with GLORYS12 currents.
- Plot real versus virtual trajectory/sampling coverage.
- Compare whether each sampling strategy intercepts the selected target signal or event.

### Stage 3 — polished demonstrator

- Wrap a single regional case in a Pluto notebook.
- Let a user choose a region, time window, release/deployment strategy, and polygon.
- Generate a concise final set of outputs:
  1. environmental context map;
  2. observed and virtual sampling paths;
  3. polygon-based time series or event index;
  4. quantitative coverage/representativeness comparison.

The key is to keep the initial repository executable with a small packaged example, while allowing GLORYS12-scale processing and more ambitious virtual experiments to run on Engaging.

