# CPR-BEAMS Workshop — MATLAB → R

R translation of the analysis workflow taught at the **[Continuous Plankton Recorder (CPR) Workshop 2025](https://website.whoi.edu/cpr-beams/cpr-workshop/workshop2025/)** (CPR-BEAMS), for exploring and analyzing CPR zooplankton and phytoplankton data.

The original workflow was written in MATLAB by **Pierre Hélaouët** (Marine Biological Association). This repository ports those scripts to R so the same steps — loading, checking, spatial trimming, time-series construction, taxon aggregation, and gap-filling in time and space — can be run in an open, reproducible R environment.

---

## Background

The CPR program tows a device that collects plankton on silk filters from ships of opportunity across the North Atlantic, Pacific, and other waters, producing long-term, high-taxonomic-resolution records of plankton abundance and ecosystem change.

This code comes out of the **CPR Workshop 2025 (CPR-BEAMS)**, a 3-day in-person event held **November 4–6, 2025** at **Woods Hole Oceanographic Institution** (Woods Hole, MA). The workshop is a research + training component of the NSF-supported **CPR-BEAMS** project and is co-sponsored by the **NOAA Cooperative Institute for the North Atlantic Region (CINAR)**.

Useful links:
- CPR-BEAMS project: https://website.whoi.edu/cpr-beams/
- Data DOI (MBA): https://doi.mba.ac.uk/data/3567

---

## Data

The analysis starts from a MATLAB data extract (`CPR_Data_CPRBeam.mat`) or the CSV extract it was built from. Step 2 loads the CSVs and writes an R-native `CPR_Data_CPRBeam.RData` that all later steps read.

The extract contains **three plankton groups** — **large zooplankton**, **small zooplankton**, and **phytoplankton** — each split into three object types:

| Prefix   | Contents |
|----------|----------|
| `Data_*` | Sample metadata + taxa abundances. First **8 columns** are metadata: `Sample`, `Latitude`, `Longitude`, `Year`, `Month`, `Day`, `Hour`, `Minute`. Remaining columns are taxa (column names are numeric CPR taxon IDs). |
| `List_*` | Taxon list. Includes the `DRI` — the year each taxon began being routinely counted / identified — plus flags such as `Copepods`. |
| `Taxo_*` | Taxon list with full taxonomic hierarchy (e.g. `Family`). |

Notes:
- **Zeros in `Data_*` are true zeros** (real absences), not missing values.
- Sample IDs encode the tow route (3 digits + 2 letters).
- Taxon columns are keyed by numeric ID; join to `List_*` / `Taxo_*` to recover names.

Expected local layout:

```
raw/CPRBeam_DataExtract/   # the *.csv extract (input to Step 2)
CPR_Data_CPRBeam.RData     # written by Step 2, read by Steps 3–9
output/                    # optional exported CSVs / figures
```

---

## Repository structure

Scripts are numbered to run in order. Each R script names the MATLAB file it translates in its header.

| R script | Translates | Purpose |
|----------|-----------|---------|
| `CPRBeam_WS_Step2_Load.R`          | `Prog_CPRBeam_WS_Step2_Load.m`          | Read all CSVs, auto-name each table, save `CPR_Data_CPRBeam.RData`. |
| `CPRBeam_WS_Step3_Check.R`         | `Prog_CPRBeam_WS_Step3_Check.m`         | Sanity checks; crude sampling map; heatmap of sampling effort (Year × Month). |
| `CPRBeam_WS_Step4_TrimToArea.R`    | `Prog_CPRBeam_WS_Step4_TrimToArea.m`    | Subset to an area: (1) rectangle, (2) Google Earth `.kml` polygon, (3) ICES shapefile areas. |
| `CPRBeam_WS_Step5_Timeseries.R`    | `Prog_CPRBeam_WS_Step5_Timeseries.m`    | Build a single-taxon (e.g. *Calanus finmarchicus*) monthly time series; force full Year × Month grid. |
| `CPRBeam_WS_Step6_AggregatedTaxa.R`| `Prog_CPRBeam_WS_Step6_AggregatedTaxa.m`| Aggregate groups (e.g. all large copepods, all Calanidae) and build their time series. |
| `CPRBeam_WS_Step7_RegulInTime.R`   | `Prog_CPRBeam_WS_Step7_RegulInTime.m`   | Fill gaps in time: linear, pchip, makima (Akima), and the "Colebrook" method. |
| `CPRBeam_WS_Step8_RegulInSpace.R`  | `Prog_CPRBeam_WS_Step8_RegulInSpace.m`  | Regularize in space: simple grid binning and IDW interpolation. |
| `CPRBeam_WS_Step9_NightDay.R`      | `Prog_CPRBeam_WS_Step9_NightDay.m`      | Solar elevation per sample; split day vs. night abundances. |

**Helper functions** (MATLAB originals in repo; R equivalents noted):

| MATLAB function | Role | R equivalent |
|-----------------|------|--------------|
| `f_CPRBeam_kml2struct.m` | Import `.kml` polygon | Replaced with `sf::st_read()` in Step 4. |
| `f_CPRBeam_pos2dist.m`   | Great-circle distance between lat/lon pairs (used by Step 8 IDW) | Ported to `f_CPRBeam_pos2dist.R`. |
| `f_CPRBeam_SolarAzEl.m`  | Solar azimuth / elevation from UTC + position (used by Step 9) | Ported to `f_CPRBeam_SolarAzEl.R`. (`oce::sunAngle()` is an equivalent alternative.) |
---

## Running the workflow

1. Place the CSV extract in `raw/CPRBeam_DataExtract/`.
2. Run **Step 2** once to generate `CPR_Data_CPRBeam.RData`.
3. Run Steps 3–9 in order (each is standalone after Step 2 — it reloads the `.RData`).

Open the project via its `.Rproj` file so the working directory is set to the project root; all paths in the scripts are relative to it.

---

## Dependencies

Core:

```r
install.packages(c("tidyverse", "sf", "sp"))
```

By step (beyond `dplyr` / `ggplot2` / `tidyr`, which are used throughout):

- **Step 3** — `sf`, `rnaturalearth`, `rnaturalearthdata`, `rnaturalearthhires`, `ggrastr` *(optional higher-quality map alternatives)*
- **Step 4** — `sf`, `sp`, `purrr`, `viridis`, `RColorBrewer`
- **Step 5–6** — `viridis`, `scales`
- **Step 7** — `signal` (pchip), `akima` (makima-style interpolation)
- **Step 9** — `patchwork` (compose the heatmap + marginal bars); `oce` *(optional, alternative solar-elevation source)*

`rnaturalearthhires` is on GitHub:

```r
# install.packages("remotes")
remotes::install_github("ropensci/rnaturalearthhires")
```

---

## Translation status

**Done**
- Steps 2–7 fully translated.
- Step 4: all three area-selection methods (rectangle, `.kml` polygon, ICES shapefiles).
- Step 7: linear, pchip, makima, and Colebrook interpolation.
- Step 8: simple spatial binning, IDW interpolation, and all comparison maps (raw vs. regularized, sampling effort). `f_CPRBeam_pos2dist` ported to R.
- Step 9: per-sample solar elevation and day vs. night abundance split. `f_CPRBeam_SolarAzEl` ported to R.

**To do**
- Use `here::here()` for paths instead of relative paths

---

## Attribution

- **Original MATLAB workflow:** Pierre Hélaouët (Marine Biological Association), CPR-BEAMS Workshop 2025.
- **R translation:** Alexandra Cabanelas (MIT–WHOI Joint Program), November 2025.

The CPR-BEAMS project is supported by the U.S. National Science Foundation and co-sponsored by NOAA CINAR.

---

## Contact

Alexandra Cabanelas — MIT–WHOI Joint Program