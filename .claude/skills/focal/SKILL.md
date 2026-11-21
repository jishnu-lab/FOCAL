---
name: focal
description: >
  Focal (Factor-Outcome Coupling for Assessment of Linkages; formerly FIREFate) codebase and
  analysis assistant. Use whenever the user mentions Focal, FOCAL, GIFs (GRNs coupled with
  Interpretable Factor-outcomes), cellular programs (CPs), SLIDE latent factors as programs,
  episodic GRNs / episodic enrichment, transition-window GRNs, TF force curves, force waves,
  regulatory phases, EpisodeDynamics, TemporalManager, SmoothedCurvesGRN, TFForceWaves,
  StateSpecificManager, in-silico TF KO on CPs, fate-bias stratification, the dictys
  dynamic.h5 -> Focal chain, matched or unmatched scRNA + scATAC inputs (dictys JOINT=0/1,
  base-GRN binlinking masks), github.com/jishnu-lab/FOCAL, or running Focal on a new dataset
  (driver scripts, sbatch jobs, parameters, figures).
---

# Focal — codebase, API and how to run it

**Focal** = *Factor-Outcome Coupling for Assessment of Linkages* (formerly FIREFate). Focal
couples mechanistic GRNs (dictys sliding-window networks, CellOracle state GRNs) with
interpretable ML (SLIDE latent factors = cellular programs, CPs) to focus dense networks onto the
regulatory components governing fate decisions; outputs are **GIFs** (GRNs coupled with
Interpretable Factor-outcomes).

## Where things live

| Item | Value |
|---|---|
| Backend repo | github.com/jishnu-lab/FOCAL (`src/focal/`), branch `main` |
| Companion notebooks | github.com/jishnu-lab/FOCAL_notebooks (git submodule at `docs/notebooks`) |
| Docs | focal-gifs.readthedocs.io (Sphinx, `docs/`; notebooks rendered from stored outputs, never executed) |
| Package | `focal` 0.1.0, Python >= 3.10, deps numpy/pandas/scipy/matplotlib/plotly/tqdm/joblib/networkx/h5py/pyyaml; extras `[scanpy]`, `[celloracle]`, `[docs]`, `[dev]` |
| Environment | a conda env with dictys (not on PyPI) and Focal installed editable (`pip install -e .`); check `import focal; focal.__file__` points into your clone |

Dataset paths for the notebooks come from `FOCAL_notebooks/temporal/datasets.yaml` via
`focal.io.DatasetPaths.find()` / `.from_yaml()` — edit the yaml, never the notebooks.

On an HPC cluster, run Focal through the scheduler (sbatch / srun), never on a login node.

## Package layout (`src/focal/`)

```
focal/
  __init__.py          # exports TemporalManager, StateSpecificManager, EpisodeDynamics,
                       #   SmoothedCurvesGRN, SmoothedCurvesChromatin, AlignTimeScales,
                       #   TFForceWaves, TFForceValidation, StateFrequency, BaseManager, ora
  base/                # BaseManager (results dict, output_dir, save/load), enrichment.py = ORA primitive
  temporal/            # Temporal module
    manager.py         #   TemporalManager, TrajectorySegments, run_episodic_enrichment/_construction
    _episodes.py       #   EpisodeDynamics + episodic figures (dotplot, TF-target heatmap, co-regulation)
    _curves.py         #   SmoothedCurvesGRN (Gaussian smoothing of dictys windows) + expression/regulation heatmaps
    _forces.py         #   calculate_force_curves(_parallel), filter_edges_by_significance_and_direction
    _source.py         #   TFForceSource (shared beta / TF-expression / force cache per branch)
    _waves.py          #   TFForceWaves + force landscapes / heatmaps / clustering
    _phases.py         #   RegulatoryPhases, ForceWavePhases, BindingPhases, order_links(_by_phase)
    _reductions.py     #   mean_force, abs_max_force, softmax_peak, get_max_points, aggregate_max_points
    _align.py          #   AlignTimeScales, plot_main_trajectory_nodes
    _chromatin.py      #   SmoothedCurvesChromatin (TF binding score / OCR count per window)
    _states.py         #   StateFrequency + state composition bars/curves
    _validation.py     #   TFForceValidation (enriched links vs size-matched random null)
  state_specific/      # State-specific module
    manager.py         #   StateSpecificManager.enrich(...)
    _enrichment.py     #   StateSpecificEnrichment, plot_tf_enrichment_bars, plot_episode_from_csvs
    _slide.py          #   hypergeom_slide_grn_score, get_slide_grn_enrichment (TF-combo ORA), build/write_enrichment_table
    _grn.py            #   combine_cluster_grn_links, load_celloracle_grn_tables, ...
  cross_prediction/    # Cross-prediction module (scaffold; R/py code in Cross_prediction/)
  backends/dictys/     # read_h5_file, reconstruct_subset/_range, stats lcpm_tf, ChromatinGRNStat
  io/                  # DatasetPaths (yaml), load_lf_gene_colors, load_tf_target_links, load_episode_enrichment, qc_reads
  utils/               # genes (get_tf_indices, check_if_gene_in_ndict), states, curves, parallel
  cli/                 # reconstruct_networks, validate_inputs, expression_to_tsv, motifs_to_homer
```

There is **no plotting namespace**: every figure lives with its subject (`from focal.temporal
import *` gives classes and plots). Managers record everything they build in `.results`.

## The three modules

Focal works on single-cell RNA and ATAC data (sc/snRNA-seq, scATAC-seq), **matched or
unmatched**, and is three modules in one package:

| Module | Product | What it does | Code |
|---|---|---|---|
| **State-specific** | Prioritized regulatory subnetworks | SLIDE discovers sparse cellular programs (CPs) separating two contrasting states; CPs are embedded in state-specific and fused cross-state GRNs (CellOracle) to surface TFs whose regulons are enriched for each program (TF-combination ORA), plus in-silico perturbation of the enriched TFs | `focal.state_specific`, `StateSpecificManager` |
| **Temporal** | Phase-resolved dynamic regulation | Reads a dictys dynamic network along a pseudotime trajectory: smoothed expression / regulation / beta curves, TF forces, force waves and regulatory phases from pseudotemporal clustering, episodic GRNs (direction-invariant forces in the top-percentile tails) and CP enrichment against them, validation against size-matched random links, chromatin binding curves | `focal.temporal`, `TemporalManager` |
| **Cross-prediction** | Fate predisposition by transfer learning | CPs learned from fate-switching perturbations (TF-KO Perturb-seq) vs controls score uncommitted cells in unperturbed data by predicted fate bias | `focal.cross_prediction` (scaffold), `Cross_prediction/` (R + notebook) |

## Matched vs unmatched ATAC (decide this first)

The data type decides how the upstream GRNs are built; everything Focal does downstream of
`dynamic.h5` / the CellOracle GRN tables is the same. Full steps in `references/grn-inputs.md`.

- **Matched** (multiome, scRNA + scATAC in the same cells; e.g. the B-cell system): dictys runs
  with `JOINT=1`; each pseudotime window has its own ATAC reads -> macs2 peaks -> wellington
  footprints -> homer motifs -> `binding`/`tssdist`/`linking`/`binlinking`, so the TF -> target
  prior changes along pseudotime. Per-window footprints also give **chromatin binding curves**
  (`SmoothedCurvesChromatin`, `BindingPhases`).
- **Unmatched, ATAC cells available** (scRNA and scATAC from different cells of the same
  system): dictys `JOINT=0`. ATAC cells are placed on the same trajectory (`traj_cell_atac.h5`,
  `coord_atac.tsv.gz`) and `dictys.dynamic.subsets_atac` assigns them to each RNA window; the
  per-window chromatin chain then runs as in the matched case, so binding curves still exist.
- **Unmatched, base GRN only** (RNA-only cells, ATAC from a reference; what the T-cell Perturb-seq
  system used): windows are cut from RNA only (`subsets_rna`); one fixed binary TF x target mask
  from an ATAC-derived base GRN (CellOracle format, collapsed over peaks) is subset to each
  window's expressed genes and written as that window's `binlinking.tsv.gz`, then
  `dictys.network.reconstruct` runs as usual. The prior is static, only the expression-fitted
  weights change along pseudotime, and there are **no chromatin binding curves**.
- State-specific GRNs (CellOracle) take either a base GRN from matched/unmatched scATAC peaks or
  CellOracle's prebuilt promoter base GRN; the T-cell system used the murine ATAC base GRN.

## Routing

| Request | Read |
|---|---|
| Run/extend episodic enrichment, episodes, forces, parameters, caveats | `references/episodic-pipeline.md` |
| API details, managers, file formats, notebooks, install | `references/codebase.md` |
| Building `dynamic.h5` / state GRNs from matched or unmatched RNA + ATAC (dictys, CellOracle, sbatch steps) | `references/grn-inputs.md` |
| Force / ORA / invariance math | `references/math-operators.md` |

## Style

Editable text in figures (`pdf.fonttype = 42`, `svg.fonttype = 'none'`). In prose say
"top-percentile tails", "nominate TFs", "downstream targets", "contrasting" (not
"discriminative"), "modules" (not "API products").

## Gotchas learned the hard way

- `TemporalManager.build_transition_window(links=...)` raises (Series vs DataFrame) — the
  transition-window figures use the notebooks' hand path (`get_beta_curves` -> `calculate_force_curves`).
- `run_episodic_construction` writes `episode_<i>.parquet`; `get_episodic_grn_subset` globs
  `episode_<i>.pkl` — build the edges-by-episode table yourself from the parquets.
- Episodic tables made before 2026-08-31 are wrong (force x wrong TF expression; episode-1
  expression reused) — regenerate, never reuse `legacy_*` folders.
- No FDR anywhere; notebooks use raw-p cutoffs (0.003 / 0.01 / 0.05 / 0.08) and >= 2 program
  targets and >= 2 downstream targets per TF for display.
- Regulator families `ZNF*` / `ZBTB*` are dropped inside `build_episode_grn` (hard-coded).
- dictys gene names must match the program's names exactly (e.g. upper-case human symbols).
