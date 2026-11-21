# Episodic GRNs and episodic enrichment — how to run and read them

Reference implementation: `FOCAL_notebooks/temporal/analysis/run_episodic_enrichment_bcell.py`.
Code: `focal.temporal._episodes.EpisodeDynamics`.

## Definitions

- **Episode** = contiguous block of `points_per_episode` sampled pseudotime points out of
  `num_points` on one branch (`TemporalManager._episode_slice(e, q) = slice((e-1)q, eq)`,
  1-based). Convention: 5 points per episode; a branch (0,x) at 20 points = 4 episodes, a full
  range (1,x) at 40 points = 8 episodes. Episodes are NOT the link-peak clusters / waves from
  the temporal force-wave analysis — those are separate (`TFForceWaves`, phases).
- **Episodic edge**: present in the episode (row sum != 0), regulator not `ZNF*`/`ZBTB*`,
  >= 3 non-zero points, one-sample t-test p < 0.05 then p < `pval_threshold` (1e-3), sign never
  flips (direction invariance).
- **Force** F = sign(beta) * exp(log10|beta| + log10 E_TF) with the TF expression taken from the
  episode's own points; episode score = mean force over the episode.
- **Episodic GRN** = edges with |mean force| >= `percentile` (98) — or signed tails via
  `select_top_activating_and_repressing_edges(98.5, 0.5)`.
- **Episodic GRN on disk**: per episode, a DataFrame with MultiIndex `(TF, Target)` and one
  column `avg_force`, holding only the thresholded (top-percentile, invariance-filtered) edges;
  companion `avg_force_episode_{i}.parquet` holds the unthresholded avg_force of all filtered
  edges. Written by `run_episodic_construction` as `episode_{i}.parquet` +
  `avg_force_episode_{i}.parquet`.
- **Enrichment** per TF: hypergeometric over-representation of program genes among the TF's
  episodic targets: `p = hypergeom.sf(k-1, N, K, n)`, `ES = k / (nK/N)`; N = targets in the
  episodic GRN, K = program genes that are active targets, n = TF out-degree, k = overlap.
  Output columns `TF, p_value, enrichment_score, genes_in_lf, genes_dwnstrm, weights`.
  No FDR correction anywhere (state this in Methods); display cutoffs raw p 0.05 (0.08 relaxed),
  >= 2 program targets and >= 2 downstream targets per TF.

## Per-episode chain (exactly what `TemporalManager.enrich_episode` does)

```python
epi = EpisodeDynamics(dictys_dynamic_object=net, output_folder=out, trajectory_range=rng,
                      num_points=P, dist=sigma, sparsity=0.01, network_type="w")
lcpm, dtime = epi.compute_expression_curves()          # once per branch
for ep in range(1, P // q + 1):
    sl = slice((ep - 1) * q, ep * q)
    beta = epi.build_episode_grn(time_slice=sl)        # first call smooths the whole network
    filtered = epi.filter_edges(n_processes=16, chunk_size=8000, pval_threshold=1e-3)
    epi.compute_tf_expression(); epi.calculate_forces(n_processes=16, chunk_size=30000)
    edges = epi.select_top_edges(98)                   # MultiIndex (TF, target), avg_force
    edges.to_parquet(f"episodic_grn_edges/episode_{ep}.parquet")
    epi.set_lf_genes(program); epi.annotate_lf_in_grn()
    epi.calculate_enrichment().to_csv(f"enrichment_episode_{ep}.csv", index=False)
```
The expensive prefix is program-independent: enrich several programs from one run. Program =
`names` column of a SLIDE `feature_list_Z*.txt` (the B cell runs drop `HLA-` genes); check
`g in net.ndict` and log the missing ones. Run it as a batch job (16 processes, generous memory).

## Figures

```python
from focal.temporal import plot_tf_episodic_enrichment_dotplot, plot_tf_target_episodic_heatmap
fig, plot_data, tfs = plot_tf_episodic_enrichment_dotplot(dfs=[pd.read_csv(...) per episode],
    episode_labels=["Ep1", ...], figsize=(6, 3.5), p_value_threshold=0.05,
    min_significance_threshold=0.05, min_dot_size=10, max_dot_size=200, cmap_name="Purples",
    tf_order=None, log_scale=True, horizontal_layout=False)
# TF x program-target avg_force by episode: build a wide csv (TF, target, Ep1..EpN) from the
# parquets yourself (get_episodic_grn_subset expects .pkl) then plot_tf_target_episodic_heatmap(csv)
```
TF order in the dotplot = peak episode then peak ES (a display rule, not a statistic);
`sort_tfs_by_gene_similarity` (Jaccard of program-target sets, Ward) is the alternative.

## Parameters used in the notebooks

| dataset | branch `trajectory_range` | P | q | dist | percentile | p | network |
|---|---|---|---|---|---|---|---|
| B cell PB / GC (Fig3_2) | (1,2) / (1,3) | 20 | 5 | 0.001 | 98 | 1e-3 | w |
| B cell GC, KO programs (Fig5) | (1,3) | 40 | 5 | 0.001 | 98 | 1e-3 | w |
| T cell | (0,1) | 24 | 4 | 0.002 | 98 | 1e-3 | w |

## Known caveats

Fixed 2026-08-31: force x wrong TF expression; episode-1 expression reused for every episode.
Still open: `build_transition_window(links=...)` raises; parquet/pkl mismatch in
`get_episodic_grn_subset`; `get_beta_curves` returns the TF x target cross product (re-filter with
`.loc[links]`); constant edges get p = 0 from the t-test; `episodic_composition()` is a stub;
chromatin smoothing propagates NaNs; `classify_tf_global_activity` crashes on zero-variance input.
