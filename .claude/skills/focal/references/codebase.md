# Focal codebase reference

## Install / import

```bash
conda activate <env with dictys>            # dictys is not on PyPI; install it via conda first
git clone https://github.com/jishnu-lab/FOCAL.git && cd FOCAL && pip install -e .
python -c "import focal, focal.temporal; print(focal.__file__)"   # must be .../FOCAL/src/focal/__init__.py
```
No `sys.path` bootstrap is needed. In sbatch scripts, call the env's python explicitly
(e.g. a `$DICTYS_PY` variable) rather than relying on `conda activate` in a batch shell.

## Entry points

```python
import dictys.net, focal
net = dictys.net.dynamic_network.from_file("dynamic.h5")
mgr = focal.TemporalManager(net, trajectory_range=(0, 1), num_points=40, dist=0.001,
                            sparsity=0.01, output_dir="episodes/")
mgr.build_episode(1, slice(0, 5), percentile=98, pval_threshold=1e-3, n_processes=16, network_type="w")
mgr.build_all_episodes(total_episodes=8, points_per_episode=5, write=True)   # episode_<i>.parquet
mgr.enrich_episode(1, slice(0, 5), lf_genes)                                   # TF, p_value, enrichment_score, genes_in_lf, genes_dwnstrm, weights
mgr.enrich_all_episodes(lf_genes, total_episodes=8, points_per_episode=5, write=True)  # enrichment_episode_<i>.csv
mgr.force_source() / mgr.waves() / mgr.phases() / mgr.validate() / mgr.curves(mode) / mgr.time_scales()
focal.temporal.TrajectorySegments  # several branches (e.g. PB (1,2), GC (1,3)) sharing settings
focal.temporal.manager.run_episodic_enrichment(...) / run_episodic_construction(...)  # fresh-process runners
```

Lower level, what the manager composes (used directly by the notebooks and drivers):
`EpisodeDynamics(dictys_dynamic_object=net, output_folder, trajectory_range, num_points, dist,
sparsity, network_type)` with the fixed chain
`compute_expression_curves()` -> `build_episode_grn(time_slice)` -> `filter_edges(n_processes,
chunk_size, pval_threshold)` -> `compute_tf_expression()` -> `calculate_forces(n_processes,
chunk_size)` -> `select_top_edges(percentile)` (or
`select_top_activating_and_repressing_edges(98.5, 0.5)`) -> `set_lf_genes(genes)` ->
`annotate_lf_in_grn()` -> `calculate_enrichment()`. Attributes: `episodic_grn_edges`
(MultiIndex (TF, target), `avg_force`, `is_in_lf`).

ORA primitive: `focal.base.enrichment.calculate_tf_episodic_enrichment(df, total_lf_genes,
total_genes_in_grn)` = `hypergeom.sf(k-1, N, K, n)`, ES = k / (n K / N); shared by temporal and
state-specific modules. `focal.ora` is the exported alias.

State-specific: `focal.StateSpecificManager().enrich(...)`, `state_specific._enrichment.
StateSpecificEnrichment` (CellOracle cluster GRNs + SLIDE LFs, TF-combination ORA via
`_slide.get_slide_grn_enrichment`), figures `plot_tf_enrichment_bars`, `plot_episode_from_csvs`.

Figures (all in `focal.temporal`): `plot_tf_episodic_enrichment_dotplot(dfs, episode_labels,
p_value_threshold, min_significance_threshold, min_targets_in_lf=2, min_targets_dwnstrm=2,
cmap_name, log_scale, tf_order, horizontal_layout, show_plot)` -> `(fig, plot_data, plotted_tfs)`;
`plot_tf_target_episodic_heatmap(csv_path, ...)` expects columns `TF, target, <one per episode>`;
`plot_tf_gene_coregulation_heatmap`, `sort_tfs_by_gene_similarity`, `plot_force_landscape`,
`plot_force_heatmap(_by_phase / _with_clustering)`, `plot_phase_binding_boxes`,
`plot_state_composition_bars/curves`, `plot_force_validation_*`.

## dictys objects Focal reads

`dynamic_network.from_file(dynamic.h5)`: `net.nids[0]` TF indices, `net.nname`, `net.ndict`
(gene -> idx), `net.traj.edges` (trajectory graph; a linear trajectory is a single edge
`[[0, 1]]`), `net.point`. `trajectory_range` = node ids of the branch; `num_points` sampled
pseudotime points; `dist` = Gaussian sigma in trajectory units (B cell 0.001, T cell 0.002;
roughly one window spacing for a new dataset); `sparsity` 0.01. `network_type`: `"w"` raw direct
effect (episodic GRNs, regulation curves), `"w_in"` normalized total effect (force waves, phases,
validation). Building `dynamic.h5` yourself needs BOTH `net_nweight` and `net_inweight`
(= normalize of `net_iweight`) passed to dictys `tofile`, and `traj_cell_rna.h5` cells sorted by
name.

## Notebooks (FOCAL_notebooks)

`temporal/analysis/`: `Fig3_2_episodic_enrichment.ipynb` (PB/GC), `Fig5_episodic_enrichment.ipynb`
(KO programs on GC, 8 episodes), `Fig4_1_LF_global_dynamics`, `Fig4_2_dynamic_validation`,
`Fig4_3_phase_clustered_links`, `Fig4_4_tf_binding_phases`, `Ext_chromatin_dynamics`,
`Ext_tpex_ex_forces`, and the driver `run_episodic_enrichment_bcell.py` (the reference
implementation of the episodic pipeline — copy it to start a new dataset).
`temporal/dynamic_grn/`, `temporal/trajectory/`, `state_specific/01..04`, `cross_prediction/`.
`NOTEBOOK_MIGRATION.md` has the OLD (firefate) -> NEW import cells; `temporal/SETUP.md` the env setup.
A new analysis notebook must be added to `temporal/index.md` or the strict docs build fails.
