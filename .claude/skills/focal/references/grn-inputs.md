# Building the GRN inputs: matched vs unmatched ATAC

Focal's temporal module reads a dictys `dynamic.h5`; the state-specific module reads CellOracle
GRN tables. This file covers producing them for each data type. Reference sbatch scripts live in
`multiome_dynamic_regulation/bash_scripts/` (`preproc/`, `chromatin/`, `network/`); their
`WORK_DIR`s are the original cluster paths — change them before use. Run every step as a batch
or array job (one task per window), not on a login node.

## Common steps (both data types)

1. **Expression**: gene x cell counts as `expression.tsv.gz` (`focal.cli.expression_to_tsv`;
   multiome peak features containing `:` are dropped). For multiome, imputation (MAGIC via
   Palantir) is used for trajectory inference only.
2. **Trajectory**: STREAM (or any method giving a graph + per-cell branch/position) ->
   `dictys.traj.trajectory.fromdist(...)` -> `traj_node.h5`; `dictys.traj.point.fromdist(...)` ->
   `traj_cell_rna.h5` (cells sorted by name), plus `coord_rna.tsv.gz`.
3. **Windows**: `dictys.dynamic.subsets_rna(...)` slides windows of `ncell` cells with
   `noverlap` overlap along every edge -> `tmp_dynamic/Subset<k>/` with `names_rna.txt`,
   `subsets.txt`, `subset_locs.h5`, `subset_edges.tsv.gz` (`preproc/2_traj_subset.sbatch`).
   Then slice expression per window and QC it (`preproc/3_slice_gene_by_cell.sbatch`,
   `4_qc_sliced_matrix.sbatch` -> `dictys.preproc.qc_reads`).
4. **Prior per window**: `tmp_dynamic/Subset<k>/binlinking.tsv.gz` (binary TF x target mask);
   how it is made is the matched/unmatched difference below.
5. **Network**: `python -m focal.cli.reconstruct_networks <start> <end> <WORK_DIR>` (wraps
   `dictys.network.reconstruct` with `fi_mask=binlinking`; `network/1_network_reconstruct.sbatch`),
   then `network/2_network_postprocess.sbatch` (normalize, `net_iweight`, `net_inweight`), then
   `dictys.network.tofile` (`network/3_network_to_file.sbatch`) -> `output/dynamic.h5`.
   `tofile` needs BOTH `net_nweight` and `net_inweight`, or the `w_in` network Focal uses for
   forces and waves is missing.

`focal.cli.validate_inputs` checks the inputs before the long jobs.

## Matched (multiome: RNA + ATAC in the same cells)

dictys `config.mk`: `JOINT=1`. Each window's ATAC cells are its RNA cells (`names_atac0.txt`
copied from `names_rna.txt`, then `preproc selects_atac`). Per window, in order:

| step | dictys call | script |
|---|---|---|
| window reads + peaks | `dictys.chromatin.macs2` (window BAM -> `reads.bam`, `peaks.bed`; dictys keeps at most 500k peaks by score) | `preproc/1_sort_bam.sbatch`, `chromatin/1_window_bam_peaks.sbatch` |
| footprints | `dictys.chromatin.wellington` -> `wellington.tsv.gz` | `chromatin/2_chromatin_footprinting.sbatch` |
| motifs | `dictys.chromatin.homer` (motifs via `focal.cli.motifs_to_homer`) -> `homer.tsv.gz` | `chromatin/3_chromatin_scan_motif.sbatch` |
| TF -> target prior | `binding` -> `tssdist` -> `linking` -> `binlinking` | `chromatin/4_chromatin_linking.sbatch` |

The prior changes window to window. The per-window `binding.tsv.gz` tables are what
`SmoothedCurvesChromatin` reads for TF binding-score / OCR-count curves (and `BindingPhases`).

## Unmatched, ATAC cells available (separate scRNA and scATAC of the same system)

dictys `config.mk`: `JOINT=0`. Provide ATAC cells on the same trajectory: `traj_cell_atac.h5`
and `coord_atac.tsv.gz` (e.g. from an RNA/ATAC co-embedding or label transfer). dictys then
runs `dynamic subsets_atac` to assign ATAC cells to each RNA window (`names_atac0.txt`) and
the chromatin chain above unchanged. Binding curves are available. Check that every window gets
enough ATAC cells; sparse windows give empty peaks/footprints.

## Unmatched, base GRN only (RNA-only cells; ATAC from a reference)

What the T-cell Perturb-seq system used (see `FOCAL_notebooks/temporal/dynamic_grn/t_cell_dictys_data.ipynb`):

1. Take an ATAC-derived base GRN in CellOracle format (peak_id, gene_short_name, one column per
   TF; e.g. CellOracle's murine scATAC base GRN).
2. Collapse over peaks (`groupby('gene_short_name').max()`), binarize (> 0), transpose to
   TF x target, keep only TFs that also appear as targets, set self-edges False, save as
   `data/<species>_base_grn.tsv.gz` (tab-separated, first column = TF, header = targets).
3. Per window, restrict the mask to that window's expressed genes and write
   `tmp_dynamic/Subset<k>/binlinking.tsv.gz` (`chromatin/subset_base_grn.sbatch`, array job).
4. Continue with the common network step. Skip all `chromatin/1-4` steps.

Consequences: the prior is identical in every window (only the fitted weights move along
pseudotime); there are no per-window binding tables, so no `SmoothedCurvesChromatin` /
`BindingPhases`. Episodic GRNs, forces, waves, phases, validation and enrichment work as usual.
`dictys.network.reconstruct` raises "Not enough target genes for specified principal components"
when a window's mask keeps too few targets — check the mask shape per window.

## State-specific GRNs (CellOracle)

Build one GRN per state with CellOracle from a base GRN (from matched or unmatched scATAC peaks,
or CellOracle's prebuilt promoter base GRN), then load them with
`focal.state_specific._grn.load_celloracle_grn_tables` and fuse with `combine_cluster_grn_links`
(concatenates the per-state `source, target, coef_mean, coef_abs` tables and marks
`strength == 1` for `coef_abs` at or above the 0.9 quantile). The paper Methods additionally
subset each state to its top 10,000 edges and keep the larger weight for shared edges — do that
before calling it if you need the same fusion.
