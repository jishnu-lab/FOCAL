# Building the GRN inputs for each data type

Focal's temporal module reads a dictys `dynamic.h5`; the state-specific module reads per-state GRN
tables. The input route (matched multiome, unmatched RNA + ATAC paired by optimal transport, or
RNA only) only changes how each window's TF x target prior (`binlinking.tsv.gz`) is made; the
windows, the network reconstruction and everything Focal does afterwards are shared. Run every
step as a batch or array job (one task per window), not on a login node.

## Shared steps

1. **Expression**: raw gene x cell counts per window as `expression.tsv.gz`
   (`focal.cli.expression_to_tsv`; multiome peak features containing `:` are dropped).
2. **Trajectory**: STREAM (or any method giving a graph plus a per-cell branch and position) ->
   `dictys.traj.trajectory.fromdist(edges, dist)` -> `traj_node.h5`;
   `dictys.traj.point.fromdist(traj, branch, dist)` -> `traj_cell_rna.h5` with the cells in
   **sorted name order** (`tofile` reads it positionally), plus `coord_rna.tsv.gz`.
3. **Windows**: `dictys.dynamic.subsets_rna(..., ncell, noverlap, dmax)` -> `tmp_dynamic/Subset<k>/`
   with `names_rna.txt`, plus `subsets.txt`, `subset_locs.h5`, `subset_edges.tsv.gz`. Typical
   `ncell` ~1,000 cells, `noverlap` 0.75-0.8 of it, `dmax` 0.004 x the maximum pseudotime.
   The windowing uses a dense cell x cell distance matrix: ~25 min and ~11 GB on one core for
   ~23k cells.
4. **Prior per window**: `tmp_dynamic/Subset<k>/binlinking.tsv.gz` (binary TF x target mask),
   made by one of the routes below.
5. **Network**: `dictys.network.reconstruct(fi_exp, fi_mask=binlinking, ..., device="cuda")`
   (`python -m focal.cli.reconstruct_networks <start> <end> <WORK_DIR>` wraps it), then per window
   `normalize` -> `net_nweight`, `indirect` -> `net_iweight`, `normalize(net_iweight)` ->
   `net_inweight`, then `dictys.network.tofile(data, tmp_dynamic, subsets.txt, dynamic.h5,
   dynamic=True, nettype="n")`. `tofile` needs both `net_nweight` and `net_inweight`, or the
   total-effect network `w_in` (used for forces and waves) is all zero. The torch build must
   include the GPU's architecture (see SKILL.md, Setup).

`focal.cli.validate_inputs` checks the inputs before the long jobs. A window whose mask keeps too
few targets makes `reconstruct` stop with "Not enough target genes"; check the mask shapes.

## Route A: matched multiome (RNA + ATAC in the same cells)

dictys `config.mk`: `JOINT=1`. Each window's ATAC cells are its RNA cells. Per window:

| step | dictys call |
|---|---|
| reads + peaks | `dictys.chromatin.macs2` (window BAM -> `reads.bam`, `peaks.bed`) |
| footprints | `dictys.chromatin.wellington` -> `wellington.tsv.gz` |
| motifs | `dictys.chromatin.homer` (motifs via `focal.cli.motifs_to_homer`) -> `homer.tsv.gz` |
| TF -> target prior | `binding` -> `tssdist` -> `linking` -> `binlinking` |

The prior changes window to window, and the per-window `binding.tsv.gz` tables feed
`SmoothedCurvesChromatin` (TF binding-score and OCR-count curves) and `BindingPhases`. Reference
batch scripts: `multiome_dynamic_regulation/bash_scripts/` (`preproc/`, `chromatin/`, `network/`).

## Route B: unmatched scRNA + scATAC, paired by optimal transport

For separate scRNA-seq and scATAC-seq of the same system (same samples or conditions). No
footprints are needed: the window's accessibility comes from its RNA cells' ATAC partners.

1. **Pairing** (moscot env, GPU). Per sample, a GW `TranslationProblem` from ATAC LSI (without
   the depth-tracking component; 30 dims) to RNA PCA (30 dims), `epsilon` 1e-3 (2e-3 on
   non-convergence), `scale_cost="mean"`. With state labels in both modalities, fused GW with a
   0/1 linear cost over allowed (ATAC state, RNA state) pairs keeps the mass within compatible
   states. Outputs per sample: `T_<sample>.npy` (ATAC x RNA), the ATAC and RNA cell names.
2. **Coupling matrix**. For each RNA cell, its column of `T` divided by its sum, top-k (k = 50)
   entries kept and renormalised -> `coupling.npz` (RNA cells x ATAC cells, rows sum to 1), with
   `atac_peaks.npz` (ATAC cells x peaks raw counts) beside it. Inferred accessibility of RNA cell
   `i` = `coupling[i] @ atac_peaks` (computed per window, never stored per cell).
3. **Peak -> gene**. Peaks within +- 250 kb of a gene's TSS, weight `exp(-d / 50 kb)`
   (`P2G.npz`, peaks x genes) from the species GTF, by gene symbol.
4. **Motif scan** (once, GPU). FIMO of every peak against a species-matched pfm database
   (HOCOMOCO v11 human / mouse, gimmemotifs format, with its `motif2factors` table), p < 1e-4 ->
   `peak_by_motif_scores.npz`, `peak_names.csv.gz`, `motif_ids.csv.gz`, `motif2tfs.json`.
   A CPU fallback (gimmemotifs scanner, FPR 0.01) gives the same format.
5. **Mask per window**. Window accessibility = `mean(coupling[window cells]) @ atac_peaks`, keep
   the 30,000 most accessible peaks that carry a motif hit; write `wellington.tsv.gz`
   (region x `<TF>_<motif>` = accessibility; regions as `chr:start:end`) and `homer.tsv.gz`
   (same shape, motif score), keep TF columns whose TF is expressed in the window, then
   `binding(cuth=0, cutw=0)` -> `tssdist` (TSS bed of the species) -> `linking` ->
   `binlinking(n=20)`.

The prior changes along pseudotime because the window's ATAC partners change. Gene and TF names
must match between the expression, the motif table and the TSS bed (mouse symbols are mixed case).

## Route C: RNA only (base GRN)

1. Take an ATAC-derived base GRN in CellOracle format (peak_id, gene_short_name, one column per
   TF; e.g. CellOracle's mouse scATAC atlas base GRN).
2. Collapse over peaks (`groupby('gene_short_name').max()`), binarise (> 0), transpose to
   TF x target, keep TFs that are also targets, set self-edges False, match names to the dataset.
3. Per window, restrict the mask to the window's expressed genes -> `binlinking.tsv.gz`.

The prior is identical in every window (only the fitted weights move along pseudotime) and there
are no binding curves; forces, waves, phases, episodic GRNs and enrichment work as usual.

## State-specific GRNs

Built per state with CellOracle from a base GRN (from scATAC peaks or CellOracle's promoter base
GRN), loaded with `focal.state_specific._grn.load_celloracle_grn_tables` and fused with
`combine_cluster_grn_links` (per-state `source, target, coef_mean, coef_abs` tables;
`strength == 1` for `coef_abs` at or above the 0.9 quantile; the paper Methods keep each state's
top 10,000 edges and the larger weight for shared edges).

To do: give the state-specific module routes A and B as well, i.e. a state's prior from the
state's own (or OT-paired) accessibility through the motif scan + dictys linking chain, instead
of a reference base GRN.
