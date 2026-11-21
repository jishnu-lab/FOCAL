---
name: focal
description: >
  Run Focal (Factor-Outcome Coupling for Assessment of Linkages) on single-cell RNA + ATAC data:
  build per-window dynamic GRNs along a pseudotime trajectory with dictys and turn them into TF
  forces, attractor states and episodic GRNs. Use whenever the user mentions Focal / FOCAL,
  dynamic or per-window GRNs, TF forces / force curves, attractor states along a trajectory,
  episodic GRNs or episodic enrichment, the dictys dynamic.h5 -> Focal chain, unmatched scRNA +
  scATAC paired by optimal transport (moscot GW pairing, coupling-weighted peaks, peak-to-gene),
  matched multiome inputs, github.com/jishnu-lab/FOCAL, or running Focal on a new dataset
  (driver scripts, batch jobs, parameters, figures).
---

# Focal: dynamic GRNs and TF forces from single-cell RNA + ATAC

Focal reads a dictys dynamic network (one GRN per sliding window of cells along a pseudotime
trajectory) and turns it into **TF -> target forces** along each fate branch, the **attractor
states** of the branch's cell-state composition, and **episodic GRNs** (forces that stay in the
top-percentile tails over a stretch of pseudotime), which can be enriched for gene programs.
This skill is the recipe for running it end to end on a new dataset; the CD8 T cell run below is
the worked example.

Code: github.com/jishnu-lab/FOCAL (`src/focal/`); notebooks: github.com/jishnu-lab/FOCAL_notebooks
(submodule at `docs/notebooks`); docs: focal-gifs.readthedocs.io.

## Setup

- One conda env with **dictys** (not on PyPI) and Focal installed editable
  (`pip install -e .`); `import focal; focal.__file__` must point into your clone.
- **The torch build must include your GPU's architecture.** dictys reconstructs each window with
  pyro on the GPU. The conda-forge `pytorch 2.3.1 cuda118` build that a dictys env typically pins
  stops at `sm_89`: on Hopper GPUs (H100 / H200 / GH200, `sm_90`) it sits at 0 % GPU utilisation
  on one CPU core. Use the pip wheel instead (`pip install torch==2.4.1 --index-url
  https://download.pytorch.org/whl/cu121`, NumPy-2 safe) and check that
  `torch.cuda.get_arch_list()` contains your `sm_XX` before submitting window jobs.
- Optimal-transport pairing needs an env with moscot + GPU JAX; run the solves on a large GPU
  (an H200 solves 5-7k x 3-6k cell GW problems in minutes; on CPU they do not finish).
- Run everything as batch jobs (one array task per window where possible), never on a login node.

## 1. Choose the input route (decide this first)

| Data | Route | Per-window TF -> target prior |
|---|---|---|
| **Matched multiome** (RNA + ATAC in the same cells) | dictys `JOINT=1`: per-window reads -> peaks -> footprints -> motifs | changes along pseudotime; binding curves available |
| **Unmatched scRNA + scATAC** of the same system | **OT pairing** (section 2): each RNA cell gets ATAC partners from a GW coupling; window accessibility = coupling-weighted peak counts; motif scan + dictys linking | changes along pseudotime |
| **RNA only** | a fixed ATAC-derived base GRN (e.g. CellOracle's), subset to each window's expressed genes | static; only the fitted weights move |

Details of every route and the file formats: `references/grn-inputs.md`.

## 2. Unmatched RNA + ATAC: the OT-pairing route

**a. Pair ATAC cells to RNA cells (moscot GW, per sample).** For every sample (e.g. condition x
timepoint) solve a Gromov-Wasserstein `TranslationProblem` from ATAC (`obsm['X_lsi']`, dropping
the LSI component that tracks sequencing depth) to RNA (`obsm['X_pca']`), so mass never crosses
samples:

```python
from moscot.problems.cross_modality import TranslationProblem
tp = TranslationProblem(adata_src=atac_s, adata_tgt=rna_s).prepare(
    src_attr={"attr": "obsm", "key": "X_lsi"}, tgt_attr={"attr": "obsm", "key": "X_pca"})
tp = tp.solve(epsilon=1e-3, alpha=1.0, scale_cost="mean")       # retry at epsilon 2e-3 if not converged
(_, prob), = tp.problems.items()
T = prob.solution.transport_matrix                                # ATAC cells x RNA cells
```

When both modalities carry state labels, use **fused GW** (`alpha < 1`) with a linear cost of 0
for allowed (ATAC state, RNA state) pairs and 1 otherwise, so cells pair only within compatible
states (the "constrained" coupling). Save `T` with its source and target cell names per sample.

**b. Coupling -> per-cell accessibility.** Normalise each RNA cell's column of `T` to a
distribution over ATAC cells and keep its top-k partners (k = 50, renormalised): a sparse
RNA x ATAC matrix `C`. A cell's inferred accessibility is `C[i] @ A`, with `A` the ATAC cell x peak
counts; do not materialise the full product (it is nearly dense). Report the mass kept by top-k
per sample as QC.

**c. Peak -> gene.** Distance-decay links from each peak to the genes whose TSS lies within
+- 250 kb, weight `exp(-d / 50 kb)` (a peaks x genes matrix). Inside the per-window chain, dictys
`tssdist` + `linking` do the peak -> target linking.

**d. Motif scan (once).** Score every peak against a species-matched motif database (e.g.
HOCOMOCO v11) with FIMO on the GPU (tangermeme / torch FIMO, p < 1e-4); output a sparse
peak x motif matrix and a motif -> TF-symbol table. About 1.5 min for 134k peaks x 356 motifs on
one GH200.

**e. Windows.** STREAM (or any tree with a per-cell branch and position) -> dictys
`trajectory.fromdist` -> `traj_node.h5`, `point.fromdist` -> `traj_cell_rna.h5` (cells in sorted
name order), then `dictys.dynamic.subsets_rna(ncell, noverlap, dmax)` -> `tmp_dynamic/Subset<k>/`.
It works on a dense cell x cell distance matrix, so ~23k cells take ~25 min and ~11 GB on one
core: give the job the time.

**f. Per-window mask.** For each window: raw counts of its cells (`expression.tsv.gz`, genes
detected in >= 10 cells); accessibility = (mean of `C` over the window's cells) @ `A`, top 30,000
peaks; write it as dictys' `wellington.tsv.gz` (region x TF_motif, value = accessibility) next to
`homer.tsv.gz` (motif scores) and run the native chain `binding -> tssdist -> linking ->
binlinking(n=20)` -> `binlinking.tsv.gz` (the TF x target mask).

**g. Per-window network.** `dictys.network.reconstruct` on the GPU, then `normalize` ->
`net_nweight`, `indirect` -> `net_iweight`, `normalize(net_iweight)` -> `net_inweight`, then
`dictys.network.tofile(..., dynamic=True)` -> `output/dynamic.h5`. `tofile` needs both
`net_nweight` and `net_inweight`, or the total-effect network `w_in` is all zero. Windows are
independent: spread them over every GPU available; a lock file per window lets several jobs
share one list without computing a window twice.

## 3. Forces and attractor states per fate branch

For each branch (root -> leaf of the tree), order the windows on the path by distance from the
root and take the state composition of each window's cells. Forces use Focal's definition

    F_link(w) = sign(beta) * exp(log10|beta| + log10 TF_expr),   TF_expr = log2(CPM + 1) in the window

on the total-effect network (`net_iweight`), over the top 20,000 links by max |F|. Smooth along
window position (Gaussian, sigma 1.5 windows) and compare

- the **total force** (sum of |F| over the links), drawn black, with
- the **composition speed** v = 0.5 * sum_states |dp/dx|, drawn red; its local minima are the
  **attractor states** (quasi-stationary composition).

Score each branch: Spearman of total force, |d total force / dx| and link-peak density against v;
distance from each v minimum to the nearest force minimum and peak; the dominant state at each
attractor; whether the last attractor is the branch's end state; agreement with an expected
state order. Colour states as in the dataset's own UMAP.

Episodic GRNs and program enrichment on the same `dynamic.h5`: `references/episodic-pipeline.md`.

## Worked example: CD8 T cells in acute vs chronic LCMV

Unmatched longitudinal scRNA-seq and scATAC-seq of CD8 T cells after LCMV Armstrong (acute,
memory) and Clone 13 (chronic, exhaustion) infection, days 8, 15 and 30.

| Step | Setting / size |
|---|---|
| Trajectory | joint non-naive Arm + Cl13 STREAM2 tree, 23,146 cells, root in the day-8 pool; branches root -> **Mem** and root -> **Exh-Term** |
| Pairing | state-constrained fused GW per sample (6 samples); top-50 partners keep a median 0.79-1.00 of each cell's mass; 35,121 ATAC cells, 133,795 peaks |
| Motifs | HOCOMOCO v11 mouse, 356 motifs, 5.4 M hits (torch FIMO on a GPU) |
| Windows | 116 windows of 1,042 cells (overlap 782) |
| Masks | ~190 TFs x ~10,000 targets per window, density ~10 % |
| Networks | reconstruct spread over GH200, H200, A100 and A40 GPUs |
| Read-out | per branch: composition bars, black total force, red composition speed and its minima; a branch metrics table |

## 4. To do: state-specific GRNs from the same inputs

The state-specific module builds one GRN per state with CellOracle from a base GRN. Implement the
same matched / unmatched handling there: for unmatched data, a state's accessibility = (mean of
`C` over the state's RNA cells) @ `A`, then the motif scan + dictys linking chain (or CellOracle
on the state's accessible peaks) gives a state-specific prior from the data's own ATAC rather
than a reference base GRN.

## References

| Need | Read |
|---|---|
| Building `dynamic.h5` for each input route, file formats | `references/grn-inputs.md` |
| Episodic GRNs, episodic enrichment, parameters | `references/episodic-pipeline.md` |
| API, managers, dictys objects, notebooks, install | `references/codebase.md` |
| Force / ORA / invariance math | `references/math-operators.md` |

Figures: keep text editable (`pdf.fonttype = 42`, `svg.fonttype = 'none'`).
