# Tests for `focal.state_specific`

Coverage of the SLIDE–GRN enrichment path that `StateSpecificEnrichment.run_enrichment`
drives: `src/focal/state_specific/_enrichment.py`, `_slide.py` and `_grn.py`.
Plotting is deliberately out of scope — these tests cover the maths and the
control flow.

## Running

```bash
conda activate fftemporal
pytest tests/state_specific -q            # ~13 s
pytest tests/state_specific -q -rx        # also list the xfail reasons
```

`ISSUES_demo.ipynb` runs on the same `fftemporal` kernel. `03_lf_enrich_backend.ipynb`
needs CellOracle, so it runs on `celloracle_env` instead.

No celloracle needed. `load_celloracle_grn_tables` is the only function in the
module that imports it, and it is not exercised here; the tests set
`GRN_links_after_fit` / `GRN_TFs` / `slide_features` directly, which are exactly
the three attributes `run_enrichment` guards on.

## Layout

| file | subject |
| --- | --- |
| `conftest.py` | hand-built GRN link tables, flat edge tables, and `cc_dict` record helpers |
| `test_hypergeom.py` | `hypergeom_slide_grn_score` — agreement with scipy, guards, and four classes of unvalidated input |
| `test_grn.py` | `combine_cluster_grn_links`, `grn_edges_from_combined_links`, `filter_network_scores_for_clusters` |
| `test_slide_enrichment.py` | `get_slide_grn_enrichment` — combination/condition enumeration, the co-regulation requirement, parallel-edge dedup |
| `test_enrichment_table.py` | `build_enrichment_table` / `write_enrichment_table` — the order-1 vs order-2 asymmetry |
| `test_run_enrichment.py` | `StateSpecificEnrichment` end to end: fusions, file outputs, post-processing, `load_results` round trip |
| `ISSUES.md` | the eight problems, each with the test that reproduces it |
| `ISSUES_demo.ipynb` | the same eight problems run live — each one's description next to the code breaking it. Executed, outputs saved. Needs no CellOracle object; runs in seconds |
| `03_lf_enrich_backend.ipynb` | the B-cell latent-factor enrichment driven through `StateSpecificEnrichment` instead of notebook-local functions |

## Conventions

Inherited from `tests/temporal`:

* A test named `test_currently_...` pins behaviour that is wrong but load-bearing;
  the matching `xfail(strict=True)` test asserts what the behaviour should be.
* `xfail(strict=True)` means a fix turns the test into an `XPASS` failure. That
  is the prompt to delete the `currently` test and drop the marker.
* Every expected value is either derived on paper (and stated in a comment) or
  computed by an independent implementation (`scipy.stats.hypergeom`,
  `numpy.quantile`) — never by calling the function under test.

## What the suite found

Eight problems, written up in `ISSUES.md`. The two that matter most:

1. **Asking for "both TFs strongly regulate this gene" always returns nothing.**
   When scoring TF pairs, the strength patterns `(1, 1)` and `(0, 0)` can never
   match, because a check that should only count distinct TFs also demands
   distinct strength values. `(1, 1)` is in the default pattern list for order-2
   runs, so the most obvious question you would ask of a TF pair silently
   returns an empty answer.
2. **One thin cluster kills the whole run.** If a cluster has fewer usable TFs
   than the order you asked for, `run_enrichment` stops with
   `ValueError: Columns must be same length as key` and everything computed so
   far is lost. The message points at pandas rather than the real cause.

## Not covered

| area | why not |
| --- | --- |
| all plotting functions | out of scope by request; they are ~45% of `_enrichment.py` |
| `load_celloracle_grn_tables` | needs celloracle and a real `.celloracle.oracle` file |
| `load_slide_features_from_files` / `_combined` | file-format parsing; worth a follow-up suite with fixture TSVs |
| `StateSpecificManager` (`manager.py`) | thin wrapper, untested |

The SLIDE feature loaders are the most valuable next target: they parse
experiment TSVs and populate `positive_corr_names` / `negative_corr_names`,
which drive the colour proportions in the published figures.
