# Problems found by the `state_specific` tests

Eight problems, each reproduced by a test in this folder. Nothing in the library
was changed.

Each one is pinned by two tests:

* `test_currently_...` — records what the code does **today**, so nobody changes
  it by accident.
* `test_should_...` — asserts what the code **should** to do. It is marked
  `xfail`, so it fails today on purpose. Once someone fixes the bug, that test
  starts passing, pytest reports it as `XPASS`, and that is the signal to delete
  the `currently` test.

Ranked by how much damage each one can do.

---

## 1. Asking for "both TFs strongly regulate this gene" always returns nothing

**Affects:** any run with `order_fr_tfcomb` of 2 or more. Runs with order 1 are fine.

### What happens

When you score a *pair* of TFs, you also pick a strength pattern: should TF1's
edge be strong or weak, and TF2's? There are four patterns: `(0,0)`, `(0,1)`,
`(1,0)`, `(1,1)`.

Two of them — `(1,1)` and `(0,0)` — never return a single gene, no matter what
the data says.

Build a network where TFA and TFB both strongly regulate G1, G2, G3 and G4. Ask
for pattern `(1,1)`. The right answer is all four genes. You get none.

### Why

A gene is kept only if this is true:

```python
filtered.groupby("target")[["source", "strength"]].nunique().eq(len(filter_df)).all(axis=1)
```

`len(filter_df)` is 2 here, because there are two TFs. The check then demands
**two different `source` values** — correct, the gene must be hit by both TFs —
and **two different `strength` values**.

That second demand is the bug. Under pattern `(1,1)` every surviving row has
strength 1, so there is only *one* distinct strength value, not two. The gene is
thrown away.

The strength check looks unnecessary: the merge just above already forced each TF
to the strength the pattern asked for.

### Why it matters

`build_enrichment_table` uses `[(1,1), (0,1), (1,0)]` as its default for order-2
runs. So the default list includes a pattern that can never return anything — and
it is the pattern most people would care about, two TFs jointly and strongly
driving the same gene.

Order 1 is unaffected: there is one TF, `len(filter_df)` is 1, one strength value,
and the check passes.

**Tests:** `test_slide_enrichment.py::TestHomogeneousStrengthConditions`

---

## 2. One thin cluster kills the whole run

**Affects:** any run where some cluster has fewer TFs than the order you asked for.

### What happens

`run_enrichment` stops with:

```
ValueError: Columns must be same length as key
```

Everything already computed for earlier cluster fusions is lost, because nothing
catches the error.

### Why

For each cluster fusion the code builds a list of TF combinations:

```python
itertools.combinations(sorted(total_source), ord_tf)
```


If that cluster only has one candidate TF but you asked for pairs, the list is
empty. No results get recorded, so the results list stays `[]`. The table builder
then tries to split an empty list into `score` and `p_value` columns, and pandas
refuses.

The usable-TF list is `(SLIDE genes in the graph + their regulators) ∩ GRN_TFs`,
which runs dry easily — a sparse cluster, a short `GRN_TFs` list, or just asking
for `order_fr_tfcomb = [1, 2]` when one cluster only supports order 1.

(A gene is a candidate when it's a known TF (listed in GRN_TFs) and either:
- it directly regulates a SLIDE gene in that cluster's thresholded network, or
- it is itself a SLIDE gene appearing in that network)
### Why it matters

This is a crash, not a wrong number, and the message points at pandas rather than
at the real cause, so it is hard to diagnose. A sensible outcome would be to skip
that fusion and carry on.

**Tests:** `test_run_enrichment.py::TestEmptyCombinationCrash`

---

## 3. The downstream-gene list is wrapped in an extra list when nothing overlaps

**Affects:** order-2 output files. Order-1 files are unaffected.

### What happens

When a TF has downstream genes but none of them are SLIDE genes, the code records
the gene list one level too deep:

```
recorded:  [['D0', 'D1', 'D2', 'D3', 'D4']]      <- a list holding one list
expected:  ['D0', 'D1', 'D2', 'D3', 'D4']
```

So `len(dwnstrm)` reports **1** when the real answer is **5**.

### Why

Every other branch appends a plain list. This one wraps it:

```python
cc_dict[...].append([tf_comb, condition, (0, 1), ([], [dwn_list]), case, (...)])
#                                                      ^^^^^^^^^^
```

### Why it matters

At order 1 these rows get dropped anyway (their score is 0 and the filter needs a
score above 0), so you never see it.

At order 2 failing rows are **kept** — the code zeroes their score instead of
removing them — so the malformed list is written into `enriched_df_*.csv`.
Anything reading that column back gets a nested list, and the "more than 2
downstream genes" filter is counting the wrong thing.

**Tests:** `test_slide_enrichment.py::TestNoOverlapBranch`

---

## 4. Order-2 tables are sorted by scores they no longer contain

**Affects:** the row order of order-2 output files. The numbers are right.

### What happens

Rows are sorted by score, and *then* the failing rows have their scores set to 0.
So a row showing `score = 0.0` can sit above a row showing `score = 2.0`.

### Why

In `_records_to_enrichment_df` the `sort_values(by="score", ...)` call runs before
the mask that zeroes the failing rows.

### Why it matters

Only readability — but "sorted by score" is the natural assumption when reading
one of these files, and it is not true.

**Tests:** `test_enrichment_table.py::TestOrderTwoZeroesInsteadOfDropping::test_currently_sorted_by_pre_zeroing_score`

---

## 5. Setting `slide_starting_genes` too low silently empties the results

**Affects:** any run where a TF has more downstream genes than the SLIDE universe size.

### What happens

`slide_starting_genes` is the size of the gene universe for the hypergeometric
test. If a TF combination has more downstream genes than that number, the test is
being asked something impossible — more successes in the population than the
population holds.

scipy returns `nan` for the p-value. The score still comes back as a normal-looking
number.

### Why it matters

Nothing checks for it, and nothing warns. The row then meets `p_value < 0.05`,
which is `False` for `nan`, so it is quietly dropped (order 1) or zeroed
(order 2).

The visible symptom is an empty or near-empty results table with no error — which
reads like "no enrichment found" rather than "the universe size is wrong".

Worth a check at the top of the run: `slide_starting_genes` should be at least as
large as the biggest downstream-gene count in the network.

**Tests:** `test_hypergeom.py::TestInconsistentPopulation`

---

## 6. Impossible overlap counts are accepted without complaint

**Affects:** nothing today; it is a missing guard rather than an active bug.

### What happens

Two inputs that cannot occur are accepted and produce confident results —
`p_value = 0.0` and a large positive score:

* more overlapping genes than there are SLIDE genes (`s > S`)
* more overlapping genes than there are downstream genes (`s > p`)

The function does guard the other degenerate cases (`s <= 0`, `S <= 0`, `P <= 0`),
so the gap is inconsistent rather than deliberate.

### Why it matters

Only reachable if something upstream is already wrong. But if that ever happens,
the result looks like the strongest possible enrichment instead of an error.

**Tests:** `test_hypergeom.py::TestImpossibleSuccessCounts`

---

## 7. Score and p-value contradict each other when there are no downstream genes

**Affects:** nothing today; reachable only through bad input.

### What happens

With zero downstream genes but a positive overlap:

* the score is `0.0` — "no enrichment"
* the p-value is `0.0` — "as significant as possible"

### Why

The score has a divide-by-zero guard that returns 0. The p-value has no matching
guard, so scipy computes it anyway.

**Tests:** `test_hypergeom.py::TestZeroPopulationSuccesses`

---

## 8. Cluster names have to be numbers

**Affects:** `filter_network_scores_for_clusters` only. Not reached by `run_enrichment`.

### What happens

```
ValueError: invalid literal for int() with base 10: 'ActB_1'
```

### Why

The function converts each cluster id with `int(cid)`. Everywhere else accepts
any label: `clusters_of_interest` is typed `list[Any]`, and
`combine_cluster_grn_links` works with whatever keys the oracle has.

### Why it matters

The B cell oracle happens to use `'0'`–`'11'`, which survive `int()`. A named
clustering — `'ActB_1'`, `'GC_1'` — would not.

**Tests:** `test_grn.py::TestFilterNetworkScoresForClusters`

---

## Checked and working

These were tested and behave correctly. The tests stay so they cannot break
quietly.

| behaviour | test |
| --- | --- |
| Edge-strength cutoff is taken across all fused clusters together, not per cluster, and an edge exactly on the cutoff counts as strong | `test_grn.py::TestCombineClusterGrnLinks` |
| The same TF→gene edge coming from two clusters is kept as two edges, not merged | `test_grn.py::TestGrnEdgesFromCombinedLinks` |
| Deduplication keeps the edge with the largest *magnitude*, so a strong negative edge beats a weaker positive one | `test_slide_enrichment.py::TestParallelEdgeDeduplication` |
| A TF with no edges contributes nothing, instead of letting its partner's genes through | `test_slide_enrichment.py::TestCommonTargetRequirement` |
| Order 1 removes failing rows; order 2 keeps them with the score zeroed | `test_enrichment_table.py` |
| Score and p-value match `scipy.stats.hypergeom` exactly | `test_hypergeom.py::TestAgreesWithScipy` |
| SLIDE genes absent from the network still count toward the universe size, lowering scores | `test_slide_enrichment.py::TestSingleTFScoring` |
