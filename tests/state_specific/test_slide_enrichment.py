"""Combinatorial core -- :func:`focal.state_specific._slide.get_slide_grn_enrichment`.

This is where ``run_enrichment`` does its real work: for every TF combination and
every binary edge-strength condition, find the targets co-regulated by the whole
combination under that condition, deduplicate parallel edges, and score the
overlap with the SLIDE gene set.

Expected overlaps are worked out by hand from the fixtures in ``conftest.py``;
expected scores come from ``scipy.stats.hypergeom`` directly.

``xfail(strict=True)`` marks a defect; see ``ISSUES.md``.
"""

import itertools
import math

import pandas as pd
import pytest
from scipy.stats import hypergeom

from focal.state_specific._slide import get_slide_grn_enrichment

from conftest import SLIDE_UNIVERSE

CF = ("c0",)


def run(edges, ord_tf, slide, sources, universe=SLIDE_UNIVERSE, case="slide"):
    """Call the function under test and return its record list."""
    cc: dict = {}
    get_slide_grn_enrichment(
        edges, cc, CF, ord_tf, slide, universe, sources, case, show_progress=False
    )
    return cc[CF][ord_tf]


def by_condition(records):
    """``{condition: record}`` for a single-TF-combination run."""
    return {tuple(r[1]): r for r in records}


def expected(universe, dwngene, n_slide, common):
    """Independent hypergeometric reference -- scipy, not the library."""
    score = math.log2((common / n_slide) / (dwngene / universe))
    pval = float(1.0 - hypergeom.cdf(common - 1, universe, dwngene, n_slide))
    return score, pval


class TestRecordLayout:
    def test_six_slot_record_structure(self, simple_edges, slide_features):
        rec = run(simple_edges, 1, slide_features, {"TFA"})[0]
        tf_comb, condition, es, genes, case, dfs = rec
        assert isinstance(tf_comb, tuple) and isinstance(condition, tuple)
        assert isinstance(es, tuple) and len(es) == 2
        assert isinstance(genes, tuple) and len(genes) == 2
        assert case == "slide"
        assert isinstance(dfs, tuple) and len(dfs) == 2
        assert all(isinstance(d, pd.DataFrame) for d in dfs)

    def test_case_label_is_passed_through(self, simple_edges, slide_features):
        for label in ("slide", "net", "rnd"):
            rec = run(simple_edges, 1, slide_features, {"TFA"}, case=label)[0]
            assert rec[4] == label

    def test_enumerates_every_combination_and_condition(self, simple_edges, slide_features):
        # 2 TFs choose 1 = 2 combinations, 2^1 = 2 conditions -> 4 records
        assert len(run(simple_edges, 1, slide_features, {"TFA", "TFB"})) == 4
        # 2 choose 2 = 1 combination, 2^2 = 4 conditions -> 4 records
        assert len(run(simple_edges, 2, slide_features, {"TFA", "TFB"})) == 4

    def test_tf_combinations_are_sorted(self, simple_edges, slide_features):
        recs = run(simple_edges, 2, slide_features, {"TFB", "TFA"})
        assert all(r[0] == ("TFA", "TFB") for r in recs)

    def test_conditions_are_the_full_binary_product(self, simple_edges, slide_features):
        recs = run(simple_edges, 2, slide_features, {"TFA", "TFB"})
        assert {tuple(r[1]) for r in recs} == set(itertools.product([0, 1], repeat=2))

    def test_cc_dict_is_filled_in_place_and_returned(self, simple_edges, slide_features):
        cc: dict = {}
        out = get_slide_grn_enrichment(
            simple_edges, cc, CF, 1, slide_features, SLIDE_UNIVERSE,
            {"TFA"}, "slide", show_progress=False,
        )
        assert out is cc and CF in cc and 1 in cc[CF]


class TestSingleTFScoring:
    """ord_tf == 1: hand-checkable overlaps on ``simple_edges``."""

    def test_strong_condition_selects_only_strength_one_edges(
        self, simple_edges, slide_features
    ):
        # TFA strong targets: G1, G2.  SLIDE = {G2, G3, G9} -> overlap {G2}.
        rec = by_condition(run(simple_edges, 1, slide_features, {"TFA"}))[(1,)]
        assert sorted(rec[3][1]) == ["G1", "G2"]
        assert sorted(rec[3][0]) == ["G2"]

    def test_weak_condition_selects_only_strength_zero_edges(
        self, simple_edges, slide_features
    ):
        # TFA weak targets: G3, which is a SLIDE gene.
        rec = by_condition(run(simple_edges, 1, slide_features, {"TFA"}))[(0,)]
        assert sorted(rec[3][1]) == ["G3"]
        assert sorted(rec[3][0]) == ["G3"]

    def test_score_matches_scipy_for_strong_condition(self, simple_edges, slide_features):
        rec = by_condition(run(simple_edges, 1, slide_features, {"TFA"}))[(1,)]
        # dwngene = 2 (G1, G2), |SLIDE| = 3, common = 1 (G2)
        assert rec[2] == pytest.approx(expected(SLIDE_UNIVERSE, 2, 3, 1))

    def test_score_matches_scipy_for_a_fully_overlapping_tf(
        self, simple_edges, slide_features
    ):
        # TFB strong targets G2, G3 -- both SLIDE genes.
        rec = by_condition(run(simple_edges, 1, slide_features, {"TFB"}))[(1,)]
        assert sorted(rec[3][0]) == ["G2", "G3"]
        assert rec[2] == pytest.approx(expected(SLIDE_UNIVERSE, 2, 3, 2))

    def test_slide_genes_outside_the_graph_still_count_in_the_sample_size(
        self, simple_edges
    ):
        # G9 is not a GRN node. It cannot be an overlap, but it does enlarge the
        # SLIDE sample size S, which lowers the score. This pins that the caller
        # is responsible for intersecting with the graph, not this function.
        with_g9 = by_condition(run(simple_edges, 1, {"G2", "G3", "G9"}, {"TFB"}))[(1,)]
        without = by_condition(run(simple_edges, 1, {"G2", "G3"}, {"TFB"}))[(1,)]
        assert with_g9[2][0] < without[2][0]


class TestCommonTargetRequirement:
    """A TF combination only scores targets regulated by *every* member."""

    def test_order_two_requires_every_tf_to_regulate_the_target(
        self, simple_edges, slide_features
    ):
        # Targets shared by TFA and TFB are {G2, G3}; G1 (TFA only) and
        # G4 (TFB only) must never appear in any order-2 record.
        for rec in run(simple_edges, 2, slide_features, {"TFA", "TFB"}):
            assert "G1" not in rec[3][1] and "G4" not in rec[3][1]

    def test_tf_with_no_edges_yields_no_targets(self, simple_edges, slide_features):
        # TFC is absent from edges_df entirely. Pairing it with TFA must produce
        # empty records rather than silently falling back to TFA's own targets.
        for rec in run(simple_edges, 2, slide_features, {"TFA", "TFC"}):
            assert rec[3][0] == [] and rec[2] == (0, 1)

    def test_combination_with_no_shared_targets_is_neutral(self, slide_features):
        disjoint = pd.DataFrame(
            [
                {"source": "TFA", "target": "G1", "key": 0, "weight": 1.0, "strength": 1},
                {"source": "TFB", "target": "G2", "key": 0, "weight": 1.0, "strength": 1},
            ]
        )
        for rec in run(disjoint, 2, slide_features, {"TFA", "TFB"}):
            assert rec[2] == (0, 1)
            assert rec[3] == ([], [])


class TestHomogeneousStrengthConditions:
    """ISSUES.md #1 -- the highest-severity finding in this suite.

    ``get_slide_grn_enrichment`` keeps a target only if::

        filtered.groupby("target")[["source", "strength"]].nunique()
                .eq(len(filter_df)).all(axis=1)

    The ``source`` half is right: a target must be hit by all ``ord_tf`` TFs.
    The ``strength`` half compares the number of *distinct* strength values
    against ``ord_tf``.  Under a homogeneous condition such as ``(1, 1)`` every
    surviving row has strength 1, so that count is 1, never 2 -- and the target
    is discarded no matter how strong the evidence.

    Net effect: for ``ord_tf >= 2`` only *mixed* conditions can ever return a
    result.  ``build_enrichment_table`` nevertheless defaults to
    ``[(1, 1), (0, 1), (1, 0)]``, so the condition most likely to be of interest
    -- both TFs strongly co-regulating the target -- is structurally dead.

    ``ord_tf == 1`` is unaffected: ``len(filter_df) == 1`` and a single strength
    value gives ``nunique == 1``, which matches.
    """

    def test_currently_homogeneous_condition_returns_nothing(
        self, co_regulated_edges, slide_features
    ):
        # TFA and TFB both strongly regulate G1..G4; G2 and G3 are SLIDE genes.
        # The correct answer for (1, 1) is four downstream targets, two common.
        rec = by_condition(run(co_regulated_edges, 2, slide_features, {"TFA", "TFB"}))
        assert rec[(1, 1)][3] == ([], [])
        assert rec[(1, 1)][2] == (0, 1)

    def test_currently_all_order_two_conditions_are_empty_on_this_network(
        self, co_regulated_edges, slide_features
    ):
        for rec in run(co_regulated_edges, 2, slide_features, {"TFA", "TFB"}):
            assert rec[3] == ([], [])

    @pytest.mark.xfail(
        strict=True,
        reason="ISSUES.md #1: nunique(strength) is compared against ord_tf",
    )
    def test_should_find_co_regulated_targets_under_condition_1_1(
        self, co_regulated_edges, slide_features
    ):
        rec = by_condition(run(co_regulated_edges, 2, slide_features, {"TFA", "TFB"}))
        assert sorted(rec[(1, 1)][3][1]) == ["G1", "G2", "G3", "G4"]
        assert sorted(rec[(1, 1)][3][0]) == ["G2", "G3"]

    @pytest.mark.xfail(
        strict=True, reason="ISSUES.md #1: homogeneous (0, 0) is equally dead"
    )
    def test_should_find_co_regulated_targets_under_condition_0_0(self, slide_features):
        weak = pd.DataFrame(
            [
                {"source": tf, "target": g, "key": 0, "weight": 0.1, "strength": 0}
                for tf in ("TFA", "TFB")
                for g in ("G2", "G3")
            ]
        )
        rec = by_condition(run(weak, 2, slide_features, {"TFA", "TFB"}))
        assert sorted(rec[(0, 0)][3][1]) == ["G2", "G3"]

    def test_mixed_condition_does_work(self, slide_features):
        # The companion case: TFA strong and TFB weak on the same two targets
        # passes the same check, which is what isolates the defect to the
        # homogeneous case rather than to the check as a whole.
        mixed = pd.DataFrame(
            [
                {"source": "TFA", "target": g, "key": 0, "weight": 1.0, "strength": 1}
                for g in ("G2", "G3")
            ]
            + [
                {"source": "TFB", "target": g, "key": 0, "weight": 0.1, "strength": 0}
                for g in ("G2", "G3")
            ]
        )
        rec = by_condition(run(mixed, 2, slide_features, {"TFA", "TFB"}))
        assert sorted(rec[(1, 0)][3][1]) == ["G2", "G3"]
        assert sorted(rec[(1, 0)][3][0]) == ["G2", "G3"]

    def test_order_one_is_unaffected(self, simple_edges, slide_features):
        # Guards the claim that the defect is specific to ord_tf >= 2.
        rec = by_condition(run(simple_edges, 1, slide_features, {"TFB"}))[(1,)]
        assert sorted(rec[3][1]) == ["G2", "G3"]


class TestParallelEdgeDeduplication:
    def test_keeps_the_largest_absolute_weight(self, duplicate_edges):
        rec = by_condition(run(duplicate_edges, 1, {"G1", "G2"}, {"TFA"}))[(1,)]
        deduped = rec[5][1]
        g1 = deduped[deduped["target"] == "G1"]
        assert len(g1) == 1
        assert g1["weight"].iloc[0] == pytest.approx(-9.0)

    def test_dedup_is_by_absolute_value_not_signed_maximum(self, duplicate_edges):
        # +3.0 is the largest signed weight; -9.0 is the largest magnitude.
        rec = by_condition(run(duplicate_edges, 1, {"G1", "G2"}, {"TFA"}))[(1,)]
        weights = rec[5][1].set_index("target")["weight"]
        assert weights["G1"] == pytest.approx(-9.0)

    def test_downstream_count_uses_deduplicated_targets(self, duplicate_edges):
        # G1 appears three times but must be counted once.
        rec = by_condition(run(duplicate_edges, 1, {"G1", "G2"}, {"TFA"}))[(1,)]
        assert sorted(rec[3][1]) == ["G1", "G2"]

    def test_both_filtered_and_deduplicated_tables_are_returned(self, duplicate_edges):
        rec = by_condition(run(duplicate_edges, 1, {"G1", "G2"}, {"TFA"}))[(1,)]
        filtered, deduped = rec[5]
        assert len(filtered) == 4      # all rows matching the condition
        assert len(deduped) == 2       # one per (source, target)


class TestNoOverlapBranch:
    """ISSUES.md #3 -- the ``common == 0`` branch nests ``dwn_list``."""

    @pytest.fixture
    def five_targets(self):
        return pd.DataFrame(
            [
                {"source": "TFX", "target": f"D{i}", "key": 0, "weight": 1.0, "strength": 1}
                for i in range(5)
            ]
        )

    def test_currently_wraps_downstream_list_in_an_extra_list(self, five_targets):
        # No SLIDE gene is downstream, so the common == 0 branch runs. It appends
        # ``([], [dwn_list])`` where every other branch appends a flat list.
        rec = by_condition(run(five_targets, 1, {"ZZZ"}, {"TFX"}))[(1,)]
        assert rec[3][1] == [["D0", "D1", "D2", "D3", "D4"]]
        assert len(rec[3][1]) == 1          # reports 1 downstream gene...
        assert len(rec[3][1][0]) == 5       # ...when there are 5

    @pytest.mark.xfail(
        strict=True, reason="ISSUES.md #3: dwn_list is double-wrapped when common == 0"
    )
    def test_should_return_a_flat_downstream_list(self, five_targets):
        rec = by_condition(run(five_targets, 1, {"ZZZ"}, {"TFX"}))[(1,)]
        assert sorted(rec[3][1]) == ["D0", "D1", "D2", "D3", "D4"]

    def test_successful_branch_is_flat(self, five_targets):
        # The contrast case: with an overlap, dwnstrm is a flat list.
        rec = by_condition(run(five_targets, 1, {"D0"}, {"TFX"}))[(1,)]
        assert sorted(rec[3][1]) == ["D0", "D1", "D2", "D3", "D4"]

    def test_neutral_score_when_no_slide_gene_is_downstream(self, five_targets):
        rec = by_condition(run(five_targets, 1, {"ZZZ"}, {"TFX"}))[(1,)]
        assert rec[2] == (0, 1)


class TestEmptyInputs:
    def test_no_matching_source_produces_no_records(self, simple_edges, slide_features):
        assert run(simple_edges, 1, slide_features, {"TF_ABSENT"}) == []

    def test_empty_source_set_produces_no_records(self, simple_edges, slide_features):
        assert run(simple_edges, 1, slide_features, set()) == []

    def test_empty_slide_set_gives_neutral_scores(self, simple_edges):
        for rec in run(simple_edges, 1, set(), {"TFA", "TFB"}):
            assert rec[2] == (0, 1)
