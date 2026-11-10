"""Record-to-table conversion -- ``build_enrichment_table`` / ``write_enrichment_table``.

``build_enrichment_table`` applies the same significance mask at both orders but
*acts* on it differently: order 1 drops failing rows, order 2 keeps them and
zeroes their statistics.  These tests pin that asymmetry and the default
condition filters, using hand-made ``cc_dict`` records from ``conftest.py``.
"""

import pandas as pd
import pytest

from focal.state_specific._slide import build_enrichment_table, write_enrichment_table

from conftest import make_record

CF = ("c0",)

#: One record that passes every order-1 criterion, plus one failing each
#: criterion in isolation. Criteria: p < 0.05, len(dwnstrm) > 2,
#: len(common) > 1, score > 0, case in {slide, net, rnd}.
PASS = make_record(("A",), (1,), 2.0, 0.01, ["g1", "g2"], ["d1", "d2", "d3"])
FAIL_P = make_record(("B",), (1,), 2.0, 0.20, ["g1", "g2"], ["d1", "d2", "d3"])
FAIL_DWN = make_record(("C",), (1,), 2.0, 0.01, ["g1", "g2"], ["d1", "d2"])
FAIL_COMMON = make_record(("D",), (1,), 2.0, 0.01, ["g1"], ["d1", "d2", "d3"])
FAIL_SCORE = make_record(("E",), (1,), -1.0, 0.01, ["g1", "g2"], ["d1", "d2", "d3"])
WRONG_COND = make_record(("F",), (0,), 9.0, 0.01, ["g1", "g2"], ["d1", "d2", "d3"])
WRONG_CASE = make_record(
    ("G",), (1,), 3.0, 0.01, ["g1", "g2"], ["d1", "d2", "d3"], case="other"
)

ALL_ORDER_1 = [PASS, FAIL_P, FAIL_DWN, FAIL_COMMON, FAIL_SCORE, WRONG_COND, WRONG_CASE]


def cc(records, order):
    return {CF: {order: records}}


def tfs(df):
    return [t[0] for t in df["TF"]]


class TestOrderOneDropsFailingRows:
    def test_only_the_fully_passing_record_survives(self):
        out = build_enrichment_table(cc(ALL_ORDER_1, 1), CF, 1)
        assert tfs(out) == ["A"]

    @pytest.mark.parametrize(
        "record,label",
        [
            (FAIL_P, "p_value >= 0.05"),
            (FAIL_DWN, "len(dwnstrm) <= 2"),
            (FAIL_COMMON, "len(common) <= 1"),
            (FAIL_SCORE, "score <= 0"),
            (WRONG_CASE, "case not in {slide, net, rnd}"),
        ],
    )
    def test_each_criterion_excludes_on_its_own(self, record, label):
        out = build_enrichment_table(cc([PASS, record], 1), CF, 1)
        assert tfs(out) == ["A"], f"record failing {label} should be dropped"

    def test_boundary_p_value_is_exclusive(self):
        # p < 0.05, so exactly 0.05 must fail.
        edge = make_record(("H",), (1,), 2.0, 0.05, ["g1", "g2"], ["d1", "d2", "d3"])
        assert build_enrichment_table(cc([edge], 1), CF, 1).empty

    def test_boundary_counts_are_exclusive(self):
        # len(dwnstrm) > 2 and len(common) > 1 -> 3 and 2 are the first passes.
        ok = make_record(("I",), (1,), 2.0, 0.01, ["g1", "g2"], ["d1", "d2", "d3"])
        assert tfs(build_enrichment_table(cc([ok], 1), CF, 1)) == ["I"]

    def test_index_is_reset(self):
        out = build_enrichment_table(cc(ALL_ORDER_1, 1), CF, 1)
        assert list(out.index) == list(range(len(out)))

    def test_default_condition_filter_is_strong_only(self):
        # WRONG_COND has condition (0,) and a very high score; the default
        # filter_conditions=[(1,)] must exclude it before the mask even runs.
        out = build_enrichment_table(cc([WRONG_COND], 1), CF, 1)
        assert out.empty

    def test_explicit_condition_filter_overrides_the_default(self):
        out = build_enrichment_table(cc([WRONG_COND], 1), CF, 1, filter_conditions=[(0,)])
        assert tfs(out) == ["F"]

    def test_none_of_the_above_yields_empty_frame_with_columns(self):
        out = build_enrichment_table(cc([FAIL_SCORE], 1), CF, 1)
        assert out.empty
        assert {"TF", "condition", "case", "score", "p_value", "common", "dwnstrm"} <= set(
            out.columns
        )


class TestOrderTwoZeroesInsteadOfDropping:
    PASS_2 = make_record(("A", "B"), (1, 1), 2.0, 0.01, ["g1", "g2"], ["d1", "d2", "d3"])
    FAIL_2 = make_record(("C", "D"), (1, 1), 2.0, 0.20, ["g1", "g2"], ["d1", "d2", "d3"])

    def test_failing_rows_are_retained(self):
        out = build_enrichment_table(cc([self.PASS_2, self.FAIL_2], 2), CF, 2)
        assert len(out) == 2

    def test_failing_rows_have_statistics_neutralised(self):
        out = build_enrichment_table(cc([self.PASS_2, self.FAIL_2], 2), CF, 2)
        failed = out[out["TF"] == ("C", "D")].iloc[0]
        assert failed["score"] == 0.0 and failed["p_value"] == 1.0

    def test_passing_rows_keep_their_statistics(self):
        out = build_enrichment_table(cc([self.PASS_2, self.FAIL_2], 2), CF, 2)
        passed = out[out["TF"] == ("A", "B")].iloc[0]
        assert passed["score"] == pytest.approx(2.0)
        assert passed["p_value"] == pytest.approx(0.01)

    def test_default_conditions_include_the_structurally_dead_one(self):
        # ISSUES.md #1: (1, 1) is in the default set but can never be produced
        # by get_slide_grn_enrichment. Pinned here so the default and the defect
        # stay linked -- fixing #1 is what makes this default meaningful.
        recs = [
            make_record(("A", "B"), c, 1.0, 0.01, ["g1", "g2"], ["d1", "d2", "d3"])
            for c in [(1, 1), (0, 1), (1, 0), (0, 0)]
        ]
        out = build_enrichment_table(cc(recs, 2), CF, 2)
        assert {tuple(c) for c in out["condition"]} == {(1, 1), (0, 1), (1, 0)}

    def test_currently_sorted_by_pre_zeroing_score(self):
        # ISSUES.md #4: sort_values runs before the mask is applied, so a row
        # that ends up with score 0 can still sort above a row that kept a
        # positive score. Output order does not reflect final scores.
        high_but_failing = make_record(
            ("E", "F"), (0, 1), 9.0, 0.20, ["g1", "g2"], ["d1", "d2", "d3"]
        )
        out = build_enrichment_table(cc([self.PASS_2, high_but_failing], 2), CF, 2)
        assert list(out["TF"]) == [("E", "F"), ("A", "B")]
        assert out.iloc[0]["score"] == 0.0
        assert out.iloc[1]["score"] == pytest.approx(2.0)

    @pytest.mark.xfail(strict=True, reason="ISSUES.md #4: sorted before zeroing")
    def test_should_be_sorted_by_final_score(self):
        high_but_failing = make_record(
            ("E", "F"), (0, 1), 9.0, 0.20, ["g1", "g2"], ["d1", "d2", "d3"]
        )
        out = build_enrichment_table(cc([self.PASS_2, high_but_failing], 2), CF, 2)
        assert list(out["score"]) == sorted(out["score"], reverse=True)


class TestSharedBehaviour:
    def test_es_and_genes_are_split_into_named_columns(self):
        out = build_enrichment_table(cc([PASS], 1), CF, 1)
        assert "ES" not in out.columns and "Genes" not in out.columns
        row = out.iloc[0]
        assert row["score"] == pytest.approx(2.0)
        assert row["p_value"] == pytest.approx(0.01)
        assert row["common"] == ["g1", "g2"]
        assert row["dwnstrm"] == ["d1", "d2", "d3"]

    def test_order_one_output_is_sorted_by_score_descending(self):
        recs = [
            make_record((n,), (1,), s, 0.01, ["g1", "g2"], ["d1", "d2", "d3"])
            for n, s in [("low", 1.0), ("high", 5.0), ("mid", 3.0)]
        ]
        out = build_enrichment_table(cc(recs, 1), CF, 1)
        assert tfs(out) == ["high", "mid", "low"]

    @pytest.mark.parametrize("order", [0, 3, -1])
    def test_unsupported_order_raises(self, order):
        with pytest.raises(ValueError, match="must be 1 or 2"):
            build_enrichment_table(cc([PASS], order), CF, order)


class TestWriteEnrichmentTable:
    def test_creates_parent_directories(self, tmp_path):
        out = tmp_path / "a" / "b" / "enriched.csv"
        write_enrichment_table(build_enrichment_table(cc([PASS], 1), CF, 1), out)
        assert out.exists()

    def test_embedded_dataframes_are_json_serialised(self, tmp_path):
        rec = make_record(("A",), (1,), 2.0, 0.01, ["g1", "g2"], ["d1", "d2", "d3"])
        rec[5] = pd.DataFrame({"source": ["TFA"], "target": ["G1"]})
        out = tmp_path / "e.csv"
        write_enrichment_table(build_enrichment_table(cc([rec], 1), CF, 1), out)
        text = out.read_text()
        assert "TFA" in text and "source" in text

    def test_serialisation_can_be_disabled(self, tmp_path):
        out = tmp_path / "raw.csv"
        write_enrichment_table(
            build_enrichment_table(cc([PASS], 1), CF, 1), out, serialize_dfs=False
        )
        assert out.exists()

    def test_round_trips_through_pandas(self, tmp_path):
        out = tmp_path / "e.csv"
        table = build_enrichment_table(cc([PASS], 1), CF, 1)
        write_enrichment_table(table, out)
        back = pd.read_csv(out)
        assert len(back) == len(table)
        assert back.iloc[0]["score"] == pytest.approx(2.0)

    def test_does_not_mutate_the_input_frame(self, tmp_path):
        table = build_enrichment_table(cc([PASS], 1), CF, 1)
        before = table["dfs"].iloc[0]
        write_enrichment_table(table, tmp_path / "e.csv")
        assert table["dfs"].iloc[0] is before
