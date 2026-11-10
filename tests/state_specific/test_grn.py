"""State-specific GRN construction -- :mod:`focal.state_specific._grn`.

Covers cluster fusion, quantile edge thresholding, the MultiDiGraph edge view,
and network-score subsetting.  Threshold expectations are derived from the
pooled coefficient list in ``conftest.py`` and cross-checked against
``numpy.quantile``, never against the function under test.
"""

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from focal.state_specific._grn import (
    combine_cluster_grn_links,
    filter_network_scores_for_clusters,
    grn_edges_from_combined_links,
)

from conftest import MEDIAN_COEF_ABS, POOLED_COEF_ABS


class TestCombineClusterGrnLinks:
    def test_concatenates_every_cluster_and_tags_provenance(self, mock_links):
        combined, _ = combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5)
        assert len(combined) == len(mock_links[0]) + len(mock_links[1]) == 7
        assert sorted(combined["cluster"].unique()) == [0, 1]

    def test_threshold_is_the_requested_quantile_of_pooled_coefficients(self, mock_links):
        _, threshold = combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5)
        assert threshold == pytest.approx(np.quantile(POOLED_COEF_ABS, 0.5))
        assert threshold == pytest.approx(MEDIAN_COEF_ABS)  # 1.0, by hand

    def test_threshold_is_pooled_not_per_cluster(self, mock_links):
        # Cluster 0 alone has median |coef| 1.0; cluster 1 alone has 0.85.
        # The fused threshold must come from the 7 pooled values, not either one.
        _, fused = combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5)
        _, only0 = combine_cluster_grn_links(mock_links, (0,), quantile=0.5)
        assert fused == pytest.approx(1.0)
        assert only0 == pytest.approx(1.0)
        _, only1 = combine_cluster_grn_links(mock_links, (1,), quantile=0.5)
        assert only1 == pytest.approx(np.quantile([0.2, 1.5], 0.5))

    def test_strength_is_one_at_or_above_threshold(self, mock_links):
        combined, threshold = combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5)
        expected = (combined["coef_abs"] >= threshold).astype(int)
        pd.testing.assert_series_equal(
            combined["strength"], expected, check_names=False
        )

    def test_threshold_boundary_is_inclusive(self, mock_links):
        # TFA->G2 has coef_abs exactly 1.0 == threshold, so it must be strong.
        combined, threshold = combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5)
        boundary = combined[combined["coef_abs"] == threshold]
        assert not boundary.empty
        assert (boundary["strength"] == 1).all()

    @pytest.mark.parametrize("quantile,expected_strong", [(0.0, 7), (1.0, 1)])
    def test_extreme_quantiles(self, mock_links, quantile, expected_strong):
        # q=0 -> threshold is the minimum, every edge is strong.
        # q=1 -> threshold is the maximum, only the single largest edge is strong.
        combined, _ = combine_cluster_grn_links(mock_links, (0, 1), quantile=quantile)
        assert int(combined["strength"].sum()) == expected_strong

    def test_cluster_order_is_preserved(self, mock_links):
        combined, _ = combine_cluster_grn_links(mock_links, (1, 0), quantile=0.5)
        assert list(combined["cluster"]) == [1, 1, 0, 0, 0, 0, 0]

    def test_does_not_mutate_the_input_tables(self, mock_links):
        before = {k: v.copy() for k, v in mock_links.items()}
        combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5)
        for k, v in before.items():
            pd.testing.assert_frame_equal(mock_links[k], v)

    def test_single_cluster_fusion(self, mock_links):
        combined, _ = combine_cluster_grn_links(mock_links, (0,), quantile=0.5)
        assert len(combined) == 5
        assert set(combined["cluster"]) == {0}

    def test_histogram_is_written_when_requested(self, mock_links, tmp_path):
        out = tmp_path / "nested" / "hist.pdf"
        combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5, histogram_path=out)
        assert out.exists() and out.stat().st_size > 0


class TestGrnEdgesFromCombinedLinks:
    def test_returns_multidigraph_and_matching_edge_table(self, mock_links):
        combined, _ = combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5)
        grn, edges_df = grn_edges_from_combined_links(combined)
        assert isinstance(grn, nx.MultiDiGraph)
        assert grn.number_of_edges() == len(edges_df) == 7

    def test_parallel_edges_get_distinct_keys(self, mock_links):
        # TFA->G1 exists in both clusters and must survive as two edges.
        combined, _ = combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5)
        _, edges_df = grn_edges_from_combined_links(combined)
        parallel = edges_df[(edges_df["source"] == "TFA") & (edges_df["target"] == "G1")]
        assert len(parallel) == 2
        assert sorted(parallel["key"]) == [0, 1]
        assert sorted(parallel["cluster"]) == [0, 1]

    def test_coef_mean_becomes_weight_preserving_sign(self, mock_links):
        combined, _ = combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5)
        _, edges_df = grn_edges_from_combined_links(combined)
        negative = edges_df[(edges_df["source"] == "TFA") & (edges_df["target"] == "G3")]
        assert negative["weight"].iloc[0] == pytest.approx(-0.5)

    def test_edge_table_columns(self, mock_links):
        combined, _ = combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5)
        _, edges_df = grn_edges_from_combined_links(combined)
        assert set(edges_df.columns) == {
            "source", "target", "key", "weight", "cluster", "strength"
        }

    def test_predecessors_drive_tf_candidate_selection(self, mock_links):
        # run_enrichment calls grn.predecessors(gene) to widen the TF candidate
        # set, so direction matters: G2 is regulated by both TFs.
        combined, _ = combine_cluster_grn_links(mock_links, (0, 1), quantile=0.5)
        grn, _ = grn_edges_from_combined_links(combined)
        assert sorted(grn.predecessors("G2")) == ["TFA", "TFB"]
        assert sorted(grn.successors("TFB")) == ["G2", "G3", "G4"]


class TestFilterNetworkScoresForClusters:
    @pytest.fixture
    def scores(self):
        return pd.DataFrame(
            {"degree_centrality_out": [0.9, 0.5, 0.7, 0.2], "cluster": [0, 0, 1, 1]},
            index=["GA", "GB", "GA", "GC"],
        )

    def test_subsets_to_requested_clusters(self, scores):
        assert set(filter_network_scores_for_clusters((0,), scores).index) == {"GA", "GB"}

    def test_sorted_by_out_centrality_descending(self, scores):
        out = filter_network_scores_for_clusters((0, 1), scores)
        assert list(out["degree_centrality_out"]) == sorted(
            out["degree_centrality_out"], reverse=True
        )

    def test_duplicate_gene_keeps_highest_centrality_row(self, scores):
        # GA is in both clusters (0.9 and 0.7). Sorting happens before the
        # duplicate drop, so the 0.9 row survives.
        out = filter_network_scores_for_clusters((0, 1), scores)
        assert out.loc["GA", "degree_centrality_out"] == pytest.approx(0.9)
        assert list(out.index) == ["GA", "GB", "GC"]

    def test_accepts_numeric_strings(self, scores):
        assert not filter_network_scores_for_clusters(("0", "1"), scores).empty

    def test_currently_raises_on_non_numeric_cluster_ids(self, scores):
        # ISSUES.md #8. StateSpecificEnrichment types clusters_of_interest as
        # list[Any] and combine_cluster_grn_links accepts any hashable, but this
        # helper hard-casts with int().  CellOracle cluster names are routinely
        # strings such as "ActB_1".
        named = scores.copy()
        named["cluster"] = ["ActB_1", "ActB_1", "GC_1", "GC_1"]
        with pytest.raises(ValueError, match="invalid literal for int"):
            filter_network_scores_for_clusters(("ActB_1",), named)

    @pytest.mark.xfail(
        strict=True, reason="ISSUES.md #8: int() cast rejects string cluster names"
    )
    def test_should_accept_any_cluster_label(self, scores):
        named = scores.copy()
        named["cluster"] = ["ActB_1", "ActB_1", "GC_1", "GC_1"]
        out = filter_network_scores_for_clusters(("ActB_1",), named)
        assert set(out.index) == {"GA", "GB"}
