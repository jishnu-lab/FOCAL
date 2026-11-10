"""``StateSpecificEnrichment.run_enrichment`` and its post-processing.

``run_enrichment`` is an orchestrator: it fuses clusters, thresholds edges,
widens the TF candidate set through the graph, then delegates scoring to
``_slide``.  These tests drive it end to end on the mock GRN.

``load_grn_data`` imports celloracle, which is heavy and not a declared
dependency of the package.  Every test here therefore sets
``GRN_links_after_fit`` / ``GRN_TFs`` / ``slide_features`` directly -- exactly
the three attributes ``run_enrichment`` guards on -- so the suite exercises the
real orchestration without needing celloracle installed.
"""

import pickle

import pandas as pd
import pytest

from focal.state_specific import StateSpecificEnrichment

from conftest import MOCK_LINKS, MOCK_SLIDE_FEATURES, SLIDE_UNIVERSE


def make_sse(tmp_path, *, clusters=(0, 1), ord_clust=(1,), ord_tf=(1,), label="", **kw):
    sse = StateSpecificEnrichment(
        grn_workdir=tmp_path,
        oracle_object_name="unused.celloracle.oracle",
        feature_folder=tmp_path,
        out_path=tmp_path,
        slide_starting_genes=SLIDE_UNIVERSE,
        clusters_of_interest=list(clusters),
        order_fr_clust=list(ord_clust),
        order_fr_tfcomb=list(ord_tf),
        quantile=kw.pop("quantile", 0.5),
        experiment_label=label,
        **kw,
    )
    sse.GRN_links_after_fit = {k: v.copy() for k, v in MOCK_LINKS.items()}
    sse.GRN_TFs = ["TFA", "TFB"]
    sse.slide_features = set(MOCK_SLIDE_FEATURES)
    return sse


@pytest.fixture
def sse(tmp_path):
    return make_sse(tmp_path)


def run(obj):
    obj.run_enrichment(save_pickle=True, show_progress=False, save_strength_diagnostics=False)


class TestConstruction:
    def test_creates_output_directories_eagerly(self, tmp_path):
        make_sse(tmp_path)
        assert (tmp_path / "figures").is_dir()
        assert (tmp_path / "out_files" / "SLIDE_LF_enrichment").is_dir()

    def test_accepts_string_paths(self, tmp_path):
        sse = make_sse(tmp_path)
        assert sse.out_path == tmp_path

    def test_state_starts_empty(self, tmp_path):
        sse = StateSpecificEnrichment(
            grn_workdir=tmp_path, oracle_object_name="x", feature_folder=tmp_path,
            out_path=tmp_path, slide_starting_genes=10, clusters_of_interest=[0],
            order_fr_clust=[1], order_fr_tfcomb=[1],
        )
        assert sse.GRN_links_after_fit is None
        assert sse.slide_features is None
        assert sse.cc_dicts == {} and sse.enrichment_dfs == {}


class TestPreconditions:
    @pytest.mark.parametrize("missing", ["GRN_links_after_fit", "GRN_TFs", "slide_features"])
    def test_requires_grn_and_slide_data(self, tmp_path, missing):
        sse = make_sse(tmp_path)
        setattr(sse, missing, None)
        with pytest.raises(RuntimeError, match="load_grn_data"):
            sse.run_enrichment(show_progress=False)


class TestClusterFusions:
    @pytest.mark.parametrize(
        "clusters,orders,expected",
        [
            ([0, 1], [1], [(0,), (1,)]),
            ([0, 1], [2], [(0, 1)]),
            ([0, 1], [1, 2], [(0,), (1,), (0, 1)]),
            ([0, 1, 2], [2], [(0, 1), (0, 2), (1, 2)]),
        ],
    )
    def test_fusions_are_combinations_of_the_requested_orders(
        self, tmp_path, clusters, orders, expected
    ):
        sse = make_sse(tmp_path, clusters=clusters, ord_clust=orders)
        sse.build_cluster_fusions()
        assert sse.cluster_fusions == expected

    def test_rebuilding_does_not_accumulate(self, sse):
        sse.build_cluster_fusions()
        first = list(sse.cluster_fusions)
        sse.build_cluster_fusions()
        assert sse.cluster_fusions == first


class TestRunEnrichment:
    def test_populates_cc_dicts_and_enrichment_dfs(self, sse):
        run(sse)
        for cf in sse.cluster_fusions:
            assert 1 in sse.cc_dicts[cf]
            assert isinstance(sse.enrichment_dfs[cf][1], pd.DataFrame)

    def test_writes_one_csv_and_one_pickle_per_fusion_and_order(self, sse):
        run(sse)
        d = sse.out_path / "out_files" / "SLIDE_LF_enrichment"
        assert len(list(d.glob("enriched_df_*.csv"))) == len(sse.cluster_fusions)
        assert len(list(d.glob("cc_dict_*.pickle"))) == len(sse.cluster_fusions)

    def test_filename_stem_encodes_order_fusion_and_label(self, tmp_path):
        sse = make_sse(tmp_path, label="expt")
        run(sse)
        d = sse.out_path / "out_files" / "SLIDE_LF_enrichment"
        assert (d / "enriched_df_1_TFs_(0,)_expt.csv").exists()
        assert (d / "cc_dict_1_TFs_(0,)_expt.pickle").exists()

    def test_no_label_means_no_suffix(self, sse):
        run(sse)
        d = sse.out_path / "out_files" / "SLIDE_LF_enrichment"
        assert (d / "enriched_df_1_TFs_(0,).csv").exists()

    def test_save_pickle_false_skips_pickles(self, sse):
        sse.run_enrichment(save_pickle=False, show_progress=False, save_strength_diagnostics=False)
        d = sse.out_path / "out_files" / "SLIDE_LF_enrichment"
        assert list(d.glob("enriched_df_*.csv"))
        assert not list(d.glob("cc_dict_*.pickle"))

    def test_pickle_round_trips_to_the_in_memory_structure(self, sse):
        run(sse)
        cf = sse.cluster_fusions[0]
        path = sse.out_path / "out_files" / "SLIDE_LF_enrichment" / f"cc_dict_1_TFs_{cf}.pickle"
        with open(path, "rb") as f:
            loaded = pickle.load(f)
        assert loaded.keys() == sse.cc_dicts[cf][1].keys()

    def test_strength_diagnostics_write_a_figure(self, sse):
        sse.run_enrichment(save_pickle=False, show_progress=False, save_strength_diagnostics=True)
        assert list((sse.out_path / "figures").glob("combined_links_key_strength_*.pdf"))

    def test_multiple_tf_orders_are_all_run(self, tmp_path):
        # Fusion (0, 1) keeps both TFA and TFB as candidates, so order 2 has a
        # combination to enumerate. (Fusion (1,) alone keeps only TFB and would
        # hit the empty-combination crash -- see TestEmptyCombinationCrash.)
        sse = make_sse(tmp_path, ord_clust=(2,), ord_tf=(1, 2))
        run(sse)
        for cf in sse.cluster_fusions:
            assert set(sse.cc_dicts[cf]) == {1, 2}

    def test_tf_candidates_are_restricted_to_known_tfs(self, tmp_path):
        # run_enrichment intersects graph predecessors with GRN_TFs, so a
        # regulator absent from GRN_TFs is never scored. Fusion (0,) is used
        # because TFA is a candidate there; see TestEmptyCombinationCrash for
        # what happens when the intersection empties out.
        sse = make_sse(tmp_path, clusters=(0,))
        sse.GRN_TFs = ["TFA"]
        run(sse)
        cf = sse.cluster_fusions[0]
        seen = {r[0][0] for r in sse.cc_dicts[cf][1][cf][1]}
        assert seen == {"TFA"}

    def test_quantile_changes_edge_strengths(self, tmp_path):
        # At q=0 every edge is strong, so the weak-condition records carry no
        # downstream genes; at q=1 almost everything is weak.
        lenient = make_sse(tmp_path / "a", quantile=0.0)
        strict = make_sse(tmp_path / "b", quantile=1.0)
        run(lenient)
        run(strict)
        cf = lenient.cluster_fusions[0]
        weak_lenient = [r for r in lenient.cc_dicts[cf][1][cf][1] if r[1] == (0,)]
        weak_strict = [r for r in strict.cc_dicts[cf][1][cf][1] if r[1] == (0,)]
        assert all(r[3][1] in ([], [[]]) or not r[3][1] for r in weak_lenient)
        assert any(r[3][1] for r in weak_strict)

    def test_records_helper_navigates_the_double_index(self, sse):
        # cc_dicts[cf][ord] stores the whole cc_dict, itself keyed by cf then
        # ord -- the layout kept for pickle compatibility with the TA pipeline.
        run(sse)
        cf = sse.cluster_fusions[0]
        assert sse._records(cf, 1) is sse.cc_dicts[cf][1][cf][1]


class TestPostProcessing:
    @pytest.fixture
    def ran(self, sse):
        run(sse)
        return sse

    def test_load_results_restores_state_from_disk(self, ran):
        cf = ran.cluster_fusions[0]
        fresh = make_sse(ran.out_path)
        fresh.load_results(cf, 1)
        assert fresh.cc_dicts[cf][1].keys() == ran.cc_dicts[cf][1].keys()
        assert len(fresh.enrichment_dfs[cf][1]) == len(ran.enrichment_dfs[cf][1])

    def test_used_weight_df_concatenates_deduplicated_edges(self, ran):
        cf = ran.cluster_fusions[0]
        out = ran.build_used_weight_df(cf, 1)
        assert "case" in out.columns
        assert (ran.out_path / "out_files" / "SLIDE_LF_enrichment"
                / f"1_TFs_{cf}_used_weights.csv").exists()

    def test_used_weight_df_rows_match_the_record_tables(self, ran):
        cf = ran.cluster_fusions[0]
        out = ran.build_used_weight_df(cf, 1)
        expected = sum(len(r[5][1]) for r in ran._records(cf, 1))
        assert len(out) == expected

    def test_slide_lf_enriched_writes_its_csv(self, ran):
        cf = ran.cluster_fusions[0]
        ran.build_used_weight_df(cf, 1)
        ran.build_slide_lf_enriched(cf, 1)
        assert (ran.out_path / "out_files" / "SLIDE_LF_enrichment"
                / f"enriched_df_1_TFs_{cf}_used_weights_strong.csv").exists()

    def test_slide_lf_enriched_columns(self, ran):
        cf = ran.cluster_fusions[0]
        ran.build_used_weight_df(cf, 1)
        out = ran.build_slide_lf_enriched(cf, 1)
        assert list(out.columns[:2]) == ["source", "target"]

    def test_tf_enrichment_scores_are_the_all_strong_condition(self, ran):
        cf = ran.cluster_fusions[0]
        scores = ran.tf_enrichment_scores(cf, 1)
        strong = {r[0][0]: r[2][0] for r in ran._records(cf, 1) if r[1] == (1,)}
        assert scores == strong

    @pytest.mark.parametrize("order", [2, 3])
    def test_tf_enrichment_scores_rejects_higher_orders(self, tmp_path, order):
        sse = make_sse(tmp_path, ord_clust=(2,), ord_tf=(1, 2))
        run(sse)
        with pytest.raises(ValueError, match="only defined for ord_tf=1"):
            sse.tf_enrichment_scores(sse.cluster_fusions[0], order)


class TestEmptyCombinationCrash:
    """ISSUES.md #2 -- ``run_enrichment`` aborts when a fusion yields too few TFs.

    ``get_slide_grn_enrichment`` enumerates
    ``itertools.combinations(sorted(total_source), ord_tf)``.  When the candidate
    set is smaller than ``ord_tf`` that is empty, so no records are appended and
    ``cc_dict[cf][ord]`` stays ``[]``.  ``_records_to_enrichment_df`` then builds
    a zero-column DataFrame and the assignment

        enrichment_df[["score", "p_value"]] = pd.DataFrame(...tolist()...)

    raises ``ValueError: Columns must be same length as key``.

    Nothing catches it, so one sparse cluster fusion aborts the whole run and
    discards the fusions already computed.  The candidate set is
    ``(slide_in_graph | predecessors) & GRN_TFs``, which empties out easily on a
    sparse cluster or a restrictive TF list.
    """

    def test_currently_crashes_when_fewer_candidates_than_order(self, tmp_path):
        sse = make_sse(tmp_path, clusters=(1,), ord_tf=(2,))
        with pytest.raises(ValueError, match="Columns must be same length as key"):
            run(sse)

    def test_currently_crashes_when_no_candidate_tf_survives(self, tmp_path):
        sse = make_sse(tmp_path, clusters=(0,))
        sse.GRN_TFs = ["TF_NOT_IN_GRN"]
        with pytest.raises(ValueError, match="Columns must be same length as key"):
            run(sse)

    def test_currently_build_enrichment_table_cannot_handle_empty_records(self):
        from focal.state_specific._slide import build_enrichment_table

        for order in (1, 2):
            with pytest.raises(ValueError, match="Columns must be same length as key"):
                build_enrichment_table({("c",): {order: []}}, ("c",), order)

    @pytest.mark.xfail(
        strict=True, reason="ISSUES.md #2: empty record list is not handled"
    )
    def test_should_produce_an_empty_table_not_a_crash(self):
        from focal.state_specific._slide import build_enrichment_table

        out = build_enrichment_table({("c",): {1: []}}, ("c",), 1)
        assert out.empty

    @pytest.mark.xfail(
        strict=True, reason="ISSUES.md #2: a sparse fusion aborts the whole run"
    )
    def test_should_complete_when_one_fusion_has_too_few_tfs(self, tmp_path):
        # Fusion (0,) has two candidates and fusion (1,) has one. Requesting
        # order 2 should skip the sparse fusion, not destroy the run.
        sse = make_sse(tmp_path, clusters=(0, 1), ord_clust=(1,), ord_tf=(2,))
        run(sse)
        assert (0,) in sse.cc_dicts
