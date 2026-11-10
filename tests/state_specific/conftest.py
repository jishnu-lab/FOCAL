"""Shared fixtures for the ``focal.state_specific`` test-suite.

Everything here is hand-built mock data small enough that every expected value
can be worked out on paper, so a failing test points at the library and not at
the fixture.

Mock GRN
--------
Two clusters of CellOracle-style coefficient tables over 3 TFs and 5 targets.
``coef_abs`` values are chosen so the pooled distribution has an unambiguous
median, which is what ``combine_cluster_grn_links`` thresholds on::

    cluster 0                        cluster 1
    TFA -> G1   coef  2.0            TFA -> G1   coef  1.5
    TFA -> G2   coef  1.0            TFB -> G3   coef  0.2
    TFA -> G3   coef -0.5
    TFB -> G2   coef  3.0
    TFB -> G4   coef  0.1

pooled |coef| = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0]  -> median 1.0

so at ``quantile=0.5`` the threshold is exactly ``1.0`` and the edges with
``coef_abs >= 1.0`` (TFA->G1 twice, TFA->G2, TFB->G2) get ``strength == 1``.

``TFA->G1`` deliberately appears in *both* clusters so the MultiDiGraph parallel
edge keys and the max-|weight| dedup inside
:func:`~focal.state_specific._slide.get_slide_grn_enrichment` are both covered.

Co-regulation fixture
---------------------
:func:`co_regulated_edges` is a separate, deliberately trivial network where two
TFs *both* strongly regulate the same four targets.  It exists to make the
homogeneous-condition defect (ISSUES.md #1) unmistakable: under condition
``(1, 1)`` the correct answer is "all four targets", and anything else is a bug
in the library rather than an artefact of a sparse fixture.
"""

from __future__ import annotations

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")  # no interactive backend in CI


# --------------------------------------------------------------------------- #
# GRN link tables                                                              #
# --------------------------------------------------------------------------- #

#: cluster -> CellOracle-style edge table, as produced by load_celloracle_grn_tables
MOCK_LINKS: dict[int, pd.DataFrame] = {
    0: pd.DataFrame(
        {
            "source": ["TFA", "TFA", "TFA", "TFB", "TFB"],
            "target": ["G1", "G2", "G3", "G2", "G4"],
            "coef_mean": [2.0, 1.0, -0.5, 3.0, 0.1],
            "coef_abs": [2.0, 1.0, 0.5, 3.0, 0.1],
        }
    ),
    1: pd.DataFrame(
        {
            "source": ["TFA", "TFB"],
            "target": ["G1", "G3"],
            "coef_mean": [1.5, 0.2],
            "coef_abs": [1.5, 0.2],
        }
    ),
}

#: pooled |coef| values across both clusters, sorted -- the median is the cutoff
POOLED_COEF_ABS = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0]
MEDIAN_COEF_ABS = 1.0

#: SLIDE latent-factor genes. G9 is absent from the GRN on purpose, so that
#: ``slide_features`` is never a subset of the graph nodes.
MOCK_SLIDE_FEATURES = {"G2", "G3", "G9"}

#: SLIDE universe size passed as ``slide_starting_genes``
SLIDE_UNIVERSE = 100


@pytest.fixture
def mock_links() -> dict[int, pd.DataFrame]:
    """Per-cluster CellOracle coefficient tables (fresh copy per test)."""
    return {k: v.copy() for k, v in MOCK_LINKS.items()}


@pytest.fixture
def slide_features() -> set[str]:
    return set(MOCK_SLIDE_FEATURES)


# --------------------------------------------------------------------------- #
# Flat edge tables (what get_slide_grn_enrichment actually consumes)           #
# --------------------------------------------------------------------------- #

def make_edges(rows) -> pd.DataFrame:
    """Build an ``edges_df`` from ``(source, target, key, weight, strength)`` tuples."""
    return pd.DataFrame(
        [
            {"source": s, "target": t, "key": k, "weight": w, "strength": st}
            for s, t, k, w, st in rows
        ]
    )


@pytest.fixture
def simple_edges() -> pd.DataFrame:
    """A 2-TF network with a mix of strong and weak edges.

    ::

        TFA -> G1  strength 1      TFB -> G2  strength 1
        TFA -> G2  strength 1      TFB -> G3  strength 1
        TFA -> G3  strength 0      TFB -> G4  strength 0

    Targets shared by both TFs: {G2, G3}.
    """
    return make_edges(
        [
            ("TFA", "G1", 0, 2.0, 1),
            ("TFA", "G2", 0, 1.0, 1),
            ("TFA", "G3", 0, -0.5, 0),
            ("TFB", "G2", 0, 3.0, 1),
            ("TFB", "G3", 0, 0.2, 1),
            ("TFB", "G4", 0, 0.1, 0),
        ]
    )


@pytest.fixture
def co_regulated_edges() -> pd.DataFrame:
    """TFA and TFB both *strongly* regulate G1..G4 -- nothing weak, nothing partial.

    Under condition ``(1, 1)`` the only defensible answer is all four targets.
    """
    return make_edges(
        [
            (tf, g, 0, 1.0, 1)
            for tf in ("TFA", "TFB")
            for g in ("G1", "G2", "G3", "G4")
        ]
    )


@pytest.fixture
def duplicate_edges() -> pd.DataFrame:
    """``TFA->G1`` three times with different weights; max |weight| is -9.0."""
    df = make_edges(
        [
            ("TFA", "G1", 0, 0.5, 1),
            ("TFA", "G1", 1, -9.0, 1),
            ("TFA", "G1", 2, 3.0, 1),
            ("TFA", "G2", 0, 1.0, 1),
        ]
    )
    df["cluster"] = ["c0", "c1", "c2", "c0"]
    return df


# --------------------------------------------------------------------------- #
# cc_dict records                                                              #
# --------------------------------------------------------------------------- #

def make_record(tf, condition, score, p_value, common, dwnstrm, case="slide"):
    """One ``cc_dict`` record in the 6-slot layout the pipeline uses.

    ``[TF, condition, (score, p), (common, dwnstrm), case, (filtered, deduped)]``
    """
    return [
        tuple(tf),
        tuple(condition),
        (score, p_value),
        (list(common), list(dwnstrm)),
        case,
        (pd.DataFrame(), pd.DataFrame()),
    ]
