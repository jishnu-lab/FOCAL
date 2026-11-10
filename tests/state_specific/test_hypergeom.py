"""Maths of :func:`focal.state_specific._slide.hypergeom_slide_grn_score`.

Every expected value here is computed independently -- either on paper or with
``scipy.stats.hypergeom`` called directly -- never by invoking the function
under test.

Parameter naming follows the call site in
:func:`~focal.state_specific._slide.get_slide_grn_enrichment`::

    hypergeom_slide_grn_score(slide_starting_genes,  # P  SLIDE universe
                              dwngene,               # p  downstream targets
                              len(slide_features),   # S  LF genes
                              common)                # s  LF genes downstream

``xfail(strict=True)`` marks a defect; see ``ISSUES.md``.
"""

import math

import numpy as np
import pytest
from scipy.stats import hypergeom

from focal.state_specific._slide import hypergeom_slide_grn_score


class TestAgreesWithScipy:
    """The happy path must reproduce scipy exactly."""

    @pytest.mark.parametrize(
        "P,p,S,s",
        [
            (100, 20, 10, 5),
            (100, 50, 10, 1),      # depleted: score must be negative
            (500, 120, 40, 20),
            (50, 10, 10, 3),
        ],
    )
    def test_score_and_pvalue_match_independent_computation(self, P, p, S, s):
        score, pval = hypergeom_slide_grn_score(P, p, S, s)
        assert pval == pytest.approx(float(1.0 - hypergeom.cdf(s - 1, P, p, S)))
        assert score == pytest.approx(math.log2((s / S) / (p / P)))

    def test_enrichment_above_expectation_is_positive(self):
        # expected overlap = S*p/P = 10*20/100 = 2; observing 5 is enrichment
        score, _ = hypergeom_slide_grn_score(100, 20, 10, 5)
        assert score > 0

    def test_depletion_below_expectation_is_negative(self):
        # expected = 10*50/100 = 5; observing 1 is depletion
        score, _ = hypergeom_slide_grn_score(100, 50, 10, 1)
        assert score < 0

    def test_observed_equals_expected_is_zero(self):
        # s/S == p/P  ->  log2(1) == 0.  2/10 == 20/100
        score, _ = hypergeom_slide_grn_score(100, 20, 10, 2)
        assert score == pytest.approx(0.0)

    def test_score_is_ratio_of_rates_not_counts(self):
        # Doubling the sample while holding the rate fixed must not move the score.
        a, _ = hypergeom_slide_grn_score(100, 20, 10, 4)
        b, _ = hypergeom_slide_grn_score(100, 20, 20, 8)
        assert a == pytest.approx(b)


class TestGuards:
    """Degenerate inputs return the neutral ``(0.0, 1.0)`` rather than raising."""

    @pytest.mark.parametrize(
        "P,p,S,s",
        [
            (100, 20, 0, 0),    # empty sample
            (0, 0, 0, 0),       # empty universe
            (100, 20, 10, 0),   # no overlap
            (100, 20, 10, -1),  # negative overlap
            (-5, 20, 10, 5),    # negative universe
        ],
    )
    def test_degenerate_input_is_neutral(self, P, p, S, s):
        assert hypergeom_slide_grn_score(P, p, S, s) == (0.0, 1.0)

    def test_neutral_result_is_plain_floats(self):
        # The values flow into a DataFrame and then a CSV; numpy scalars or None
        # would change the serialised dtype.
        score, pval = hypergeom_slide_grn_score(100, 20, 10, 0)
        assert isinstance(score, float) and isinstance(pval, float)


class TestInconsistentPopulation:
    """``p > P`` -- more downstream targets than the declared SLIDE universe.

    Reachable whenever ``slide_starting_genes`` is set smaller than the number of
    downstream targets a TF combination has in the GRN.  Nothing in the call
    chain validates it.  See ISSUES.md #5.
    """

    def test_currently_returns_nan_pvalue(self):
        # Pins present behaviour: scipy yields nan for an out-of-range population.
        _, pval = hypergeom_slide_grn_score(50, 80, 10, 5)
        assert np.isnan(pval)

    def test_currently_returns_finite_score_alongside_nan_pvalue(self):
        # The score is computed from the ratio alone, so it stays finite even
        # though the p-value is undefined -- the pair is internally inconsistent.
        score, pval = hypergeom_slide_grn_score(50, 80, 10, 5)
        assert np.isfinite(score)
        assert np.isnan(pval)

    @pytest.mark.xfail(strict=True, reason="ISSUES.md #5: p > P is not validated")
    def test_should_reject_or_neutralise_impossible_population(self):
        _, pval = hypergeom_slide_grn_score(50, 80, 10, 5)
        assert not np.isnan(pval), "an undefined p-value should not reach the caller"

    def test_nan_pvalue_silently_fails_downstream_significance_filter(self):
        # Documents the consequence: NaN < 0.05 is False, so a row carrying an
        # undefined p-value is dropped (ord 1) or zeroed (ord 2) with no warning.
        _, pval = hypergeom_slide_grn_score(50, 80, 10, 5)
        assert (pval < 0.05) is False


class TestImpossibleSuccessCounts:
    """``s > p`` and ``s > S`` are arithmetically impossible but accepted."""

    def test_currently_accepts_more_hits_than_population_successes(self):
        # s=5 hits drawn from only p=3 successes cannot happen.
        score, pval = hypergeom_slide_grn_score(100, 3, 10, 5)
        assert score > 0 and pval == pytest.approx(0.0)

    def test_currently_accepts_more_hits_than_sample_size(self):
        # s=8 hits in a sample of S=5 cannot happen.
        score, pval = hypergeom_slide_grn_score(100, 20, 5, 8)
        assert score > 0 and pval == pytest.approx(0.0)

    @pytest.mark.xfail(strict=True, reason="ISSUES.md #6: no bounds validation")
    def test_should_reject_impossible_success_counts(self):
        with pytest.raises(ValueError):
            hypergeom_slide_grn_score(100, 3, 10, 5)


class TestZeroPopulationSuccesses:
    """``p == 0`` with ``s > 0``: the score guard and the p-value disagree."""

    def test_currently_score_is_zero_but_pvalue_is_zero_too(self):
        # expected = S*p/P = 0 trips the `expected > 0` guard, so score is 0.0
        # ("no enrichment"), while the p-value is 0.0 ("maximally significant").
        score, pval = hypergeom_slide_grn_score(100, 0, 10, 3)
        assert score == 0.0
        assert pval == pytest.approx(0.0)

    @pytest.mark.xfail(
        strict=True, reason="ISSUES.md #7: score guard and p-value disagree"
    )
    def test_should_be_self_consistent(self):
        score, pval = hypergeom_slide_grn_score(100, 0, 10, 3)
        assert not (score == 0.0 and pval < 0.05), (
            "a neutral score must not come with a maximally significant p-value"
        )

    def test_guard_prevents_division_by_zero(self):
        # Whatever the semantics, it must not raise.
        hypergeom_slide_grn_score(100, 0, 10, 3)
