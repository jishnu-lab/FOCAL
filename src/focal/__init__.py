"""Focal: Functional and Interpretable Regulatory Encoding of cell Fate decisions.

Focal is three modules in one package, each answering a different question about
regulation:

``focal.temporal``
    **FocalTemporal** -- how does TF regulation change *along* a trajectory?
    Transition-window and episodic GRNs, force waves, regulatory phases
    (capabilities 3, 4, dynamic 5).

``focal.state_specific``
    **FocalStateSpecific** -- what separates two *fixed* states, and what happens
    if we perturb it? State and cross-state GRNs, enrichment against them
    (capabilities 2, static 5; 1 and 6 to follow).

``focal.cross_prediction``
    **FocalCrossPrediction** -- do programs learned on one dataset *transfer* to
    stratify another? (capability 7; scaffolded, port in progress).

Each module's entry point is its manager. Supporting namespaces:
:mod:`focal.base` (shared contracts and the ORA primitive),
:mod:`focal.backends` (dictys and friends), :mod:`focal.io`,
:mod:`focal.utils`, :mod:`focal.cli`.

There is no separate plotting namespace: every figure lives in the module that
owns its subject, so ``focal.temporal`` carries the force landscapes and phase
heatmaps and ``focal.state_specific`` carries the enrichment bars.

Examples
--------
>>> import focal as ff
>>> mgr = ff.TemporalManager(dyn_obj, output_dir="./episodes")
>>> enrichment = mgr.enrich_episode(1, slice(0, 5), lf_genes=lf_blimp1)
"""

__version__ = "0.1.0"
__author__ = "Akanksha Sachan"

from focal import (
    base,
    cross_prediction,
    io,
    state_specific,
    temporal,
    utils,
)
from focal.base import BaseManager, ora
from focal.state_specific import StateSpecificManager
from focal.temporal import (
    AlignTimeScales,
    EpisodeDynamics,
    SmoothedCurvesChromatin,
    SmoothedCurvesGRN,
    StateFrequency,
    TemporalManager,
    TFForceValidation,
    TFForceWaves,
)

__all__ = [
    "AlignTimeScales",
    "BaseManager",
    "EpisodeDynamics",
    "SmoothedCurvesChromatin",
    "SmoothedCurvesGRN",
    "StateFrequency",
    "StateSpecificManager",
    "TFForceValidation",
    "TFForceWaves",
    "TemporalManager",
    "base",
    "cross_prediction",
    "io",
    "ora",
    "state_specific",
    "temporal",
    "utils",
]
