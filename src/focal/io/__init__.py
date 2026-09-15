"""Reading what Focal consumes, and locating it.

:class:`DatasetPaths` replaces the old ``analysis/config.py``: dataset locations
live in a YAML file beside the notebooks, never inside the installed package.
"""

from focal.io._paths import DatasetPaths
from focal.io._qc import qc_reads
from focal.io._readers import (
    load_episode_enrichment,
    load_lf_gene_colors,
    load_tf_target_links,
)

__all__ = [
    "DatasetPaths",
    "load_episode_enrichment",
    "load_lf_gene_colors",
    "load_tf_target_links",
    "qc_reads",
]
