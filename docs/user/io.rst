I/O
~~~
Locating the data Focal reads, and reading it.

.. module:: focal.io
.. currentmodule:: focal

Dataset locations live in a YAML file next to the notebooks, never inside the
installed package::

    from focal.io import DatasetPaths
    config = DatasetPaths.from_yaml("datasets.yaml")

.. autosummary::
    :toctree: genapi

    io.DatasetPaths
    io.load_lf_gene_colors
    io.load_tf_target_links
    io.load_episode_enrichment
    io.qc_reads
