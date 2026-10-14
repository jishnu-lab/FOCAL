Installation
============

Python **3.9+** is required. Focal runs on top of stock
`dictys <https://github.com/pinellolab/dictys>`_ 1.1.0, which comes from conda (channel
``lingfeiwang``), not PyPI, so ``pip install .`` alone gives a package that cannot be imported.

Build the conda environment shipped with the companion
`FOCAL_notebooks <https://github.com/jishnu-lab/FOCAL_notebooks>`_ repository (checked out here
as the ``docs/notebooks`` submodule), then install Focal into it:

.. code-block:: bash

   git clone --recurse-submodules https://github.com/jishnu-lab/FOCAL Focal
   conda env create -f Focal/docs/notebooks/environment.yml   # creates the `fftemporal` env
   conda activate fftemporal
   pip install -e Focal --no-deps

``--no-deps`` keeps pip from re-resolving packages conda already placed. The environment pins
exact versions and is linux-64 only. It covers the dictys-based dynamic GRN work; the custom
dictys network stats Focal needs (``lcpm_tf``, ``ChromatinGRNStat``) ship inside Focal
(:mod:`focal.backends.dictys`), so dictys itself is not patched.

Optional dependency groups (defined in ``pyproject.toml``):

.. code-block:: bash

   pip install ".[scanpy]"
   pip install ".[celloracle]"
   pip install ".[sklearn]"   # TF clustering in episodic heatmaps

Documentation tools (also used on Read the Docs):

.. code-block:: bash

   pip install ".[docs]"

Build the HTML docs locally:

.. code-block:: bash

   sphinx-build -W -b html docs docs/_build/html

Then open ``docs/_build/html/index.html`` in a browser.
