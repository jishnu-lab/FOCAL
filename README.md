# Focal

**Focal** (Factor-Outcome Coupling for Assessment of Linkages) is an open-source toolkit that combines mechanistic models of networks with interpretable machine learning to focus dense state-specific and dynamic networks onto the regulatory components that govern cell fate decisions. It works on single-cell RNA and ATAC data (sc/snRNA-seq, scATAC-seq), matched or unmatched and provides GIFs (GRNs coupled with Interpretable Factor-outcomes).

![Focal overview](docs/_static/img/fig1_091526.png)

What it does:

- **Prioritized regulatory subnetworks.** Interpretable ML discovers sparse cellular programs (CPs) that separate contrasting cell states, and embeds them in state-specific and cross-state GRNs to surface the transcription factors (TFs) whose regulons are enriched for each program, including in-silico perturbation of the enriched TFs.
- **Phase-resolved dynamic regulation.** Transition-window and episodic GRNs along pseudotime assign TF–target edges to regulatory phases inferred from pseudotemporal clustering, ordering waves of TF regulation and retaining the forces that stay invariant within each episode.
- **Fate predisposition by transfer learning.** CPs learned from fate-switching perturbations (e.g. TF knockouts) versus controls stratify uncommitted populations in unperturbed data by their predicted fate bias.

The full description of all capabilities, the methods, the user and developer API, and the rendered analysis notebooks are on Read the Docs: **[focal-gifs.readthedocs.io](https://focal-gifs.readthedocs.io/)**.

## Installation

Focal runs on top of stock [dictys](https://github.com/pinellolab/dictys) 1.1.0, which comes from
conda (channel `lingfeiwang`), not PyPI, so `pip install .` alone gives a package that cannot be
imported. Build the conda environment from the companion notebooks repository first, then install
Focal into it (linux-64 only):

```bash
git clone --recurse-submodules https://github.com/jishnu-lab/FOCAL Focal
conda env create -f Focal/docs/notebooks/environment.yml   # creates the `fftemporal` env
conda activate fftemporal
pip install -e Focal --no-deps
```

`--no-deps` keeps pip from re-resolving packages conda already placed. The environment covers the
dictys-based dynamic GRN work; the custom dictys network stats Focal needs (`lcpm_tf`,
`ChromatinGRNStat`) ship inside Focal (`focal.backends.dictys`), so dictys itself is not patched.
Optional extras: `.[scanpy]`, `.[celloracle]`, `.[sklearn]` (TF clustering in episodic heatmaps).

## Notebooks

The analysis notebooks live in a companion repository,
[**FOCAL_notebooks**](https://github.com/jishnu-lab/FOCAL_notebooks), and are
pulled in here as a git submodule at `docs/notebooks` so the documentation can render
them. Clone with them:

```bash
git clone --recurse-submodules https://github.com/jishnu-lab/FOCAL
```

If you already cloned without them:

```bash
git submodule update --init docs/notebooks
```

The submodule is optional — the package installs and the API docs build without it.

## Documentation

User guide, API reference and the rendered notebooks are built with Sphinx and hosted on
**Read the Docs**. The build is configured by `.readthedocs.yml`, which also tells Read the Docs to
check out the `docs/notebooks` submodule. Point a project at this repository at
[readthedocs.org](https://readthedocs.org/); with the slug `focal` the site lands at
[https://focal-gifs.readthedocs.io/](https://focal-gifs.readthedocs.io/).

To build the HTML docs locally:

```bash
git submodule update --init docs/notebooks   # optional; without it the API docs still build
pip install ".[docs]"
sphinx-build -b html docs docs/_build/html
```

Notebooks are rendered from their stored outputs and are never executed
(`nb_execution_mode = "off"`), so the build needs no data and no GPU.
