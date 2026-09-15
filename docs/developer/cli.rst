Command line
~~~~~~~~~~~~
One module per pipeline step, each runnable as ``python -m focal.cli.<name>``. The
SLURM scripts under ``multiome_dynamic_regulation/bash_scripts/`` invoke them this way.

.. code-block:: bash

    python -m focal.cli.expression_to_tsv MATRIX_DIR expression.tsv.gz
    python -m focal.cli.validate_inputs --dir_data data --dir_makefiles makefiles
    python -m focal.cli.reconstruct_networks START END WORK_DIR
    python -m focal.cli.motifs_to_homer IN.meme OUT.motif

.. automodule:: focal.cli.expression_to_tsv
    :members:

.. automodule:: focal.cli.validate_inputs
    :members:

.. automodule:: focal.cli.reconstruct_networks
    :members:

.. automodule:: focal.cli.motifs_to_homer
    :members:
