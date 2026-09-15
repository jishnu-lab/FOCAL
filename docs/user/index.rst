User API
########

.. module:: focal.user

The public surface of the three Focal modules. Import :mod:`focal` as::

    import focal as ff

Each module's entry point is its manager: :class:`~focal.temporal.TemporalManager`,
:class:`~focal.state_specific.StateSpecificManager`. Figures are not a separate
namespace — every plotting function is exported by the module whose results it draws.

.. toctree::
    :maxdepth: 2

    temporal
    state_specific
    cross_prediction
    io
