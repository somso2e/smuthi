Plane wave coupling
===================

.. warning:: The plane-wave coupling module is still in development, and its current functionality is experimental.

If non-spherical particles are located such that their circumscribing spheres overlap, the conventional superposition T-matrix method is not applicable. A coupling method based on a temporary plane wave expansion of the scattered field was developed [Theobald 2017] in order to allow for simulations also in such cases.

Right now, the plane wave coupling can be used if

- direct matrix inversion is selected, see :ref:`SolverSettingsAnchor`.
- :code:`m_max` is set to :code:`l_max` for all particles.


Use PVWF coupling in a Smuthi simulation
----------------------------------------

To activate the PVWF coupling feature in a Smuthi simulation, set the :code:`use_pvwf_coupling` parameter of the simulation constructor to :code:`True` and provide a suitable :math:`n_{eff}` truncation and discretization with the :code:`pvwf_coupling_neff_max` parameter and the :code:`pvwf_coupling_neff_resolution` parameter of the simulation constructor::

  simul = smuthi.simulation.Simulation( ...,
                                        use_pvwf_coupling=True,
                                        pvwf_coupling_neff_max=7,
                                        pvwf_coupling_neff_resolution=1e-2)
