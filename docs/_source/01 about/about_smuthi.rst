About Smuthi
========================

Smuthi stands for 'scattering by multiple particles in thin-film systems'.
It is a Python software that allows to solve light scattering problems involving
one ore multiple particles near or inside a system of planar layer interfaces.

.. image:: drawing.png
   :scale: 40%
   :align: center

It solves the Maxwell equations (3D wave optics) in frequency domain (one wavelength per simulation).


Simulation method
------------------
Smuthi is based on the T-matrix method for the single particle scattering and on the scattering-matrix method
for the propagation through the layered medium.
See :doc:`[Egel 2018] <../literature>` and other publications listed in the :doc:`literature section <../literature>` for a description of the method.

For non spherical particles, Smuthi calls the
`NFM-DS <https://scattport.org/index.php/programs-menu/t-matrix-codes-menu/239-nfm-ds>`_
by Doicu, Wriedt and Eremin to compute the single particle T-matrix. This is a Fortran software package written by
based on the "Null-field method with discrete sources", see :doc:`[Doicu et al. 2006] <../literature>`.

Performance critical parts of the software are implemented in CUDA. When dealing with a large number of particles, Smuthi can benefit from a substantial acceleration if a suitable (NVIDIA) GPU is available.

For CPU-only execution, other acceleration concepts (including MPI parallelization, Numba JIT compilation) are currently tested. 


Range of applications
----------------------

Smuthi can be applied to any scattering problem in frequency domain involving

  - a system of plane parallel layer interfaces separating an arbitrary number of thin metallic or dielectric layers.

  - an arbitrary number of wavelength-scale scattering particles (currently available: spheres, spheroids, finite cylinders, custom particle shapes, anisotropic spheres, layered spheroids). The particles can be metallic or dielectric and rotated to an arbitrary orientation.

  - an initial field in form of a plane wave, a beam (currently available: beam with Gaussian xy-profile) or a collection of dipole sources

Thus, the range of applications spans from scattering by a single particle on a substrate to scattering by several thousand particles inside a planarly layered medium. For a number of examplary simulations, see the :doc:`examples<../04 examples/examples>` section.


Simulation output
------------------

Smuthi can compute

  - the 3D electric and/or magnetic field, for example along a cut plane and save it in the form of ascii data files, png images or gif animations. 

  - the far field power flux of the total field, the initial field or the scattered field. 
    For plane wave excitation, it can be processed to the form of 
    :doc:`differential scattering and extinction cross sections<../03 simulation guidelines/cross_sections>`.

  - For dipole sources, the dissipated power can be computed (Purcell effect).


Current limitations
---------------------

The following issues need to be considered when applying Smuthi:

  - Particles must not intersect with each other or with layer interfaces.
  - Magnetic or anisotropic materials are currently not supported (anisotropic spheres are currently tested).
  - The method is in principle valid for a wide range of particle sizes -  
    however, the numerical validity has only been tested for particle diameters up to around one wavelength.
    For larger particles, note that the number of multipole terms in the spherical wave expansion 
    grows with the particle size. For further details, see the 
    :doc:`hints for the selection of the multipole truncation order <../03 simulation guidelines/simulation_guidelines>`.
  - Particles in a single homogeneous medium (or in free space) can be treated 
    by setting a trivial two layer system with the same refractive index.
    However, Smuthi was not designed for that use case and we believe that 
    there is better software for that case.
  - Smuthi was designed for particles on a substrate or particles near or inside a thin-film system 
    with layer thicknesses of up to a few wavelengths. 
    Simulations involving thick layers might fail or return wrong results due to numerical instability.
    Maybe a more stable algorithm for the layer system response does exist - help is welcome.
  - Smuthi does not provide error checking of user input, nor does it check if 
    numerical parameters specified by the user are sufficient for accurate 
    simulation results. It is thus required that the user develops some 
    understanding of the influence of various numerical parameters on the 
    validity of the results. 
    See the :doc:`simulation guidelines <../03 simulation guidelines/simulation_guidelines>`.
  - A consequence of using the T-matrix method is that the electric field inside the circumscribing
    sphere of a particle cannot be correctly computed, see for example `Auguié et al. (2016) <https://doi.org/10.1088/2040-8978/18/7/075007>`_. 
    In the electric field plots, the circumscribing sphere is displayed as a dashed circle around the particle
    as a reminder that there, the computed near fields cannot be trusted.
  - Particles with initersecting circumscribing spheres can lead to incorrect results. 
    The use of Smuthi is therefore limited to geometries with particles that have disjoint circumscribing spheres.
  - If particles are located near interfaces, such that the circumscribing shere of the particle intersects the 
    interface, a correct simulation result can in principle be achieved. However, special care has to be taken
    regarding the selection of the truncation of the spherical and plane wave expansion, see
    the :doc:`hints for the selection of the wavenumber truncation<../03 simulation guidelines/simulation_guidelines>`. 
  - Dipole sources must not be placed inside the circumscribing sphere of a non-spherical particle (exception: it is OK if the particle is in a different layer)

License
-------

The software is licensed under the `MIT license <https://en.wikipedia.org/wiki/MIT_License>`_.

How to cite this software
-------------------------

If you use SMUTHI for a publication, please consider to cite :doc:`[Egel et al. 2021] <../literature>`.
If you use SMUTHI with non-spherical particles for a publication, please consider to additionally cite :doc:`[Doicu et al. 2006] <../literature>`.
If you use SMUTHI with periodic boundary conditions, please consider to additionally cite :doc:`[Theobald et al. 2021] <../literature>`.



.. _ContactAnchor:

Contact
-------

Email to the author under |emailpic| or to the Smuthi mailing list under smuthi@googlegroups.com for questions, feature requests or if you would like to contribute.

.. |emailpic| image:: email.png


Acknowledgments
---------------

The following persons are/were involved in the Smuthi development: Amos Egel, Dominik Theobald, Krzysztof Czajkowski, Konstantin Ladutenko, Lorenzo Pattelli, Alexey Kuznetsov, Parker Wray.

The authors wish to thank Adrian Doicu, Thomas Wriedt and Yuri Eremin for the
`NFM-DS <https://scattport.org/index.php/programs-menu/t-matrix-codes-menu/239-nfm-ds>`_ package, a copy of which
is distributed with Smuthi.

Ilia Rasskazov, Giacomo Mazzamuto, Fabio Mangini, Refet Ali Yalcin and Johanne Heitmann Solheim have helped with useful comments, bug reports and code additions.

We thank Håkan T Johansson for making his pywigjxpf software availible through PyPi and also under Windows.

The creation of Smuthi was supervised by Uli Lemmer and Guillaume Gomard during the research project
`LAMBDA <http://gepris.dfg.de/gepris/projekt/278746617>`_, funded by the `DFG <http://www.dfg.de/>`_ 
in the priority programme `tailored disorder <http://gepris.dfg.de/gepris/projekt/255652081>`_.
