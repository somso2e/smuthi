Sphere on a substrate
=======================

.. image:: drawing.png
   :width: 25%
   :align: right

The "sphere on a substrate" is a minimal Smuthi simulation. 
It investigates scattering by a single glass sphere on a glass substrate and evaluates the total scattering cross section.

Click :download:`here <../../../examples/tutorials/01_sphere_on_substrate/dielectric_sphere_on_substrate.py>` 
to download the Python script.

The simulation output is the total scattering cross section of the sphere.

.. image:: console_screenshot.png
   :width: 75%

In general, a Smuthi simulation script contains the following building blocks:

- Definition of the optical system: the initial field, the layer system and a list of scattering particles are defined
- Definition of the simulation object: the simulation object is initialized with the ingredients of the optical system. Further numerical settings can be applied.
- Simulation start: The calculation is launched with the command `simulation.run()`
- Post processing: The results are processed into the desired output (in our example: scattering cross section).

See also the section :ref:`ProgramStructureAnchor`.


   