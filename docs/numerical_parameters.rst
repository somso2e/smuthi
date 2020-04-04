.. _NumericalParametersAnchor:

Numerical parameters
====================

We want to explain the meaning of the parameters that control Smuthi's performance regarding accuracy and runtime.

.. _MultipoleCutOffAnchor:

Multipole cut-off
-----------------
The scattering properties of each particle are represented by its T-matrix :math:`T_{plm,p'l'm'}`
where :math:`plm` and :math:`p'l'm'` are the multipole polarization, degree and order of the scattered
and incoming field, respectively, see sections 3.3 and 2.3.2 of :doc:`[Egel 2018] <literature>`.
In practice, the T-matrix is truncated at some multipole degree :math:`l_{max} \ge 1` and order
:math:`0 \le m_{max} \le l_{max}` to obtain a finite system of linear equations.

In general we can say:

 - Large particles require higher multipole orders than small particles.

 - Particles very close to each other, very close to an interface or very close to a point dipole
   source require higher multipole orders than those that stand freely.

 - Larger multipole cutoff parameters imply better accuracy, but also a quickly growing numerical effort.

Literature offers various rules of thumb for the selection of the multipole truncation in the
case of spherical particles, see for example :doc:`[Neves 2012] <literature>` 
or :doc:`[Wiscombe 1980] <literature>`.

Otherwise, you can use Smuthi's built-in automatic parameter selection feature 
to estimate a suitable multipole truncation 

.. _SommerfeldParametersAnchor:

Sommerfeld integral contour
---------------------------

Sommerfeld integrals arise in the treatment of the layer system response to the scattered field. 
Their numerical evaluation relies on an integral contour that is deflected into the complex plane in order to avoid
sharp features stemming from waveguide mode singularities (see :ref:`SommerfeldAnchor`).

.. image:: images/contour.png
   :scale: 70%
   :align: center

This deflected contour is typically described by the following three parameters:

- `neff_max`
- `neff_imag`
- `neff_resolution`

Contour truncation
~~~~~~~~~~~~~~~~~~

The contour truncation scale `neff_max` is a real number which specifies where the contour ends.
It should be larger than the refractive index of the layer in which the particle resides. The offset :math:`n_\mathrm{eff}-n`
should be chosen with regard to the distance between the particles (and point sources) to the next layer interface.
If that distance is large, the truncation scale is uncritical, whereas whereas point sources or particles whose
center is very close to a layer interface require a larger offset.
	
At a :math:`z`-distance of :math:`\Delta z`, evanescent waves with an effective refractive index of 
:math:`n_\mathrm{eff}` are damped by a factor of 

.. math:: \exp\left(2\pi\mathrm{i}\frac{\Delta z}{\lambda} \sqrt{n_\mathrm{eff}^2-n^2}\right),	

where :math:`\lambda` is the vacuum wavelength and :math:`n` is the refractive index of the medium.
	 
.. image:: images/delta_z.png
   :scale: 50%
   :align: center
	 
	 
To select a reasonable `neff_max`, we should consider that the shortest possible interaction path is *twice* the :math:`z`-distance between some particle center (or dipole position) and the next layer interface.
	 
	
.. admonition:: Uncritical example

   A layer system consists of a substrate (:math:`n=1.5`), covered with a 1000nm thick layer of titania (:math:`n=2.1`) under air (:math:`n=1`).
   A silica sphere is immersed in the middle of the titania layer. The system is illuminated with a plane wave at vacuum wavelength of 550nm.
	 
   Then, :math:`\Delta z= 2\times 500\mathrm{nm}` such that evanescent waves with :math:`n_\mathrm{eff}=2.3` are already damped by a factor of
   :math:`\exp(-2\pi \frac{1000\mathrm{nm}}{550\mathrm{nm}} \sqrt{(2.3^2-2.1^2)}) \approx 2\times 10^{-5}` when they propagate to the layer interface and back to the sphere.
   Waves beyond that effective refractive index thus can be safely neglected in the particle-layer system interaction, such that a truncation parameter of :math:`n_\mathrm{eff, max}=2.3` is reasonable.

.. admonition:: Critical example

   A layer system consists of a substrate (:math:`n=1.5`), under air (:math:`n=1`).
   A point dipole source of wavelength 550nm is located 10nm above the substrate/air interface.
	 
   Here we need to consider :math:`\Delta z= 2\times 10\mathrm{nm}` such that Then, evanescent waves with 
   :math:`n_\mathrm{eff}=2.3` are only damped by a factor of
   :math:`\exp(-2\pi \frac{20nm}{550nm} \sqrt{(2.3^2-1^2)}) \approx 0.62` when scattered by the layer interface.
   Even a truncation of :math:`n_\mathrm{eff, max}=10` would only lead to an evanescent damping of
   :math:`\exp(-2\pi \frac{20nm}{550nm} \sqrt{(10^2-1^2)}) \approx 0.1` which might still not be enough.
   
   


.. todo:: explain

Particle coupling lookup parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo:: explain resolution and interpolation order

Near field calculation parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo:: explain discretization and cut-off for the plane wave expansion for the near field calculation

Far field discretization
~~~~~~~~~~~~~~~~~~~~~~~~
 
.. todo:: explain discretization of the far field in direction space

General strategies to select numerical parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Empirical tests check the consistency of simulation results for a given numerical setting.
By "consistency" we mean for example the agreement to accurate benchmark results that can be
analytical results, results from other software or Smuthi results for a more accurate setting.

In certain cases, we can also check how accurately energy is conserved as a consistency criterion.
However, this criterion is suited only for certain numerical parameters.
In other cases, it is misleading.

.. note::
  In certain cases, even inaccurate simulations yield an accurately conserved energy. This will happen for
  example in case of a too small multipole truncation.

.. note::
   Smuthi currently supports only the evaluation of optical power in the far field 
   (and, in addition, the dissipated power of dipole sources).
   Therefore, it is only possible to check the conservation of energy in systems with no absorbing materials and no waveguiding.

Rules of thumb on the other hand can stem from heuristical reasoning or represent former experience.
They can be fit formulae to earlier results from empirical tests,
see for example :doc:`[Wiscombe 1980] <literature>` or :doc:`[Neves 2012] <literature>` for the selection of multipole truncation
or :doc:`[Egel2017] <literature>` for the truncation of Sommerfeld integrals.