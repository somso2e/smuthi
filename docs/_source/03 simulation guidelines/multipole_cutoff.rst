.. _MultipoleCutOffAnchor:

Multipole cut-off
-----------------
The scattering properties of each particle are represented by its T-matrix :math:`T_{plm,p'l'm'}`
where :math:`plm` and :math:`p'l'm'` are the multipole polarization, degree and order of the scattered
and incoming field, respectively, see sections 3.3 and 2.3.2 of :doc:`[Egel 2018] <../literature>`.
In practice, the T-matrix is truncated at some multipole degree :math:`l_{max} \ge 1` and order
:math:`0 \le m_{max} \le l_{max}` to obtain a finite system of linear equations.

Specify the cut-off parameters for each particle like this::

   large_sphere = smuthi.particles.Sphere( ...
                                           l_max=10,
                                           m_max=10,
                                           ...)
																					 
	 
   small_sphere = smuthi.particles.Sphere( ...
                                           l_max=3,
                                           m_max=3,
                                           ...)

In general, we can say:

 - Large particles require higher multipole orders than small particles.

 - Particles very close to each other, very close to an interface or very close to a point dipole
   source require higher multipole orders than those that stand freely.

 - Larger multipole cutoff parameters imply better accuracy, but also a quickly growing numerical effort.

 - When simulating flat particles near planar interfaces, the multipole truncation should be chosen with regard to the Sommerfeld integral truncation. See :doc:`[Egel et al. 2017] <../literature>`.

Literature offers various rules of thumb for the selection of the multipole truncation in the
case of spherical particles, see for example :doc:`[Neves 2012] <../literature>` 
or :doc:`[Wiscombe 1980] <../literature>`.

Otherwise, you can use Smuthi's built-in automatic parameter selection feature 
to estimate a suitable multipole truncation, see section on :ref:`AutoParamAnchor`.
