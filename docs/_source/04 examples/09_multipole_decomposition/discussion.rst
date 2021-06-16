:orphan:

Multipole decomposition
=======================

.. image:: extinction.svg
   :width: 90%
   :align: center

Smuthi allows to analyze the contribution of individual multipole moments to the overall extinction cross section.
Click :download:`here <../../../../examples/tutorials/09_multipole_decomposition/decompose_extinction_for_one_sphere.py>` 
to download an example script which demonstrates this use case. It reproduces the results from 
I.Sinev et al. "Polarization control over electric and magnetic dipole resonances of dielectric nanoparticles on metallic films" [1].
In the mentioned example you can see decomposition of dipole extinction and conversion from spherical to Cartesian coordinates. 
For this conversion we applied the following formulas 
from C. Rockstuhl et al. "An electromagnetic multipole expansion beyond the long-wavelength approximation" [2]:
	.. math::
		\mathbf{\hat{e}}_{1} = -\frac{\mathbf{\hat{x}} + i\mathbf{\hat{y}}}{\sqrt{2}}, \\
		\mathbf{\hat{e}}_{0} = \mathbf{\hat{z}}, \\
		\mathbf{\hat{e}}_{-1} = -\frac{\mathbf{\hat{x}} - i\mathbf{\hat{y}}}{\sqrt{2}}.
Hence 
	.. math::
		\mathit{p}_{x}^{\omega} = \frac{\mathit{a}_{1-1}^{\omega} - \mathit{a}_{11}^{\omega}}{\sqrt{2}}, \\
		\mathit{p}_{y}^{\omega} = \frac{\mathit{a}_{1-1}^{\omega} + \mathit{a}_{11}^{\omega}}{\sqrt{2\mathit{i}}}, \\
		\mathit{p}_{z}^{\omega} = \mathit{a}_{0}^{\omega}
and
	.. math::
		\mathit{m}_{x}^{\omega} = \frac{\mathit{a}_{1-1}^{\omega} - \mathit{a}_{11}^{\omega}}{\sqrt{2}}, \\
		\mathit{m}_{y}^{\omega} = \frac{\mathit{a}_{1-1}^{\omega} + \mathit{a}_{11}^{\omega}}{\sqrt{2\mathit{i}}}, \\
		\mathit{m}_{z}^{\omega} = \mathit{a}_{0}^{\omega}

Formulas for a general multipole also can be found in this paper.

| [1] Laser Photonics Rev. 10, No. 5, 799–806 (2016), http://dx.doi.org/10.1002/lpor.201600055
| [2] Optics Communications Volume 407, 15 January 2018, Pages 17-21, Appendix B, https://doi.org/10.1016/j.optcom.2017.08.064
