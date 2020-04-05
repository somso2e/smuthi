Physical units
--------------
Smuthi is commited to a "relative units" philosophy. 
That means, all quantities have only relative meaning.

Length units
~~~~~~~~~~~~
The user is free to select the unit in which all lengths are provided. 
Just make sure that particle sizes, layer thicknesses and wavelengths are all specified in the same unit.
Results will automatically refer to the same unit. 
For example, if you specify the wavelength in nanometers, resulting cross sections will be in square nanometers.
Besides, quantities with an inverse length dimension (wavenumbers) also implicitly refer to the selected length unit.

Field strength units
~~~~~~~~~~~~~~~~~~~~
When the electromagnetic fields are computed, their absolute value has no physical meaning.
Only relative quantities can be used for further analysis. 
For example, the scattered field strength divided by the amplitude of the initial field *does* have a physical meaning.


Power units
~~~~~~~~~~~
Also power units have no meaning as absolute values. 
To get meaningful information, power-related figures always need to be guarded in reference to other power-related figures. 
Some examples: 
 - Scattering cross section as the quotient of scattered (angular) intensity and incident (power-per-area) intensity.
 - Diffuse reflectivity as the total back scattered far field power divided by the initial Gaussian beam power.
 - Purcell factor as the dissipated power of a dipole source divided by the dissipated power of the same source in the absence of planar interfaces and scattering particles.

 
