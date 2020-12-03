Installation
=============

We **recommend to use Linux operating systems** to run Smuthi. Otherwise, Smuthi can run on Windows, too, but issues regarding dependencies or performance are more likely.

Installing Smuthi under Ubuntu (recommended)
--------------------------------------------

Prerequisites
~~~~~~~~~~~~~

`python3` with `pip`, `gfortran` and `gcc` usually are shipped with the operating system. However, Smuthi requires a Python version of 3.6 or newer. Check the installed Python version by::

  python3 --version
	
If the version is 3.5 or less, please install a newer Python version. You can have multiple Python versions installed in parallel. Depending on your  configuration, you might need to replace the command :code:`python3` in the below by the command that belongs to the newly installed Python, e.g. :code:`python3.8`.

Make sure that the Foreign Function Interface library is available (needed for pywigxjpf)::

  sudo apt-get install libffi6 libffi-dev

Installation
~~~~~~~~~~~~

To install Smuthi from PyPi, simply type::

  sudo python3 -m pip install smuthi

Alternatively, you can install it locally from source (see below section :ref:`local_install`).

Installing Smuthi under Windows
-------------------------------

Prerequisites
~~~~~~~~~~~~~

First make sure that a 64 Bit Python 3.6 or newer is installed on your computer. 
You can install for example 
`Anaconda <https://www.continuum.io/downloads>`_ 
or `WinPython <https://winpython.github.io/>`_ 
to get a full Python environment.

C compiler under Windows (for pywigxjpf, optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. note:: 
	If you skip the installation of a C compiler, Smuthi will still work, but you might not benefit from fast evaluation of Wigner3j symbols through the wigxjpf libary. This can be an issue when your simulation involves large multipole degrees.

To benefit from faster evaluation of Wigner3j symbols through the pywigxjpf package, you need a C compiler.
If you have Microsoft Visual Studio installed, `MS VC` is probably already there. Otherwise, open the Visual Studio setup and install the Visual C compiler. If you don't have Microsoft Visual Studio, see 
`the Python Wiki <https://wiki.python.org/moin/WindowsCompilers>`_ 
for further instructions.

.. _gfortranAnchor:

gfortran under Windows (for NFM-DS, optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. note:: 
	Smuthi can be installed from PyPi through binary wheels which include a compiled binary of the NFM-DS module, so if you install Smuthi with `python -m pip install smuthi`, you don't need a Fortran compiler. The following information is only relevant if you plan to install Smuthi from source.

Visit the `MinGW getting started page <http://mingw.org/wiki/Getting_Started>`_ and follow the instructions to install `gfortran`. 
Also make sure to add the bin folder of your MinGW installation to the Windows PATH variable. See `Environment Settings` section of the `MinGW getting started page <http://mingw.org/wiki/Getting_Started>`_ for instructions.

.. note::
  The MinGW version needs to fit to your Python installation. If you have 64 Bit Python, make sure to download a `Mingw-64 <https://sourceforge.net/projects/mingw-w64/>`_

Installation
~~~~~~~~~~~~

Open a command window and type::

  python -m pip install smuthi

Depending on where pip will install the package, you might need administrator rights for that.

Alternatively, install locally from source (see below section :ref:`local_install`).


.. _local_install:

Installing Smuthi from source
-----------------------------

This option allows to install a non-release version of Smuthi or to modify the source code and then run your custom version of Smuthi.

Ubuntu
~~~~~~
Clone Smuthi and install it locally by::

  git clone https://gitlab.com/AmosEgel/smuthi.git
  cd smuthi/
  sudo python3 -m pip install -e .

Windows
~~~~~~~
`Download <https://gitlab.com/AmosEgel/smuthi/tags>`_ or git clone the Smuthi project folder from the `gitlab repository <https://gitlab.com/AmosEgel/smuthi.git>`_. Open a command prompt and change directory to the Smuthi
project folder and enter::

  python -m pip install -e .


Verification
~~~~~~~~~~~~

After installation from source you can check the unit tests:

Ubuntu::

  sudo python3 -m pip install nose2
  nose2

Windows::

  python -m pip install nose2
  nose2


.. _GPUAnchor:

GPU-acceleration (optional)
---------------------------
.. note:: 
	PyCuda support is recommended if you run heavy simulations with many particles. In addition, it can speed up certain post processing steps like the evaluation of the electric field on a grid of points, e.g. when you create images of the field distribution. 
	For simple simiulations involving one particle on a substrate, you might well go without.

If you want to benefit from fast simulations on the GPU, you need:

* A CUDA-capable NVIDIA GPU
* The `NVIDIA CUDA toolkit <https://developer.nvidia.com/cuda-toolkit>`_ installed
* PyCuda installed

Under Ubuntu, install PyCuda simply by::

  sudo python3 -m pip install pycuda

Under Windows, installing PyCuda this is not as straightforward as under Linux.
There exist prebuilt binaries on `Christoph Gohlke's homepage <https://www.lfd.uci.edu/~gohlke/pythonlibs/#pycuda>`_. 
See for example `these instructions <https://www.ibm.com/developerworks/community/blogs/jfp/entry/Installing_PyCUDA_On_Anaconda_For_Windows?lang=en>`_ 
for the necessary steps to get it running. 


Troubleshooting
---------------

Windows: Unable to import the nfmds module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Try to install Smuthi from source. You will need the `gfortran` compiler, see :ref:`gfortranAnchor`.