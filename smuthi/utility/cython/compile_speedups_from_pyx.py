#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 10:46:08 2022

This code is used to make the .c files and compile the Cython code from the
.pyx file. Note: The user will use the "compile_direct_coupling_block_from_c.c"
file, NOT this file. 


@author: parkerwray
"""


from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        "cython_speedups",
        ["cython_speedups.pyx"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
        include_dirs=[numpy.get_include()]
    )
]

setup(
    name='cython_speedups',
    ext_modules=cythonize(ext_modules),
)



"""
To compile cython_speedups.pyx, do the following:
    1) conda env list
        - This list paths to enviroments. Find one that has Cython installed
    2) conda activate [env path]
        - Run python scripts with this env.
    3) python compile_speedups_from_pyx.py build_ext --inplace
        - compile the c code using python to direct the compilation process.

NOTE: The function will first compile the .pyx code to a .c code. The function 
will then compile the .c code, using a c compiler, to a .so that can be run. 
The .so file is saved in the "build" folder. This folder houses the operating system
specific executables. If you want, you can copy the .so file to the main location 
to make importing the code easier.  

"""