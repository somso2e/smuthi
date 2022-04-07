#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 10:46:08 2022

This code is used to compile the cython .c file into machine code of the users
specific operating system. This method is the prefered method by Cython for 
distributing Cython code. 
https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html

The user should use this method to compile the existing .c file that will be 
distributed with smuthi. In this way, the user does not need to have Cython 
installed. They can compile using Python's setuptools module. 


Note: It is possible to combine with file with the compile_direct_coupling_block_in_c_from_pyx.py 
file using a "try" statement and adding an additional command line argument during 
compile time (See Cython link). But, I am unsure how the code will be compiled during user download 
of Smuthi. So, I opt to use the easiest code and have the other file exist only
for users who know what they are doing..


@author: parkerwray
"""


from setuptools import Extension, setup
import numpy

setup(
    ext_modules = [Extension("cython_speedups",
			 ["cython_speedups.c"],       
			 extra_compile_args=['-fopenmp'],
        		 extra_link_args=['-fopenmp'],
        		 include_dirs=[numpy.get_include()])]
)


"""
To compile direct_coupling_block_in_c.pyx, do the following:
    1) conda env list
    2) conda activate [env path]
        - Run python scripts with this env.
    3) python compile_direct_coupling_block_in_c_from_c.py build_ext --inplace
        - compile the c code using python to direct the compilation process.

"""
