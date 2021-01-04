import sys
import subprocess

# -*- coding: utf-8 -*-
"""Test the python only t-matrix"""
import unittest
import smuthi.fields as flds
import smuthi.particles as part
import smuthi.linearsystem.tmatrix.nfmds.stlmanager as stlc
import os
import numpy as np

# We compute the T-matrix of a cylinder with 100nm height and 50nm radius

lmax = 3
mmax = 3
nrank = 10

nmed = 1
wl = 500

axsym_cylinder = part.FiniteCylinder(position=(0, 0, 0),
                                     refractive_index=2,
                                     cylinder_radius=50,
                                     cylinder_height=100,
                                     l_max=lmax,
                                     m_max=mmax,
                                     n_rank=nrank)

t_axsym = axsym_cylinder.compute_t_matrix(vacuum_wavelength=wl,
                                          n_medium=nmed)

axsym_cylinder.python_tmt = True
axsym_cylinder.nrank = 10

t_python = axsym_cylinder.compute_t_matrix(vacuum_wavelength=wl,
                                           n_medium=nmed)


def test_python_cylinder():
    pass # test_python_cylinder()  <- uncomment when Alan's code is on PyPi
    print(t_axsym[0, 0])
    print(t_python[0, 0])

    err = np.linalg.norm(t_axsym - t_python) / np.linalg.norm(t_axsym)
    print('error t-matrix:', err)
    assert err < 1e-3


if __name__ == '__main__':
    test_python_cylinder()  <- uncomment when Alan's code is on PyPi
