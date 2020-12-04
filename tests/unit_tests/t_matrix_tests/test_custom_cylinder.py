# -*- coding: utf-8 -*-
"""Test the custom particle t-matrix"""
import unittest
import smuthi.fields as flds
import smuthi.particles as part
import smuthi.linearsystem.tmatrix.nfmds.stlmanager as stlc
import os 
import numpy as np

file_path = os.path.dirname(os.path.realpath(__file__))

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

stlc.convert_stl_to_fem(file_path + "/cylinder.stl", file_path + "/cylinder.fem")

custom_cylinder = part.CustomParticle(position=(0, 0, 0), 
                                     refractive_index=2, 
                                     fem_filename=file_path + "/cylinder.fem", 
                                     scale=100, 
                                     l_max=lmax, 
                                     m_max=mmax, 
                                     n_rank=nrank)

t_axsym = axsym_cylinder.compute_t_matrix(vacuum_wavelength=wl, 
                                          n_medium=nmed)
                                          
t_custom = custom_cylinder.compute_t_matrix(vacuum_wavelength=wl, 
                                            n_medium=nmed)

def test_custom_cylinder():
    print(t_axsym[0, 0])
    print(t_custom[0, 0])

    err = np.linalg.norm(t_axsym - t_custom) / np.linalg.norm(t_axsym)
    print('error t-matrix:', err)
    assert err < 1e-2

if __name__ == '__main__':
    test_custom_cylinder()
