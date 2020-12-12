# -*- coding: utf-8 -*-
"""Test the custom particle t-matrix"""
import unittest
import smuthi.fields as flds
import smuthi.particles as part
import smuthi.linearsystem.tmatrix.nfmds.stlmanager as stlc
import os 
import numpy as np

file_path = os.path.dirname(os.path.realpath(__file__))


lmax = 3
mmax = 3
nrank = 10

nmed = 1
wl = 500

sphere = part.Sphere(position=(0, 0, 0), 
                     refractive_index=2, 
                     radius=100,
                     l_max=lmax, 
                     m_max=mmax)

custom_sphere = part.CustomParticle(position=(0, 0, 0),
                                    refractive_index=2,
                                    geometry_filename=file_path + "/sphere_fine.stl",
                                    scale=100,
                                    l_max=lmax,
                                    m_max=mmax,
                                    n_rank=nrank)

t_mie = sphere.compute_t_matrix(vacuum_wavelength=wl, n_medium=nmed)
                                          
t_custom = custom_sphere.compute_t_matrix(vacuum_wavelength=wl, n_medium=nmed)


def test_custom_sphere():
    print(t_mie[0, 0])
    print(t_custom[0, 0])

    err = np.linalg.norm(t_mie - t_custom) / np.linalg.norm(t_mie)
    print('error t-matrix:', err)
    assert err < 1e-2


if __name__ == '__main__':
    test_custom_sphere()
