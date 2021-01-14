import unittest
import numpy as np
import smuthi.linearsystem.tmatrix.nfmds.nfmds as nfmds
import smuthi.linearsystem.tmatrix.nfmds.indexconverter as nfic

vacuum_wavelength = 550
layer_refractive_index = 1.3
particle_refractive_index = 1.8 + 0.01j
half_axis_z = 100
half_axis_xy = 200
use_ds = True
n_int = 300
n_rank = 8
l_max = 4

cylinder_height = 100
cylinder_radius = 200


class TestNFMDSWrapper(unittest.TestCase):
    def test_spheroid_tmatrix_against_prototype(self):
        typegeom = 1
        surf = [half_axis_z, half_axis_xy]
        
        surf=np.array(surf)
        Mrank = n_rank
        Nmax = n_rank + Mrank * (2 * n_rank - Mrank + 1)
        nparam = 1
        """
        t_s = taxs.tmatrix_spheroid(vacuum_wavelength=vacuum_wavelength,
                                    layer_refractive_index=layer_refractive_index,
                                    particle_refractive_index=particle_refractive_index,
                                    semi_axis_c=half_axis_z, semi_axis_a=half_axis_xy, use_ds=use_ds,
                                    nint=n_int, nrank=n_rank, l_max=4, m_max=4)
        """                            
        tnfmds = nfmds.taxsym(surf, Nmax, typegeom=typegeom, nparam=nparam,
                              wavelength=vacuum_wavelength, ind_refrel=particle_refractive_index/layer_refractive_index,
                              nrank=n_rank,nint=n_int,ind_refmed=layer_refractive_index,ds=use_ds)
        t_s = nfic.nfmds_to_smuthi_matrix(tnfmds,l_max=l_max)

        t00 = -0.416048522578639 + 0.462839918856895j
        self.assertAlmostEqual(t_s[0, 0], t00, places=5)

        t4210 = -4.663469643281004e-04 - 3.215630661547245e-04j
        self.assertAlmostEqual(t_s[42, 10], t4210, places=5)

    def test_layered_sphere_tmatrix_against_prototype(self):
        ind_host = 4+0j
        ind_shell = 2+0j
        R=100
        H=10
        wavelength = 847.5
        ref = 10.35853 # reference value from Kostya
        area = np.pi*(R+H)**2
        surf=np.array([[R+H,R+H],[R,R]])
        zpart=np.array([0,0])
        nparam=np.array([1,1])
        nsurf=np.array([2,2])
        ind_ref=np.array([ind_shell,ind_host])
        nrankpmax=5
        l_max=2

        k=2*np.pi/wavelength
        tnfmds = nfmds.tlay(k,ind_ref,surf,nrankpmax)
        diagt = np.diag(tnfmds)
        n=1
        qscat = 2*np.pi/k**2*(2*n+1)*np.sum(np.abs(diagt)**2)/area/3 #orientation averaged scattering efficiency
        err=np.abs(qscat-ref)/ref
        self.assertAlmostEqual(err, 0.0, places=3)

    def test_cylinder_tmatrix_against_prototype(self):
        typegeom = 2    
        nparam = 3
        surf = [cylinder_height / 2, cylinder_radius]
        
        surf=np.array(surf)
        Mrank = n_rank
        Nmax = n_rank + Mrank * (2 * n_rank - Mrank + 1)
        tnfmds = nfmds.taxsym(surf, Nmax, typegeom=typegeom, nparam=nparam,
                              wavelength=vacuum_wavelength, ind_refrel=particle_refractive_index/layer_refractive_index,
                              nrank=n_rank,nint=n_int,ind_refmed=layer_refractive_index,ds=use_ds)
        """
        t_c = taxs.tmatrix_cylinder(vacuum_wavelength=vacuum_wavelength,
                                    layer_refractive_index=layer_refractive_index,
                                    particle_refractive_index=particle_refractive_index,
                                    cylinder_height=cylinder_height, cylinder_radius=cylinder_radius,
                                    use_ds=use_ds, nint=n_int, nrank=n_rank, l_max=4, m_max=4)
        """
        t_c = nfic.nfmds_to_smuthi_matrix(tnfmds,l_max=l_max)

        t00 = -0.119828956584570 + 0.282351044628953j
        self.assertAlmostEqual(t_c[0, 0], t00, places=5)

        t4210 = -0.001017849863151 - 0.000754036833086j
        self.assertAlmostEqual(t_c[42, 10], t4210, places=5)


if __name__ == '__main__':
    unittest.main()
