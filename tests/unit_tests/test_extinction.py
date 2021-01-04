import unittest
import numpy as np
import smuthi
import smuthi.simulation
import smuthi.initial_field
import smuthi.layers
import smuthi.particles
import smuthi.postprocessing.far_field as ff

class TestExtinction(unittest.TestCase):
    def test_lower_zlimit(self):
        radius = 100
        l_max=5

        two_layers = smuthi.layers.LayerSystem(thicknesses=[0, 0],
                                       refractive_indices=[1.52, 1])

        # Scattering particle
        first_sphere = smuthi.particles.Sphere(position=[0, 0, radius],
                                 refractive_index=4 + 0.5j,
                                 radius=radius,
                                 l_max=l_max)

        second_sphere = smuthi.particles.Sphere(position=[0, 0, 3*radius],
                                 refractive_index=4 + 0.5j,
                                 radius=radius,
                                 l_max=l_max)

        # list of all scattering particles (only one in this case)
        two_spheres = [first_sphere, second_sphere]

        # Initial field
        plane_wave = smuthi.initial_field.PlaneWave(vacuum_wavelength=950,
                                            polar_angle=np.pi,    # from top
                                            azimuthal_angle=0,
                                            polarization=0)       # 0=TE 1=TM

        # Initialize and run simulation
        simulation = smuthi.simulation.Simulation(layer_system=two_layers,
                                            particle_list=two_spheres,
                                            initial_field=plane_wave)
        simulation.run()

        general_transmitted_ecs = ff.extinction_cross_section(initial_field=plane_wave,
                                                              particle_list=two_spheres,
                                                              layer_system=two_layers,
                                                              extinction_direction='transmission')
        general_reflected_ecs = ff.extinction_cross_section(initial_field=plane_wave,
                                                            particle_list=two_spheres,
                                                            layer_system=two_layers,
                                                            extinction_direction='reflection')

        ecs_bottom=0  
        ecs_top=0                                      
        for l in range(1, l_max+1):
            for tau in range(2):
                for pol in range(2):
                    ecs_top += ff.extinction_cross_section(initial_field=plane_wave,
                                                           particle_list=two_spheres,
                                                           layer_system=two_layers,
                                                           only_l=l, only_tau=tau, only_pol=pol,
                                                           extinction_direction='reflection')
                    ecs_bottom += ff.extinction_cross_section(initial_field=plane_wave,
                                                              particle_list=two_spheres,
                                                              layer_system=two_layers,
                                                              only_l=l, only_tau=tau, only_pol=pol,
                                                              extinction_direction='transmission')

        top_difference = np.abs(ecs_top - general_reflected_ecs)
        bottom_difference = np.abs(ecs_bottom - general_transmitted_ecs)
        
        self.assertEqual(top_difference < 10**(-9), True)
        self.assertEqual(bottom_difference < 10**(-9), True)


if __name__ == '__main__':
    unittest.main()
