"""
This module contains functions to compute the direct (i.e., not layer 
mediated) particle coupling coefficients.
"""

import numpy as np
import smuthi.fields as flds
import smuthi.fields.transformations as trf
from smuthi.fields.spherical_harmonic_translations import svwf_translate
import smuthi.utility.math as sma
import smuthi.utility.memoizing as memo
import scipy.optimize
import scipy.special
import sys
import warnings


def direct_coupling_block(vacuum_wavelength, receiving_particle, emitting_particle, layer_system):
    r"""Direct particle coupling matrix :math:`W` for two particles that do not have intersecting circumscribing spheres.
       This routine is explicit.
       
       To reduce computation time, this routine relies on two internal accelerations. 
       First, in most cases the number of unique maximum multipole indicies,
       :math:`(\tau, l_{max}, m_{max})`, is much less than the number of unique particles. 
       Therefore, all calculations that depend only on multipole indicies are stored in an 
       intermediate hash table. Second, Cython acceleration is used by default to leverage 
       fast looping. If the Cython files are not supported, this routine will 
       fall back on equivalent Python looping.
       
       Cython acceleration can be between 10-1,000x faster compared to the Python 
       equivalent. Speed variability depends on the number of unique multipoles indicies,
       the size of the largest multipole order, and if particles share the same z coordinate.
       

    Args:
        vacuum_wavelength (float):                          Vacuum wavelength :math:`\lambda` (length unit)
        receiving_particle (smuthi.particles.Particle):     Particle that receives the scattered field
        emitting_particle (smuthi.particles.Particle):      Particle that emits the scattered field
        layer_system (smuthi.layers.LayerSystem):           Stratified medium in which the coupling takes place

    Returns:
        Direct coupling matrix block as numpy array.
    """
    
    # Note: Direct coupling is transpose of translation. 
    # The translation convention is emit = 1 -> recieve = 2
    # The transpose of this is emit = 2 -> recieve = 1
    lmax1 = int(receiving_particle.l_max) 
    mmax1 = int(receiving_particle.m_max)
    lmax2 = int(emitting_particle.l_max)
    mmax2 = int(emitting_particle.m_max)    

    
    # Check if particles are in the same layer.
    rS1 = receiving_particle.position
    rS2 = emitting_particle.position
    iS1 = layer_system.layer_number(rS1[2])
    iS2 = layer_system.layer_number(rS2[2])
    
    if (iS1 == iS2 or layer_system.is_degenerate()) and not emitting_particle == receiving_particle:
        # Initialize variables from abstract classes to be passed as simple 
        # data types
        omega = flds.angular_frequency(vacuum_wavelength)
        k = complex(omega * layer_system.refractive_indices[iS1])
        d = [rS1[i]-rS2[i] for i in range(3)]
        return svwf_translate(k, d, lmax1, mmax1, lmax2, mmax2, kind = 'outgoing to regular', threaded = False)

    else:
        blocksize1 = int(flds.blocksize(lmax1, mmax1)) # Rows of translation matrix
        blocksize2 = int(flds.blocksize(lmax2, mmax2)) # Cols of translation matrix
        # Note: Direct coupling returns zeros such that you recieve 0 on multiplication (no coupling to self..)
        # Note: Smuthi assumes that if you try to calculate direct coupling of two particles across an interface
        # that the function returns zero instead of an exception! This allows for lazy functions that do not 
        # check for the validity of direct coupling before calling. 
        return np.zeros((blocksize1, blocksize2), dtype = np.complex128, order = 'C')
    
    # To Do: lazy evaluation of direct coupling could be replaced with an explicit check?
    #raise Exception('''
        #                You are attempting to couple particles directly that exist in different layers. 
        #                Direct coupling (2D or 3D) is only sensible for particles in the same layer.
        #                You need to first deal with the effect of the layer interface!
        #                ''')


        dx = float(rS1[0] - rS2[0])
        dy = float(rS1[1] - rS2[1])
        dz = float(rS1[2] - rS2[2])

        #Note: It is allways better to use the hash tables so just do it. 
        if dz == 0:
            a5leg_array, b5leg_array = ab5_coefficient_and_legendre_hash_table(lmax1, lmax2, mmax1, mmax2)
            w = direct_coupling_block_2D_from_hash_table(blocksize1, blocksize2,
                                                        w, sph, k, dx, dy,dz,
                                                        lmax1, lmax2, mmax1, mmax2,
                                                        a5leg_array, b5leg_array,                              
                                                        threaded = False)
        else:
            a5_array, b5_array = ab5_coefficient_hash_table(lmax1, lmax2, mmax1, mmax2)
            w = direct_coupling_block_3D_from_hash_table(blocksize1, blocksize2,
                                                        w, sph, k, dx, dy,dz,
                                                        lmax1, lmax2, mmax1, mmax2,
                                                        a5_array, b5_array,                              
                                                        threaded = False)           
        
    return w



def direct_coupling_matrix(vacuum_wavelength, particle_list, layer_system):
    """Return the direct particle coupling matrix W for a particle collection in a layered medium.

    Args:
        vacuum_wavelength (float):                                  Wavelength in length unit
        particle_list (list of smuthi.particles.Particle obejcts:   Scattering particles
        layer_system (smuthi.layers.LayerSystem):                   The stratified medium
   
    Returns:
        Ensemble coupling matrix as numpy array.
    """
    # indices
    blocksizes = [flds.blocksize(particle.l_max, particle.m_max)
                  for particle in particle_list]

    # initialize result
    w = np.zeros((sum(blocksizes), sum(blocksizes)), dtype=complex)

    for s1, particle1 in enumerate(particle_list):
        idx1 = np.array(range(sum(blocksizes[:s1]), sum(blocksizes[:s1+1])))
        for s2, particle2 in enumerate(particle_list):
            idx2 = range(sum(blocksizes[:s2]), sum(blocksizes[:s2+1]))
            w[idx1[:, None], idx2] = direct_coupling_block(vacuum_wavelength, particle1, particle2, layer_system)

    return w


try: 
    from smuthi.utility.cython.cython_speedups import _ab5_coefficient_and_legendre_hash_table
    @memo.Memoize
    def ab5_coefficient_and_legendre_hash_table(lmax1, lmax2, mmax1, mmax2):
        return _ab5_coefficient_and_legendre_hash_table(int(lmax1), int(lmax2), int(mmax1), int(mmax2))
except:
    sys.stdout.write("Cython acceleration could not be loaded. \nFalling back on Python equivalents...")
    sys.stdout.flush()
    @memo.Memoize
    def ab5_coefficient_and_legendre_hash_table(lmax1, lmax2, mmax1, mmax2):
        r"""Creates a hash table of the elements
        :math:`a5(l_1,m_1,l_2,m_2,l_d)*P_{l_d}^{m_1-m_2}(0)` and :math:`a5(l_1,m_1,l_2,m_2,l_d)*P_{l_d}^{m_1-m_2}(0)`
        found in appendix B of [Egel 2018 diss],where a5 and b5 are the coefficients used in the evaluation of the SVWF translation
        operator and :math:`P_l^m(\cos\theta)` are the normalized associated Legendre functions. This hash table is usefull in 
        reducing computation time when calculating the coupling between two particles that exist in the same layer 
        and have the same z coordinate.


        Args:
            lmax1 (int):           Largest polar quantum number of the recieving particle
            lmax2  (int):          Largest polar quantum number  of the emitting particle
            mmax1  (int):          Largest azimuthal quantum number of the recieving particle
            mmax2  (int):          Largest azimuthal quantum number of the emitting particle

    
        Returns:
            a5leg_array (ndarray):          hash table of the elements in :math:`A` not dependent on :math:`/phi` or :math:`kd`
            b5leg_array (ndarray):          hash table of the elements in :math:`B` not dependent on :math:`/phi` or :math:`kd`          

        """
        
        a5leg_array = []
        b5leg_array = []
        #ab5leg_key = []
        
        legendre_2D = sma.legendre_normalized(0, 1, lmax1+lmax2+1)[0][:,:,0]
        for im1, m1 in enumerate(range(-mmax1, mmax1 + 1)):
            for im2, m2 in enumerate(range(-mmax2, mmax2 + 1)):
                for il1, l1 in enumerate(range(max(1, abs(m1)), lmax1 + 1)):
                    for il2, l2 in enumerate(range(max(1, abs(m2)), lmax2 + 1)):
                        for ild, ld in enumerate(range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1)): 
                            
                            a5, b5 = trf.ab5_coefficients(l2, m2, l1, m1, ld)
                            leg =  legendre_2D[ld][abs(m1 - m2)]
                           
                            a5leg = a5*leg
                            b5leg = b5*leg
                            #ab5leg_key.append([l1, l2, m1, m2, ld])
                            a5leg_array.append(a5leg)
                            b5leg_array.append(b5leg)
                          
    
        return a5leg_array, b5leg_array


try:
    from smuthi.utility.cython.cython_speedups  import _ab5_coefficient_hash_table
    @memo.Memoize
    def ab5_coefficient_hash_table(lmax1, lmax2, mmax1, mmax2):
        return _ab5_coefficient_hash_table(int(lmax1), int(lmax2), int(mmax1), int(mmax2))  
except:
    @memo.Memoize
    def ab5_coefficient_hash_table(lmax1, lmax2, mmax1, mmax2):
        r"""Creates a hash table of the elements
        :math:`a5(l_1,m_1,l_2,m_2,l_d)` and :math:`a5(l_1,m_1,l_2,m_2,l_d)`
        found in appendix B of [Egel 2018 diss],where a5 and b5 are the coefficients used in the evaluation of the SVWF translation
        operator. This hash table is usefull in reducing computation time when calculating the coupling between two particles
        that exist in the same layer but do not have the same z coordinate.


        Args:
            lmax1 (int):           Largest polar quantum number of the recieving particle
            lmax2  (int):          Largest polar quantum number  of the emitting particle
            mmax1  (int):          Largest azimuthal quantum number of the recieving particle
            mmax2  (int):          Largest azimuthal quantum number of the emitting particle

    
        Returns:
            a5_array (ndarray):          Hash table of the elements in :math:`A` not dependent on :math:`/theta`, :math:`/phi` or :math:`kd`
            b5_array (ndarray):          Hash table of the elements in :math:`B` not dependent on :math:`/theta`, :math:`/phi` or :math:`kd`          

        """        
        a5_array = []
        b5_array = []
        #ab5_key = []
        for im1, m1 in enumerate(range(-mmax1, mmax1 + 1)):
            for im2, m2 in enumerate(range(-mmax2, mmax2 + 1)):
                for il1, l1 in enumerate(range(max(1, abs(m1)), lmax1 + 1)):
                    for il2, l2 in enumerate(range(max(1, abs(m2)), lmax2 + 1)):
                        for ild, ld in enumerate(range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1)): 
                            a5, b5 = trf.ab5_coefficients(l2, m2, l1, m1, ld)
                            #ab5_key.append([l1, l2, m1, m2, ld])
                            a5_array.append(a5)
                            b5_array.append(b5)
                          
        return a5_array, b5_array


try:
    from smuthi.utility.cython.cython_speedups import direct_coupling_block_2D_from_hash_table 
except:
    def direct_coupling_block_2D_from_hash_table(blocksize1,blocksize2,
                                        w,sph,k,dx, dy, dz,
                                        lmax1, lmax2, mmax1, mmax2,
                                        a5leg_hash_table,b5leg_hash_table,
                                        threaded = False):

        r"""Subroutine to calculate the direct coupling between two particles 
            that are in the same layer and have the same z coordinate. This 
            subroutine is called internally by direct_coupling_block. 
            If Cython is enabled then this subroutine is Cython accelerated.
            Otherwise, the Python equivalent subroutine is used.


        Args:
            blocksize1 (int):           Number of columns in the direct coupling block
            blocksize2 (int):           Number of rows in the direct coupling block
            w  (ndarray):               Zero initialized direct coupling matrix block as numpy array
            sph  (ndarray):             Zero initialized Hankel function matrix as numpy array
            k (complex):                Wavenumber in the shared media
            dx (float):                 x-coordinate of the distance between the two particles
            dy  (float):                y-coordinate of the distance between the two particles
            dz  (float):                z-coordinate of the distance between the two particles
            lmax1 (int):           Largest polar quantum number of the recieving particle
            lmax2  (int):          Largest polar quantum number  of the emitting particle
            mmax1  (int):          Largest azimuthal quantum number of the recieving particle
            mmax2  (int):          Largest azimuthal quantum number of the emitting particle
            a5leg_array (ndarray):          Hash table of the elements in :math:`A` not dependent on :math:`/phi` or :math:`kd`
            b5leg_array (ndarray):          Hash table of the elements in :math:`B` not dependent on :math:`/phi` or :math:`kd`          
            threaded (bool):          Flag to enable multithreading (valid only for Cython where the Gil can be released). Currently hard coded to False.  
            
        Returns:
            w  (ndarray):          Direct coupling matrix block as numpy array.
        """


        # calculate distance vector
        d = np.sqrt(dx**2 + dy**2 + dz**2)
        phi = np.arctan2(dy, dx)
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = sma.spherical_hankel(int(n), complex(k * d))
        
        # the particle coupling operator is the transpose of the SVWF translation operator
        # therefore, (l1,m1) and (l2,m2) are interchanged:
        idx_ab5leg = 0
        for m1 in range(-mmax1, mmax1 + 1):
            for m2 in range(-mmax2, mmax2 + 1):
                eimph_val = np.exp(1j * (m2 - m1) * phi)
                for l1 in range(max(1, abs(m1)), lmax1 + 1):
                    for l2 in range(max(1, abs(m2)), lmax2 + 1):
                        A, B = complex(0), complex(0)
                        for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1):  # if ld<abs(m1-m2) then P=0

                            A += a5leg_hash_table[idx_ab5leg]*sph[ld] 
                            B += b5leg_hash_table[idx_ab5leg]*sph[ld]
                            idx_ab5leg += 1
                            
                        A, B = eimph_val * A, eimph_val * B
                        for tau1 in range(2):
                            n1 = flds.multi_to_single_index(tau1, l1, m1, lmax1, mmax1)
                            for tau2 in range(2):
                                n2 = flds.multi_to_single_index(tau2, l2, m2, lmax2, mmax2)
                                if tau1 == tau2:
                                    w[n1, n2] = A
                                else:
                                    w[n1, n2] = B
    
        return w
    

try:
    from smuthi.utility.cython.cython_speedups import direct_coupling_block_3D_from_hash_table
except:  
    def direct_coupling_block_3D_from_hash_table(blocksize1,blocksize2,
                                        w,sph,k,dx, dy, dz,
                                        lmax1, lmax2, mmax1, mmax2,
                                        a5_hash_table, b5_hash_table,
                                        threaded = False):
        
        r"""Subroutine to calculate the direct coupling between two particles 
            that are in the same layer and do not have the same z coordinate. This 
            subroutine is called internally by direct_coupling_block. 
            If Cython is enabled then this subroutine is Cython accelerated.
            Otherwise, the Python equivalent subroutine is used.


        Args:
            blocksize1 (int):           Number of columns in the direct coupling block
            blocksize2 (int):           Number of rows in the direct coupling block
            w  (ndarray):               Zero initialized direct coupling matrix block as numpy array
            sph  (ndarray):             Zero initialized Hankel function matrix as numpy array
            k (complex):                Wavenumber in the shared media
            dx (float):                 x-coordinate of the distance between the two particles
            dy  (float):                y-coordinate of the distance between the two particles
            dz  (float):                z-coordinate of the distance between the two particles
            lmax1 (int):           Largest polar quantum number of the recieving particle
            lmax2  (int):          Largest polar quantum number  of the emitting particle
            mmax1  (int):          Largest azimuthal quantum number of the recieving particle
            mmax2  (int):          Largest azimuthal quantum number of the emitting particle
            a5_array (ndarray):          Hash table of the elements in :math:`A` not dependent on :math:`/theta`, :math:`/phi` or :math:`kd`
            b5_array (ndarray):          Hash table of the elements in :math:`B` not dependent on :math:`/theta`, :math:`/phi` or :math:`kd`          
            threaded (bool):          Flag to enable multithreading (valid only for Cython where the Gil can be released). Currently hard coded to False.  
            
        Returns:
            w  (ndarray):          Direct coupling matrix block as numpy array.
        """        

        d = np.sqrt(dx**2 + dy**2 + dz**2)
        cos_theta = dz / d
        sin_theta = np.sqrt(dx**2 + dy**2) / d
        phi = np.arctan2(dy, dx)
    
        # spherical functions
        for n in range(lmax1 + lmax2 + 1):
            sph[n] = sma.spherical_hankel(int(n), complex(k * d))
        legendre = sma.legendre_normalized(cos_theta, sin_theta, lmax1 + lmax2)[0][:,:,0]
        
        # the particle coupling operator is the transpose of the SVWF translation operator
        # therefore, (l1,m1) and (l2,m2) are interchanged:
        idx_ab5 = 0
        for m1 in range(-mmax1, mmax1 + 1):
            for m2 in range(-mmax2, mmax2 + 1):
                eimph_val = np.exp(1j * (m2 - m1) * phi)
                for l1 in range(max(1, abs(m1)), lmax1 + 1):
                    for l2 in range(max(1, abs(m2)), lmax2 + 1):
                        A, B = complex(0), complex(0)
                        for ld in range(max(abs(l1 - l2), abs(m1 - m2)), l1 + l2 + 1):  # if ld<abs(m1-m2) then P=0
                            A += a5_hash_table[idx_ab5]*sph[ld]*legendre[ld][abs(m1 - m2)] 
                            B += b5_hash_table[idx_ab5]*sph[ld]*legendre[ld][abs(m1 - m2)] 
                            idx_ab5 += 1
                            
                        A, B = eimph_val * A, eimph_val * B
                        for tau1 in range(2):
                            n1 = flds.multi_to_single_index(tau1, l1, m1, lmax1, mmax1)
                            for tau2 in range(2):
                                n2 = flds.multi_to_single_index(tau2, l2, m2, lmax2, mmax2)
                                if tau1 == tau2:
                                    w[n1, n2] = A
                                else:
                                    w[n1, n2] = B
    
        return w



###############################################################################
#                  PVWF coupling - experimental!                              #                               
###############################################################################
"""The following code section contains functions to compute the particle 
coupling via a PVWF expansion. This allows in principle to treat particles
with intersecting circumscribing spheres, see 
Theobald et al.: "Plane-wave coupling formalism for T-matrix simulations of 
light scattering by nonspherical particles", Phys Rev A, 2018"""


def spheroids_closest_points(ab_halfaxis1, c_halfaxis1, center1, orientation1, ab_halfaxis2, c_halfaxis2, center2, 
                             orientation2):
    """ Computation of the two closest points of two adjacent spheroids.
    For details, see: Stephen B. Pope, Algorithms for Ellipsoids, Sibley School of Mechanical & Aerospace Engineering, 
    Cornell University, Ithaca, New York, February 2008
    
    Args:
        ab_halfaxis1 (float):        Half axis orthogonal to symmetry axis of spheroid 1
        c_halfaxis1 (float):         Half axis parallel to symmetry axis of spheroid 1
        center1 (numpy.array):       Center coordinates of spheroid 1
        orientation1 (numpy.array):  Orientation angles of spheroid 1
        ab_halfaxis2 (float):        Half axis orthogonal to symmetry axis of spheroid 2
        c_halfaxis2 (float):         Half axis parallel to symmetry axis of spheroid 2
        center2 (numpy.array):       Center coordinates of spheroid 2
        orientation2 (numpy.array):  Orientation angles of spheroid 2
        
    Retruns:
        Tuple containing:
          - closest point on first particle (numpy.array)
          - closest point on second particle (numpy.array)
          - first rotation Euler angle alpha (float)
          - second rotation Euler angle beta (float)
    """
    
    def rotation_matrix(ang):
        rot_mat = (np.array([[np.cos(ang[0]) * np.cos(ang[1]), -np.sin(ang[0]), np.cos(ang[0]) * np.sin(ang[1])],
                             [np.sin(ang[0]) * np.cos(ang[1]), np.cos(ang[0]), np.sin(ang[0]) * np.sin(ang[1])],
                             [-np.sin(ang[1]), 0, np.cos(ang[1])]]))
        return rot_mat
    
    rot_matrix_1 = rotation_matrix(orientation1)
    rot_matrix_2 = rotation_matrix(orientation2)
        
    a1, a2 = ab_halfaxis1, ab_halfaxis2
    c1, c2 = c_halfaxis1, c_halfaxis2
    ctr1, ctr2 = np.array(center1), np.array(center2)
    
    eigenvalue_matrix_1 = np.array([[1 / a1 ** 2, 0, 0], [0, 1 / a1 ** 2, 0], [0, 0, 1 / c1 ** 2]])
    eigenvalue_matrix_2 = np.array([[1 / a2 ** 2, 0, 0], [0, 1 / a2 ** 2, 0], [0, 0, 1 / c2 ** 2]])
    
    E1 = np.dot(rot_matrix_1, np.dot(eigenvalue_matrix_1, np.transpose(rot_matrix_1)))
    E2 = np.dot(rot_matrix_2, np.dot(eigenvalue_matrix_2, np.transpose(rot_matrix_2)))
    S = np.matrix.getH(np.linalg.cholesky(E1))
    
    # transformation of spheroid E1 into the unit-sphere with its center at origin / same transformation on E2
    # E1_prime = np.dot(np.transpose(np.linalg.inv(S)), np.dot(E1, np.linalg.inv(S)))
    # ctr1_prime = np.zeros_like(ctr1)
    E2_prime = np.dot(np.transpose(np.linalg.inv(S)), np.dot(E2, np.linalg.inv(S)))
    ctr2_prime = -(np.dot(S, (ctr1 - ctr2)))  
    E2_prime_L = np.linalg.cholesky(E2_prime)
        
    H = np.dot(np.linalg.inv(E2_prime_L), np.transpose(np.linalg.inv(E2_prime_L)))
    p = np.array([0, 0, 0])
    f = np.dot(np.transpose(ctr2_prime - p), np.transpose(np.linalg.inv(E2_prime_L)))
    
    def minimization_fun(y_vec):
        fun = 0.5 * np.dot(np.dot(np.transpose(y_vec), H), y_vec) + np.dot(f, y_vec)
        return fun

    def constraint_fun(x):
        eq_constraint = (x[0] ** 2 + x[1] ** 2 + x[2] ** 2) ** 0.5 - 1
        return eq_constraint
    
    bnds = ((-1, 1), (-1, 1), (-1, 1))
    length_constraints = {'type' : 'eq', 'fun' : constraint_fun}
    
    flag = False
    while not flag:
        x0 = -1 + np.dot((1 + 1), np.random.rand(3))
        optimization_result = scipy.optimize.minimize(minimization_fun, x0, method='SLSQP', bounds=bnds,
                                                      constraints=length_constraints, tol=None, callback=None, options=None)
        x_vec = np.transpose(np.dot(np.transpose(np.linalg.inv(E2_prime_L)), optimization_result['x'])
                             + np.transpose(ctr2_prime))
        if optimization_result['success'] == True:
            if np.linalg.norm(x_vec) <= 1:
                raise ValueError("particle error: particles intersect")
            elif np.linalg.norm(x_vec) < np.linalg.norm(ctr2_prime):
                flag = True
            else:
                print('wrong minimum ...')
        else:
            print('No minimum found ...')

    p2_prime = x_vec
    p2 = np.dot(np.linalg.inv(S), p2_prime) + ctr1
    
    E1_L = np.linalg.cholesky(E1)
    H = np.dot(np.linalg.inv(E1_L), np.transpose(np.linalg.inv(E1_L)))
    p = p2
    f = np.dot(np.transpose(ctr1 - p), np.transpose(np.linalg.inv(E1_L)))
       
    flag = False
    while not flag:
        x0 = -1 + np.dot((1 + 1), np.random.rand(3))
        optimization_result2 = scipy.optimize.minimize(minimization_fun, x0, method='SLSQP', bounds=bnds,
                                                      constraints=length_constraints, tol=None, callback=None, options=None)
        p1 = np.transpose(np.dot(np.transpose(np.linalg.inv(E1_L)), optimization_result2['x']) + np.transpose(ctr1))
        if optimization_result2['success'] == True:
            if np.linalg.norm(p1 - p) < np.linalg.norm(ctr1 - p):
                flag = True
            else:
                print('wrong minimum ...')
        else:
            print('No minimum found ...')
    
    p1p2 = p2 - p1
    azimuth = np.arctan2(p1p2[1], p1p2[0])
    elevation = np.arctan2(p1p2[2], (p1p2[0] ** 2 + p1p2[1] ** 2) ** 0.5)

    if p1p2[2] < 0:
        beta = (np.pi / 2) + elevation
    else:
        beta = (-np.pi / 2) + elevation
    alpha = -azimuth
              
    return p1, p2, alpha, beta


def get_separating_plane(receiving_particle, emitting_particle):
    """Estimate the orientation of a plane that separates two particles.
    Such a plane always exists for convex particles.

    For some special cases (pairs of spheroids, pairs of non-rotated cylinders)
    a separating plane is constructed.

    For all other cases, the orientation is chosen along the center-to-center
    vector. Note that this choice is not guaranteed
    to yield a correct PVWF coupling result.

        Args:
        receiving_particle (smuthi.particles.Particle):        Receiving particle
        emitting_particle (smuthi.particles.Particle):         Emitting particle

    Retruns:
        Tuple containing Euler angles for the active rotation into a frame such
        that the emitting particle is above some z-plane and the receiving particle
        is below that plane:
          - first rotation Euler angle alpha (float)
          - second rotation Euler angle beta (float)
    """

    relpos = np.array(emitting_particle.position) - np.array(receiving_particle.position)

    if type(receiving_particle).__name__ == 'Spheroid' \
            and type(emitting_particle).__name__ == 'Spheroid':
        _, _, alpha, beta = spheroids_closest_points(
            emitting_particle.semi_axis_a, emitting_particle.semi_axis_c,
            emitting_particle.position, emitting_particle.euler_angles,
            receiving_particle.semi_axis_a, receiving_particle.semi_axis_c,
            receiving_particle.position, receiving_particle.euler_angles)
        return alpha, beta

    elif (type(receiving_particle).__name__ == 'FiniteCylinder'
            and type(emitting_particle).__name__ == 'FiniteCylinder'
            and np.linalg.norm(receiving_particle.euler_angles) == 0
            and np.linalg.norm(emitting_particle.euler_angles) == 0):
        in_plane_distance2 = relpos[0]**2 + relpos[1]**2
        sum_radii2 = (emitting_particle.cylinder_radius + receiving_particle.cylinder_radius)**2
        if in_plane_distance2 > sum_radii2:
            alpha = - np.arctan2(relpos[1], relpos[0])
            beta = - np.pi / 2
            return alpha, beta
        if relpos[2] > (emitting_particle.cylinder_height + receiving_particle.cylinder_height) / 2.0:
            return 0, 0
        if -relpos[2] > (emitting_particle.cylinder_height + receiving_particle.cylinder_height) / 2.0:
            return 0, np.pi

    warnings.warn("Unable to estimate suitable separating plane. Using center-to-center vector for PVWF coupling.")
    alpha = - np.arctan2(relpos[1], relpos[0])
    beta = - np.arccos(relpos[2] / np.linalg.norm(relpos))
    return alpha, beta


def direct_coupling_block_pvwf_mediated(vacuum_wavelength, receiving_particle, emitting_particle, layer_system,
                                        k_parallel, alpha=None, beta=None):
    """Direct particle coupling matrix :math:`W` for two particles (via plane vector wave functions).
    For details, see:
    Dominik Theobald et al., Phys. Rev. A 96, 033822, DOI: 10.1103/PhysRevA.96.033822 or arXiv:1708.04808

    The plane wave coupling is performed in a rotated coordinate system,
    which must be chosen such that both particles can be separated by a plane
    that is parallel to the xy-plane (such that the emitting particle is
    entirely above that plane and the receiving particle is entirely below that
    plane).

    Two angles (alpha and beta) are required to specify the active rotation into
    that coordinate system, i.e., the rotation which rotates the particle
    locations such that the abovementioned condition is fulfilled.

    If the angle arguments are omitted, Smuthi tries to estimate a suitable
    separating plane.

    Args:
        vacuum_wavelength (float):                          Vacuum wavelength :math:`\lambda` (length unit)
        receiving_particle (smuthi.particles.Particle):     Particle that receives the scattered field
        emitting_particle (smuthi.particles.Particle):      Particle that emits the scattered field
        layer_system (smuthi.layers.LayerSystem):           Stratified medium in which the coupling takes place
        k_parallel (numpy.array):                           In-plane wavenumber for plane wave expansion
        alpha (float):                                      First Euler angle, rotation around z-axis, in rad
        beta (float):                                       Second Euler angle, rotation around y'-axis in rad

    Returns:
        Direct coupling matrix block (numpy array).
    """    
    lmax1 = receiving_particle.l_max
    mmax1 = receiving_particle.m_max
    assert lmax1 == mmax1, 'PVWF coupling requires lmax == mmax for each particle.'
    lmax2 = emitting_particle.l_max
    mmax2 = emitting_particle.m_max
    assert lmax2 == mmax2, 'PVWF coupling requires lmax == mmax for each particle.'
    lmax = max([lmax1, lmax2])
    blocksize1 = flds.blocksize(lmax1, mmax1)
    blocksize2 = flds.blocksize(lmax2, mmax2)

    # initialize result
    w = np.zeros((blocksize1, blocksize2), dtype=complex)

    # Check if particles are in the same layer.
    rS1 = receiving_particle.position
    rS2 = emitting_particle.position
    iS1 = layer_system.layer_number(rS1[2])
    iS2 = layer_system.layer_number(rS2[2])
    if (iS1 != iS2) or emitting_particle == receiving_particle:
        return w

    n_medium = layer_system.refractive_indices[layer_system.layer_number(receiving_particle.position[2])]
    
    if alpha is None or beta is None:
        alpha, beta = get_separating_plane(receiving_particle, emitting_particle)

    # positions
    r1 = np.array(receiving_particle.position)
    r2 = np.array(emitting_particle.position)
    r21_lab = r1 - r2  # laboratory coordinate system
    
    # distance vector in rotated coordinate system
    r21_rot = np.dot(np.dot([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [- np.sin(beta), 0, np.cos(beta)]],
                           [[np.cos(alpha), - np.sin(alpha), 0], [np.sin(alpha), np.cos(alpha), 0], [0, 0, 1]]), 
                    r21_lab)
    rho21 = (r21_rot[0] ** 2 + r21_rot[1] ** 2) ** 0.5 
    phi21 = np.arctan2(r21_rot[1], r21_rot[0])
    z21 = r21_rot[2]
    
    # wavenumbers
    omega = flds.angular_frequency(vacuum_wavelength)
    k = omega * n_medium
    kz = flds.k_z(k_parallel=k_parallel, vacuum_wavelength=vacuum_wavelength, refractive_index=n_medium)
    if z21 < 0:
        kz_var = -kz
    else:
        kz_var = kz
        
    # Bessel lookup 
    bessel_list = []
    for dm in range(mmax1 + mmax2 + 1):
        bessel_list.append(scipy.special.jn(dm, k_parallel * rho21))
    
    # legendre function lookups
    ct = kz_var / k
    st = k_parallel / k
    _, pilm_list, taulm_list = sma.legendre_normalized(ct, st, lmax)
    
    # prefactor
    const_arr = k_parallel / (kz * k) * np.exp(1j * (kz_var * z21))
                        
    for m1 in range(-mmax1, mmax1 + 1):
        for m2 in range(-mmax2, mmax2 + 1):
            jmm_eimphi_bessel = 4 * 1j ** abs(m2 - m1) * np.exp(1j * phi21 * (m2 - m1)) * bessel_list[abs(m2 - m1)]
            prefactor = const_arr * jmm_eimphi_bessel
            for l1 in range(max(1, abs(m1)), lmax1 + 1):
                for l2 in range(max(1, abs(m2)), lmax2 + 1):
                    for tau1 in range(2):
                        n1 = flds.multi_to_single_index(tau1, l1, m1, lmax1, mmax1)
                        for tau2 in range(2):
                            n2 = flds.multi_to_single_index(tau2, l2, m2, lmax2, mmax2)
                            for pol in range(2):
                                B_dag = trf.transformation_coefficients_vwf(tau1, l1, m1, pol, pilm_list=pilm_list,
                                                                        taulm_list=taulm_list, dagger=True)
                                B = trf.transformation_coefficients_vwf(tau2, l2, m2, pol, pilm_list=pilm_list,
                                                                            taulm_list=taulm_list, dagger=False)
                                integrand = prefactor * B * B_dag
                                w[n1, n2] += np.trapz(integrand, k_parallel) 
                                
    rot_mat_1 = trf.block_rotation_matrix_D_svwf(lmax1, mmax1, 0, beta, alpha)
    rot_mat_2 = trf.block_rotation_matrix_D_svwf(lmax2, mmax2, -alpha, -beta, 0)
    
    return np.dot(np.dot(np.transpose(rot_mat_1), w), np.transpose(rot_mat_2))
