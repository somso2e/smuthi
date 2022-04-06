#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 23:07:49 2022

@author: parkerwray
"""



import numpy as np
import smuthi.fields as flds
import smuthi.fields.transformations as trf
import smuthi.utility.math as sma
import sys
import cmath

from numba.core.typing import cffi_utils
from pywigxjpf_ffi import ffi, lib
import pywigxjpf_ffi
cffi_utils.register_module(pywigxjpf_ffi)
nb_wig3jj = pywigxjpf_ffi.lib.wig3jj

lib.wig_table_init(100,9)
lib.wig_temp_init(100)





def is_same(val1, val2):
    if (cmath.isclose(val1, val2) and cmath.isclose(np.real(val1), np.real(val2)) 
    and cmath.isclose(np.imag(val1), np.imag(val2))):
        return 0
    else:
        print()
        print(str(val1))
        print(str(val2))
        print()
        return 1
    
    
    
#%%   
##############################################################################
#  Make hard lookup table of ab5 and Legendre for a particle order pair 
#  This is for 2D coupling (Cython uses the GIL)
##############################################################################
# """
#     This lookup table is considered "hard" because it is dependent on the order 
#     of which particle is the emitter and which is the reciever. This lookup 
#     table is particularly good when all particles have the same order. 
# """


from smuthi.utility.cython.cython_speedups import hard_ab5_coefficient_and_legendre_lookup
def cab5_coefficient_and_legendre_lookup(lmax1, lmax2, mmax1, mmax2):
    return hard_ab5_coefficient_and_legendre_lookup(int(lmax1), int(lmax2), int(mmax1), int(mmax2))


def ab5_coefficient_and_legendre_lookup(lmax1, lmax2, mmax1, mmax2):
    r"""Creates a lookup table of the elements
    :math:`a5(l_1,m_1,l_2,m_2,l_d)*P_{l_d}^{m_1-m_2}(0)` and :math:`a5(l_1,m_1,l_2,m_2,l_d)*P_{l_d}^{m_1-m_2}(0)`
    found in appendix B of [Egel 2018 diss],where a5 and b5 are the coefficients used in the evaluation of the SVWF translation
    operator and :math:`P_l^m(\cos\theta)` are the normalized associated Legendre functions. This lookup table is usefull in 
    reducing computation time when calculating the coupling between two particles that exist in the same layer 
    and have the same z coordinate.


    Args:
        lmax1 (int):           Largest polar quantum number of the recieving particle
        lmax2  (int):          Largest polar quantum number  of the emitting particle
        mmax1  (int):          Largest azimuthal quantum number of the recieving particle
        mmax2  (int):          Largest azimuthal quantum number of the emitting particle


    Returns:
        a5leg_array (ndarray):          Lookup table of the elements in :math:`A` not dependent on :math:`/phi` or :math:`kd`
        b5leg_array (ndarray):          Lookup table of the elements in :math:`B` not dependent on :math:`/phi` or :math:`kd`          

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
                        a5leg_array.append(a5leg)
                        b5leg_array.append(b5leg)
                      

    return a5leg_array, b5leg_array

def test1(lmax1, lmax2, mmax1, mmax2):
    a5leg_array, b5leg_array = ab5_coefficient_and_legendre_lookup(lmax1, lmax2, mmax1, mmax2)
    ca5leg_array, cb5leg_array = cab5_coefficient_and_legendre_lookup(lmax1, lmax2, mmax1, mmax2)
    counter = 0
    for i in range(len(a5leg_array)):
        counter += is_same(a5leg_array[i], ca5leg_array[i])
        counter += is_same(b5leg_array[i], cb5leg_array[i])
    print(f'Number of errors is {counter} out of {2*len(a5leg_array)}')
    return a5leg_array, ca5leg_array

p,c = test1(8,8,8,8)


#%%
##############################################################################
#  Make hard lookup table of ab5  for a particle order pair 
#  This is for 3D coupling (Cython uses the GIL)
##############################################################################
# """
#     This lookup table is considered "hard" because it is dependent on the order 
#     of which particle is the emitter and which is the reciever. This lookup 
#     table is particularly good when all particles have the same order. 
# """


from smuthi.utility.cython.cython_speedups  import hard_ab5_coefficient_lookup
def cab5_coefficient_lookup(lmax1, lmax2, mmax1, mmax2):
    return hard_ab5_coefficient_lookup(int(lmax1), int(lmax2), int(mmax1), int(mmax2))  


def ab5_coefficient_lookup(lmax1, lmax2, mmax1, mmax2):
    r"""Creates a lookup table of the elements
    :math:`a5(l_1,m_1,l_2,m_2,l_d)` and :math:`a5(l_1,m_1,l_2,m_2,l_d)`
    found in appendix B of [Egel 2018 diss],where a5 and b5 are the coefficients used in the evaluation of the SVWF translation
    operator. This lookup table is usefull in reducing computation time when calculating the coupling between two particles
    that exist in the same layer but do not have the same z coordinate.


    Args:
        lmax1 (int):           Largest polar quantum number of the recieving particle
        lmax2  (int):          Largest polar quantum number  of the emitting particle
        mmax1  (int):          Largest azimuthal quantum number of the recieving particle
        mmax2  (int):          Largest azimuthal quantum number of the emitting particle


    Returns:
        a5_array (ndarray):          Lookup table of the elements in :math:`A` not dependent on :math:`/theta`, :math:`/phi` or :math:`kd`
        b5_array (ndarray):          Lookup table of the elements in :math:`B` not dependent on :math:`/theta`, :math:`/phi` or :math:`kd`          

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

def test1(lmax1, lmax2, mmax1, mmax2):
    a5leg_array, b5leg_array = ab5_coefficient_lookup(lmax1, lmax2, mmax1, mmax2)
    ca5leg_array, cb5leg_array = cab5_coefficient_lookup(lmax1, lmax2, mmax1, mmax2)
    counter = 0
    for i in range(len(a5leg_array)):
        counter += is_same(a5leg_array[i], ca5leg_array[i])
        counter += is_same(b5leg_array[i], cb5leg_array[i])
    print(f'Number of errors is {counter} out of {2*len(a5leg_array)}')
    return a5leg_array, ca5leg_array

p,c = test1(8,8,8,8)

#%%
##############################################################################
#  2D Direct coupling block using hard lookup table (Cython releases the Gil)
##############################################################################

try:
    from smuthi.utility.cython.cython_speedups import direct_coupling_block_2D_lookup 
except:
    def direct_coupling_block_2D_lookup(blocksize1,blocksize2,
                                        w,sph,k,dx, dy, dz,
                                        lmax1, lmax2, mmax1, mmax2,
                                        a5leg_lookup,b5leg_lookup,
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
            a5leg_array (ndarray):          Lookup table of the elements in :math:`A` not dependent on :math:`/phi` or :math:`kd`
            b5leg_array (ndarray):          Lookup table of the elements in :math:`B` not dependent on :math:`/phi` or :math:`kd`          
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

                            A += a5leg_lookup[idx_ab5leg]*sph[ld] 
                            B += b5leg_lookup[idx_ab5leg]*sph[ld]
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
    
    
##############################################################################
#  3D Direct coupling block using hard lookup table (Cython releases the Gil)
##############################################################################  
    
try:
    from smuthi.utility.cython.cython_speedups import direct_coupling_block_3D_lookup
except:  
    def direct_coupling_block_3D_lookup(blocksize1,blocksize2,
                                        w,sph,k,dx, dy, dz,
                                        lmax1, lmax2, mmax1, mmax2,
                                        a5_lookup, b5_lookup,
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
            a5_array (ndarray):          Lookup table of the elements in :math:`A` not dependent on :math:`/theta`, :math:`/phi` or :math:`kd`
            b5_array (ndarray):          Lookup table of the elements in :math:`B` not dependent on :math:`/theta`, :math:`/phi` or :math:`kd`          
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
                            A += a5_lookup[idx_ab5]*sph[ld]*legendre[ld][abs(m1 - m2)] 
                            B += b5_lookup[idx_ab5]*sph[ld]*legendre[ld][abs(m1 - m2)] 
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




































