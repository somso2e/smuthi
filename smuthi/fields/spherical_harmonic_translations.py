"""Functions for the translation of spherical vector wave functions."""

import smuthi.fields as flds
import smuthi.utility.math as sma
import smuthi.utility.memoizing as memo
import smuthi.fields.transformations as trf
import numpy as np
import sys

##############################################################################
#                       Highest Level Translation Method 
#                   Used to directely interface with Smuthi. 
#          All other funcitons in this file are internal to this method.
############################################################################## 
# Note: The direct coupling block is just an outgoing to incoming translation


def svwf_translate(k, d, lmax1, mmax1, lmax2, mmax2, kind, threaded = False):
    
    # Cython has a non-intuitive convention to dealing with strings/chars between python2 and python3
    # To avoid this, convert the readable string flag to a binary int flag. This prevents compatibility problems.
    # Find more at https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
    if kind == 'outgoing to regular':
        ckind = 0 
    elif kind == 'outgoing to outgoing':
        ckind = 1
    else:
        raise Exception('''
                        The Spherical harmonic translations only support the following conversions:
                                            "outgoing to regular"
                                                    and 
                                            "outgoing to outgoing"
                        You specified a kind that was not one of those two strings.
                        ''')
        
    # If the translation kind is valid, precalculate necessary inputs for translation.
    # Note: Translations go from 1 -> 2
    blocksize1 = int(flds.blocksize(lmax1, mmax1)) # Rows of translation matrix
    blocksize2 = int(flds.blocksize(lmax2, mmax2)) # Cols of translation matrix
    w = np.zeros((blocksize1, blocksize2), dtype = np.complex128, order = 'C')
    sph = np.zeros((lmax1+lmax2+1), dtype=np.complex128, order = 'C')
    dx = float(d[0])
    dy = float(d[1])
    dz = float(d[2])
    
    #Note: It is allways better to use the hash tables so just do it. 
    if dz == 0:
        # Note: The hash table is memoized. So if this lmax and mmax combo was used 
        # prior (such as in the solution of the linear system) it will not be re-calculated
        a5leg_array, b5leg_array =  ab5_coefficient_and_legendre_hash_table(lmax1, lmax2, mmax1, mmax2)
        return svwf_translate_2D_from_hash_table(blocksize1, blocksize2,
                                                    w, sph, k, dx, dy,dz,
                                                    lmax1, lmax2, mmax1, mmax2,
                                                    a5leg_array, b5leg_array,                              
                                                    ckind, threaded)
    else:
        # Note: The hash table is memoized. So if this lmax and mmax combo was used 
        # prior (such as in the solution of the linear system) it will not be re-calculated
        a5_array, b5_array =  ab5_coefficient_hash_table(lmax1, lmax2, mmax1, mmax2)
        return svwf_translate_3D_from_hash_table(blocksize1, blocksize2,
                                                    w, sph, k, dx, dy,dz,
                                                    lmax1, lmax2, mmax1, mmax2,
                                                    a5_array, b5_array,                              
                                                    ckind, threaded)           
    



##############################################################################
#  Make  hash table of ab5  for a particle order pair 
#  This is for 3D coupling (Cython uses the GIL)
##############################################################################
# """
# This hash table is particularly good when all particles have the same order. 
# """

try:
    from smuthi.utility.cython.cython_speedups  import _ab5_coefficient_hash_table
    @memo.Memoize
    def  ab5_coefficient_hash_table(lmax1, lmax2, mmax1, mmax2):
        return _ab5_coefficient_hash_table(int(lmax1), int(lmax2), int(mmax1), int(mmax2))  
except:
    @memo.Memoize
    def  ab5_coefficient_hash_table(lmax1, lmax2, mmax1, mmax2):
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

##############################################################################
#         3D translation and direct coupling block using hash table 
#                          (Cython releases the Gil)
############################################################################## 


# Call to lower level translation funciton. Use Try to enable Cython if able. 
try:
    from smuthi.utility.cython.cython_speedups import svwf_translate_3D_from_hash_table
except:  
    def svwf_translate_3D_from_hash_table(blocksize1,blocksize2,
                                        w,sph,k,dx, dy, dz,
                                        lmax1, lmax2, mmax1, mmax2,
                                        a5_hash_table, b5_hash_table,
                                        kind, threaded):
        
        r"""Subroutine to calculate the translation between two locations
            that are in the same layer and do not have the same z coordinate. This 
            subroutine is called internally when no Cython accelorations are found.
            This subroutine is the Python equivalent of the Cython accelorated code.


        Args:
            blocksize1 (int):           Number of columns in the direct coupling block
            blocksize2 (int):           Number of rows in the direct coupling block
            w  (ndarray):               Zero initialized direct coupling matrix block as numpy array
            sph  (ndarray):             Hankel or Bessel function matrix as numpy array
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
    
        legendre = sma.legendre_normalized(cos_theta, sin_theta, lmax1 + lmax2)[0][:,:,0]
        
        if kind == 0: # 0 = outgoing to incoming wave, 1 = outgoing to outgoing wave. 
        #NOTE: an long flag is used instead of a string because Cython 
        # has troublesome conventions about datatyping strings from python v2 to v3
        # See more on this at https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
            for n in range(lmax1 + lmax2 + 1):
                sph[n] = sma.spherical_hankel(n, k*d) 
        else:
            for n in range(lmax1 + lmax2 + 1):
                sph[n] = sma.spherical_bessel(n, k*d) 
                
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



##############################################################################
#  Make  hash table of ab5 and Legendre for a particle order pair 
#  This is for 2D coupling (Cython uses the GIL)
##############################################################################
# """
# This hash table is particularly good when all particles have the same order. 
# """

try: 
    from smuthi.utility.cython.cython_speedups import _ab5_coefficient_and_legendre_hash_table
    @memo.Memoize
    def ab5_coefficient_and_legendre_hash_table(lmax1, lmax2, mmax1, mmax2):
        return _ab5_coefficient_and_legendre_hash_table(int(lmax1), int(lmax2), int(mmax1), int(mmax2))
except:
    sys.stdout.write(
"""
Cython acceleration could not be loaded.
Falling back on Python equivalents... 
"""
        )
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





##############################################################################
#  2D Direct coupling block using  hash table (Cython releases the Gil)
##############################################################################

try:
    from smuthi.utility.cython.cython_speedups import svwf_translate_2D_from_hash_table
except:  
    def svwf_translate_2D_from_hash_table(blocksize1,blocksize2,
                                    w,sph,k,dx, dy, dz,
                                    lmax1, lmax2, mmax1, mmax2,
                                    a5leg_hash_table,b5leg_hash_table,
                                    kind, threaded):

        r"""Subroutine to calculate the direct coupling between two particles 
            that are in the same layer and have the same z coordinate. This 
            subroutine is called internally by direct_coupling_block. 
            If Cython is enabled then this subroutine is Cython accelerated.
            Otherwise, the Python equivalent subroutine is used.
    
    
        Args:
            blocksize1 (int):           Number of columns in the direct coupling block
            blocksize2 (int):           Number of rows in the direct coupling block
            w  (ndarray):               Zero initialized direct coupling matrix block as numpy array
            sph  (ndarray):             Hankel or Bessel function matrix as numpy array
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
        
        if kind == 0: # 0 = outgoing to incoming wave, 1 = outgoing to outgoing wave. 
        #NOTE: an long flag is used instead of a string because Cython 
        # has troublesome conventions about datatyping strings from python v2 to v3
        # See more on this at https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
            for n in range(lmax1 + lmax2 + 1):
                sph[n] = sma.spherical_hankel(n, k*d) 
        else:
            for n in range(lmax1 + lmax2 + 1):
                sph[n] = sma.spherical_bessel(n, k*d) 
                
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
    
    