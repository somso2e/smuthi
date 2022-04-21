"""This module contains functions to compute Ewald lattice sums."""

import smuthi.periodicboundaries.ewald_helper as pbeh
import smuthi.utility.math as sf
from numba import njit
import numpy as np



@njit()
def D_LM(L, M, k, k0t, a1, a2, eta, c=np.zeros(3, np.float64)):
    """ Computation of the Ewald sum to account for the coupling between one particle
        and a periodic particle arrangement.
    Args:
        L (int):                multipole degree
        M (int):                multipole order
        k (complex):            wavenumber
        k0t (numpy.ndarray):    complex in-plane wave vector in Carthesian coordinates 
        a1 (numpy.ndarray):     lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):     lattice vector 2 in Carthesian coordinates      
        eta (float):            Ewald sum separation parameter
        c (numpy.ndarray):      displacment vector between emitting and receiving particle
    Returns:
        D_LM (complex):         Ewald sum D1_LM + D2_LM + D3_00
    """        
    A = np.linalg.norm(np.cross(a1, a2))
    
    if np.linalg.norm(c) == 0:
        D1_LM_val = D1_LM(L, M, k, k0t, a1, a2, eta, A)
        D2_LM_val = D2_LM(L, M, k, k0t, a1, a2, eta)                
        if L == 0:
            D3_00_val = D3_00(k, eta)
            return D1_LM_val + D2_LM_val + D3_00_val
        else:
            return D1_LM_val + D2_LM_val
    else:
        # distance from receiving to emitting particle (opposing to the common distance vector)
        D1_LM_val = D1_LM_ij(L, M, -c, k, k0t, a1, a2, eta, A)
        D2_LM_val = D2_LM_ij(L, M, -c, k, k0t, a1, a2, eta)
        return D1_LM_val + D2_LM_val  
    
    
@njit()
def D1_LM(L, M, k, k0t, a1, a2, eta, A, c=np.zeros(3)):
    """ Ewald sum's reciprocal space summand to account for the coupling between one particle
        with it's own periodic arrangement. Beutel, arXiv:2004.08098,  Appendix D, (D2a).
    Args:
        L (int):                multipole degree
        M (int):                multipole order
        k (complex):            wavenumber
        k0t (numpy.ndarray):    complex in-plane wave vector in Carthesian coordinates 
        a1 (numpy.ndarray):     lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):     lattice vector 2 in Carthesian coordinates
        eta (float):            Ewald sum separation parameter 
        A (float):              area of the unit cell
        c (numpy.ndarray):      displacement vector betwenn particle S0j and S0j' in Carthesian coordinates
    Returns:
        D1_LM (complex):        reciprocal space summand of Ewald sum D_LM
    """
    if (L - abs(M)) % 2:
        return 0
        
    b1, b2 = pbeh.reciprocal_lattice_vec(a1, a2)      
    root = (2 * L + 1) ** 0.5  \
            * pbeh.numba_factorial(np.array([L - M], np.int32))[0, 0] ** 0.5 * pbeh.numba_factorial(np.array([L - M], np.int32))[0, 1] ** 0.5 \
            * pbeh.numba_factorial(np.array([L + M], np.int32))[0, 0] ** 0.5 * pbeh.numba_factorial(np.array([L + M], np.int32))[0, 1] ** 0.5
    prefac =  (1j) ** M * root / (A * k * (2 * k) ** L)
    
    x=0        
    n1, n2 = pbeh.n1n2_indices(x)
    D1_LM = compute_sum_D1LM(L, M, k, k0t, n1, n2, b1, b2, eta, c)
    
    flag = 0
    while flag == 0:
        D1_LM_prime = D1_LM        
        x += 1
        n1, n2 = pbeh.n1n2_indices(x)
        D1_LM += compute_sum_D1LM(L, M, k, k0t, n1, n2, b1, b2, eta, c)  
        if not D1_LM_prime == complex(0):
            if np.linalg.norm(np.array([D1_LM - D1_LM_prime])) / np.linalg.norm(np.array([D1_LM_prime])) < 1e-5:
                flag = 1
    
    return prefac * D1_LM


@njit()
def compute_sum_D1LM(L, M, k, k0t, n1, n2, b1, b2, eta, c=np.zeros(3)):
    """ Helper function for D1_LM(). """    
    shape = n1.shape[0]    
    gamma_g = np.zeros(shape, np.complex128)
    kgt_r = np.zeros(shape, np.float64)
    kgt_phi = np.zeros(shape, np.float64)
    
    Bgt = n1.reshape(shape, 1) * b1.reshape(1, 2) + n2.reshape(shape, 1) * b2.reshape(1, 2)
    kgt = k0t + Bgt
    for idx in range(shape):
        kgt_r[idx] = np.linalg.norm(kgt[idx, :])
        kgt_phi[idx] = np.arctan2(kgt[idx, 1], kgt[idx, 0])
    
    gamma_g = np.sqrt(k ** 2 - kgt_r ** 2)  
    gamma_g[np.where(gamma_g == complex(0))[0]] += 1e-10j                 
    eiMphi = np.exp(1j * M * kgt_phi)
    emnikgtc = np.exp(-1j * np.dot(kgt, c[:2])) 
    
    n = np.arange(int((L - abs(M)) / 2) + 1)
    denom = pbeh.numba_factorial(n)[:, 0] * pbeh.numba_factorial(n)[:, 1] \
            * pbeh.numba_factorial(int((L + M) / 2) - n)[:, 0] * pbeh.numba_factorial(int((L + M) / 2) - n)[:, 1] \
            * pbeh.numba_factorial(int((L - M) / 2) - n)[:, 0] * pbeh.numba_factorial(int((L - M) / 2) - n)[:, 1]
    nom = gamma_g.reshape(shape, 1) ** (2 * n - 1).reshape(1, n.shape[0]) * kgt_r.reshape(shape, 1) ** (L - 2 * n).reshape(1, n.shape[0])
    gamma_fun = pbeh.upper_incomplete_gamma_fun(n[-1], -gamma_g ** 2 / (4 * eta ** 2))  
    return np.sum(emnikgtc * eiMphi * np.sum(gamma_fun * nom / denom, axis=1))


@njit()
def D2_LM(L, M, k, k0t, a1, a2, eta):
    """ Ewald sum's real space summand to account for the coupling between one particle
        with it's own periodic arrangement. Beutel, arXiv:2004.08098,  Appendix D, (D2b).
    Args:
        L (int):                multipole degree
        M (int):                multipole order
        k (complex):            wavenumber
        k0t (numpy.ndarray):    complex in-plane wave vector in Carthesian coordinates 
        a1 (numpy.ndarray):     lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):     lattice vector 2 in Carthesian coordinates
        eta (float):            Ewald sum separation parameter      
    Returns:
        D2_LM (complex):        real space summand of Ewald sum D_LM
    """
    if (L - abs(M)) % 2:
        return 0
    
    frac = -1j * (-1) ** ((L + M) / 2) / (2 ** (L + 1) * np.pi \
            * pbeh.numba_factorial(np.array([(L - M) / 2], np.int32))[0, 0] * pbeh.numba_factorial(np.array([(L - M) / 2], np.int32))[0, 1] \
            * pbeh.numba_factorial(np.array([(L + M) / 2], np.int32))[0, 0] * pbeh.numba_factorial(np.array([(L + M) / 2], np.int32))[0, 1]) 
    root = (2 * L + 1) ** 0.5  \
            * pbeh.numba_factorial(np.array([L - M], np.int32))[0, 0] ** 0.5 * pbeh.numba_factorial(np.array([L - M], np.int32))[0, 1] ** 0.5 \
            * pbeh.numba_factorial(np.array([L + M], np.int32))[0, 0] ** 0.5 * pbeh.numba_factorial(np.array([L + M], np.int32))[0, 1] ** 0.5
    prefac = frac * root
      
    x=1
    n1, n2 = pbeh.n1n2_indices(x)
    D2_LM = compute_sum_D2LM(L, M, k, k0t, n1, n2, a1, a2, eta)
    
    flag = 0
    while flag == 0:
        D2_LM_prime = D2_LM       
        x += 1
        n1, n2 = pbeh.n1n2_indices(x)
        D2_LM += compute_sum_D2LM(L, M, k, k0t, n1, n2, a1, a2, eta) 
        if not D2_LM_prime == complex(0):
            if np.linalg.norm(np.array([D2_LM - D2_LM_prime])) / np.linalg.norm(np.array([D2_LM_prime])) < 1e-5:
                flag = 1
       
    return prefac * D2_LM


@njit()
def compute_sum_D2LM(L, M, k, k0t, n1, n2, a1, a2, eta):
    """ Helper function for D2_LM(). """ 
    shape = n1.shape[0] 
    Rn_r = np.zeros(shape, np.float64)
    Rn_phi = np.zeros(shape, np.float64)
    
    Rn = n1.reshape(shape, 1) * a1[:2].reshape(1, 2) + n2.reshape(shape, 1) * a2[:2].reshape(1, 2)
    for idx in range(shape):
        Rn_r[idx] = np.linalg.norm(Rn[idx, :])
        Rn_phi[idx] = np.arctan2(Rn[idx, 1], Rn[idx, 0])       
    k0tRn = np.dot(k0t, np.transpose(Rn))

    int_val = (k ** 2 / 4) ** (L + 0.5) * pbeh.int_recursion(L, eta, k, Rn_r) 
    return  np.sum(np.exp(1j * (k0tRn + M * (Rn_phi + np.pi))) / k * (2 * Rn_r / k) ** L * int_val)


@njit()
def D3_00(k, eta):
    """ Central lattice point correction for the Ewald sum to account for the coupling between
        one particle and it's own periodic arrangement. Beutel, arXiv:2004.08098,  Appendix D, (D2c).
    Args:
        k (complex):        wavenumber
        eta (float):        Ewald sum separation parameter      
    Returns:
        D3_00 (complex):    central lattice point correction of Ewald sum D_LM
    """
    x = np.array([-k ** 2 / (4 * eta ** 2)])
    return 1 / (4 * np.pi) * pbeh.upper_incomplete_gamma_fun(1, x)[0, -1]


@njit()
def D1_LM_ij(L, M, c, k, k0t, a1, a2, eta, A):
    """ Ewald sum's reciprocal space summand to account for the coupling between one particle (S_0i)
        and a different particle's (S_0j) periodic arrangement. Kambe, Z. Naturforschung 23a, (3.17).
    Args:        
        L (int):                multipole degree
        M (int):                multipole order
        c (numpy.ndarray):      displacement vector betwenn particle S0j and S0j' in Carthesian coordinates
        k (complex):            wavenumber
        k0t (numpy.ndarray):    complex in-plane wave vector in Carthesian coordinates 
        a1 (numpy.ndarray):     lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):     lattice vector 2 in Carthesian coordinates
        eta (float):            Ewald sum separation parameter  
        A (float):              area of one unit cell
    Returns:
        D1LM (complex):        reciprocal space summand of Ewald sum D_LM
    """   
    if c[2] == 0:
        return D1_LM(L, M, k, k0t, a1, a2, eta, A, c)
    else:  
        b1, b2 = pbeh.reciprocal_lattice_vec(a1, a2)
        root = (2 * L + 1) ** 0.5  \
            * pbeh.numba_factorial(np.array([L - M], np.int32))[0, 0] ** 0.5 * pbeh.numba_factorial(np.array([L - M], np.int32))[0, 1] ** 0.5 \
            * pbeh.numba_factorial(np.array([L + M], np.int32))[0, 0] ** 0.5 * pbeh.numba_factorial(np.array([L + M], np.int32))[0, 1] ** 0.5
        prefac = (-1j) ** M * root / ((-2) ** L * A * k ** 2)
        
        x=0        
        n1, n2 = pbeh.n1n2_indices(x)
        D1LM = compute_sum_D1LM_ij(L, M, k, k0t, n1, n2, b1, b2, eta, A, c)
                
        flag = 0
        while flag == 0:
            D1_LM_prime = D1LM        
            x += 1
            n1, n2 = pbeh.n1n2_indices(x)
            D1LM += compute_sum_D1LM_ij(L, M, k, k0t, n1, n2, b1, b2, eta, A, c)  
            if not D1_LM_prime == complex(0):
                if np.linalg.norm(np.array([D1LM - D1_LM_prime])) / np.linalg.norm(np.array([D1_LM_prime])) < 1e-5:
                    flag = 1
        return prefac * D1LM
    
    
@njit()
def compute_sum_D1LM_ij(L, M, k, k0t, n1, n2, b1, b2, eta, A, c):
    """ Helper function for D1_LM_ij(). """ 
    shape = n1.shape[0]
    gamma_g = np.zeros(shape, np.complex128)
    kgt_r = np.zeros(shape, np.float64)
    kgt_phi = np.zeros(shape, np.float64)
    
    Bgt = n1.reshape(shape, 1) * b1.reshape(1, 2) + n2.reshape(shape, 1) * b2.reshape(1, 2)
    kgt = k0t + Bgt
    for idx in range(shape):
        kgt_r[idx] = np.linalg.norm(kgt[idx, :])
        kgt_phi[idx] = np.arctan2(kgt[idx, 1], kgt[idx, 0])
    
    gamma_g = np.sqrt(k ** 2 - kgt_r ** 2)      
    gamma_g[np.where(gamma_g == complex(0))[0]] += 1e-10j            
    eiMphi = np.exp(1j * M * kgt_phi)
    emnikgtc = np.exp(-1j * np.dot(kgt, c[:2]))   
    
    n = np.arange(0, L - abs(M) + 1)
    inner_sum = np.zeros((shape, n.shape[0]), np.complex128)
    for nn in n:    
        s = np.arange(nn, min(L - abs(M), 2 * nn) + 1)
        if not (L - abs(M)) % 2:
            s = s[np.argwhere(s % 2 == 0).flatten()]
        else:
            s = s[np.argwhere(s % 2).flatten()]
            
        if not s.shape[0] == 0:    
            denom = pbeh.numba_factorial(2 * nn - s)[:, 0] * pbeh.numba_factorial(2 * nn - s)[:, 1] \
                    * pbeh.numba_factorial(s - nn)[:, 0] * pbeh.numba_factorial(s - nn)[:, 1] \
                    * pbeh.numba_factorial(np.asarray((L + abs(M) - s) / 2, np.int32))[:, 0] \
                    * pbeh.numba_factorial(np.asarray((L + abs(M) - s) / 2, np.int32))[:, 1] \
                    * pbeh.numba_factorial(np.asarray((L - abs(M) - s) / 2, np.int32))[:, 0] \
                    * pbeh.numba_factorial(np.asarray((L - abs(M) - s) / 2, np.int32))[:, 1]
            nom = (-k * c[2]) ** (2 * nn - s) * (kgt_r / k).reshape(shape, 1) ** (L - s)
            inner_sum[:, nn] = np.sum(nom / denom, axis=1)
            
    delta = pbeh.delta_n(n[-1], gamma_g, c[2], eta)
    return np.sum(emnikgtc * eiMphi * np.sum((gamma_g / k).reshape(shape, 1) ** (2 * n - 1) * delta * inner_sum, axis=1))    


@njit()
def D2_LM_ij(L, M, c, k, k0t, a1, a2, eta):
    """ Ewald sum's real space summand to account for the coupling between one particle (S_0i)
        and a different particle's (S_0j) periodic arrangement. Kambe, Z. Naturforschung 22a, (46b).
    Args:
        L (int):                multipole degree
        M (int):                multipole order
        c (numpy.ndarray):      displacement vector betwenn particle S0j and S0j' in Carthesian coordinates
        k (complex):            wavenumber
        k0t (numpy.ndarray):    complex in-plane wave vector in Carthesian coordinates 
        a1 (numpy.ndarray):     lattice vector 1 in Carthesian coordinates
        a2 (numpy.ndarray):     lattice vector 2 in Carthesian coordinates
        eta (float):            Ewald sum separation parameter 
    Returns:
        D2_LM (complex):        real space summand of Ewald sum D_LM
    """
    if c[2] == 0 and (L - abs(M)) % 2:
        return 0
  
    prefac = -1j * np.sqrt(2 / np.pi)
      
    x=0
    n1, n2 = pbeh.n1n2_indices(x)
    D2_LM = compute_sum_D2LM_ij(L, M, k, k0t, n1, n2, a1, a2, eta, c)
    
    flag = 0
    while flag == 0:
        D2_LM_prime = D2_LM       
        x += 1
        n1, n2 = pbeh.n1n2_indices(x)
        D2_LM += compute_sum_D2LM_ij(L, M, k, k0t, n1, n2, a1, a2, eta, c)       
        if not D2_LM_prime == complex(0):
            if np.linalg.norm(np.array([D2_LM - D2_LM_prime])) / np.linalg.norm(np.array([D2_LM_prime])) < 1e-5:
                flag = 1
    return prefac * D2_LM


@njit()
def compute_sum_D2LM_ij(L, M, k, k0t, n1, n2, a1, a2, eta, c):
    """ Helper function for D2_LM_jjp(). """ 
    shape = n1.shape[0]
    mnRnc_r = np.zeros(shape, np.float64)
    mnRnc_theta = np.zeros(shape, np.float64)
    mnRnc_phi = np.zeros(shape, np.float64)
    
    Rn = n1.reshape(shape, 1) * a1[:2].reshape(1, 2) + n2.reshape(shape, 1) * a2[:2].reshape(1, 2)
    mnRnc = -(n1.reshape(shape, 1) * a1.reshape(1, 3) + n2.reshape(shape, 1) * a2.reshape(1, 3) + c)
    for idx in range(shape):
        mnRnc_r[idx] = np.linalg.norm(mnRnc[idx, :])
        mnRnc_theta[idx] = np.arctan2(np.linalg.norm(mnRnc[idx, :2]), mnRnc[idx, 2])
        mnRnc_phi[idx] = np.arctan2(mnRnc[idx, 1], mnRnc[idx, 0])
             
    eik0tRn = np.exp(1j * np.dot(k0t, np.transpose(Rn)))
    Y_LM = sf.legendre_normalized_numbed(np.cos(mnRnc_theta), np.sin(mnRnc_theta),
                               max(1, L))[0][L, abs(M)] * np.exp(1j * M * mnRnc_phi) / pbeh.normalization_spherical_harmonics(M)
    
    int_val = (1 / 2) ** (L + 1.5) * pbeh.int_recursion(L, eta, k, mnRnc_r)
    return  np.sum(eik0tRn * (k * mnRnc_r) ** L * Y_LM * int_val)