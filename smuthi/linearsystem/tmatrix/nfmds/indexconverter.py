import numpy as np
import smuthi.fields as flds
from numba import jit


@jit(nopython=True, cache=True)
def single_index_to_multi_nfmds(index, Nrank, Mrank):
    """Converts single index to (tau,l,m) tuple in NFMDS convention.

    Args:
        index (int):    single index in NFMDS convention    
        Nrank (int):    NFMDS Nrank parameter
        Mrank (int):    NFMDS Mrank parameter       

    Returns:
        tau (int):      SVWF polarization (0 for spherical TE, 1 for spherical TM)
        l (int):        SVWF degree
        m (int):        SVWF order    
    """      
    
    index = index + 1
    nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
    if index > nmax:
        tau = 1
        index -= nmax
    else:
        tau = 0

    for m in range(Mrank + 1):
        N0 = (Nrank + (abs(m + 1) - 1) * (2 * Nrank - abs(m + 1) + 2))
        if index <= N0:
            break
    if m > 0:
        N0 = (Nrank + (m - 1) * (2 * Nrank - m + 2))
        index -= N0
    if index > (Nrank - abs(m) + 1):
        index -= (Nrank - abs(m) + 1)
        m = -m
    l = index + (abs(m) - 1) * (abs(m) > 0)
    return tau, l, m


def multi_index_to_single_nfmds(tau, l, m, Nrank, Mrank):
    """Converts a (tau,l,m) index to single index in NFMDS convention.

    Args:
        tau (int):      SVWF polarization (0 for spherical TE, 1 for spherical TM)
        l (int):        SVWF degree
        m (int):        SVWF order        
        Nrank (int):    NFMDS Nrank parameter
        Mrank (int):    NFMDS Mrank parameter              

    Returns:
        index (int):    single index in NFMDS convention 
    """      
    nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
    if m == 0:
        index = l
    else:
        N0 = Nrank + (abs(m) - 1) * (2 * Nrank - abs(m) + 2)

        if m < 0:
            N0 = N0 + Nrank - abs(m) + 1
        index = N0 + l - abs(m) - 1
    if tau:
        index = index + nmax
    return index - 1  # cause python numbers from 0 and fortran from 1


@jit(nopython=True, cache=True)
def nfmds_to_smuthi_matrix(T, Nrank=None, Mrank=None, l_max=None, m_max=None):
    """Converts a T-matrix obtained with NFMDS to SMUTHI compatible format.

    Args:
        T (array):      T-matrix in NFMDS convention
        Nrank (int):    NFMDS Nrank parameter
        Mrank (int):    NFMDS Mrank parameter
        l_max (int):    Maximal multipole degree used for the spherical wave expansion of incoming and
                        scattered field
        m_max (int):    Maximal multipole order used for the spherical wave expansion of incoming and
                        scattered field
       
    Returns:
        Tsm (array):    T-matrix in SMUTHI convention
    """

    if Nrank is None:
        nmax = T.shape[0] / 2
        Nrank = int(-1 + np.sqrt(1 + nmax))
        Mrank = int(Nrank)

    if Mrank is None:
        Mrank = Nrank

    if l_max is None:
        l_max = Nrank

    if m_max is None:
        m_max = l_max

    Tsm = np.zeros((flds.blocksize(l_max, m_max), flds.blocksize(l_max, m_max)), dtype=np.complex128)
    num_cols = range(T.shape[0])
    num_rows = range(T.shape[1])
    for ii in num_cols:
        tau, l_ii, m_ii = single_index_to_multi_nfmds(ii, Nrank, Mrank)
        if l_ii > l_max or abs(m_ii) > m_max:
            continue
        ii_sm = flds.multi_to_single_index(tau, l_ii, m_ii, l_max, m_max)
        for jj in num_rows:
            tau, l_jj, m_jj = single_index_to_multi_nfmds(jj, Nrank, Mrank)
            if l_jj > l_max or abs(m_jj) > m_max:
                continue
            jj_sm = flds.multi_to_single_index(tau, l_jj, m_jj, l_max, m_max)
            Tsm[ii_sm, jj_sm] = T[ii, jj]
    return Tsm


@jit(nopython=True, cache=True)
def python_to_smuthi_matrix(T, Nrank, Mrank=None, l_max=None, m_max=None):
    """Converts a T-matrix obtained with Alan's code to SMUTHI compatible format.

    Args:
        T (array):      T-matrix in NFMDS convention
        Nrank (int):    Alan's lmax parameter
        Mrank (int):    Alan's lmax parameter
        l_max (int):    Maximal multipole degree used for the spherical wave expansion of incoming and
                        scattered field
        m_max (int):    Maximal multipole order used for the spherical wave expansion of incoming and
                        scattered field

    Returns:
        Tsm (array):    T-matrix in SMUTHI convention
    """

    if Mrank is None:
        Mrank = Nrank

    if l_max is None:
        l_max = Nrank

    if m_max is None:
        m_max = l_max

    Tsm = np.zeros((flds.blocksize(l_max, m_max), flds.blocksize(l_max, m_max)), dtype=np.complex128)

    mlim = min(m_max, Mrank)
    llim = min(l_max, Nrank)

    for tau1 in range(2):
        for m1 in range(-mlim, mlim + 1):
            for l1 in range(max(1, abs(m1)), llim + 1):
                n1_py = flds.multi_to_single_index(tau1, l1, m1, Nrank, Mrank)
                n1_sm = flds.multi_to_single_index(tau1, l1, m1, l_max, m_max)
                for tau2 in range(2):
                    for m2 in range(-mlim, mlim + 1):
                        for l2 in range(max(1, abs(m2)), llim + 1):
                            n2_py = flds.multi_to_single_index(tau2, l2, m2, Nrank, Mrank)
                            n2_sm = flds.multi_to_single_index(tau2, l2, m2, l_max, m_max)
                            Tsm[n1_sm, n2_sm] = T[n1_py, n2_py]

    return Tsm
