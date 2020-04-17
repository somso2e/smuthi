import numpy as np
import smuthi.fields as flds


def single_index_to_multi_nfmds(index, Nrank, Mrank):
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


def multi_index_to_single_smuthi(tau, l, m, Nrank, Mrank, l_max=None, m_max=None):
    tau_blocksize = Mrank * (Mrank + 2) + (Nrank - Mrank) * (2 * Mrank + 1)
    if l_max is None:
        l_max=Nrank
        m_max=Mrank
    else:
        m_max=l_max # for now
    tau_blocksizesm = m_max * (m_max + 2) + (l_max - m_max) * (2 * m_max + 1)
    n = tau * tau_blocksizesm
    if (l - 1) <= Mrank:
        n += (l - 1) * (l - 1 + 2)
    else:
        n += Mrank * (Mrank + 2) + (l - 1 - Mrank) * (2 * Mrank + 1)
    n += m + min(l, Mrank)
    return n


def nfmds_to_smuthi_index(ii, Nrank, Mrank):
    tau, l, m = single_index_to_multi_nfmds(ii, Nrank, Mrank)
    ii_sm = multi_index_to_single_smuthi(tau, l, m, Nrank, Mrank)
    return ii_sm


def nfmds_to_smuthi_matrix(T, Nrank=None, Mrank=None,l_max=None,m_max=None):
    if Mrank is None:
        Mrank = Nrank
    if Nrank is None:
        nmax = T.shape[0] / 2
        Nrank = int(-1 + np.sqrt(1 + nmax))
        Mrank = int(Nrank)
    if l_max is not None:
        m_max=l_max # for now
        tau_blocksizesm = m_max * (m_max + 2) + (l_max - m_max) * (2 * m_max + 1)
    else:
        tau_blocksizesm = Mrank * (Mrank + 2) + (Nrank - Mrank) * (2 * Mrank + 1)
        l_max=Nrank
    Tsm = np.zeros((2*tau_blocksizesm,2*tau_blocksizesm),dtype=complex)
    for ii in range(T.shape[0]):
        tau, l_ii, m_ii = single_index_to_multi_nfmds(ii, Nrank, Mrank)
        ii_sm = multi_index_to_single_smuthi(tau, l_ii, m_ii, Nrank, Mrank, l_max)
        for jj in range(T.shape[1]):
            tau, l_jj, m_jj = single_index_to_multi_nfmds(jj, Nrank, Mrank)
            jj_sm = multi_index_to_single_smuthi(tau, l_jj, m_jj, Nrank, Mrank, l_max)
            if not(l_ii>l_max or l_jj>l_max):
                Tsm[ii_sm, jj_sm] = T[ii, jj]
    return Tsm

def taxsym_read_tmatrix(t_nfmds, l_max, m_max):
    n_rank = l_max
    m_rank = m_max
    t_matrix = np.zeros((flds.blocksize(l_max, m_max), flds.blocksize(l_max, m_max)), dtype=complex)

    for m in range(-m_max, m_max + 1):
        n_max_nfmds = n_rank - max(1, abs(m)) + 1
        for tau1 in range(2):
            for l1 in range(max(1, abs(m)), l_max + 1):
                n1 = flds.multi_to_single_index(tau=tau1, l=l1, m=m, l_max=l_max, m_max=m_max)
                l1_nfmds = l1 - max(1, abs(m))
                n1_nfmds = 2 * n_rank * abs(m) + tau1 * n_max_nfmds + l1_nfmds
                for tau2 in range(2):
                    for l2 in range(max(1, abs(m)), l_max + 1):
                        n2 = flds.multi_to_single_index(tau=tau2, l=l2, m=m, l_max=l_max, m_max=m_max)
                        l2_nfmds = l2 - max(1, abs(m))
                        n2_nfmds = tau2 * n_max_nfmds + l2_nfmds
                        if abs(m) <= m_rank:
                            if m >= 0:
                                t_matrix[n1, n2] = t_nfmds[n1_nfmds][n2_nfmds]
                            else:
                                t_matrix[n1, n2] = t_nfmds[n1_nfmds][n2_nfmds] * (-1) ** (tau1 + tau2)

    return t_matrix
