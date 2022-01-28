"""Test Ewald sum with reference values provided by Dominik Beutel."""

import numpy as np
import smuthi.periodicboundaries.ewald_lattice_sums as pbe


test_case = {'1':{'k': complex(1.2),
                  'k0t': np.array([0, 0], dtype=float),
                  'a1': np.array([1, 0, 0], dtype=float),
                  'a2': np.array([0, 1, 0], dtype=float),
                  'c': np.array([0, 0, 0], dtype=float),
                  'line0': 9},
             '2':{'k': 2 + 0.1j,
                  'k0t': np.array([0.1, -0.15], dtype=float),
                  'a1': np.array([1, 0, 0], dtype=float),
                  'a2': np.array([0, 1.1, 0], dtype=float),
                  'c': np.array([0, 0, 0], dtype=float),
                  'line0': 33},
             '3':{'k': complex(1.2),
                  'k0t': np.array([0, 0], dtype=float),
                  'a1': np.array([1, 0, 0], dtype=float),
                  'a2': np.array([0, 1, 0], dtype=float),
                  'c': np.array([0.4, -0.3, 0], dtype=float),
                  'line0': 58},
             '4':{'k': 2 + 0.1j,
                  'k0t': np.array([0.1, -1.5], dtype=float),
                  'a1': np.array([1, 0, 0], dtype=float),
                  'a2': np.array([0, 1.1, 0], dtype=float),
                  'c': np.array([0.4, -0.3, 0], dtype=float),
                  'line0': 84}}




def testEwaldSum(case):
    
    k = test_case[str(case)]['k']
    k0t = test_case[str(case)]['k0t']
    a1 = test_case[str(case)]['a1']
    a2 = test_case[str(case)]['a2']
    c = test_case[str(case)]['c']
    line0 = test_case[str(case)]['line0']
    
    A = np.linalg.norm(np.cross(a1, a2)) 
    eta = pbe.separation_parameter_eta(k, k0t, a1, a2, magM=1)

    d1LM = np.zeros([15, 2], dtype=complex)
    d2LM = np.zeros([15, 2], dtype=complex)

    ''' import reference data '''
    filename = 'test_ewaldsum_reference_values.txt'  
    with open(filename, 'r') as n_file:
        lines = n_file.readlines()  
      
    for idx, line in enumerate(lines[line0:line0+15]):
        line = line.replace('i', 'j')
        split_line = line.split()
        d1LM[idx, 1] = complex(split_line[2])
        d2LM[idx, 1] = complex(split_line[3])


    ''' compute DLM '''
    idx = 0
    for L in [0, 1, 2, 3, 4]:
        for M in range(-L, L+1):
            if not (L - abs(M)) % 2:
                if np.linalg.norm(c) == 0:
                    d1LM[idx, 0] = pbe.D1_LM(L, M, k, k0t, a1, a2, eta, A)
                    d2LM[idx, 0] = pbe.D2_LM(L, M, k, k0t, a1, a2, eta)
                    idx += 1
                else:
                    d1LM[idx, 0] = pbe.D1_LM_ij(L, M, c, k, k0t, a1, a2, eta, A)
                    d2LM[idx, 0] = pbe.D2_LM_ij(L, M, c, k, k0t, a1, a2, eta)
                    idx += 1
    
           
    np.testing.assert_allclose(d1LM[:, 0], d1LM[:, 1], rtol=1e-4, atol=1e-12)
    np.testing.assert_allclose(d2LM[:, 0], d2LM[:, 1], rtol=1e-4, atol=1e-12)
            
        
if __name__ == '__main__':
    testEwaldSum(case=1)
    # testEwaldSum(case=2)
    testEwaldSum(case=3)
    # testEwaldSum(case=4)