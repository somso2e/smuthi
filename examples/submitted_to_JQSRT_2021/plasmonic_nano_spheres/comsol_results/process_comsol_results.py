import numpy as np

comsol_data = np.loadtxt("2Au_spheres_on_silicon_wl_dsweep_2nm_gap_data.txt", comments="%")
comsol_data2 = np.loadtxt("2Au_spheres_on_silicon_wl_dsweep_2nm_gap_intermediate_steps_data.txt", comments="%")
comsol_data = np.append(comsol_data, comsol_data2, axis=0)
rows = comsol_data[:, 1] > 3  # skip 2nm gap data
comsol_data = comsol_data[rows, :]
np.savetxt("comsol_results.dat", np.transpose([comsol_data[:, 0] / 1000, comsol_data[:, 1], comsol_data[:, 7]]), header="wavelength distance enhancement")