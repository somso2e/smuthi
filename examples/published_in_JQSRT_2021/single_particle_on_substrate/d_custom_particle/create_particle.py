import numpy as np
import trimesh

mesh = trimesh.load("sphere_fine.stl")

#mesh.show()

mesh.vertices -= mesh.center_mass

for vertex in mesh.vertices:
    x = vertex[0]
    y = vertex[1]
    z = vertex[2]
    r = np.linalg.norm(vertex)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    fac = (1 - 0.15 * np.cos(2 * theta) + 0.15 * np.sin(theta)**2 * np.cos(-theta) * np.sin(3 * phi) + 0.15 * np.sin(theta)**2 * np.cos(theta) * np.sin(4 * phi))
    x = x * fac
    y = y * fac
    z = z * fac
    vertex[0], vertex[1], vertex[2] = x, y, z

mesh.vertices -= mesh.center_mass
vol = mesh.volume
target_vol = 0.4**3 * np.pi * 4 / 3
mesh.vertices *= (target_vol / vol) ** (1 / 3)

min_z = 1000
for vertex in mesh.vertices:
    if vertex[2] < min_z:
        min_z = vertex[2]

print("volume:", mesh.volume)
print("min z:", min_z)
#mesh.show()
mesh.export("freeform_particle_vol%.3f_minz%.3f.stl"%(mesh.volume, min_z), "stl_ascii")