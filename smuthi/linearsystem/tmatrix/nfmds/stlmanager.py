import numpy as np


def convert_stl_to_fem(stlname, femname):
    """Reads an STL file generated with GMSH and converts it to FEM (NFMDS geometry file) format.

    Args:
        stlname(string):    Path to STL file to be read
        femname(string):    Path to FEM file to be written

    """

    surfaces = []
    with open(stlname) as f:
        for x in f:
            if "solid" in x and "end" not in x:
                nsurf = int(x.split(" ")[-1])
                print(nsurf)
                if nsurf > 1:
                    surf = {
                        'nsurf': nsurf - 1,
                        'baris': np.array(baris),
                        'areas': np.array(areas),
                        'normals': np.array(normals)
                    }
                    surfaces.append(surf)
                baricenter = 0
                baris = []
                areas = []
                triangle = []
                normals = []
            if "facet" in x and "end" not in x:
                normal = np.array([float(i) for i in x.split(" ")[2:]])
                normals.append(normal)
                if np.any(np.abs(baricenter) > 0) and bool(triangle):
                    ab = np.array(triangle[1] - triangle[0])
                    ac = np.array(triangle[2] - triangle[0])
                    areas.append(np.linalg.norm(np.cross(ab, ac)) / 2)
                    baris.append(baricenter)
                baricenter = 0
                triangle = []
            if "vertex" in x:
                vertex = np.array([float(i) for i in x.split(" ")[5:]])
                baricenter += vertex / 3
                triangle.append(vertex)
    surf = {
        'nsurf': nsurf - 1,
        'baris': np.array(baris),
        'areas': np.array(areas),
        'normals': np.array(normals)
    }
    surfaces.append(surf)

    fid = open(femname, 'w+')
    print('%7i' % 1, file=fid)
    print('%7i' % len(areas), file=fid)
    ii0 = 1
    for surface in surfaces:
        areas = surface['areas']
        F = surface['normals']
        P = surface['baris']
        for ii in range(len(areas)):
            print('%7i%17.7E%17.7E%17.7E%17.7E%17.7E%17.7E%17.7E' % (
            ii0 + ii, P[ii, 0], P[ii, 1], P[ii, 2], F[ii, 0], F[ii, 1], F[ii, 2], areas[ii]), file=fid)
        ii0 += len(areas)
    fid.close()
