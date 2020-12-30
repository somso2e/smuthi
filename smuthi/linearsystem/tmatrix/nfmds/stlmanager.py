import numpy as np


def readstl(stlname):
    """
    Reads surface information from STL file
    Args:
        stlname (string): name of STL file

    Returns:
        A list of dictionaries with information about faces of scatterer geometry.
    """
    surfaces = []
    with open(stlname) as f:
        nsurf = 0
        for line in f:
            if "solid" in line and "endsolid" not in line:
                nsurf += 1
                if nsurf > 1:
                    surf = {'nsurf': nsurf - 1,
                            'baris': np.array(baris),
                            'areas': np.array(areas),
                            'normals': np.array(normals)}
                    surfaces.append(surf)
                baricenter = 0
                baris = []
                areas = []
                triangle = []
                normals = []
            if "facet" in line and "end" not in line:
                normal = np.array([float(i) for i in line.split(" ")[2:]])
                normals.append(normal)
                if np.any(np.abs(baricenter) > 0) and bool(triangle):
                    ab = np.array(triangle[1] - triangle[0])
                    ac = np.array(triangle[2] - triangle[0])
                    areas.append(np.linalg.norm(np.cross(ab, ac)) / 2)
                    baris.append(baricenter)
                baricenter = 0
                triangle = []
            if "vertex" in line:
                vertex = np.array([float(i) for i in line.split()[1:]])
                baricenter += vertex / 3
                triangle.append(vertex)

    surf = {'nsurf': nsurf - 1,
            'baris': np.array(baris),
            'areas': np.array(areas),
            'normals': np.array(normals)}
    surfaces.append(surf)
    return surfaces


def writefem(femname,surfaces):
    """
    Writes information about particle geometry to FEM file.
    Args:
        femname (string):    name of FEM file
        surfaces  (list):    information about faces of scatterer geometry
    """    
    fid = open(femname, 'w+')
    lenareas = 0
    print('%7i' % len(surfaces), file=fid)
    ii0 = 1
    for surface in surfaces:
        areas = surface['areas']
        F = surface['normals']
        P = surface['baris']
        print('%7i' % len(areas), file=fid)
        for ii in range(len(areas)):
            print('%7i%17.7E%17.7E%17.7E%17.7E%17.7E%17.7E%17.7E' % (
            ii0 + ii, P[ii, 0], P[ii, 1], P[ii, 2], F[ii, 0], F[ii, 1], F[ii, 2], areas[ii]), file=fid)
        ii0 += len(areas)
    fid.close()


def convert_stl_to_fem(stlname, femname):
    """
    Converts STL to FEM file
    Args:
        stlname (string):    name of STL file
        femname (string):    name of FEM file
    """    
    try:
        surfaces=readstl(stlname)
    except(UnicodeDecodeError):
        raise Exception("Binary STL files are not supported.")
    writefem(femname,surfaces)
