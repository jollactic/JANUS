#!/usr/bin/env python3
import numpy as np
import argparse
import sys
import copy
from scipy.spatial import ConvexHull
import ase
from ase.io import Trajectory, read, write


def find_facets(atoms, facet=0, mirror=0, offset=[0.0, 0.0, 0.0], angle=0):
    pos = atoms.get_positions()
    # FIND SURFACE
    Hull = ConvexHull(pos)

    atoms_out = copy.deepcopy(atoms)

    # MERGE FACETS
    N = len(Hull.simplices)
    facets_add = np.zeros(N)
    facets_matrix = np.zeros([N, N])
    i = -1
    for tmp1 in Hull.equations:
        i = i + 1
        j = -1
        for tmp2 in Hull.equations:
            j = j + 1
            if(np.dot(tmp1[0:3], tmp2[0:3]) > 0.95 and np.dot(tmp1[0:3], tmp2[0:3]) < 1.05 and j >= i):
                if(facets_add[j] == 0):
                    facets_add[j] = 1
                    facets_matrix[i][j] = 1
    facets = []
    normals = []
    nfacets = 0
    for i in range(N):
        if(np.linalg.norm(facets_matrix[i][0:-1]) > 0.5):
            nfacets = nfacets + 1
            facets_list = []
            tmp_normal = np.zeros(3)
            for j in range(N):
                if(facets_matrix[i][j] == 1):
                    facets_list.append(Hull.simplices[j])
                    tmp_normal = tmp_normal + Hull.equations[j][0:3]
            tmp_normal = tmp_normal / np.linalg.norm(tmp_normal)
            facets.append(facets_list)
            normals.append(tmp_normal)

    # /MERGE FACETS

    # ALIGN NORMAL TO X
    a = normals[facet]
    if(abs(np.dot(a, [1, 0, 0])) < 0.95):
        b = np.cross(a, [1, 0, 0])
        b = b / np.linalg.norm(b)
        c = np.cross(a, b)
    else:
        b = np.cross(a, [0, 1, 0])
        b = b / np.linalg.norm(b)
        c = np.cross(a, b)

    N = [a, b, c]
    # N=np.transpose(N)
    N = np.linalg.inv(N)
    pos = np.dot(pos, N)
    # /ALIGN NORMAL TO X

    # CALCULATE AREA
    area_points = []
    nc = 0
    center = [0., 0., 0.]
    for tri in facets[facet]:
        for p in tri:
            nc = nc + 1
            center = center + pos[p]
            area_points.append(pos[p])
            area_points.append(pos[p] + [1., 0., 0.])
    HullFace = ConvexHull(np.asarray(area_points))
    face_area = HullFace.volume
    # /CALCULATE AREA

    # /SHIFT SURFACE AND ROTATE ABOUT X
    center = center / nc
    for j in range(atoms.get_number_of_atoms()):
        pos[j] = pos[j] - center

    cos_theta = np.cos(angle * np.pi / 180.)
    sin_theta = np.sin(angle * np.pi / 180.)
    Rmat = [[1.0, 0, 0], [0.0, cos_theta, sin_theta],
            [0.0, -sin_theta, cos_theta]]
    pos = np.matmul(pos, Rmat)

    for j in range(atoms.get_number_of_atoms()):
        pos[j] = pos[j] - [offset]
        if(mirror == 1):
            pos[j][0] = -pos[j][0]
    # /SHIFT SURFACE AND ROTATE ABOUT X

    atoms_out = copy.deepcopy(atoms)
    atoms_out.set_positions(pos)

    return atoms_out, nfacets, face_area


def janus(atomsA, atomsB, rx, ryz):
    traj = Trajectory('JANUS.traj', 'w')

    tmp, NA, AA = find_facets(atomsA)
    tmp, NB, AB = find_facets(atomsB, mirror=1)

    # FIND LARGEST AREA OF PARTICLE A
    AArec = 0
    for i in range(NA):
        atomsAout, tmp, AA = find_facets(atomsA, facet=i)
        if(AA > AArec):
            AArec = AA
            largest_AA = i

    # FIND LARGEST AREA OF PARTICLE B
    ABrec = 0
    for i in range(NB):
        atomsBout, tmp, AB = find_facets(atomsB, facet=i)
        if(AB > ABrec):
            ABrec = AB
            largest_AB = i

    print(
        "--- Area " +
        str(largest_AA) +
        " : " +
        str(AA) +
        " --- Area " +
        str(largest_AB) +
        " : " +
        str(AB))

    for i in range(12):
        ang = i * 15
        off = [rx, 0., 0.]
        atomsAout, tmp, AA = find_facets(
            atomsA, facet=largest_AA, angle=ang, offset=off)
        atomsBout, tmp, AB = find_facets(atomsB, facet=largest_AB, mirror=1)
        atomsBout.extend(atomsAout)
        traj.write(atomsBout)
        for j in range(12):
            ang_disp = j * 30
            off = [
                rx,
                ryz *
                np.cos(
                    ang_disp *
                    np.pi /
                    180.),
                ryz *
                np.sin(
                    ang_disp *
                    np.pi /
                    180.)]
            atomsAout, tmp, AA = find_facets(
                atomsA, facet=largest_AA, angle=ang, offset=off)
            atomsBout, tmp, AB = find_facets(
                atomsB, facet=largest_AB, mirror=1)
            atomsBout.extend(atomsAout)
            traj.write(atomsBout)


# MAIN
def main():

    struct = sys.argv[1]
    atomsA = read(struct)

    struct = sys.argv[2]
    atomsB = read(struct)

    rx = float(sys.argv[3])  # SEPARATION OF PARTICLES ALONG X
    ryz = float(sys.argv[4])  # MAGNITUDE OF DISPLACMET IN THE YZ-PLANE

    janus(atomsA, atomsB, rx, ryz)


if __name__ == '__main__':
    main()
