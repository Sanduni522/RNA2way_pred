import numpy as np
import glob
import sys

def get_coordinates_pdb(filename, allowed_atoms=None):
    """
    Get coordinates from the first chain in a pdb file
    and return a vectorset with all the coordinates.
    """
    V = []
    f = open(filename)
    lines = f.readlines()
    f.close()
    dict = {}
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        curr_atom_name = line[12:16].strip()
        curr_resi = int(line[22:26])
        if allowed_atoms is not None:
            if curr_resi not in allowed_atoms:
                continue
        if curr_atom_name.startswith("H"):
            continue
        if curr_atom_name.startswith("O2'"):
            continue
        if curr_atom_name.startswith("P"):
            continue
        if curr_atom_name.startswith("OP"):
            continue
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]
        #if curr_resi not in dict:
        #    dict[curr_resi] = []
       # dict[curr_resi].append([curr_atom_name, float(x), float(y), float(z)])
       # print(curr_atom_name)
        V.append(np.asarray([x, y, z], dtype=float))
    V = np.asarray(V)
    # print filename, resi_set, len(resi_set)
    return V


def get_beads_pdb(filename, allowed_atoms=None):
    """
    Get coordinates from the first chain in a pdb file
    and return a vectorset with all the coordinates.
    """
    V = []
    f = open(filename)
    lines = f.readlines()
    f.close()
    dict = {}
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        curr_atom_name = line[12:16].strip()
        curr_resi = int(line[22:26])
        if allowed_atoms is not None:
            if curr_atom_name not in allowed_atoms:
                continue
        if curr_atom_name.startswith("H"):
            continue
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]
        if curr_resi not in dict:
            dict[curr_resi] = []
        dict[curr_resi].append([curr_atom_name, float(x), float(y), float(z)])

    keys = sorted(dict.keys())
    for k in keys:
        data = dict[k]
        avg = [0.0, 0.0, 0.0]
        for d in data:
            avg[0] += d[1]
            avg[1] += d[2]
            avg[2] += d[3]
        avg = np.array(avg) / len(data)
        V.append(avg)
    V = np.asarray(V)
    # print filename, resi_set, len(resi_set)
    return V


def kabsch_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD
    """
    P = rotate(P, Q)
    return rmsd(P, Q)


def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters:
    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix

    Returns:
    U -- Rotation matrix

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U


def rotate(P, Q):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm
    """
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return P


def centroid(X):
    """
    Calculate the centroid from a vectorset X
    """
    C = sum(X) / len(X)
    return C


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
     #   print(v,w)
        rmsd += sum([(v[i] - w[i]) ** 2.0 for i in range(D)])
    return np.sqrt(rmsd / N)


def compute_rmsd_from_coords(Q, P):
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc
    rmsd = round(kabsch_rmsd(P, Q), 2)
    return rmsd


