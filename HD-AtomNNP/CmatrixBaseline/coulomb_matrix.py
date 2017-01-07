from __future__ import print_function
import numpy as np
#from ase.io import read

# toy example of reometry dead
# just replace with your reader
#geometry = read('/data/work/ani/NEB/water.xyz')

# provide 2D array of coordinates
#xyz = geometry.positions

# provide list of chemical symbols
#atoms = geometry.get_chemical_symbols()


def periodicfunc(element):
    """
    Helper function to output atomic number for each element in the periodic table
    """
    
    # replace with your path if necessary
    f = open("CoulombMatrix/pt.txt")
    atomicnum = [line.split()[1] for line in f if line.split()[0] == element]
    f.close()
    return int(atomicnum[0])


def CMatrix(xyzmatrix, atomlist, dim, sort=True):
    """
    This function takes in an xyz, number of atoms in the biggest molecule and computes sorted coulomb Matrix 
    """

    xyzheader = len(atomlist)
    
    cij=np.zeros((dim,dim))
    chargearray = np.zeros((xyzheader,1))
    
    chargearray = [periodicfunc(symbol)  for symbol in atomlist]
    
    for i in range(0,xyzheader):
        for j in range(i,xyzheader):
            if i == j:
                cij[i,j] = 0.5*chargearray[i]**2.4
            else:
                dist = np.linalg.norm(xyzmatrix[i,:] - xyzmatrix[j,:])
                cij[j,i] = cij[i,j] = chargearray[i]*chargearray[j]/dist
    
    if sort==True:
        summation = np.array([sum(x**2) for x in cij])
        sorted_mat = cij[np.argsort(summation)[::-1,],:]
        return sorted_mat.ravel()
    
    else: 
        return cij.ravel()


def GenCMatData(xyz,typ,N):

    M = xyz.shape[0]
    cmats = np.empty([xyz.shape[0], N * N + 1], dtype=np.float32)

    for i in range(0,M):
        cmats[i,0:N*N] = CMatrix(xyz[i],typ,N)

    return cmats


def GenCMatData2(xyz,typ,N):

    M = xyz.shape[0]
    cmats = np.empty([xyz.shape[0], N * N], dtype=np.float32)

    for i in range(0,M):
        cmats[i,0:N*N] = CMatrix(xyz[i],typ,N)

    return cmats

def computerISE(typ):
    H= -0.500607632585
    C= -37.8302333826
    N= -54.5680045287
    O= -75.0362229210

    Eshf = 0.0
    for i in typ:
        if i == 'H':
            Eshf += H
        elif i == 'C':
            Eshf += C
        elif i == 'N':
            Eshf += N
        elif i == 'O':
            Eshf += O

    return Eshf


#c = CMatrix(xyz, atoms, 4)
#print(c)
