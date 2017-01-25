__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import graphtools as gt
import numpy as np
import matplotlib.pyplot as plt
import time as tm
from scipy import stats as st
from os import listdir

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def generate_angular_data():
    dmat = np.zeros([crds.shape[0],int((Na*(Na+1))/2)],np.float)


    return

dtdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_01/testdata/'

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk'
cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

files = listdir(dtdir)

_timeloop = tm.time()
#for i in files:
    # Read NC data
xyz,typ,Eact_W,readf = gt.readncdat(dtdir + files[0])

print (xyz[0,1]-xyz[0,1])

dist = gt.generatedmats(xyz,len(typ))

    #print("Distances: \n",dist)

_timeloop2 = (tm.time() - _timeloop)
print('Computation complete. Time: ' + "{:.4f}".format(_timeloop2)  + 'ms')