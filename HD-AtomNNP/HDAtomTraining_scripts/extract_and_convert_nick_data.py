import numpy as np
import hdnntools as hdn
from pyNeuroChem import cachegenerator as cg
import re

def remove_bad_data(file, data):
    f = open(file,'r').read()

    reg = re.compile(' +?(\d+?) +?\S+? +?\S+ +?\S+? +?\d+?\.\d+?\s*?\n')
    m = reg.findall(f)

    index = []
    for i in m:
        index.append(int(i)-1)

    return np.delete(data,index)

#xyz_file = '/home/jujuman/dataset-qm9/data-qm9R.npy'
#eng_file = '/home/jujuman/dataset-qm9/data-qm9U0.npy'
#spc_file = '/home/jujuman/dataset-qm9/data-qm9S.npy'
#atn_file = '/home/jujuman/dataset-qm9/data-qm9Z.npy'

xyz_file = '/home/jujuman/Research/QM-7TEST/dataset-qm7/data-qm7R.npy'
eng_file = '/home/jujuman/Research/QM-7TEST/dataset-qm7/data-qm7T.npy'
spc_file = '/home/jujuman/Research/QM-7TEST/dataset-qm7/data-qm7Z.npy'

xyz = np.load(xyz_file)
eng = np.load(eng_file)[0]
spc = np.load(spc_file)
#atn = np.load(atn_file)

store_dir = "/home/jujuman/Research/QM-7TEST/tester/"
saef   = "/home/jujuman/Research/QM-7TEST/tester/sae_6-31gd.dat"

data_index = np.array(range(eng.shape[0]))

#unc_file = '/home/jujuman/dataset-qm9/uncharacterized.txt'
#data_index = remove_bad_data(unc_file, data_index)
np.random.shuffle(data_index)
print(data_index.shape)

listt = data_index[:int(0.8*len(data_index))]
listv = data_index[int(0.8*len(data_index)):int(0.9*len(data_index))]
listte = data_index[int(0.9*len(data_index)):]

cachet = cg('_train', saef, store_dir)
cachev = cg('_valid', saef, store_dir)

eng = eng / hdn.hatokcal

print('max: ', eng.max(), ' min: ', eng.min())

for n,i in enumerate(listt):
    print(n)
    x = xyz[i]
    e = eng[i]
    z = spc[i]
    #z = atn[i]

    z = z[~((z == 0))]
    Na = z.shape[0]

    xyz_t = np.array(x[0:Na], dtype = np.float32, order='C').reshape(1,Na,3)
    spc_t = [hdn.convertatomicnumber(i) for i in z]

    #print(xyz_t)
    #print(spc_t)

    #if 'F' not in spc_t:
    cachet.insertdata(xyz_t, np.array([e],dtype=np.float64), list(spc_t))

for i in listv:
    x = xyz[i]
    e = eng[i]
    z = spc[i]
    #z = atn[i]

    z = z[~((z == 0))]
    Na = z.shape[0]

    xyz_t = np.array(x[0:Na], dtype = np.float32, order='C').reshape(1,Na,3)
    #spc_t = s[0:Na]
    spc_t = [hdn.convertatomicnumber(i) for i in z]

    #if 'F' not in spc_t:
    cachev.insertdata(xyz_t, np.array([e],dtype=np.float64), list(spc_t))

cachet.makemetadata()
cachev.makemetadata()

import pyanitools as pyt

path = "/home/jujuman/Scratch/Research/QM-7TEST/QM7-test-ho.h5"
dpack = pyt.datapacker(path)

print('Storing test data...')
for i in listte:
    x = xyz[i]
    e = eng[i]
    z = spc[i]
    #z = atn[i]

    z = z[~((z == 0))]
    Na = z.shape[0]
    spc_t = [hdn.convertatomicnumber(i) for i in z]

    xyz_t = np.array(x[0:Na], dtype = np.float32, order='C').reshape(1,Na*3)
    #spc_t = s[0:Na]
    spc_t = np.array([hdn.convertatomicnumber(i) for i in z])

    # Prepare and store the data
    dpack.store_data('QM9' + "/mol" + str(i), xyz_t, np.array(e), spc_t)

dpack.cleanup()