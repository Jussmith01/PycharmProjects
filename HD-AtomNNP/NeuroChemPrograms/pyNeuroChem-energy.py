__author__ = 'jujuman'

# Import pyNeuroChem
import pyNeuroChem as pync
import numpy as np
import hdnntools as gt
import time

import matplotlib.pyplot as plt

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Dropbox/ChemSciencePaper.AER/networks/ANI-c08e-ntwk/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

dtdir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/Poster-GTC-May-2017/Timings/'
#xyz,typ,Eact,tmp    = gt.readncdat(dtdir + 'gdb11_s01-1_test.dat',np.float32)
#xyz,typ,Na = gt.readxyz(dtdir + 'benzene.xyz')

#print (xyz)

# Construct pyNeuroChem class
mol = pync.molecule(cnstfile, saefile, nnfdir, 0)
#con = pync.conformers(cnstfile, saefile, nnfdir, 0)

params = mol.getntwkparams(0,3)

print(params["weights"].shape)

print(params["weights"])
print(params["biases"])

# Set the conformers in NeuroChem

'''
mol.setMolecule(coords=xyz[0],types=typ)
con.setConformers(confs=xyz,types=typ)

# Print some data from the NeuroChem
print( 'Number of Atoms Loaded: ' + str(mol.getNumAtoms()) )
print( 'Number of Confs Loaded: ' + str(mol.getNumConfs()) )

#print ('Optimizing...')
#O = nc.optimize(conv=0.000001,max_iter=250)
#print(O)

#print(len(typ))
#for i in range(0,len(typ)):
#    print (typ[i] + ' ' + "{:.10f}".format(O[i,0]) + ' ' + "{:.10f}".format(O[i,1]) + ' ' + "{:.10f}".format(O[i,2]))


# Compute Energies of Conformations
print('Calculating energies and forces...')

E1 = mol.energy()
F1 = mol.force()

start_time = time.time()
E1 = mol.energy()
molE_t = time.time() - start_time
print('[ANI-MOL Energy time:', molE_t, 'seconds]')

start_time = time.time()
F1 = mol.force()
molF_t = time.time() - start_time
print('[ANI-MOL Force time: ', molF_t, 'seconds]')

start_time = time.time()
E2 = con.energy()
conE_t = time.time() - start_time
print('[ANI-CON Energy time:', conE_t, 'seconds]')

start_time = time.time()
F2 = con.force()
conF_t = time.time() - start_time
print('[ANI-CON Force time: ', conF_t, 'seconds]')

print('\nDELTA E |MOL-CON|:')
print('Error: ', abs(E1-E2).sum())
print('dtime: ', (conE_t - molE_t)*1000.0,'ms')

print('\nDELTA F sum(|MOL-CON|):')
print('Error: ', (abs(F1-F2)/float(F1.size)).sum())
print('Max Delta: ', (F1-F2).max())
print('dtime: ', (conF_t - molF_t)*1000.0,'ms')

print('\nCON E+F time: ', (conF_t+conE_t))
print('MOL E+F time: ', (molE_t+molF_t))
print('\nCON E+F time: ', (conF_t+conE_t)*1000.0,'ms')
print('MOL E+F time: ', (molE_t+molF_t)*1000.0,'ms')
print('E+F dtime: ', ((conF_t+conE_t) - (molE_t+molF_t))*1000.0,'ms')

AE = nc.aenergies()

print('\nAE:')
print(AE)
print(typ)

AV1 = nc.activations(atom_idx=0, layer_idx=2, molec_idx=0)

print('\nAV1:')
print(AV1)

AV2 = nc.activations(atom_idx=0, layer_idx=2, molec_idx=1)

print('\nAV2:')
print(AV2)

AV3 = AV1 - AV2

print('\nAV3:')
print(AV3)

IDX = np.arange(0,AV3.shape[0],1,dtype=float) + 1
plt.plot (IDX,AV1,marker='o', color='red',  label='ANI-X',  linewidth=2)
plt.plot (IDX,AV2,marker='o', color='blue',  label='ANI-X',  linewidth=2)

plt.title("C10H20 - ANI vs DFT")

plt.ylabel('E cmp (kcal/mol)')
plt.xlabel('E act (kcal/mol)')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=16)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

plt.show()
'''