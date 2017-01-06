from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
import random
import numpy as np
import re

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl

import pyNeuroChem as pync

def formatsmilesfile(file):
    ifile = open(file, 'r')
    contents = ifile.read()
    ifile.close()    

    p = re.compile('([^\s]*).*\n')
    smiles = p.findall(contents)
    
    ofile = open(file, 'w')
    for mol in smiles:    
        ofile.write(mol + '\n')
    ofile.close()

#-------- Parameters -----------

R = 0.3
fpf = 'gdb11_s10' #Filename prefix
#wdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_10/' #working directory
smfile = '/home/jujuman/Research/ANN-Test-Data/GDB-11-B3LYP-6-31gd/smiledata/gdb11_size06.smi' # Smiles file
At = ['C', 'O', 'N'] # Hydrogens added after check
P = 1.0

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/smallAEV_testing/train_384-256-128-64-1_c08e/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

#------- End Parameters ---------

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 0)

#fix the file
formatsmilesfile(smfile)

#molecules = Chem.SmilesMolSupplier('/home/jujuman/Research/ANN-Test-Data/GDB-11/gdb11_size02.smi', nameColumn=0)
molecules = Chem.SmilesMolSupplier(smfile, nameColumn=0)
Nmol = 0
NDat = 0

c_sp1 = []
c_sp2 = []
c_sp3 = []
n_sp2 = []
n_sp3 = []

for m in molecules:
    if m is None: continue

    typecheck = False
    for a in m.GetAtoms():
        sym = str(a.GetSymbol())
        count = 0

        for i in At:
            if i is sym:
                count = 1

        if count is 0:
            typecheck = True

    if typecheck is False and random.random() < P:

        m = Chem.AddHs(m) # Add Hydrogens
        AllChem.EmbedMolecule(m) # Embed in 3D Space
        AllChem.UFFOptimizeMolecule(m) # Classical Optimization

        #---------------------------------------------

        print('Molecule ', str(Nmol) ,': ', Chem.MolToSmiles(m))

        print('Number of Atoms: ', m.GetNumAtoms())
        print('Number of Bonds: ', m.GetNumBonds())
        print('Number of Conformers: ', m.GetNumConformers())
        
        if m.GetNumConformers() > 1:
            print('MORE THAN ONE CONFORMER!')
            exit(1)

        xyz = np.ndarray(shape=(m.GetNumAtoms(),3), dtype=np.float32)
        typ = []
        hbr = []

        for i in range (0,m.GetNumAtoms()):
            pos = m.GetConformer().GetAtomPosition(i)
            sym = m.GetAtomWithIdx(i).GetSymbol()
            hyb = m.GetAtomWithIdx(i).GetHybridization()
            #print (' ' + str(sym) + ' ' + str(hyb) + ' ' + "{:.5f}".format(pos.x) + ' ' + "{:.5f}".format(pos.y) + ' ' + "{:.5f}".format(pos.z))
            xyz[i,0] = pos.x
            xyz[i,1] = pos.y
            xyz[i,2] = pos.z
            typ.append(sym)
            hbr.append(str(hyb))

        print (hbr)
        print (typ)
        print (xyz)

        nc.setMolecule(coords=xyz,types=typ)

        O = nc.optimize(conv=0.00001,max_iter=250)
        print(O)

        E1 = nc.energy()
        print('Energy:  ' + str(E1))

        AE1 = np.copy(nc.aenergies(sae=False))
        print('Aenergy: ' + str(AE1))

        for i in range (0,m.GetNumAtoms()):
            if typ[i] == 'C' and hbr[i] == 'SP':
                print ('SP  Carbon Found! : ' + str(AE1[i]))
                c_sp1.append(AE1[i])

            if typ[i] == 'C' and hbr[i] == 'SP2':
                print ('SP2 Carbon Found! : ' + str(AE1[i]))
                c_sp2.append(AE1[i])

            if typ[i] == 'C' and hbr[i] == 'SP3':
                print ('SP3 Carbon Found! : ' + str(AE1[i]))
                c_sp3.append(AE1[i])

            if typ[i] == 'N' and hbr[i] == 'SP2':
                print ('SP3 Carbon Found! : ' + str(AE1[i]))
                n_sp2.append(AE1[i])

            if typ[i] == 'N' and hbr[i] == 'SP3':
                print ('SP3 Carbon Found! : ' + str(AE1[i]))
                n_sp3.append(AE1[i])

        Nmol += 1 #increment counter


print ("|---------Computations Complete----------|")
print (c_sp1)
print (c_sp2)
print (c_sp3)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 12}

plt.rc('font', **font)

fig, axes = plt.subplots(nrows=1, ncols=1)

axes.set_title("Comparison of atomic energy from ANI-c08e(384) for SP, SP2, and SP3 carbons in GDB-04")
axes.set_ylabel('Energy count')
axes.set_xlabel('Energy (Ha)')

axes.hist(c_sp1, 500, color='red',normed=True, label='C SP',linewidth=2,alpha=0.6)
axes.hist(c_sp2, 50, color='blue',normed=True, label='C SP2',linewidth=2,alpha=0.6)
axes.hist(c_sp3, 50, color='orange',normed=True, label='C SP3',linewidth=2,alpha=0.6)
axes.hist(n_sp2, 50, color='green',normed=True, label='N SP2',linewidth=2,alpha=0.6)
axes.hist(n_sp3, 50, color='black',normed=True, label='N SP3',linewidth=2,alpha=0.6)

plt.legend(bbox_to_anchor=(0.6, 0.98), loc=2, borderaxespad=0., fontsize=14)

# -----
# PLOT
# -----
plt.show()