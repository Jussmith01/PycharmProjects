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
smfile = '/home/jujuman/Research/ANN-Test-Data/GDB-11-B3LYP-6-31gd/smiledata/gdb11_size05.smi' # Smiles file
At = ['C', 'O', 'N'] # Hydrogens added after check
P = 1.0

# Set required files for pyNeuroChem
wkdir    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/smallAEV_testing/train_384-256-128-64-1_c08e/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = wkdir + 'sae_6-31gd.dat'
nnfdir   = wkdir + 'networks/'

#------- End Parameters ---------

# Construct pyNeuroChem class
nc = pync.pyNeuroChem(cnstfile, saefile, nnfdir, 1)

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

AV_T1 = np.zeros((1,64),dtype=np.float32)
AV_T2 = np.zeros((1,64),dtype=np.float32)
AV_T3 = np.zeros((1,64),dtype=np.float32)

print (AV_T1)

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

        #O = nc.optimize(conv=0.00001,max_iter=250)
        #print(O)

        E1 = nc.energy()
        print('Energy:  ' + str(E1))

        AE1 = np.copy(nc.aenergies(sae=False))
        print('Aenergy: ' + str(AE1))

        for i in range (0,m.GetNumAtoms()):

            if typ[i] == 'H':
               nei = m.GetAtomWithIdx(i).GetNeighbors()
               for j in nei:
                   bond = m.GetBondBetweenAtoms(i,j.GetIdx())
                   if j.GetSymbol() == 'C':
                       #if str(bond.GetBondType()) == "DOUBLE":
                        if str(j.GetHybridization()) == "SP":
                            AV = nc.activations(atom_idx=i, layer_idx=2)
                            AV_T1 = np.vstack([AV_T1,AV])
                            print ('Carbon = oxygen double bond: ' + str(bond.GetBondType()))

            if typ[i] == 'H':
               nei = m.GetAtomWithIdx(i).GetNeighbors()
               for j in nei:
                   bond = m.GetBondBetweenAtoms(i,j.GetIdx())
                   if j.GetSymbol() == 'C':
                       if str(j.GetHybridization()) == "SP2":
                            AV = nc.activations(atom_idx=i, layer_idx=2)
                            AV_T2 = np.vstack([AV_T2,AV])
                            print ('Carbon = carbon double bond: ' + str(bond.GetBondType()))

            if typ[i] == 'C':
               nei = m.GetAtomWithIdx(i).GetNeighbors()
               for j in nei:
                   bond = m.GetBondBetweenAtoms(i,j.GetIdx())
                   if j.GetSymbol() == 'C':
                       if str(j.GetHybridization()) == "SP3":
                            AV = nc.activations(atom_idx=i, layer_idx=2)
                            AV_T3 = np.vstack([AV_T3,AV])
                            print ('Carbon = carbon double bond: ' + str(bond.GetBondType()))

            if typ[i] == 'C' and hbr[i] == 'SP':
                print ('SP  Carbon Found! : ' + str(AE1[i]))
                c_sp1.append(AE1[i])

                #AV = nc.activations(atom_idx=i, layer_idx=2)
                #for j in range(0,AV_T1.shape[0]):
                    #if AV[j] < AV_T1[j]:
                        #AV_T1[j] = AV[j]

            if typ[i] == 'C' and hbr[i] == 'SP2':
                print ('SP2 Carbon Found! : ' + str(AE1[i]))
                c_sp2.append(AE1[i])

                #AV = nc.activations(atom_idx=i, layer_idx=2)
                #for j in range(0,AV_T2.shape[0]):
                #    if AV[j] < AV_T2[j]:
                #        AV_T2[j] = AV[j]

            if typ[i] == 'C' and hbr[i] == 'SP3':
                print ('SP3 Carbon Found! : ' + str(AE1[i]))
                c_sp3.append(AE1[i])

                #AV = nc.activations(atom_idx=i, layer_idx=2)
                #for j in range(0,AV_T3.shape[0]):
                #    if AV[j] < AV_T3[j]:
                #        AV_T3[j] = AV[j]

            if typ[i] == 'N' and hbr[i] == 'SP2':
                print ('SP3 Nitrogen Found! : ' + str(AE1[i]))
                n_sp2.append(AE1[i])

            if typ[i] == 'N' and hbr[i] == 'SP3':
                print ('SP3 Nitrogen Found! : ' + str(AE1[i]))
                n_sp3.append(AE1[i])

        Nmol += 1 #increment counter


print ("|---------Computations Complete----------|")

AV_T1 = AV_T1[1:]
AV_T2 = AV_T2[1:]
AV_T3 = AV_T3[1:]

AV_T1_std = np.zeros(64,dtype=np.float32)
AV_T2_std = np.zeros(64,dtype=np.float32)
AV_T3_std = np.zeros(64,dtype=np.float32)

for i in range(0,AV_T1_std.shape[0]):
    AV_T1_std[i] = AV_T1[:,i].std()
    AV_T2_std[i] = AV_T2[:,i].std()
    AV_T3_std[i] = AV_T3[:,i].std()

IDX = np.arange(0,AV_T1.shape[1],1,dtype=float) + 1
plt.plot (IDX,AV_T1_std,marker='o', color='red',  label='H-C(SP)',  linewidth=2)
plt.plot (IDX,AV_T2_std,':',marker='o', color='blue',  label='H-C(SP2)',  linewidth=2)
plt.plot (IDX,AV_T3_std,'--',marker='o', color='green',  label='H-C(SP3)',  linewidth=2)

plt.title("Common bonding features (std. dev.) of last hidden layer for hydrogen atoms in GDB-05")

plt.ylabel('Activation')
plt.xlabel('Element')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=16)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

plt.show()

AV_T1_mean = np.zeros(64,dtype=np.float32)
AV_T2_mean = np.zeros(64,dtype=np.float32)
AV_T3_mean = np.zeros(64,dtype=np.float32)

for i in range(0,AV_T1_std.shape[0]):
    AV_T1_mean[i] = AV_T1[:,i].mean()
    AV_T2_mean[i] = AV_T2[:,i].mean()
    AV_T3_mean[i] = AV_T3[:,i].mean()

IDX = np.arange(0,AV_T1.shape[1],1,dtype=float) + 1
plt.plot (IDX,AV_T1_mean,marker='o', color='red',  label='H-C(SP)',  linewidth=2)
plt.plot (IDX,AV_T2_mean,':',marker='o', color='blue',  label='H-C(SP2)',  linewidth=2)
plt.plot (IDX,AV_T3_mean,'--',marker='o', color='green',  label='H-C(SP3)',  linewidth=2)

plt.title("Common bonding features (mean) of last hidden layer for hydrogen atoms in GDB-05")

plt.ylabel('Activation')
plt.xlabel('Element')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=16)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

plt.show()

AV_T1_min = np.zeros(64,dtype=np.float32)
AV_T2_min = np.zeros(64,dtype=np.float32)
AV_T3_min = np.zeros(64,dtype=np.float32)

for i in range(0,AV_T1_std.shape[0]):
    AV_T1_min[i] = AV_T1[:,i].min()
    AV_T2_min[i] = AV_T2[:,i].min()
    AV_T3_min[i] = AV_T3[:,i].min()

IDX = np.arange(0,AV_T1.shape[1],1,dtype=float) + 1
plt.plot (IDX,AV_T1_min,marker='o', color='red',  label='H-C(SP)',  linewidth=2)
plt.plot (IDX,AV_T2_min,':',marker='o', color='blue',  label='H-C(SP2)',  linewidth=2)
plt.plot (IDX,AV_T3_min,'--',marker='o', color='green',  label='H-C(SP3)',  linewidth=2)

plt.title("Common bonding features (min) of last hidden layer for hydrogen atoms in GDB-05")

plt.ylabel('Activation')
plt.xlabel('Element')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=16)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

plt.show()

print (c_sp1)
print (c_sp2)
print (c_sp3)

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 12}

plt.rc('font', **font)

fig, axes = plt.subplots(nrows=1, ncols=1)

axes.set_title("Comparison of atomic energy from ANI-c08e(384) for SP, SP2, and SP3 carbons in GDB-05")
axes.set_ylabel('Energy count')
axes.set_xlabel('Energy (Ha)')

axes.hist(c_sp1, 500, color='red'  ,normed=True, label='C SP',linewidth=2,alpha=0.6)
axes.hist(c_sp2, 50, color='blue'  ,normed=True, label='C SP2',linewidth=2,alpha=0.6)
axes.hist(c_sp3, 50, color='orange',normed=True, label='C SP3',linewidth=2,alpha=0.6)
axes.hist(n_sp2, 50, color='green' ,normed=True, label='N SP2',linewidth=2,alpha=0.6)
axes.hist(n_sp3, 50, color='black' ,normed=True, label='N SP3',linewidth=2,alpha=0.6)

plt.legend(bbox_to_anchor=(0.6, 0.98), loc=2, borderaxespad=0., fontsize=14)

# -----
# PLOT
# -----
plt.show()