# Import pyNeuroChem
from __future__ import print_function

# Neuro Chem
from ase_interface import ANI
import pyNeuroChem as pync

import hdnntools as gt
import numpy as np
import matplotlib.pyplot as plt
import time as tm
from scipy import stats as st
import time

import hdnntools as hdt

from rdkit import Chem
from rdkit.Chem import AllChem

# ASE
import  ase
from ase.io import read, write
from ase.optimize import BFGS, LBFGS
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.units import Bohr
from ase.calculators.calculator import Calculator, all_changes
from ase import Atoms

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

#def make_atoms

#--------------Parameters------------------
smfile = '/home/jujuman/Research/RawGDB11Database/gdb11_size09.smi' # Smiles file

wkdir = '/home/jujuman/Research/CrossValidation/GDB-09-Retrain/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile = wkdir + 'sae_6-31gd.dat'

At = ['C', 'O', 'N'] # Hydrogens added after check

Nnc = 5

#-------------------------------------------
#nnfdir   = wkdir + 'cv_c08e_ntw_' + str(0) + '/networks/'

# Construct pyNeuroChem classes
nc =  [pync.molecule(cnstfile, saefile, wkdir + 'cv_c08e_ntw_' + str(l) + '/networks/', 0) for l in range(Nnc)]

molecules = Chem.SmilesMolSupplier(smfile, nameColumn=0)

total_mol = 0
total_bad = 0

#mols = [molecules[i] for i in range(217855,217865)]
f = open('/home/jujuman/Research/CrossValidation/GDB-09-High-sdev/gdb-09-1.0sdev_fix.dat','w')
for k,m in enumerate(molecules):
    if m is None: continue

    typecount = 0

    #print (Chem.MolToSmiles(m))

    typecheck = False
    for a in m.GetAtoms():
        sym = str(a.GetSymbol())
        count = 0

        for i in At:
            if i is sym:
                count = 1

        if count is 0:
            typecheck = True

    if typecheck is False:
        total_mol = total_mol + 1
        m = Chem.AddHs(m) # Add Hydrogens
        embed = AllChem.EmbedMolecule(m, useRandomCoords=True)
        if embed is 0: # Embed in 3D Space was successful
            check = AllChem.MMFFOptimizeMolecule(m, maxIters=1000)  # Classical Optimization

            xyz = np.zeros((m.GetNumAtoms(),3),dtype=np.float32)
            spc = []

            Na = m.GetNumAtoms()
            for i in range (0,Na):
                pos = m.GetConformer().GetAtomPosition(i)
                sym = m.GetAtomWithIdx(i).GetSymbol()

                spc.append(sym)
                xyz[i, 0] = pos.x
                xyz[i, 1] = pos.y
                xyz[i, 2] = pos.z

            mol = Atoms(symbols=spc, positions=xyz)
            mol.set_calculator(ANI(False))
            mol.calc.setnc(nc[0])

            dyn = LBFGS(mol,logfile='logfile.txt')
            dyn.run(fmax=0.001,steps=10000)
            conv = True if dyn.get_number_of_steps() == 10000 else False

            xyz = np.array(mol.get_positions(),dtype=np.float32).reshape(xyz.shape[0],3)

            #if conv:
            #    print('Failed to converge!!!')

            energies = np.zeros((Nnc),dtype=np.float64)

            N = 0
            for comp in nc:
                comp.setMolecule(coords=xyz, types=list(spc))
                energies[N] = comp.energy()[0]
                N = N + 1

            if np.std(hdt.hatokcal*energies) > 1.0:
                total_bad = total_bad + 1
                perc = int(100.0 * total_bad/float(total_mol))
                output = '  ' + str(k) + ' ' + str(total_bad) + '/' + str(total_mol) + ' ' + str(perc) + '% (' + str(Na) + ') : stps=' + str(dyn.get_number_of_steps()) + ' : ' + str(energies) + ' : std(kcal/mol)=' + str(np.std(hdt.hatokcal*energies)) + ' : ' + Chem.MolToSmiles(m)
                f.write(output+'\n')
                print(output)

print('Total Molecs: ', total_mol)
print('Total Bad:    ', total_bad)
print('Percent Bad:  ', int(100.0 * total_bad/float(total_mol)), '%')
f.close()
#print('End...')

