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

wkdir = '/home/jujuman/Research/CrossValidation/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile = wkdir + 'sae_6-31gd.dat'

At = ['C', 'O', 'N'] # Hydrogens added after check

#-------------------------------------------
#nnfdir   = wkdir + 'cv_c08e_ntw_' + str(0) + '/networks/'

# Construct pyNeuroChem classes
#nc = pync.molecule(cnstfile, saefile, nnfdir, 0)

nc =  [pync.molecule(cnstfile, saefile, wkdir + 'cv_c08e_ntw_' + str(l) + '/networks/', 0) for l in range(5)]

molecules = Chem.SmilesMolSupplier(smfile, nameColumn=0)

#mols = [molecules[i] for i in range(0,1000)]
#print(mols)
f = open('/home/jujuman/Research/CrossValidation/GDB-09-High-sdev/gdb-09-0.5sdev.dat','w')
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
        m = Chem.AddHs(m) # Add Hydrogens
        AllChem.EmbedMolecule(m) # Embed in 3D Space
        check = AllChem.MMFFOptimizeMolecule(m, maxIters=1000)  # Classical Optimization

        #print('CLASSICAL OPT:', check)

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

        #print(xyz)
        #print(spc)

        #print('Optimize...')
        #start_time = time.time()
        dyn = LBFGS(mol,logfile='logfile.txt')
        dyn.run(fmax=0.001,steps=5000)
        conv = True if dyn.get_number_of_steps() == 10000 else False
        #print('[ANI Total time:', time.time() - start_time, 'seconds]')

        xyz = np.array(mol.get_positions(),dtype=np.float32).reshape(xyz.shape[0],3)
        xyz2 = np.array(mol.get_positions(),dtype=np.float32).reshape(1,xyz.shape[0],3)

        if conv:
            print('Failed to converge!!!')

        #print('    SPECIES: ', spc)
        energies = np.zeros((5),dtype=np.float64)
        energies2 = np.zeros((5),dtype=np.float64)

        N = 0
        for comp in nc:
            #print('Comp:', l)
            #comp = pync.molecule(cnstfile, saefile, wkdir + 'cv_c08e_ntw_' + str(l) + '/networks/', 0)
            comp.setMolecule(coords=xyz, types=list(spc))
            energies[N] = comp.energy()[0]
            #forces = comp.force().flatten()

            #print('FORCE: ', np.abs(forces - forces2).sum()/(3*Na))
            N = N + 1

        if np.std(hdt.hatokcal*energies) > 0.5:
            output = '  ' + str(k) + ' (' + str(Na) + ') : stps=' + str(dyn.get_number_of_steps()) + ' : ' + str(energies) + ' : std(kcal/mol)=' + str(np.std(hdt.hatokcal*energies)) + ' : ' + Chem.MolToSmiles(m)
            f.write(output+'\n')
            print(output)
        #print('              ',energies)


f.close()
        #print('End...')

