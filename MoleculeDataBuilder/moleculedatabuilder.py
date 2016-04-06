from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import re

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

R = 0.1
fpf = 'gdb11_s06' #Filename prefix
wdir = '/home/jujuman/Research/ANN-Test-Data/GDB-11-B3LYP-6-31gd/dnntsgdb11_06/' #working directory
smfile = '/home/jujuman/Research/ANN-Test-Data/GDB-11-B3LYP-6-31gd/smiledata/gdb11_size06.smi' # Smiles file
At = ['C', 'O', 'N'] # Hydrogens added after check

TSS='1000' # Training Set Size
VSS='250' # Validation set size
LOT='UB3LYP/6-31g*' # High level of theory
rdm='uniform' #Random dist
type='moldyn'
stsize='2000'
Ntraj='50'
MaxKE='200000'
SCF='Tight'

#------- End Parameters ---------

#fix the file
formatsmilesfile(smfile)

#molecules = Chem.SmilesMolSupplier('/home/jujuman/Research/ANN-Test-Data/GDB-11/gdb11_size02.smi', nameColumn=0)
molecules = Chem.SmilesMolSupplier(smfile, nameColumn=0)
Nmol = 0

#mdcrd = open(wdir + 'molecules.xyz' , 'w')

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

    if typecheck is False:

        f = open(wdir + fpf + '-' + str(Nmol) + '.ipt' , 'w')

        #---------- Write Input Variables ------------
        dfname=fpf + '-' + str(Nmol) + '_train.dat'
        vdfname=fpf + '-' + str(Nmol) + '_valid.dat'

        f.write ('TSS=' + TSS + ' \n')
        f.write ('VSS=' + VSS + ' \n')
        f.write ('LOT=' + LOT + ' \n')
        f.write ('rdm=' + rdm + '\n')
        f.write ('type=' + type + '\n')
        f.write ('stsize=' + stsize + '\n')
        f.write ('Ntraj=' + Ntraj + '\n')
        f.write ('MaxKE=' + MaxKE + '\n')
        f.write ('SCF=' + SCF + '\n')
        f.write ('dfname=' + dfname + ' \n')
        f.write ('vdfname=' + vdfname + ' \n')
        f.write ('optimize=1 \n')

        #---------------------------------------------

        m = Chem.AddHs(m) # Add Hydrogens
        AllChem.EmbedMolecule(m) # Embed in 3D Space
        AllChem.UFFOptimizeMolecule(m) # Classical Optimization

        print('Molecule ', str(Nmol) ,': ', Chem.MolToSmiles(m))

        #print('#Number of Atoms: ', m.GetNumAtoms())
        #print('#Number of Bonds: ', m.GetNumBonds())
        #print('#Number of Conformers: ', m.GetNumConformers())
        
        if m.GetNumConformers() > 1:
            print('MORE THAN ONE CONFORMER!')
            exit(1)

        f.write ('\n')
        f.write('#Smiles: ' + Chem.MolToSmiles(m))
        f.write ('\n\n')
        f.write ('$coordinates\n')

        mdcrd = open(wdir + 'molecule-' + str(Nmol) + '.xyz' , 'w')
        mdcrd.write('\n' + str(m.GetNumAtoms()) + '\n')

        for i in range (0,m.GetNumAtoms()):
            pos = m.GetConformer().GetAtomPosition(i)
            sym = m.GetAtomWithIdx(i).GetSymbol()
            f.write (' ' + str(sym) + ' ' + str(sym) + ' ' + "{:.5f}".format(pos.x) + ' ' + "{:.5f}".format(pos.y) + ' ' + "{:.5f}".format(pos.z) + ' ' + "{:.3f}".format(R) + '\n')
            mdcrd.write (str(sym) + ' ' + "{:.5f}".format(pos.x) + ' ' + "{:.5f}".format(pos.y) + ' ' + "{:.5f}".format(pos.z) + '\n')

        mdcrd.close()
        f.write ('&\n\n')

        f.write ('$connectivity\n')
        for b in m.GetBonds():
           f.write (' ' + str(b.GetBeginAtomIdx() + 1) + ' ' + str(b.GetEndAtomIdx() + 1) + '\n')
        f.write ('&\n\n')

        f.write ('$randrange\n')
        f.write (' NONE\n')
        f.write ('&\n\n')

        Nmol += 1 #increment counter
    else:
        print('Not Using Structure with Smiles: ', Chem.MolToSmiles(m))
