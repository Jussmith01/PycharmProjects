from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem

#-------- Parameters -----------

R = 0.3
At = ['C', 'O', 'N', 'H']

#------- End Parameters ---------

molecules = Chem.SmilesMolSupplier('/home/jujuman/Research/GDB-11/gdb11_size02.smi', nameColumn=0)
for m in molecules:
    if m is None: continue
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m)
    AllChem.UFFOptimizeMolecule(m)

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
        print('#Number of Atoms: ', m.GetNumAtoms())
        print('#Number of Bonds: ', m.GetNumBonds())
        print('#Number of Conformers: ', m.GetNumConformers())

        if m.GetNumConformers() > 1:
            print('MORE THAN ONE CONFORMER!')
            exit(1)

        print('#Smiles: ', Chem.MolToSmiles(m))
        print ('\n')
        print ('$coordinates')
        for i in range (0,m.GetNumAtoms()):
            pos = m.GetConformer().GetAtomPosition(i)
            sym = m.GetAtomWithIdx(i).GetSymbol()
            print (' ', sym, ' ', sym, ' ', pos.x, ' ', pos.y, ' ', pos.z, ' ', R)
        print ('&\n')

        print ('$connectivity')
        for b in m.GetBonds():
           print (' ', b.GetBeginAtomIdx() + 1, ' ', b.GetEndAtomIdx() + 1)
        print ('&\n')

    else:
        print('Not Using Structure with Smiles: ', Chem.MolToSmiles(m))