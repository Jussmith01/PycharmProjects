from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.SDMolSupplier('/home/jujuman/Downloads/Structure2D_CID_46195494.sdf')[0]
#mol = Chem.SDMolSupplier('/home/jujuman/Downloads/Structure2D_CID_123133709.sdf')[0]

print('Atoms: ',mol.GetNumAtoms())

print('Adding Hs...')
mol = Chem.AddHs(mol)

print('Atoms: ',mol.GetNumAtoms())

#print('Computing 2D coords...')
#AllChem.Compute2DCoords(mol)

#print('Write SDF...')
#w = Chem.SDWriter('mol1.sdf')
#w.write(mol)

print('Embedding in 3D...')
AllChem.EmbedMolecule(mol)

print('Write SDF...')
w = Chem.SDWriter('mol2.sdf')
w.write(mol)

#print('Optimizing...')
#AllChem.UFFOptimizeMolecule(mol)

#print('Write SDF...')
#w = Chem.SDWriter('mol3.sdf')
#w.write(mol)

print('Finished.')