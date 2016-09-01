import sys
sys.path.append("/usr/local/lib/python2.7/site-packages/")
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Draw

m = Chem.MolFromSmiles('NCCc1c[nH]cn1')
mp_A = Chem.MolFromSmiles('[NH3+]CCc1c[nH]cn1')
mp_B = Chem.MolFromSmiles('NCCc1c[nH]c[nH+]1')
mpp = Chem.MolFromSmiles('[NH3+]CCc1c[nH]c[nH+]1')
ref_A = Chem.MolFromSmiles('CCN')
ref_B = Chem.MolFromSmiles('c1c[nH]cn1')
ref_A_fps = AllChem.GetMorganFingerprint(ref_A,3,fromAtoms=None)
ref_B_fps = AllChem.GetMorganFingerprint(ref_B,3,fromAtoms=None)

tmp=AllChem.Compute2DCoords(m)

Draw.MolToFile(tmp,'~/cdk2_mol1.o.png')

def find_ref(reactant,product,rxn):

    charged_atoms = Chem.MolFromSmarts('[+,-]')
    # find which atom is deprotonated
    charged_atoms_in_reactant = reactant.GetSubstructMatches(charged_atoms)
    charged_atoms_in_product = product.GetSubstructMatches(charged_atoms)
    deprotonated_atom = [x[0] for x in charged_atoms_in_reactant if x not in charged_atoms_in_product]

    # fingerprint of product around deprotonated atom
    product = AllChem.GetMorganFingerprint(product,3,fromAtoms=deprotonated_atom)

    similarity_to_ref_A = DataStructs.DiceSimilarity(product,ref_A_fps)
    similarity_to_ref_B = DataStructs.DiceSimilarity(product,ref_B_fps)

    print (deprotonated_atom, similarity_to_ref_A, similarity_to_ref_B)

    if similarity_to_ref_A > similarity_to_ref_B:
        print ("the best reference for "+rxn+" is ref_A")
    else:
        print ("the best reference for "+rxn+" is ref_B")

    return;

find_ref(mp_A,m,"m+_A -> m + H+")
find_ref(mp_B,m,"m+_B -> m + H+")
find_ref(mpp,mp_A,"m++ -> mp_A + H+")
find_ref(mpp,mp_B,"m++ -> mp_B + H+")