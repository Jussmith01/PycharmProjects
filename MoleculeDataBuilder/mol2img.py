from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from PIL import Image
'''
dir = '/home/jujuman/Dropbox/Research/ChemSciencePaper/TestCases/C10H20Isomers/'

molnames = [ ['P-Menthane',        'Structure2D_CID_7459.sdf']
             ,['N-Butylcyclohexane','Structure2D_CID_15506.sdf']
             ,['T-Butylcyclohexane','Structure2D_CID_18508.sdf']
             ,['Pentylcyclopentane','Structure2D_CID_19540.sdf']
             ,['Trans-2-Decene',    'Structure2D_CID_5364559.sdf']
             ,['Trans-4-Decene',    'Structure2D_CID_5364458.sdf']
             ,['Trans-3-Decene',    'Structure2D_CID_5362724.sdf']
             ,['Trans-5-Decene',    'Structure2D_CID_637820.sdf']
             ,['Cis-4-Decene',      'Structure2D_CID_5364504.sdf']
             ,['Cis-5-Decene',      'Structure2D_CID_5364449.sdf']
             ,['Cis-2-Decene',      'Structure2D_CID_5364451.sdf']
             ,['Cis-3-Decene',      'Structure2D_CID_5362723.sdf']
             ,['Dec-1-ene',         'Structure2D_CID_13381.sdf']]
'''
dir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/Figures/Results/EtotCorrMulti/Remake/'

molnames = [  ['17',  '[H]OC(C(=NN([H])[H])N([H])C([H])=C([H])[H])(C([H])([H])[H])C([H])([H])[H]']
	     ,['66',  '[H]ON=C1C([H])([H])C(=C([H])[H])C(O[H])(C([H])([H])[H])C1([H])[H]']
	     ,['67',  '[H]C(=O)C1=C([H])C(=O)C(C([H])([H])[H])(C([H])([H])[H])N1[H]']
	     ,['78',  '[H]ON1C(=O)C(C([H])([H])[H])(C([H])([H])C([H])=C([H])[H])C1([H])[H]']
	     ,['94',  '[H]N=C(OC([H])([H])[H])N([H])C([H])(C([H])([H])[H])C([H])([H])N([H])C([H])([H])[H]']
	     ,['112', '[H]N=C1N([H])N([H])C([H])([H])C(O[H])(C([H])([H])C([H])([H])[H])C1([H])[H]']
	     ,['119', '[H]OC1(C([H])=C([H])[H])C(=O)C([H])=C(N([H])[H])C1([H])[H]']
	     ,['125', '[H]N=C([H])N(N=C([H])C([H])([H])[H])C([H])([H])C([H])(O[H])C([H])([H])[H]']
	     ,['131', '[H]N=C([H])N(N([H])C(=O)C([H])([H])[H])C([H])(C([H])([H])[H])C([H])([H])[H]']]

lst = []
for i in molnames:
    mol = Chem.MolFromSmiles(i[1])
    tmp = AllChem.Compute2DCoords(mol)
    lst.append(mol)
    Draw.MolToFile(mol,dir + 'f_line_' + i[0] + '.png',size=(200,200))

leg = []
for i in range(0,len(molnames)):
    leg.append( str(i) + ') ' + molnames[i][0])

print(leg)

img=Draw.MolsToGridImage(lst,molsPerRow=3,subImgSize=(200,200),legends=[x for x in leg])
img = img.resize((2048, 2048), Image.ANTIALIAS)
img.save(dir + 'molecules.png', 'JPEG', quality=1000)
