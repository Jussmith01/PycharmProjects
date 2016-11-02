from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from PIL import Image
#from cairo import ImageSurface

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

dir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/Retinol/'

molnames = [ ['P-Menthane',        'Retinol.mol'] ]


lst = []
for i in molnames:
    suppl = Chem.SDMolSupplier(dir+'/'+i[1])
    ms = [x for x in suppl if x is not None]
    for m in ms: tmp=AllChem.Compute2DCoords(m)
    lst.append(ms[0])

    #Draw.MolToFile(ms[0],dir + i[0] + '/' + i[0] + '.png')

leg = []
for i in range(0,len(molnames)):
    leg.append( str(i) + ') ' + molnames[i][0])

print(leg)

Draw.MolToFile(lst[0],'images/cdk2_mol1.o.png')

#img=Draw.MolsToGridImage(lst,molsPerRow=3,subImgSize=(200,200),legends=[x for x in leg],)
#img = img.resize((4096, 4096), Image.ANTIALIAS)
#img.save(dir + 'C10H20_isomers.png', 'JPEG', quality=100)