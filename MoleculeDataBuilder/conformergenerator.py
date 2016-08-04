from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
import random
import numpy as np
import re

def generateconformations(m, n):
    m = Chem.AddHs(m)
    print ("number of atoms:" + str(m.GetNumAtoms()))
    ids=AllChem.EmbedMultipleConfs(m, numConfs=n)

    for id in ids:
        AllChem.UFFOptimizeMolecule(m, confId=id)

    # EmbedMultipleConfs returns a Boost-wrapped type which
    # cannot be pickled. Convert it to a Python list, which can.
    return m, list(ids)

#-------- Parameters -----------

dir = "/home/jujuman/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/PeptideCases/"
file = "pp_case_1.mol"

R = 0.3
fpf = 'peptide_conformers_01' #Filename prefix
wdir = '/home/jujuman/Research/ANN-Test-Data/GDB-11-W98XD-6-31gd/PeptideCases/' #working directory

TSS = 1
LOT='WB97X/6-31g*' # High level of theory
rdm='uniform' #Random dist
type='nmrandom'
Temp='400.0'
SCF='Tight'
MEM='2048'

#------- End Parameters ---------

m = Chem.MolFromMolFile(dir + file)

if m is None:
    print("Error: cannot open file!")
    exit(1)

m,ids = generateconformations(m,8)

rmslist = []
AllChem.AlignMolConformers(m, RMSlist=rmslist)
for i in rmslist:
    print(i)

print( m.GetNumConformers() )

confcrd = open(wdir+'conformers.xyz' , 'w')

for i in range(0,m.GetNumConformers()):
    print ("-------Conformer " + str(i) + "------")
    c = m.GetConformer(i)

    #---------- Write Input Variables ------------
    dfname=fpf + '-' + str(i) + '_train.dat'
    vdfname=fpf + '-' + str(i) + '_valid.dat'
    edfname=fpf + '-' + str(i) + '_test.dat'

    f = open(wdir + fpf + '-' + str(i) + '.ipt' , 'w')

    DOF = (3 * m.GetNumAtoms() - 6)
    f.write ('TSS=' + str(int(TSS * DOF)) + ' \n')
    f.write ('VSS=0 \n')
    f.write ('ESS=0 \n')

    f.write ('LOT=' + LOT + ' \n')
    f.write ('rdm=' + rdm + '\n')
    f.write ('type=' + type + '\n')
    f.write ('Temp=' + Temp + '\n')
    f.write ('mem=' + MEM + '\n')
    f.write ('SCF=' + SCF + '\n')
    f.write ('dfname=' + dfname + ' \n')
    f.write ('vdfname=' + vdfname + ' \n')
    f.write ('edfname=' + edfname + ' \n')
    f.write ('optimize=1 \n')
    f.write ('frequency=1 \n')

    f.write ('\n')

    if i is 0:
        f.write('#Conformer RMS from 0: 0.0')
    else:
        f.write('#Conformer RMS from 0: ' + str( rmslist[i-1]) )

    f.write ('\n\n')
    f.write ('$coordinates\n')

    confcrd.write("\n")
    confcrd.write(str(m.GetNumAtoms()) + "\n")

    for j in range (0,m.GetNumAtoms()):
        pos = c.GetAtomPosition(j)
        typ = m.GetAtomWithIdx(j).GetSymbol()

        f.write (' ' + str(typ) + ' ' + str(typ) + ' ' + "{:.5f}".format(pos.x) + ' ' + "{:.5f}".format(pos.y) + ' ' + "{:.5f}".format(pos.z) + '\n')
        confcrd.write(typ + ' ' + str(pos.x) + ' ' + str(pos.y)+ ' ' + str(pos.z) + '\n')

    f.write ('&\n\n')

    f.write ('$connectivity\n')
    f.write ('   NONE\n')
    f.write ('&\n\n')

    f.write ('$normalmodes\n')
    f.write (' NEED TO COMPUTE\n')
    f.write ('&\n\n')