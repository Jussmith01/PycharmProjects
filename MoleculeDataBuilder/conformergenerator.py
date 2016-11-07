from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import random
import numpy as np
import re

def generateconformations(m, n):
    m = Chem.AddHs(m)

    tors = m.GetSubstructMatches(Chem.MolFromSmarts('O=C[NH](*)'),)
    #print(tors)

    print ( "number of atoms:" + str(m.GetNumAtoms()) )
    ids=AllChem.EmbedMultipleConfs(m, numConfs=n, enforceChirality=True)

    print ("Initial optimization")
    AllChem.UFFOptimizeMoleculeConfs(m, maxIters=500)

    print ("Applying Restraint")
    for id in ids:
        #fp = AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94')
        #ff = AllChem.MMFFGetMoleculeForceField(m, fp, confId = id)

        #for t in tors:
        #    if m.GetAtomWithIdx(t[3]).GetSymbol() == 'H':
         #       before = AllChem.GetDihedralDeg(m.GetConformer(id),
         #           t[0], t[1], t[2], t[3])
         #       ff.UFFAddTorsionConstraint(t[0], t[1], t[2], t[3],
         #           True, 180-before, 180-before, 9999)

        if AllChem.UFFOptimizeMolecule(m, maxIters=500) != 0:
            if AllChem.UFFOptimizeMolecule(m, maxIters=500) != 0:
                print ("Optimization Failed!")
            else:
                 print ("Optimization Sucsess (2)!")
        else:
            print ("Optimization Sucsess (1)!")
        #energy_value = ff.CalcEnergy()
        #after = AllChem.GetDihedralDeg(m.GetConformer(id),
        #    t[0], t[1], t[2], t[3])
        #m.SetProp('DIHEDRAL', '{0:.2f},{1:.2f}'.format \
        #    (before, after))
        #m.SetProp('ENERGY', '{0:.2f}'.format(energy_value))

    # EmbedMultipleConfs returns a Boost-wrapped type which
    # cannot be pickled. Convert it to a Python list, which can.
    return m, list(ids)

def molconformergenerator (m,fpf,wdir,Nconf,idx):
    TSS = 1
    LOT='WB97X/6-31g*' # High level of theory
    rdm='uniform' #Random dist
    type='nmrandom'
    Temp='400.0'
    SCF='Tight'
    MEM='2048'

    dpf = fpf + '-'

    #------- End Parameters ---------

    if m is None:
        print("Error: cannot open file!")
        exit(1)

    m,ids = generateconformations(m,Nconf)

    rmslist = []
    AllChem.AlignMolConformers(m, RMSlist=rmslist)
    for i in rmslist:
        print(i)

    print( m.GetNumConformers() )

    confcrd = open(wdir+ dpf + str(idx) + '.xyz' , 'w')

    for i in range(0,m.GetNumConformers()):
        print ("-------Conformer " + str(i) + "------")
        c = m.GetConformer(i)
        tmp = AllChem.Compute2DCoords(m)
        Draw.MolToFile(m,wdir + 'c_line_' + str(i) + '.png',size=(200,200))

        #---------- Write Input Variables ------------
        dfname=dpf + str(idx+i) + '_train.dat'
        vdfname=dpf + str(idx+i) + '_valid.dat'
        edfname=dpf + str(idx+i) + '_test.dat'

        filename = wdir + dpf + str(idx+i) + '.ipt'
        f = open(filename , 'w')
        if not f:
            print('Cannot open file: ' + filename)

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
            f.write('#Conformer RMS from 0: 0.0 Mole: ' + fpf)
        else:
            f.write('#Conformer RMS from 0: ' + str( rmslist[i-1]) + ' Mole: ' + fpf )

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

def gendipeptidelist(AAlist):
    fasta = []
    nlist = []
    for i in AAlist:
        for j in AAlist:
            fasta.append( ">1AKI:A|PDBID|CHAIN|SEQUENCE\n" + i + j )
            nlist.append("dipeptide-" + i + j)

    return fasta,nlist
#-------- Parameters -----------

name = 'Retinol'
wdir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/' + name + '/TEST_REMAKE/' #working directory
fpf = name #Filename prefix
file = name + ".mol"

m = Chem.MolFromMolFile(wdir + file)

Nconf = 8
counter = 0

molconformergenerator (m,fpf,wdir,Nconf,counter)

