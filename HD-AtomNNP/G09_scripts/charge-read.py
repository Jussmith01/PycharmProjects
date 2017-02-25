import hdnntools as hd
import pygau09tools as gau
import numpy as np
import pyNeuroChem as pync

dir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/TestCases/AEChargeTesting/xyz/'

num = [0, 50, 142, 3284, 6453, 7123, 8239, 10291, 11203, 12353]

cf = open(dir+'charge_data_C.dat','w')
hf = open(dir+'charge_data_H.dat','w')
of = open(dir+'charge_data_O.dat','w')
nf = open(dir+'charge_data_N.dat','w')

for i in num:
    file = 'mol-gdb08_opt-'+str(i)+'.xyz.out'
    fxyz = 'mol-gdb08_opt-'+str(i)+'.xyz'

    print('Working on file: ', fxyz)

    c_mul, s_mul = gau.read_charge(dir+file, 'mulliken')
    c_esp, s_esp = gau.read_charge(dir+file, 'esp')

    xyz,typ,Na = hd.readxyz2(dir+fxyz)

    # Set required files for pyNeuroChem
    anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk'
    cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
    saefile  = anipath + '/sae_6-31gd.dat'
    nnfdir   = anipath + '/networks/'

    # Construct pyNeuroChem class
    nc = pync.molecule(cnstfile, saefile, nnfdir, 0)

    # Set the conformers in NeuroChem
    nc.setMolecule(coords=xyz[0],types=list(typ))

    # ANI - Energy
    E = nc.energy()

    # Atomic energy return
    AE = np.copy(nc.aenergies(sae=False))

    for m,e,a,s in zip(c_mul, c_esp, list(AE), typ):
        if s == 'C':
            cf.write("{:.6f}".format(a) + ' ' + m + ' ' + e + '\n')
        elif s == 'H':
            hf.write("{:.6f}".format(a) + ' ' + m + ' ' + e + '\n')
        elif s == 'O':
            of.write("{:.6f}".format(a) + ' ' + m + ' ' + e + '\n')
        elif s == 'N':
            nf.write("{:.6f}".format(a) + ' ' + m + ' ' + e + '\n')
