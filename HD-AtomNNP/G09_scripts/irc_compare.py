import pygau09tools as g09
import hdnntools as hdn
import numpy as np
import matplotlib.pyplot as plt
import pyNeuroChem as pync

#dir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/scans_double_bond_migration/'
#dir = '/home/jujuman/Dropbox/IRC_DBondMig/rxn1/'
dir = '/home/jujuman/Dropbox/IRC_DBondMig/original_rxn/'

filef = dir + 'IRC_fwd.log'
fileb = dir + 'IRC_bck.log'

dataf = g09.read_irc(filef)
datab = g09.read_irc(fileb)

xyz = np.concatenate([np.flipud(datab[1]),dataf[1]])
eng = hdn.hatokcal * np.concatenate([np.flipud(datab[0]),dataf[0]])

print (xyz.shape)

hdn.writexyzfile(dir + 'scan.xyz',xyz,dataf[2])

# Set required files for pyNeuroChem

# Set required files for pyNeuroChem
anipath  = '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk'
cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

# Construct pyNeuroChem class
nc = pync.conformers(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=dataf[2])

ani1 = hdn.hatokcal * nc.energy()

anipath    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/ANI-c08e-ntwk_retrain/'

cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'

# Construct pyNeuroChem class
nc = pync.conformers(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=dataf[2])

ani2 = hdn.hatokcal * nc.energy()


anipath    = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/ANI-c08e-ntwk_newtrain/'
cnstfile = anipath + '/rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = anipath + '/sae_6-31gd.dat'
nnfdir   = anipath + '/networks/'


# Construct pyNeuroChem class
nc = pync.conformers(cnstfile, saefile, nnfdir, 0)

# Set the conformers in NeuroChem
nc.setConformers(confs=xyz,types=dataf[2])

ani3 = hdn.hatokcal * nc.energy()

rmse1 = hdn.calculaterootmeansqrerror(ani1,eng)
rmse2 = hdn.calculaterootmeansqrerror(ani2,eng)
rmse3 = hdn.calculaterootmeansqrerror(ani3,eng)

xv = xyz[:,3,:]-xyz[0,3,:]
x = np.linalg.norm(xv, axis=1)

eng = eng - eng.min()
ani1 = ani1 - ani1.min()
ani2 = ani2 - ani2.min()
ani3 = ani3 - ani3.min()

plt.plot(x,eng, 'r-',  color='black', label='DFT',linewidth=3)
plt.plot(x,ani1,'r-.', color='blue', label= 'ANI Original Er: ' + "{:.2f}".format(rmse1),linewidth=3)
#plt.plot(x,ani2,'r--', color='red', label=  'ANI RT Er: ' + "{:.2f}".format(rmse2),linewidth=3)
plt.plot(x,ani3,'r--', color='red', label='ANI New Er: ' + "{:.2f}".format(rmse3),linewidth=3)

#--------------Parameters------------------
wkdir = '/home/jujuman/Research/CrossValidation/'
cnstfile = wkdir + 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile = wkdir + 'sae_6-31gd.dat'
#-------------------------------------------

# Build networks
nc =  [pync.conformers(cnstfile, saefile, wkdir + 'cv_c08e_ntw_' + str(l) + '/networks/', 0) for l in range(5)]
for i,nn in enumerate(nc):
    # Set the conformers in NeuroChem
    nn.setConformers(confs=xyz, types=dataf[2])

    Ecv = hdn.hatokcal * nn.energy()
    Ecv = Ecv - Ecv.min()
    if i is 0:
        plt.plot(x,Ecv,'r--', color='orange',label='ANI CV',linewidth=1)
    else:
        plt.plot(x, Ecv, 'r--', color='orange', linewidth=1)

plt.title("Reaction barrier for double bond migration\n(CV = cross validation; Err in $kcal*{mol}^-1$)")
#plt.xlabel('Conformation Pair (Count 49)')
plt.xlabel('Reaction Coordinate ($\AA$)')
plt.ylabel('Calculated $\Delta E$ ($kcal*{mol}^-1$)')
plt.legend(bbox_to_anchor=(0.65, 0.975), loc=2, borderaxespad=0.)

# -----
# PLOT
# -----
plt.show()
