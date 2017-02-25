import pygau09tools as g09
import hdnntools as hdn
import numpy as np

wkdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/scans_double_bond_migration/'

filef = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/scans_double_bond_migration/IRC_fwd.log'
fileb = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/scans_double_bond_migration/IRC_bck.log'

dataf = g09.read_irc(filef)
datab = g09.read_irc(fileb)

xyz = np.concatenate([np.flipud(datab[1]),dataf[1]])

for i,x in enumerate(xyz):
    hdn.write_rcdb_input(x,dataf[2],i,wkdir,'double_B_mig',50,'wb97x/6-31g*','1000.0',opt='0')