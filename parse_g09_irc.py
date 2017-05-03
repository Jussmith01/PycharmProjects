import numpy as np
import pygau09tools as pyg
import hdnntools as hdt

n = '10'

fpf = 'benz_rxn_' + str(n)
wkdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/benz_dbm_rxns_full/irc_dbm_benz_'+str(n)+'/'
#dtdir = '/home/jujuman/Dropbox/IRC_DBondMig/Benzene_rxn/rxn_ben'+str(n)+'/'
dtdir = '/home/jujuman/Dropbox/IRC_DBondMig/Benzene_rxn2/rxn_ben'+str(n)+'/IRC.log'

#en, cd, ty  = pyg.get_irc_data(dtdir+'IRC_fwd.log',dtdir+'IRC_bck.log',dtdir+'saddle_ts.log')
en, cd, ty, Rc = pyg.read_irc(dtdir)

np.savez(wkdir+'reaction_coordinate.npz',x=Rc)

Na = len(ty)
Nc = en.shape[0]

print (cd.shape,' ',en.shape)

#for i,x in enumerate(cd):
    #hdt.write_rcdb_input(x,ty,i,wkdir,fpf,5,'wb97x/6-31g*','600.0',opt='0')

hdt.writexyzfile(wkdir+'irc.xyz',cd,ty)

f = open(wkdir+'irc.dat','w')
f.write("comment\n")
f.write(str(Nc)+'\n')
f.write(str(Na)+',')
for j in ty: f.write(j+',')
f.write('\n')
mol = 0
for l,i in enumerate(cd):

    for j in i:
        for k in j:
            f.write(str(k)+',')
        f.write('')
    f.write(str(en[l]) + ',' + '\n')
    mol += 1
f.close()