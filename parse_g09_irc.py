import numpy as np
import pygau09tools as pyg

en1, cd1, ty  = pyg.read_irc('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/scans_double_bond_migration/IRC_fwd.log')
en2, cd2, ty2 = pyg.read_irc('/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_rxns/scans_double_bond_migration/IRC_bck.log')

en = np.concatenate([np.fliplr(en2),en1])
print(en)

cd = cd1 + cd2

Na = len(ty)
Nc = en.shape[0]

print (cd.shape,' ',en.shape)

import matplotlib.pyplot as plt

plt.plot (en, color='blue',  label='SRC+ANI',  linewidth=2)

#plt.plot (x, morse(x,1.0,1.9,1.53), color='grey',  label='SRC+ANI',  linewidth=2)

plt.title("C-C Dissociation")

plt.ylabel('E (kcal/mol)')
plt.xlabel('Distance $\AA$')
plt.legend(bbox_to_anchor=(0.05, 0.95), loc=2, borderaxespad=0.,fontsize=16)

plt.show()

f = open("irc.out",'w')
f.write("comment\n")
f.write(str(Nc-2)+'\n')
f.write(str(Na)+',')
for j in ty: f.write(j+',')
f.write('\n')
mol = 0
for l,i in enumerate(cd):

    for j in i:
        for k in j:
            f.write(str(k)+',')
        f.write('')
    f.write(str(en[l]) + ',')
    mol += 1
f.close()