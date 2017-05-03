# Import pyNeuroChem
import pyNeuroChem as pync
import numpy as np
import hdnntools as hdt
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import pygau09tools as pyg
from itertools import chain
import os

def plot_irc_data(axes, file, title, ntwl, cnstfile, saefile, dir, trained):
    Eact, xyz, typ, Rc = pyg.read_irc(file)
    Rc = Rc[:,1]

    # Shift reference to reactant
    #Eact = Eact[::-1]
    Eact = hdt.hatokcal * (Eact - Eact[0])

    # Plot reference results
    axes.plot (Rc,Eact, color='black',  linewidth=3)

    # Plot ANI results
    color = cm.rainbow(np.linspace(0, 1, len(ntwl)))
    terr = np.zeros(len(ntwl))
    derr = np.zeros(len(ntwl))
    berr = np.zeros(len(ntwl))
    for i, (nt,c) in enumerate(zip(ntwl,color)):
        ncr = pync.conformers(dir + cnstfile, dir + saefile, rcdir + nt[0] + 'networks/', 0)

        # Set the conformers in NeuroChem
        ncr.setConformers(confs=xyz, types=list(typ))

        # Compute Energies of Conformations
        E1 = ncr.energy()

        # Shift ANI E to reactant
        E1 = hdt.hatokcal * (E1 - E1[0])

        # Calculate error
        errn = hdt.calculaterootmeansqrerror(E1,Eact)

        terr[i] = errn
        derr[i] = np.abs(np.abs((E1[0] - E1[-1])) - np.abs((Eact[0] - Eact[-1])))
        berr[i] = np.abs(E1.max() - Eact.max())

        # Plot
        axes.plot(Rc,E1, 'r--', color=c, label="["+nt[1]+"]: "+"{:.2f}".format(errn), linewidth=2)
        #axes.plot([Rc['x'][:,1].min(),Rc['x'][:,1].max()],[E1[-1],E1[-1]], 'r--', color=c)
        #axes.plot([Rc['x'][:,1].min(),Rc['x'][:,1].max()],[E1[0],E1[0]], 'r--', color=c)

    axes.set_xlim([Rc.min(), Rc.max()])
    axes.legend(loc="upper left",fontsize=8)
    if trained:
        axes.set_title(title,color='green',fontdict={'weight':'bold'})
    else:
        axes.set_title(title, color='red', fontdict={'weight': 'bold'})
    return terr,derr,berr

main_dir = '/home/jujuman/Dropbox/IRC_DBondMig/Benzene_rxn2/'

# Set required files for pyNeuroChem
rcdir  = '/home/jujuman/Research/ANI-DATASET/RXN1_TNET/training/'
cnstfile = 'rHCNO-4.6A_16-3.1A_a4-8.params'
saefile  = 'sae_6-31gd.dat'

ntwl = [('ANI-c08f-ntwk/', 'N'),
        #('rxn1/ani_benz_rxn_ntwk/', '1'),
        #('rxn2/ani_benz_rxn_ntwk/', '1,2'),
        #('rxn-1-2-5-6/ani_benz_rxn_ntwk/', '1,2,5,6'),
        ('rxn1to6/ani_benz_rxn_ntwk/','1-6'),
        ]

t_list = ['1-1',
          '2-1',
          '2-5',
          '2-6',
          '3-1',
          '4-1',]

# THE PROGRAM!
sub_dir = list(chain.from_iterable([[main_dir+d+'/'+i+'/IRC.log' for i in os.listdir(main_dir+d)] for d in os.listdir(main_dir) if '.docx' not in d]))

#print(lambda x:x.split('/')[-2])
sub_dir.sort(key=lambda x:(int(x.split('/')[-2].split('-')[0]), int(x.split('/')[-2].split('-')[1])))

print(len(sub_dir))
X = int(np.ceil(np.sqrt(len(sub_dir))))
Y = int(np.round(np.sqrt(len(sub_dir))))
f, axarr = plt.subplots(Y, X,sharey=False,sharex=False)

terrs = []
derrs = []
berrs = []
for i,d in enumerate(sub_dir):
    print(d)
    t_te, t_de, t_be = plot_irc_data(axarr[int(np.floor(i/X)), i % X], d, d.split('/')[-2], ntwl, cnstfile, saefile, rcdir, d.split('/')[-2] in t_list)
    terrs.append(t_te)
    derrs.append(t_de)
    berrs.append(t_be)

terrs = np.vstack(terrs)
derrs = np.vstack(derrs)
berrs = np.vstack(berrs)

cts = "Average RMSE:        "
cds = "Average R/P Delta E: "
cbs = "Barrier height:      "
for i, (ct, cd, cb) in enumerate(zip(terrs.T, derrs.T, berrs.T)):
    cts += "[" + ntwl[i][1] + "] " + "{:.2f}".format(np.mean(ct)) + " "
    cds += "[" + ntwl[i][1] + "] " + "{:.2f}".format(np.mean(cd)) + " "
    cbs += "[" + ntwl[i][1] + "] " + "{:.2f}".format(np.mean(cb)) + " "
    print("[" + ntwl[i][1] + "] " + "{:.2f}".format(np.mean(ct)) + " ")
    print("[" + ntwl[i][1] + "] " + "{:.2f}".format(np.mean(cd)) + " ")
    print("[" + ntwl[i][1] + "] " + "{:.2f}".format(np.mean(cb)) + " ")

plt.suptitle(str(len(sub_dir))+" double bond migration reactions (x axis=$R_c$;y-axis=relative E [kcal/mol])\n"+cts+"\n"+cds+"\n"+cbs,fontsize=12)

plt.show()
#for d in
#print(sub_dir)

#en, cd, ty, Rc = pyg.read_irc(dtdir)