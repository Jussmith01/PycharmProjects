__author__ = 'jujuman'

# Import pyNeuroChem
import hdnntools as gt
import numpy as np
import time as tm
from scipy import stats as st
from os import listdir

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def convertatomicnumber(X):
    if X == 'H':
        return 1
    elif X == 'C':
        return 6
    elif X == 'N':
        return 7
    elif X == 'O':
        return 8

def compute_sae(spec):
    sae = 0.0

    for s in spec:
        if s == 'H':
            sae += -0.500607632585
        elif s == 'C':
            sae += -37.8302333826
        elif s == 'N':
            sae += -54.5680045287
        elif s == 'O':
            sae += -75.0362229210
        else: 
            print('Error, unknown type: ', s)
            exit(1)
    return sae

def generate_angular_data(xyz,spc,Na,fraction=1.0):

    Nm = xyz.shape[0]
    Ns = int(Nm * fraction)

    Nd = int(((Na-1)*((Na-1)-1))/2)
    data = np.empty([Ns*Na,Nd,3],np.float)
    spec = np.empty([Ns*Na,Nd,3],np.int)

    #print(data.shape)
    randindices = np.random.permutation(Nm)[:int(Ns)]

    for idx,mi in enumerate(range(0,int(Ns))):
        m = xyz[randindices[mi]]
        #print("Data(",idx,",",Na,",",mi,"): ")
        for i in range(m.shape[0]):
            #print("    -(",i,")")
            count = 0
            for j in range(m.shape[0]):
                if i != j:
                    v1 = m[i] - m[j]
                    Dij = np.linalg.norm(v1)
                    for k in range(j+1,m.shape[0]):
                        if i != k:
                            v2 = m[i] - m[k]
                            Dik = np.linalg.norm(v2)

                            #if Dik < 3.1 and Dij < 3.1:
                            Ajik = angle_between(v1, v2)
                            data[i+idx*Na,count] = [Dij,Dik,Ajik]
                            spec[i+idx*Na,count] = [convertatomicnumber(spc[i])
                                                    ,convertatomicnumber(spc[j])
                                                    ,convertatomicnumber(spc[k])]
                            #print( '  -DAT: ', data[idx, count], ' SPC: ',spec[idx, count])
                            count += 1

    #print (data.shape)
    return data,spec

#def get_member_count(array):
#    np.array_equal(t, np.array([1, 6, 1]))

#dtdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_begdb/begdb-h2oclusters/h2o_cluster/inputs/data/'

pref = ['01',
        '02',
        '03',
        '04',
        '05',
        '06',
        '07',
        '08',]
#pref = ['01']

for p in pref:
    #dtdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_' + p + '/testdata/'
    dtdir = '/home/jujuman/Python/PycharmProjects/HD-AtomNNP/Data-ANI-1-Figures/data/minimized/GDB-' + p + '/'

    files = listdir(dtdir)
    _timeloop = tm.time()

    data = np.empty((0,3),dtype=np.float)
    spec = np.empty((0,3),dtype=np.float)
    enrg = np.empty((0),dtype=np.float)

    for en,i in enumerate(files):
        #if 'gdb11' in i and '_test.dat' in i:
        if 'gdb11' in i:
            # Read NC data
            print("Reading file (", en, " of ", len(files) ,"): ", i)
            xyz,spc,t_enrg = gt.readncdat(dtdir + i)

            sae = compute_sae(spc)
            t_enrg = t_enrg - sae

            t_data,t_spec = generate_angular_data(xyz,spc,len(spc),1.0)

            t_data = t_data.reshape((t_data.shape[0]*t_data.shape[1],3))
            t_spec = t_spec.reshape((t_spec.shape[0]*t_spec.shape[1],3))

            #print('data shape: ', t_data.shape)
            #print('spec shape: ', t_spec.shape)

            #for l in t_spec:
            #    print(l)

            if t_data.shape[0] > 0:
                data = np.vstack([data,t_data])

            if t_spec.shape[0] > 0:
                spec = np.vstack([spec,t_spec])

            if t_enrg.shape[0] > 0:
                enrg = np.concatenate([enrg,t_enrg])

    _timeloop2 = (tm.time() - _timeloop)
    print('Computation complete. Time: ' + "{:.4f}".format(_timeloop2)  + 'ms')

    print('Storing data...')
    np.savez('data/minimized/angular/GDB-' + p + '_data.npz',data)
    np.savez('data/minimized/angular/GDB-' + p + '_spec.npz',spec)
    np.savez('data/minimized/energy/GDB-' + p + '_enrg.npz' ,enrg)

'''
n_atoms = spec.shape[0]
print (n_atoms)
fraction = 0.1
randindices = np.random.permutation(n_atoms)[:int(n_atoms*fraction)]

print('randint: ',randindices)

plot_data = []
for i in randindices:
    a = data[i]
    t = spec[i]
    if nparray_compare(t,8,1,1):
        plot_data.append(a)

print ('Plot data: ',plot_data)

pld = np.vstack(plot_data)

_timeloop2 = (tm.time() - _timeloop)
print('Computation complete. Time: ' + "{:.4f}".format(_timeloop2)  + 'ms')
'''

# Creating subplots and axes dynamically
#fig, axes = plt.subplots()

# Labelling xaxes and title
#axes.set_title("pH-"+ str("{0:.2f}".format(pH)))
#axes.set_xlabel('$d H108@N\epsilon2-Y115@O$')
#axes.set_yticks(range(3))
#axes.set_xticks(range(6))

#plt.subplots_adjust(wspace=0, hspace=0)

# Deleting labels for other sharedy axes
#plt.setp([a.get_yticklabels() for a in axes[1:]], visible=False)

# Setting the number of yticklabels on first strip of plot
#string = (np.array(['$HIP-N\delta1$ $N\epsilon2$', '$HID-N\delta1$ ', '$HIE-N\epsilon2$']))
# setting how many yticks a plot can have for first f
#axes[0].set_yticks(range(3))
# Passing the string as label
#axes[0].set_yticklabels(string, rotation=45)

# Formatted pH value
#pH_f = np.float("{0:.2f}".format(pH))
#import matplotlib

# pH - 4.00
#x = 0.5 * (pld[:,0] + pld[:,1])
#y = pld[:,2]

#plt.scatter(x,y)
#plt.show()

#xmin, xmax = 0.0, 3.1
#ymin, ymax = 0.0, 3.14159

#xx, yy, f = kde_estimate(x, xmin, xmax, y, ymin, ymax)

#axes.set_xlim(xmin, xmax)
#axes.set_ylim(ymin, ymax)
# Contourf plot
#cfset = axes.contourf(xx, yy, f, 10, cmap='Reds')
#plt.show()
