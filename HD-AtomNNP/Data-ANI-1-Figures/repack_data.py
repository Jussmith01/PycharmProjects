import numpy as np
import pandas as pd

path = "/home/jujuman/Toshiba460GB/Research/ANI-1-DATA-PAPER-FILES/data/data-ani-1.h5"

# opening file
store = pd.HDFStore(path, complib='blosc',complevel=0)
print(store)



# HDFStore iterates over the names of its contents:
'''
for x in store.get_node(""):
    print("Name:", x._v_name)

    xyz = store.select(x._v_name + "/" + x.xyz._v_name, 'molecule == 1 & conformer == 1')
    energy = store.select(x._v_name + "/" + x.energy._v_name, 'molecule == 1  & conformer == 1')
    species = store.select(x._v_name + "/" + x.species._v_name, 'molecule == 1')

    print (np.asarray(xyz))
    print (np.asarray(energy).flatten())
    print (np.asarray(species).flatten())

    for i in x._v_children:
        print(i)

        # getting an item from the hdf5 file:
        z = store.select(x._v_name + "/" + i,'molecule > 0')

        #print(np.asarray(z))

        # print converting to numpy array because numpy has pretty printing
        maxlen = min(len(z), 10)
        print("Shape:", z.shape)
        print("Value:\n", z[:maxlen])
'''
store.close()