import numpy as np
import pandas as pd
import os

##Note: this code will make and delete the file store.h5 in the directory it runs in.

alphabet = "abcdefghijklmnopqrstuvwxyz"
n_arrs = 100

features = {"xyz": (3,),
            "species": (1,),
            "en": (1,)}
n_size = 10

print("Crude size estimate: {:5.5g} MB".format(4 * n_arrs * n_size / 2 * 5 / 2 ** 20))


def randomstring():
    return "".join(np.random.choice([char for char in alphabet], size=6))


def random_size_array(n, restshape):
    return np.random.normal(
        size=(np.random.randint(1, n),) + restshape).astype('float32')


def print_mb(fname):
    print("Object {} is {:5.5g} MB".format(fname, os.stat(fname).st_size / 2 ** 20))


# Dictionary of dictionary of dataframes of varying sizes:
vals = {randomstring():
            {j:
                 pd.DataFrame(random_size_array(n_size, features[j])) for j in features}
        for i in range(n_arrs)}

path = 'store.h5'
if os.path.exists(path):
    os.remove(path)
# open an HDF5 for compressed storage.
# Note that if the path exists, it will open whatever is there.
store = pd.HDFStore(path, complib='zlib', complevel=9)

# Add tables to hdf5 while tracking their size
size = 0
for k1, v1 in vals.items():
    for k2, v2 in v1.items():
        print(k1 + '/' + k2)
        store.put(k1 + '/' + k2, v2, complevel=9)
        # Track size
        s = v2.memory_usage().sum()
        size += s

print("Size of values in memory: {:5.5g} MB".format(size / 2 ** 20))
store.close()
print_mb(path)

# Note that the h5 file can be zipped for better compression...

# opening file
store = pd.HDFStore(path, complib='zlib', complevel=9)

# HDFStore iterates over the names of its contents:
for x in list(store)[:5]:
    print("Name:", x)

    # getting an item from the hdf5 file:
    z = store.select(x)

    # print converting to numpy array because numpy has pretty printing
    maxlen = min(len(z), 5)
    print("Value:", z[:maxlen])
store.close()