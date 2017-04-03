import h5py
import numpy as np
import platform

PY_VERSION = int(platform.python_version().split('.')[0]) > 3

class datapacker(object):
    def __init__(self, store_file, mode='w-', complib='gzip', complevel=6):
        """Wrapper to store arrays within HFD5 file
        """
        # opening file
        self.store = h5py.File(store_file, mode=mode)
        self.clib = complib
        self.clev = complevel

    def store_data(self, store_loc, **kwargs):
        """Put arrays to store
        """
        g = self.store.create_group(store_loc)
        for k, v, in kwargs.items():
            #print(type(v[0]))

            if type(v[0]) is np.str_ and PY_VERSION:
                v = [a.encode('utf8') for a in v]

            g.create_dataset(k, data=v, compression=self.clib, compression_opts=self.clev)

    def cleanup(self):
        """Wrapper to close HDF5 file
        """
        self.store.close()


class anidataloader(object):
    def __init__(self, store_file):
        self.store = h5py.File(store_file)

    def __iter__(self):
        for g in self.store.values():
            yield dict((d, g[d].value) for d in g)

    getnextdata = __iter__

    def get_node_list(self):
        data = [g for g in self.store]

        for d in data:
            if type(d[0]) is bytes and PY_VERSION:
                d = d.astype(str)

        return data

    def size(self):
        return len(self.get_node_list())

    def cleanup(self):
        self.store.close()

