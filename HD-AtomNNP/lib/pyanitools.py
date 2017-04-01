import h5py


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
            print(k, v)
            g.create_dataset(k, data=v, compression=self.clib, compression_opts=self.clev)
            print("post")

    def cleanup(self):
        """Wrapper to close HDF5 file
        """
        self.store.close()


class anidataloader(object):
    def __init__(self, store_file):
        self.store = h5py.File(store_file)

    def __iter__(self):
        for g in self.store.values():
            print(g)
            yield dict((d, g[d].value) for d in g)

    getnextdata = __iter__

    def get_node_list(self):
        return [g for g in self.store]

    def size(self):
        return len(self.get_node_list())
    
    def cleanup(self):
        self.store.close()
