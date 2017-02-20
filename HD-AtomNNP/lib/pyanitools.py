import numpy as np
import pandas as pd
from pyNeuroChem import cachegenerator as cg

# Extract and stores a molecules conformer data and eneriges
class ANIDataExtractor:
    def __init__(self, store, subdir, dfname):
        store_loc = subdir + "/" + dfname
        self.extract_species_attrb(store,store_loc)
        self.extract_from_df(store.select(store_loc))

    def extract_from_df(self, df_data):
        self.coords = np.asarray(df_data['coordinates'],dtype=np.float32)
        self.coords = self.coords.reshape(self.coords.shape[0],self.Na,3)
        self.energy = np.asarray(df_data['energy']).flatten()

    def extract_species_attrb(self,store,store_loc):
        self.species = store.get_storer(store_loc).attrs.species
        self.Na      = self.species.shape[0]

# Create a list of ANIDataExtractors
def aniloadhdf5(file):
    # opening file
    store = pd.HDFStore(file, complib='blosc', complevel=8)
    # HDFStore iterates over the names of its contents:
    data = []
    for x in store.get_node(""):
        print("Name:", x._v_name)
        for i in x._v_children:
            data.append(ANIDataExtractor(store, x._v_name, i))

    store.close()
    return data

# Create a list of ANIDataExtractors
def anigeneratecache(file):
    # opening file
    store = pd.HDFStore(file, complib='blosc', complevel=8)
    # HDFStore iterates over the names of its contents:
    cache = cg('_train', '/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk/sae_6-31gd.dat')

    for x in store.get_node(""):
        print("Name:", x._v_name)
        for i in x._v_children:
            data = ANIDataExtractor(store, x._v_name, i)
            cache.insertdata(data.coords,data.energy,list(data.species))

    cache.makemetadata()
    store.close()


path = "/home/jujuman/Research/test_data2.h5"

anigeneratecache(path)