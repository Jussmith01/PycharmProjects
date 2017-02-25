import numpy as np
import pandas as pd

# Extract and stores a molecules conformer data and eneriges
class anidataextractor:
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

# Load the entire data set
class anidataloader:
    def __init__(self, store_file, complib = 'blosc', complevel = 8):
        # opening file
        self.store = pd.HDFStore(store_file, complib=complib, complevel=complevel)

    def splitload(self, N):
        self.crd_list = []
        self.eng_list = []
        self.spc_list = []
        for x in self.store.get_node(""):
            print("Name:", x._v_name)
            for i in x._v_children:
                ae = anidataextractor(self.store, x._v_name, i)
                self.crd_list.append(np.array_split(ae.coords, N))
                self.eng_list.append(np.array_split(ae.energy, N))
                self.spc_list.append(ae.species)

    def getdata(self,idx=0,dl=[0]):

        if max(dl) >= len(self.crd_list):
            raise (IndexError('Index given is outside of the array range.'))

        return [np.concatenate([self.crd_list[idx][j] for j in dl]),
                np.concatenate([self.eng_list[idx][j] for j in dl]),
                self.spc_list[idx]]

    def totalload(self):
        self.crd_list = []
        self.eng_list = []
        self.spc_list = []
        for x in self.store.get_node(""):
            print("Name:", x._v_name)
            for i in x._v_children:
                ae = anidataextractor(self.store, x._v_name, i)
                self.crd_list.append([ae.coords])
                self.eng_list.append([ae.energy])
                self.spc_list.append(ae.species)

    def size(self):
        return len(self.spc_list)

    def cleanup(self):
        self.store.close()