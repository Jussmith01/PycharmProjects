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
        self.Na      = np.array(self.species).shape[0]

# Load an ANI data set stored as a pandas dataframe/HDF5 file format
class anidataloader:

    #----------------------------------------
    # Construct the class and open the store
    #----------------------------------------
    def __init__(self, store_file, complib = 'blosc', complevel = 8):
        # opening file
        self.store = pd.HDFStore(store_file, complib=complib, complevel=complevel)

    #----------------------------------------
    #----- Load and split the dataset -------
    #----------------------------------------
    def splitload(self, N):
        self.crd_list = []
        self.eng_list = []
        self.spc_list = []
        for x in self.store.get_node(""):
            for i in x._v_children:
                ae = anidataextractor(self.store, x._v_name, i)
                self.crd_list.append(np.array_split(ae.coords, N))
                self.eng_list.append(np.array_split(ae.energy, N))
                self.spc_list.append(ae.species)

    #------------------------------------------
    #------------ Get loaded data -------------
    # Arg 1: idx - file index
    # Arg 2: dl  - list of partition indicies
    #------------------------------------------
    '''
    def getdata(self,idx=0,dl=[0]):

        if max(dl) >= len(self.crd_list):
            raise (IndexError('Index given is outside of the array range.'))

        return [np.concatenate([self.crd_list[idx][j] for j in dl]),
                np.concatenate([self.eng_list[idx][j] for j in dl]),
                self.spc_list[idx]]
    '''

    def getdata(self,dl=[0]):

        if max(dl) >= len(self.crd_list):
            raise (IndexError('Index given is outside of the array range.'))

        for i in range(self.size()):
            yield [np.concatenate([self.crd_list[i][j] for j in dl]),
                   np.concatenate([self.eng_list[i][j] for j in dl]),
                   self.spc_list[i]]

    #-----------------------------------------
    #-------- Returns entire data set --------
    #-----------------------------------------
    def getalldata(self):
        return [np.concatenate([self.crd_list]),
                np.concatenate([self.eng_list]),
                self.spc_list]

    #-----------------------------------------
    #----- Loads data set without split ------
    #-----------------------------------------
    def totalload(self):
        self.crd_list = []
        self.eng_list = []
        self.spc_list = []
        for x in self.store.get_node(""):
            child = [str(i) for i in x._v_children]
            child = sorted(child, key=lambda x: int(x.split('mol')[1].split('.')[0]))
            for i in child:
                print(x._v_name)
                ae = anidataextractor(self.store, x._v_name, i)
                self.crd_list.append([ae.coords])
                self.eng_list.append([ae.energy])
                self.spc_list.append(ae.species)

    # -----------------------------------------
    # ----- Loads data set without split ------
    # -----------------------------------------
    def load_node(self, node):
        self.crd_list = []
        self.eng_list = []
        self.spc_list = []
        x = self.store.get_node(node)
        child = [str(i) for i in x._v_children]
        #print(child)
        child = sorted(child, key=lambda x: int(x.split('mol')[1].split('.')[0]))
        for i in child:
            #print(x._v_name)
            ae = anidataextractor(self.store, x._v_name, i)
            self.crd_list.append([ae.coords])
            self.eng_list.append([ae.energy])
            self.spc_list.append(ae.species)

    #--------------------------------------------
    #---------- Returns the Node list -----------
    #--------------------------------------------
    def get_node_list(self):
        return [x._v_name for x in self.store.get_node("")]

    #--------------------------------------------
    #----- The number of data files loaded ------
    #--------------------------------------------
    def size(self):
        return len(self.spc_list)

    #--------------------------------------------
    #----------- Close the store file -----------
    #--------------------------------------------
    def cleanup(self):
        self.store.close()

class datapacker:
    # ----------------------------------------
    # Construct the class and open the store
    # ----------------------------------------
    def __init__(self, store_file, complib='blosc', complevel=8):
        # opening file
        self.store = pd.HDFStore(store_file, complib=complib, complevel=complevel)
        self.clib = complib
        self.clev = complevel

    # ----------------------------------------
    #            Store the data
    # ----------------------------------------
    def store_data(self, store_loc, xyz, energy, species):
        cols = [["x" + str(z), "y" + str(z), "z" + str(z)] for z in list(range(int(xyz.shape[1]/3)))]
        cols = [item for sublist in cols for item in sublist]
        cols = [('coordinates',l) for l in cols]
        cols = pd.MultiIndex.from_tuples(cols)  # Notice these are un-named

        df_xyz = pd.DataFrame(xyz, columns=cols)
        df_xyz['energy'] = energy
        df_xyz.index.name = 'conformer'

        self.store.put(store_loc, df_xyz, complib=self.clib, complevel=self.clev, format='table')
        self.store.get_storer(store_loc).attrs.species = species

    # --------------------------------------------
    # ----------- Close the store file -----------
    # --------------------------------------------
    def cleanup(self):
        self.store.close()