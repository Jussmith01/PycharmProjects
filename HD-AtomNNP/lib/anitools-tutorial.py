# Import pyanitools
import  pyanitools as pyt

# path the the store file
store_file = '/home/jujuman/Research/ANI-DATASET/rxn_db_mig.h5'

# Declare the loader, opens the store
loader = pyt.anidataloader(store_file)

# Load the entire store into memory
loader.totalload()

# Loop over store data
for i in range(loader.size()):
    data = loader.getdata(i)
    print(data[0])

# Closes the store file
loader.cleanup()