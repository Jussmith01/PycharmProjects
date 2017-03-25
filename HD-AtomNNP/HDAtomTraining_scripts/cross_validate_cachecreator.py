import numpy as np
import pyanitools as pyt
from pyNeuroChem import cachegenerator as cg

saef   = "/home/jujuman/Dropbox/ChemSciencePaper.AER/ANI-c08e-ntwk/sae_6-31gd.dat"
h5file = "/home/jujuman/Research/ANI-DATASET/ani_data_c08e.h5"

store_dir = "/home/jujuman/Research/CrossValidation/cache-c08e-"

adl = pyt.anidataloader(h5file)

adl.split_load(10)

train_idx = [[2, 3, 4, 5, 6, 7, 8, 9],
             [0, 1, 4, 5, 6, 7, 8, 9],
             [0, 1, 2, 3, 6, 7, 8, 9],
             [0, 1, 2, 3, 4, 5, 8, 9],
             [0, 1, 2, 3, 4, 5, 6, 7]
             ]

valid_idx = [[0]
            ,[2]
            ,[4]
            ,[6]
            ,[8]
             ]

r = 0
for t,v in zip(train_idx, valid_idx):
    print("Working on index: ", r)

    cachet = cg('_train', saef, store_dir + str(r) + '/')
    cachev = cg('_valid', saef, store_dir + str(r) + '/')
    for i in range(0, adl.size()):
        print("Working on : ", i, ' from set ', r)

        t_data = adl.getdata(i,t)
        v_data = adl.getdata(i,v)

        #cn = 0
        #for x,y in zip(v_data[0],v_data[1]):
        #    print('Element ',cn,': ',x,'\n',y)

        #print(t_data[0].shape, ' : ', t_data[1].shape, ' : ', t_data[2].shape)
        #print(v_data[0].shape, ' : ', v_data[1].shape, ' : ', v_data[2].shape)

        if t_data[0].shape[0] != t_data[1].shape[0]:
            print('T DATA SIZE MISSMATCH!')
            exit(1)

        if t_data[0].shape[0] != t_data[1].shape[0]:
            print('V DATA SIZE MISSMATCH!')
            exit(1)

        ttest = np.array(t_data[0], dtype = np.float32, order='C')
        vtest = np.array(v_data[0], dtype = np.float32, order='C')

        cachet.insertdata(ttest, t_data[1], list(t_data[2]))
        cachev.insertdata(vtest, v_data[1], list(v_data[2]))

    print('Making meta files...')
    cachet.makemetadata()
    cachev.makemetadata()
    r = r + 1

adl.cleanup()
