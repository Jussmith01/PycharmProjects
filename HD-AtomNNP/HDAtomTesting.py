__author__ = 'jujuman'

import libpyNeuroChem as nc

cnstfile = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/train_08-a3.1A_r4.5_dn2/rHCNO-4.5A_32-3.1A_a8-8.params'
saefile = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/sae_6-31gd.dat'
nnfdir = '/home/jujuman/Research/trainingcases/wB97X-631gd-train-highgarden/train_08-a3.1A_r4.5_dn2/networks/'

test = nc.pyNeuroChem(cnstfile,saefile,nnfdir,0)

list = [1.0,2.0,3.0]

test.testfunction(list,3)