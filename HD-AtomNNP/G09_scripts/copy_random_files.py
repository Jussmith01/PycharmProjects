from shutil import copyfile
import random as rn
import os

srcdir = "/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnntsgdb11_08/testdata/"
dstdir = "/home/jujuman/Dropbox/GDB-08_DATA/"

files = os.listdir(srcdir)
files = [i for i in files if '_train' in i]

rn.shuffle(files)

for i in files[0:99]:
    copyfile(srcdir+i, dstdir+i)

