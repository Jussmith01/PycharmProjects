import os
import re
from shutil import copyfile

S1 = 'N'
S2 = 'O'

dir = '/home/jujuman/Dropbox/ChemSciencePaper.AER/networks/ANI-CHNOSF-1_cv/cv_c08e_ntw_4/networks/'

files = os.listdir(dir)

SF1 = [ f for f in files if '-'+S1+'.' in f and 'bparam' not in f and 'wparam' not in f]
SF2 = [ f for f in files if '-'+S2+'.' in f and 'bparam' not in f and 'wparam' not in f]

SF1.sort()
SF2.sort()

print(SF1)
print(SF2)

for f1,f2 in zip(SF1,SF2):
    print('1:',dir+f1,'->', dir+f1+'_tmp')
    copyfile(dir + f1, dir + f1 + '_tmp')
    print('2:',dir+f2,'->', dir+f1)
    copyfile(dir + f2, dir + f1)
    print('3:',dir+f1+'_tmp','->', dir+f2)
    copyfile(dir + f1 + '_tmp', dir + f2)
    print('4:(del)',dir+f1+'_tmp')
    os.remove(dir + f1 + '_tmp')
    #copyfile(dir+f1, dir+f2)
    #copyfile(dir+f2, dir+f1)


#copyfile(src, dst)