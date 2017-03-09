import numpy as np

min_dist = 1.50
tolerance = 100
box_size = [25.0, 25.0, 25.0]
types = ['H','C','N','O']
part = [0.6,0.2,0.05,0.15]

#---------------------Program-----------------
def writexyzfile (fn,xyz,typ):
    f = open(fn, 'w')
    N = len(typ)
    print('N ATOMS: ',typ)
    for m in xyz:
        f.write(str(N)+'\n')
        f.write('      comment\n')
        #print(m)
        for i in range(N):
            x = m[i,0]
            y = m[i,1]
            z = m[i,2]
            f.write(typ[i] + ' ' + "{:.7f}".format(x) + ' ' + "{:.7f}".format(y) + ' ' + "{:.7f}".format(z) + '\n')
        f.write('\n')
    f.close()

pad = 0.25 * min_dist
vol = box_size[0] * box_size[1] * box_size[2]
atom_crd = np.zeros((int(vol*min_dist*1.1),3),dtype=np.float32)
atom_spc = []

box_size = np.array(box_size)

atm_cnt = 0
rho = 0.0
tol = tolerance
while tol is not 0:
    crds = np.random.rand((3)) * (box_size - 2.0*pad) + pad
    spec = np.random.choice(types, 1, p=part)[0]

    fit = True
    for i in range (atm_cnt):
        if np.linalg.norm(atom_crd[i] - crds) < min_dist:
            fit = False
            break

    if fit:
        atom_crd[atm_cnt] = crds
        atom_spc.append(spec)
        rho = atm_cnt / vol
        atm_cnt = atm_cnt + 1
        tol = tolerance
    else:
        tol = tol - 1

    print(' Fit =', '{:>5}'.format(str(fit)),' Rho =', "{:.4f}".format(rho), ' tol =', '{:>3}'.format(tol), ' Z =', spec, ' R =', crds, ' Na =',atm_cnt)

print(len(atom_spc), ' : ', atom_crd.shape, ' : ', atm_cnt)
writexyzfile('hotmess.xyz', atom_crd[:atm_cnt].reshape(1,atm_cnt,3), atom_spc[:atm_cnt])
