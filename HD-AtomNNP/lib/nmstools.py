import numpy as np
import random as rn

mDynetoMet = 1.0e-5 * 1.0e-3 * 1.0e10
Kb = 1.38064852e-23
MtoA = 1.0e10

class nmsgenerator():
    def __init__(self,xyz,nmo,fcc,T):
        self.xyz = xyz
        self.nmo = nmo
        self.fcc = fcc
        self.T = T
        self.Na = xyz.shape[0]
        self.Nf = nmo.shape[0]
        print(self.Nf)

    def genrandomstruct(self):
        rdt = np.random.random(self.Nf+1)
        rdt[0] = 0.0
        norm = np.random.random(1)[0]
        rdt = norm*np.sort(rdt)
        rdt[self.Nf] = norm
        print(rdt)

        accnm = np.zeros((self.xyz.shape), dtype=np.float32)

        for i in range(self.Nf):
            Ki = mDynetoMet * self.fcc[i]
            ci = rdt[i+1]-rdt[i]
            Sn = -1.0 if np.random.binomial(1,0.5,1) else 1.0
            Ri = Sn * MtoA * np.sqrt((3.0 * ci * float(self.Na) * Kb * self.T)/Ki)

            for mode in self.nmo:
                accnm = accnm + Ri * mode

            print(accnm)

        rxyz = np.array((1),dtype=np.float32)
        return rxyz