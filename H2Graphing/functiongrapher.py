__author__ = 'jujuman'

from matplotlib import pyplot as plt
import numpy as np

#El = float(5.0e-4)
#Tl = 1.0e-13
#Th = 1.0e-12

#Eh = np.linspace(5.5e-4, 9.0e-2, 70, endpoint=True)
#dE = (Eh/El - 1.0)
#print(dE)
#C = Tl/Th - 1.0
#print(C)
#Z = dE*C
#print(Z)
#P = np.exp(Z)
#print(P)
#plt.plot(dE, P)

#plt.show()
#--------------------------------
#           Functions
#--------------------------------

# ------------------------------------------
#               Radial Function
# ------------------------------------------
def radialfunction(X,eta,Rc):
    F = np.exp(-eta*(X-0.25)**2.0) * (0.5 * (np.cos((np.pi * X)/Rc) + 1.0))
    return F

# ------------------------------------------
#               Angular Function
# ------------------------------------------
def angularfunction(T,zeta,lam):
    F = (2.0**(1.0-zeta)) * ((1.0 + lam * np.cos(T))**zeta)
    return F


# ------------------------------------------
# Calculate The Steps for a Radial Dataset
# ------------------------------------------
def computeradialdataset(x1,x2,pts,eta,Rc,plt,scolor,slabel):

    X = np.linspace(x1, x2, pts, endpoint=True)
    F = radialfunction(X,eta,Rc)
    plt.plot(X, F, label=slabel, color=scolor)

# ------------------------------------------
# Calculate The Steps for an angular Dataset
# ------------------------------------------
def computeangulardataset(t1,t2,pts,zeta,lam,plt,scolor,slabel):

    T = np.linspace(t1, t2, pts, endpoint=True)
    F = angularfunction(T,zeta,lam)
    plt.plot(T, F, label=slabel, color=scolor)

# ------------------------------------------
# Calculate The Steps for an angular Dataset
# ------------------------------------------
def printdatatofile(f,title,X,N):
    f.write(title + ' = [')
    for i in range(0,N):
        if i < N-1:
            s = "{:.7e}".format(X[i]) + ','
        else:
            s = "{:.7e}".format(X[i])
        f.write(s)
    f.write(']\n')

#--------------------------------
#          Parameters
#--------------------------------
#File name
pf = 'test.params' # Output filename

Nrr = 16
Nat = 4
Nar = 8
Nzt = 8

Rc = 6.0
Atyp = '[H,C,O,N]'

#--------------------------------
#           Program
#--------------------------------
Nrt = Nrr * Nat
EtaR = np.zeros(Nrr)
for i in range(0,Nrr):
    stddev = float(Nrr)/4.5
    step = 3.5 * np.exp(-(float(i-0.0)**2)/(2.0*(stddev)**2.0))
    computeradialdataset(0.25, Rc, 1000, step, 6.0, plt, 'blue', 'eta = 0.1')
    EtaR[i] = step

plt.show()

EtaA = np.zeros(Nar)
Zeta = np.zeros(Nzt)

Nat = Nar * Nzt
for i in range(0,Nzt):
<<<<<<< HEAD
    step = float(i)+1.0
=======
    step = float(i+1.0)
>>>>>>> b7f74bfd09e4302c8a6c4a8b562ccc9dd9795342
    computeangulardataset(0.0,2.0*np.pi,1000,step,1.0,plt, 'red', 'eta = 0.1')
    Zeta[i] = step

plt.show()

for i in range(0,Nar):
    stddev = float(Nar)/4.5
    step = 3.5 * np.exp(-(float(i-0.0)**2)/(2.0*(stddev)**2.0))
    computeradialdataset(0.25, Rc, 1000, step, 6.0, plt, 'blue', 'eta = 0.1')
    EtaA[i] = step

plt.show()

Nt = Nrt + Nat
print('Total N Size: ',Nt)

# Open File
f = open(pf,'w')

#Write data to parameters file
f.write('Rc = ' + "{:.7e}".format(Rc) + '\n')
printdatatofile(f,'EtaR',EtaR,Nrr)
printdatatofile(f,'EtaA',EtaA,Nar)
printdatatofile(f,'Zeta',Zeta,Nzt)
f.write('Atyp = ' + Atyp + '\n')
