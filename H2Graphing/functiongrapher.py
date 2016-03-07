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
def radialfunction(X,eta,Rc,Rs):
    F = np.exp(-eta*(X-Rs)**2.0) * (0.5 * (np.cos((np.pi * X)/Rc) + 1.0))
    return F

# ------------------------------------------
#               Angular Function
# ------------------------------------------
def angularfunction(T,zeta,lam,Ts):
    F = (2.0**(1.0-zeta)) * ((1.0 + lam * np.cos(T-Ts))**zeta)
    return F


# ------------------------------------------
# Calculate The Steps for a Radial Dataset
# ------------------------------------------
def computeradialdataset(x1,x2,pts,eta,Rc,Rs,plt,scolor,slabel):

    X = np.linspace(x1, x2, pts, endpoint=True)
    F = radialfunction(X,eta,Rc,Rs)
    plt.plot(X, F, label=slabel, color=scolor)

# ------------------------------------------
# Calculate The Steps for an angular Dataset
# ------------------------------------------
def computeangulardataset(t1,t2,pts,zeta,lam,Ts,plt,scolor,slabel):

    T = np.linspace(t1, t2, pts, endpoint=True)
    F = angularfunction(T,zeta,lam,Ts)
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
#File nam
pf = 'rHCNO-12-a4-8.params' # Output filename

Nrr = 8
Na = 4
Nar = 4
Nzt = 8

Rc = 6.0
Atyp = '[H,C,N,O]'
EtaR = 8.0
EtaA = 4.0
Zeta = 32.0

#--------------------------------
#           Program
#--------------------------------
Nrt = Nrr * Na
ShfR = np.zeros(Nrr)

#for i in range(0,Nrr):
#    stddev = float(Nrr)/4.5
#    step = 3.5 * np.exp(-(float(i-0.0)**2)/(2.0*(stddev)**2.0))
#    computeradialdataset(0.5, Rc, 1000, step, Rc,0.25, plt, 'blue', 'eta = 0.1')
#    EtaR[i] = step

#plt.show()

#Now instead of multiple etaR we use multiple shifts with a single large EtaR
for i in range(0,Nrr):
    stepsize = Rc / float(Nrr+1.0)
    #step = (i**1.5)*0.1+0.5
    step = i * stepsize + 0.5
    computeradialdataset(0.5, Rc, 1000, EtaR, Rc,step, plt, 'blue', 'eta = '+str(EtaR))
    ShfR[i] = step

plt.show()

ShfZ = np.zeros(Nzt)

Nat = Nar * Nzt

#for i in range(0,Nzt):
#    step = float(i+1.0)
#    computeangulardataset(0.0,2.0*np.pi,1000,step,1.0,0.0,plt, 'red', 'eta = 0.1')
#    Zeta[i] = step

#plt.show()

for i in range(0,Nzt):
    stepsize = (2.0 * np.pi) / (float(Nzt))
    step = i*stepsize
    computeangulardataset(0.0,2.0*np.pi,1000,Zeta,1.0,step,plt, 'red', 'zeta = ' + str(Zeta))
    ShfZ[i] = step

plt.show()

#for i in range(0,Nar):
#    stddev = float(Nar)/4.5
#    step = 3.5 * np.exp(-(float(i-0.0)**2)/(2.0*(stddev)**2.0))
#    computeradialdataset(0.25, Rc, 1000, step, 6.0,0.25, plt, 'blue', 'eta = 0.1')
#    EtaA[i] = step

#plt.show()

ShfA = np.zeros(Nar)

for i in range(0,Nar):
    stepsize = Rc / float(Nar+1.0)
    #step = (i**1.5)*0.1+0.5
    step = i * stepsize + 0.5
    computeradialdataset(0.5, Rc, 1000, EtaA, Rc,step, plt, 'blue', 'eta = '+str(EtaR))
    ShfA[i] = step

plt.show()

Nt = Nat * (Na*(Na+1)/2) + Nrt
print('Total N Size: ',int(Nt))

# Open File
f = open(pf,'w')

#Write data to parameters file
f.write('Rc = ' + "{:.7e}".format(Rc) + '\n')
f.write('EtaR = ' + "{:.7e}".format(EtaR) + '\n')
printdatatofile(f,'ShfR',ShfR,Nrr)
f.write('Zeta = ' + "{:.7e}".format(Zeta) + '\n')
printdatatofile(f,'ShfZ',ShfZ,Nzt)
f.write('EtaA = ' + "{:.7e}".format(EtaA) + '\n')
printdatatofile(f,'ShfA',ShfA,Nar)
f.write('Atyp = ' + Atyp + '\n')

f.close()
