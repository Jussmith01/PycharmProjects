__author__ = 'jujuman'

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import pyplot as plt
import numpy as np
import scipy.interpolate

#--------------------------------
#           Functions
#--------------------------------

# ------------------------------------------
#               Radial Function
# ------------------------------------------
def radialfunction(X,eta,Rc,Rs):
    F = np.sqrt(np.exp(-eta*(X-Rs)**2.0) * (0.5 * (np.cos((np.pi * X)/Rc) + 1.0)))
    return F

# ------------------------------------------
#               Radial Function
# ------------------------------------------
def radialfunction2(X,eta,Rc,Rs):

    #F = (1.0/(np.exp(-1.0))) * np.exp(-(1.0/(1.0-((X-Rs)/Rc)**2)))

    F=np.zeros(X.shape[0])
    c = (Rc-Rs)
    for i in range(0,X.shape[0]):
        if X[i] < Rc:
            F[i] = (np.exp(1.0) * np.exp(-(1.0/(1.0-((X[i]-Rs)/(Rc-Rs))**(2.0)))))**(eta)
        else :
            F[i]=0.0

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
    plt.plot(X, F, label=slabel, color=scolor, linewidth=2)

# ------------------------------------------
# Calculate The Steps for a Radial Dataset
# ------------------------------------------
def computeradial2dataset(x1,x2,pts,eta,Rc,Rs,plt,scolor,slabel):

    X = np.linspace(x1, x2, pts, endpoint=True)
    F = radialfunction2(X,eta,Rc,Rs)
    plt.plot(X, F, label=slabel, color=scolor, linewidth=2)

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

def show2dcontangulargraph (eta,zeta,Rs,Rc,M):
    N = 1000000
    x, y = 12.0 * np.random.random((2, N)) - 6.0

    print(x)

    R = np.sqrt(x**2 + y**2)
    T = np.arctan2(x,y)

    z = angularfunction(T,zeta,1.0,0.0) * radialfunction(R,eta,Rc,Rs)**2

    for i in range(1,M):
        Ts = float(i) * (2.0*np.pi/float(M))
        print(Ts)
        if Ts < np.pi:
            zt = angularfunction(T,zeta,1.0,Ts) * radialfunction(R,eta,Rc,Rs)**2
        else:
            zt = angularfunction(T,zeta,1.0,Ts) * radialfunction2(R,eta,Rc,Rs)**2

        z = z + zt

    # Set up a regular grid of interpolation points
    xi, yi = np.linspace(x.min(), y.max(), 300), np.linspace(x.min(), y.max(), 300)
    xi, yi = np.meshgrid(xi, yi)

    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')

    plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
           extent=[x.min(), x.max(), y.min(), y.max()])

    plt.colorbar()
    plt.show()

#--------------------------------
#          Parameters
#--------------------------------
#File nam
pf = 'rHCNO-24-a1-4.params' # Output filename

Nrr = 24
Na = 4
Nar = 1
Nzt = 2

Rc = 6.0
Atyp = '[H,C,N,O]'
EtaR = 4.0
EtaA1 = 0.1
EtaA2 = 0.5
Zeta = 64.0

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
    step = i * stepsize + 0.25
    computeradialdataset(0.5, Rc, 1000, EtaR, Rc,step, plt, 'blue', 'eta = '+ str(EtaR))
    ShfR[i] = step


plt.title('Radial Environment Functions (REF)')
#plt.title('SCAN: Formic Acid Energy vs. H-O-H Angle')
plt.ylabel('REF Output')
plt.xlabel('Angstroms')
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
    computeangulardataset(0.5,2.0*np.pi,1000,Zeta,1.0,step,plt, 'red', 'zeta = ' + str(Zeta))
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
    computeradial2dataset(0.5, Rc, 1000, EtaA2, Rc,step, plt, 'blue', 'eta = '+ str(EtaA2))
    computeradialdataset(0.5, Rc, 1000, EtaA1, Rc,step, plt, 'blue', 'eta = '+ str(EtaA1))
    ShfA[i] = step

plt.show()

show2dcontangulargraph(EtaA2,Zeta,0.5,Rc,Nzt)

Nt = Nat * (Na*(Na+1)/2) + Nrt
print('Total N Size: ',int(Nt))

# Open File
f = open(pf,'w')

#Write data to parameters file
f.write('Rc = ' + "{:.4}".format(Rc) + '\n')
f.write('EtaR = ' + "{:.4}".format(EtaR) + '\n')
printdatatofile(f,'ShfR',ShfR,Nrr)
f.write('Zeta = ' + "{:.4}".format(Zeta) + '\n')
printdatatofile(f,'ShfZ',ShfZ,Nzt)
f.write('EtaA = ' + "{:.4}".format(EtaA2) + '\n')
printdatatofile(f,'ShfA',ShfA,Nar)
f.write('Atyp = ' + Atyp + '\n')

f.close()
