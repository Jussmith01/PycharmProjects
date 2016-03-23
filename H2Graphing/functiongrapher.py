__author__ = 'Justin Smith'

from matplotlib import pyplot as plt
import numpy as np
import scipy.interpolate

#--------------------------------
#           Functions
#--------------------------------

# ------------------------------------------
#          Radial Function Cos w/ Sqrt
# ------------------------------------------
def radialfunctionsqrt(X,eta,Rc,Rs):
    F = np.sqrt(np.exp(-eta*(X-Rs)**2.0) * (0.5 * (np.cos((np.pi * X)/Rc) + 1.0)))
    return F

# ------------------------------------------
#          Radial Function Cos
# ------------------------------------------
def radialfunctioncos(X,eta,Rc,Rs):
    F = np.exp(-eta*(X-Rs)**2.0) * (0.5 * (np.cos((np.pi * X)/Rc) + 1.0))
    return F

# ------------------------------------------
#               Radial Function
# ------------------------------------------
def radialfunction2(X,eta,Rc,Rs):

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
    F = radialfunctioncos(X,eta,Rc,Rs)
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

# ------------------------------------------
#         Simple Addition Function
# ------------------------------------------
def add (x,y):
    return x+y

# ----------------------------------------------------
# Show a 2d Contour Plot of the Angular Env Functions
# ----------------------------------------------------
def show2dcontangulargraph (ShfA,ShfZ,eta,zeta,Rc,func):
    N = 1000000
    x, y = 12.0 * np.random.random((2, N)) - 6.0

    print(x)

    R = np.sqrt(x**2 + y**2)
    T = np.arctan2(x,y)

    z = np.zeros(N)

    for i in ShfZ:
        for j in ShfA:
            print( 'ShfZ: ' + str(i) + ' ShfA: ' + str(j) )
            zt = angularfunction(T,zeta,1.0,i) * radialfunctionsqrt(R,eta,Rc,j) * radialfunctionsqrt(R,eta,Rc,j)

            for k in range(1,z.shape[0]):
                z[k] = func(z[k],zt[k])

    # Set up a regular grid of interpolation points
    xi, yi = np.linspace(x.min(), y.max(), 300), np.linspace(x.min(), y.max(), 300)
    xi, yi = np.meshgrid(xi, yi)

    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')

    plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
           extent=[x.min(), x.max(), y.min(), y.max()])

    plt.colorbar()
    plt.show()

# ----------------------------------------------------
# Show a 2d Contour Plot of the Radial Env Functions
# ----------------------------------------------------
def show2dcontradialgraph (ShfR,eta,Rc,func):
    N = 1000000
    x, y = 14.0 * np.random.random((2, N)) - 7.0

    print(x)

    R = np.sqrt(x**2 + y**2)
    T = np.arctan2(x,y)

    z = np.zeros(N)

    for j in ShfR:
        print( 'ShfZ: ' + str(i) + ' ShfA: ' + str(j) )
        zt = radialfunctioncos(R,eta,Rc,j)

        for k in range(1,z.shape[0]):
            z[k] = func(z[k],zt[k])

    # Set up a regular grid of interpolation points
    xi, yi = np.linspace(x.min(), y.max(), 300), np.linspace(x.min(), y.max(), 300)
    xi, yi = np.meshgrid(xi, yi)

    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')

    plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
           extent=[x.min(), x.max(), y.min(), y.max()])

    plt.colorbar()
    plt.show()

# ****************************************************
#--------------------------------
#         Set Parameters
#--------------------------------
#File nam
pf = 'rHCNO-4-a1-4.params' # Output filename

Nrr = 4
Na = 4
Nar = 1
Nzt = 4

Rc = 4.0
Atyp = '[H,C,N,O]'
EtaR = 4.0
EtaA1 = 0.0001
Zeta = 6.0

# ****************************************************

#--------------------------------
#         Main Program
#    (Build Env Params File)
#--------------------------------
Nrt = Nrr * Na
ShfR = np.zeros(Nrr)

#Now instead of multiple etaR we use multiple shifts with a single large EtaR
for i in range(0,Nrr):
    stepsize = Rc / float(Nrr+1.0)
    step = i * stepsize + 0.50
    computeradialdataset(0.5, Rc, 1000, EtaR, Rc,step, plt, 'blue', 'eta = '+ str(EtaR))
    ShfR[i] = step

plt.title('Radial Environment Functions (REF)')
plt.ylabel('REF Output')
plt.xlabel('Angstroms')
plt.show()

#Uncomment for pretty contour plots of the radial environments using a sum and then max function
#show2dcontradialgraph(ShfR,EtaR,Rc,add)
#show2dcontradialgraph(ShfR,EtaR,Rc,max)

ShfZ = np.zeros(Nzt)

Nat = Nar * (Na*(Na+1)/2) * Nzt

for i in range(0,Nzt):
    stepsize = (2.0 * np.pi) / (float(Nzt))
    step = i*stepsize
    computeangulardataset(0.5,2.0*np.pi,1000,Zeta,1.0,step,plt, 'red', 'zeta = ' + str(Zeta))
    ShfZ[i] = step

plt.show()


ShfA = np.zeros(Nar)

for i in range(0,Nar):
    stepsize = Rc / float(Nar+1.0)
    step = i * stepsize + 0.5
    computeradialdataset(0.5, Rc, 1000, EtaA1, Rc,step, plt, 'blue', 'eta = '+ str(EtaA1))
    ShfA[i] = step

plt.show()

#Uncomment for pretty contour plots of the angular environments using a sum and then max function
#show2dcontangulargraph(ShfA,ShfZ,EtaA1,Zeta,Rc,add)
#show2dcontangulargraph(ShfA,ShfZ,EtaA1,Zeta,Rc,max)

Nt = Nat + Nrt
print('Total Environmental Vector Size: ',int(Nt))

# Open File
f = open(pf,'w')

#Write data to parameters file
f.write('Rc = ' + "{:.4}".format(Rc) + '\n')
f.write('EtaR = ' + "{:.4}".format(EtaR) + '\n')
printdatatofile(f,'ShfR',ShfR,Nrr)
f.write('Zeta = ' + "{:.4}".format(Zeta) + '\n')
printdatatofile(f,'ShfZ',ShfZ,Nzt)
f.write('EtaA = ' + "{:.4}".format(EtaA1) + '\n')
printdatatofile(f,'ShfA',ShfA,Nar)
f.write('Atyp = ' + Atyp + '\n')

f.close()
