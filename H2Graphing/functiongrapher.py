__author__ = 'Justin Smith'

from matplotlib import pyplot as plt
import numpy as np
import scipy.interpolate

#--------------------------------
#           Functions
#--------------------------------

# ------------------------------------------
#          Radial Function Cos
# ------------------------------------------
def radialfunctioncos(X,eta,Rc,Rs):
    F = np.exp(-eta*(X-Rs)**2.0) * (0.5 * (np.cos((np.pi * X)/Rc) + 1.0))
    return F

# ------------------------------------------
#          Radial Function Cos
# ------------------------------------------
def angularradialfunctioncos(X,eta,Rc,Rs):
    F = np.sqrt(np.exp(-eta*(X-Rs)**2.0) * (0.5 * (np.cos((np.pi * X)/Rc) + 1.0)))
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
def computeangularradialdataset(x1,x2,pts,eta,Rc,Rs,plt,scolor,slabel):

    X = np.linspace(x1, x2, pts, endpoint=True)
    F = angularradialfunctioncos(X,X,eta,Rc,Rs)
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
def show2dcontangulargraph (ShfA,ShfZ,eta,zeta,Rc,func,title):
    N = 1000000
    x, y = 12.0 * np.random.random((2, N)) - 6.0

    #print(x)

    R = np.sqrt(x**2 + y**2)
    T = np.arctan2(x,y)

    z = np.zeros(N)

    for i in ShfZ:
        for j in ShfA:
            print( 'ShfZ: ' + str(i) + ' ShfA: ' + str(j) )
            zt = angularfunction(T,zeta,1.0,i) * angularradialfunctioncos(R,eta,Rc,j) * angularradialfunctioncos(R,eta,Rc,j)
            #zt = angularradialfunctioncos(R1,R2,eta,Rc,j)

            for k in range(1,z.shape[0]):
                z[k] = func(z[k],zt[k])
                #print(z[k])

    # Set up a regular grid of interpolation points
    xi, yi = np.linspace(x.min(), y.max(), 600), np.linspace(x.min(), y.max(), 600)
    xi, yi = np.meshgrid(xi, yi)

    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')

    plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
           extent=[x.min(), y.max(), x.min(), y.max()])

    plt.title(title)
    plt.ylabel('Angstroms')
    plt.xlabel('Angstroms')

    plt.colorbar()
    plt.show()

# ----------------------------------------------------
# Show a 2d Contour Plot of the Radial Env Functions
# ----------------------------------------------------
def show2dcontradialgraph (ShfR,eta,Rc,func,title):
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

    plt.title(title)
    plt.ylabel('Angstroms')
    plt.xlabel('Angstroms')

    plt.colorbar()
    plt.show()

# ****************************************************
#--------------------------------
#         Set Parameters
#--------------------------------
#File nam
pf = 'rHCNO-32-a1-32.params' # Output filename


Nrr = 32
Na = 4
Nar = 1
Nzt = 32

Rc = 6.0
Atyp = '[H,C,O,N]'
EtaR = 64.0
EtaA1 = 0.0001
Zeta = 64.0


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
#show2dcontradialgraph(ShfR,EtaR,Rc,add,'Sum Radial Output')
#show2dcontradialgraph(ShfR,EtaR,Rc,max,'Max Radial Output')

ShfZ = np.zeros(Nzt)

Nat = Nar * (Na*(Na+1)/2) * Nzt

for i in range(0,Nzt):
    stepsize = (2.0 * np.pi) / (float(Nzt))
    step = i*stepsize
    computeangulardataset(0.0,2.0*np.pi,1000,Zeta,1.0,step,plt, 'red', 'zeta = ' + str(Zeta))
    ShfZ[i] = step

plt.title('Angular (Only Angular) Environment Functions (AAEF)')
plt.ylabel('AAEF Output')
plt.xlabel('Radians')
plt.show()


ShfA = np.zeros(Nar)

for i in range(0,Nar):
    stepsize = Rc / float(Nar+1.0)
    step = (i * stepsize + 0.5)
    computeradialdataset(0.5, Rc, 1000, EtaA1, Rc,step, plt, 'blue', 'eta = '+ str(EtaA1))
    ShfA[i] = step

plt.title('Angular (Only Radial) Environment Functions (AREF)')
plt.ylabel('AREF Output')
plt.xlabel('Angstroms')
plt.show()

#Uncomment for pretty contour plots of the angular environments using a sum and then max function
#show2dcontangulargraph(ShfA,ShfZ,EtaA1,Zeta,Rc,add,'Sum Angular Output')
#show2dcontangulargraph(ShfA,ShfZ,EtaA1,Zeta,Rc,max,'Max Angular Output')

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
