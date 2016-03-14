__author__ = 'jujuman'

"""
Animate by File
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# -----------------------
# readfile into np array
# -----------------------
def getfltsfromfile(file, cols):
    # Open and read the data file
    infile = open(file, 'r')

    infile_s = []

    for line in infile:
        row = line.strip().split(" ")
        infile_s.append(row)

    # Truncate and convert to numpy array
    nparray = np.array(infile_s)
    data = nparray[:-1, cols]
    data = np.array(data, dtype=float)
    return data

#------------ INPUT PARAMETERS ----------------
N = 500 #Number of elements per graph

filename = 'demo.mp4'

DPI = 100
h_in_inches = 14
w_in_inches = 14

#---------- END INPUT PARAMETERS --------------

# Set up formatting for the movie files
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)

fig, ax = plt.subplots()

fig.set_size_inches(w_in_inches, h_in_inches)

# ------------
# Load Data
# ------------
user = os.environ['USER']
#user = 'jujuman'
dir = '/Research/ANN-Test-Data/GDB-11/train3/'

file = 'graph.dat'

#data1 = getfltsfromfile('/home/' + user + dir + file, [0])
#data1 = data1 * 1.0 + -180.0
data2 = getfltsfromfile('/home/' + user + dir + file, [1])
data3 = getfltsfromfile('/home/' + user + dir + file, [2])

M = (data3.shape[0]+1)/N - 1
print ('Number of Fames:', M)

axis_x = data2[0:N]

target_y = data2[0:N]
actual_y = data3[0:N]

#line, = ax.plot(axis_x, target_y,color='blue',linewidth=4)
#line, = ax.plot(axis_x, actual_y,color='red',linewidth=4)

scat = ax.scatter(axis_x, target_y,color='blue',linewidth=4,label='B3LYP/6-31g* vs B3LYP/6-31g*')
scat = ax.scatter(axis_x, actual_y,color='red',linewidth=4,label='AtomicNN vs B3LYP/6-31g*')

diff = target_y.min() - target_y.max()
ax.set_ylim([target_y.min() + diff  * 0.05,target_y.max() - diff  * 2.00])
ax.set_xlim([axis_x.min(),axis_x.max()])

plt.title('Scatter Plot of $H_2 O$ During Training')
plt.xlabel('Energy (Hartrees)')
plt.ylabel('Energy (Hartrees)')
plt.legend(bbox_to_anchor=(0.1, 0.95), loc=2, borderaxespad=0.)

def animateplot(i): # Animates a line plot
    print ('Processing Frame ', int(i+1), ' of ', int(M))
    next_y = data3[i*N:N+i*N]
    line.set_ydata(next_y)  # update the data
    return line,

def animatescat(i): # Animates a scatter plot
    print ('Processing Frame ', int(i+1), ' of ', int(M))
    next_y = data3[i*N:N+i*N]
    data = np.hstack((axis_x[:N,np.newaxis], next_y[:N, np.newaxis]))
    scat.set_offsets(data)
    return scat,

ani = animation.FuncAnimation(fig, animatescat, np.arange(0, M, 1),
                              interval=50, blit=True,repeat=True)
ani.save(filename,dpi = DPI)

#plt.show()