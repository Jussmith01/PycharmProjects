#!/usr/bin/env python
"""
Draw a graph with matplotlib, color edges.
You must have matplotlib>=87.7 for this to work.
"""
__author__ = "Justin Smith"
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import matplotlib.animation as animation


# -----------------------
# readfile into np array
# -----------------------
def getfltsfromfile(file, row, set, N):
    # Open and read the data file
    infile = open(file, 'r')

    infile_s = []

    count = 0
    data = np.array(0)
    for line in infile:
        if count == row:
            sdat = line.strip().split(",")
            infile_s.append(sdat)

            # Truncate and convert to numpy array
            nparray = np.array(infile_s)
            data = nparray[0, set * N:set * N + N]
            data = np.array(data, dtype=float)

        count += 1

    return data


def normalize(array):
    sum = 0.0
    array = np.abs(array)
    for i in array:
        sum += i
    array = array / sum
    array = np.sqrt(array)
    return array

def makenodes(G, nodes, N):
    fixed_positions = dict()
    idx = 0
    for i in range(0, N):
        for j in range(0, nodes[i]):
            mult = nodes[0] / float(nodes[i])
            shift = mult * (nodes[i] - 1) / 2
            fixed_positions[idx] = (i, j * mult - shift)
            idx += 1
    fixed_nodes = fixed_positions.keys()
    pos = nx.spring_layout(G, pos=fixed_positions, fixed=fixed_nodes)

    return pos

def getnodesinfo(node_colors, nodes, a, N):
    fixed_positions = dict()
    idx = 0
    for i in range(0, N):
        for j in range(0, nodes[i]):
            node_colors.append([1.0, 1.0 - a[i][j], 1.0 - a[i][j]])

def makeedges(G, nodes, N):
    ls = 0
    Ne = 0
    for i in range(0, N - 1):
        for j in range(0, nodes[i]):
            for k in range(0, nodes[i + 1]):
                G.add_edge(j + ls, k + ls + nodes[i])

        Ne += nodes[i] * nodes[i + 1]
        ls += nodes[i]


def getedgeinfo(edge_list, edge_colors, nodes, a, N, thres):
    ls = 0
    Ne = 0
    for i in range(0, N - 1):
        for j in range(0, nodes[i]):
            for k in range(0, nodes[i + 1]):

                if (a[i][j] < thres or a[i + 1][k] < thres):
                    cshift = 0.0
                else:
                    cshift = (a[i][j] + a[i + 1][k]) / 2.0
                    edge_list.append([j + ls, k + ls + nodes[i]])
                    edge_colors.append([1.0, 1.0 - cshift, 1.0 - cshift])

        Ne += nodes[i] * nodes[i + 1]
        ls += nodes[i]

fig = plt.gcf()

Writer = animation.writers['ffmpeg']
writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=3600)

font = {'family': 'Bitstream Vera Sans',
        'weight': 'normal',
        'size': 14}

plt.rc('font', **font)

# *********PARAMETERS************

N = 5
Frames = 198

nodes = [56, 8, 4, 4, 1]

user = 'jujuman'
#dir = '/Research/ANN-Test-Data/GDB-11/train2/'
dir = '/PycharmProjects/H2Graphing/'

file = 'netfileANN-H.nnf.dat'

datafile = '/home/' + user + dir + file

movieout = 'Hnet.mp4'
DPI = 100
h_in_inches = 14
w_in_inches = 14

# *******************************

fig.set_size_inches(w_in_inches, h_in_inches)

G = nx.Graph()

makeedges(G,nodes,N)
pos = makenodes(G,nodes,N)

def update(n):
    print('Frame: ' + str(n))

    a = dict()
    for i in range(0, N):
        a[i] = getfltsfromfile(datafile, i, n+2, nodes[i])

    b = dict()
    for i in range(0, N):
        b[i] = getfltsfromfile(datafile, i, n, nodes[i])

    for i in range(0, N):
        b[i] = normalize(b[i]-a[i])

    edge_colors = []
    edge_list = []
    node_colors = []

    getedgeinfo(edge_list, edge_colors, nodes, b, N,0.05)
    getnodesinfo(node_colors, nodes, b, N)

    nx.draw(G, pos, node_color=node_colors, edge_color=edge_colors, edgelist=edge_list)
    #nx.draw_networkx_nodes(G, pos, node_color=node_colors)
    #nx.draw_networkx_edges(G, pos, edge_color=edge_colors, edgelist=edge_list)

#nx.draw(G)

anim = animation.FuncAnimation(fig, update, frames=Frames, interval=60, blit=False, repeat=True)
anim.save(movieout,dpi = DPI)
# plt.xlim([-5,9])
plt.show()
