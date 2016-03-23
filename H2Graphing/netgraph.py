#!/usr/bin/env python
"""
Draw a graph with matplotlib, color edges.
You must have matplotlib>=87.7 for this to work.
"""
__author__ = "Justin Smith"
import matplotlib.pyplot as plt
import networkx as nx

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

#*********PARAMETERS************

N = 5
nodes = [32,8,4,4,1]

#*******************************

G=nx.Graph()

ls = 0
Ne = 0
for i in range(0,N-1):
    print('i: '+str(i) + ' ls: ' + str(ls))
    for j in range(ls,nodes[i]+ls):
        for k in range(nodes[i]+ls,nodes[i+1]+nodes[i]+ls):
            print('Edge(i: ' + str(i) + ' ) ' + str(j) + ' ' + str(k))
            G.add_edge(j,k)

    Ne += nodes[i] * nodes[i+1]
    ls += nodes[i]

fixed_positions = dict()
colors = []
idx = 0
for i in range(0,N):
    shift = nodes[i]/2
    for j in range(0,nodes[i]):
        #color = '#%02x%02x%02x' % (int((255/nodes[i]) * j), 0, 0)
        colors.append([j/float(nodes[i]),0.0,0.0,1.0])
        fixed_positions[idx] = (i,j-shift)
        idx += 1

print(G.number_of_edges())
print(Ne)

#fixed_positions = {1:(0,0.5),2:(0,1.5),3:(0,2.5),4:(1,0),5:(1,1),6:(1,2),7:(1,3),8:(2,1.5)}#dict with two of the positions set
print(fixed_positions)
print(colors)
fixed_nodes = fixed_positions.keys()
pos = nx.spring_layout(G,pos=fixed_positions, fixed = fixed_nodes)

nx.draw_networkx_nodes(G,pos,node_color=colors)
nx.draw_networkx_edges(G,pos,node_color=[[0.0,1.0,0.0,1.0]])

        #plt.scatter(i, j, color=color,s=100)

#plt.xlim([-5,9])
plt.show()