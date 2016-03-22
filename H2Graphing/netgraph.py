#!/usr/bin/env python
"""
Draw a graph with matplotlib, color edges.
You must have matplotlib>=87.7 for this to work.
"""
__author__ = "Justin Smith"
import matplotlib.pyplot as plt

font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

nodes = [12,6,6,1]

for i in range(0,4):
    for j in range(0,nodes[i]):
        color = '#%02x%02x%02x' % (int((255/nodes[i]) * j), 0, 0)
        plt.scatter(i, j, color=color,linewidth=5)

plt.show()