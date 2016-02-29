__author__ = 'jujuman'

import pylab
import numpy as np

b1 = 0.0
b2 = 1.0
b3 = -1.0
c1 = 0.01
c2 = 1.0
c3 = 2.0

x = np.linspace(-5,5,1000) # 100 linearly spaced numbers
y1 = np.exp( -0.5*((1/c1)*(x-b1))**2 ) # computing the values of sin(x)/x
y2 = np.exp( -0.5*((1/c2)*(x-b2))**2 ) # computing the values of sin(x)/x
y3 = np.exp( -0.5*((1/c3)*(x-b3))**2 ) # computing the values of sin(x)/x

# compose plot
pylab.plot(x,y1) # sin(x)/x
pylab.plot(x,y2) # sin(x)/x
pylab.plot(x,y3) # sin(x)/x
pylab.show() # show the plot