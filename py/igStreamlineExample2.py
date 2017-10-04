import numpy
from scipy.integrate import odeint
from matplotlib import pylab
from igStreamVectorField import StreamVectorField
from igCartesianGrid import CartesianGrid

"""
Stream function y^2 + sin(x)^2
v = d psi ^ dz
"""

#  create the grid
ns = (31, 21, 1)
ls = (2*numpy.pi, 2.0, 1.0)
origin = (0., -1.0, 0.0)
cart = CartesianGrid(ns, ls, origin)
grid = cart.getUnstructuredGrid()

# create the stream function
def streamFunc(x):
    return numpy.sin(x[0])**2 + x[1]**2

vecField = StreamVectorField(grid)
vecField.setStreamFunction(streamFunc)

x = numpy.array([numpy.pi/2. + 0.01, 0., 0.])
ts = numpy.linspace(0., 1., 11) #numpy.linspace(0., 1000., 10001)
sol = odeint(vecField, x, ts)


pylab.plot(sol[:, 0], sol[:, 1], 'b-')

pylab.show()
