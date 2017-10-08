import numpy
import math
from scipy.integrate import odeint
from matplotlib import pylab
from igStreamVectorField import StreamVectorField
from igCubedSphereElv import CubedSphereElv
from igGridGeometry import GridGeometry
from igPolyLineWriter import PolyLineWriter

"""
Stream function y^2 + sin(x)^2
v = d psi ^ dz
"""

#  create the grid
n = 10
radius, maxRelElv = 1.0, 1.0
midRadius = radius * (1.0 + 0.5*maxRelElv)
cart = CubedSphereElv(numCellsPerTile=n, numElvs=1, radius=radius, maxRelElv=1.0)
grid = cart.getUnstructuredGrid()
geom = GridGeometry(grid)

angle = 0.0 # 0.3 #numpy.pi/6.
cos_angle = numpy.cos(angle)
sin_angle = numpy.sin(angle)

# stream function
def streamFuncExact(xyz, *args):
    x, y, z = xyz
    rhoSq = x**2 + y**2
    r = numpy.sqrt(rhoSq + z**2)
    rho = numpy.sqrt(rhoSq)
    the = math.atan2(z, rho)
    lam = math.atan2(y, x)
    # apply rotation to the coordinates
    lamp = cos_angle * lam - sin_angle * the
    thep = sin_angle * lam + cos_angle * the
    return numpy.sin(lamp)**2 + thep**2

def XYZFromLamThe(lam, the):
    rho = midRadius * numpy.cos(the)
    z = midRadius * numpy.sin(the)
    x = rho * numpy.cos(lam)
    y = rho * numpy.sin(lam)
    return x, y, z


velocityFace = StreamVectorField(grid)
velocityFace.setStreamFunction(streamFuncExact)


# time steps
ts = numpy.linspace(0., 5.0, 101) #10.0, 101)

# vary initial conditions
theStart = 0.1
index = 0
for lamStart in numpy.linspace(1., 1.4, 3):

    # initial condition
    x0, y0, z0 = XYZFromLamThe(lamStart, theStart)
    x = numpy.array([x0, y0, z0]) 

    # solve
    sol = odeint(velocityFace, x, ts)

    writer = PolyLineWriter(sol)
    writer.save('trajectory{}.vtk'.format(index))

    index += 1



