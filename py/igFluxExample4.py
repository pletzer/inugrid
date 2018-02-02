from igFluxCalculator import FluxCalculator
from igCubedSphere import CubedSphere
from igPiecewiseLinearLine import PiecewiseLinearLine
from igNodalFunctionWriter import NodalFunctionWriter
import numpy
import math

"""
Example showing how to compute the flux across stream function grad lambda, 
ie the 2-form is d lambda ^ dr.
"""

angle = 0.0 # numpy.pi/6.2324325
cos_angle = numpy.cos(angle)
sin_angle = numpy.sin(angle)


def psi(x, y):
    """
    Stream function
    x ~ longitude
    y ~ latitude
    """
    xp = cos_angle*x - sin_angle*y
    yp = sin_angle*x + cos_angle*y
    return numpy.sin(xp)**2 + numpy.sin(yp)**2

def integralFunction(xa, ya, xb, yb):
    """
    Compute flux attached to edge 
    x ~ longitude
    y ~ latitude
    a: starting point
    b: ending point
    """
    return psi(xb, yb) - psi(xa, ya)

def streamFuncExact(xyz, *args):
    x, y, z = xyz
    rho = math.sqrt(x**2 + y**2)
    the = math.atan2(z, rho)
    lam = math.atan2(y, x)
    return psi(lam, the)

n = 20
cs = CubedSphere(n)
cs.save('cs.vtk')
grid = cs.getUnstructuredGrid()
cs.save('cubedSphereGrid.vtk')

fc = FluxCalculator(grid, integralFunction)

# number of contour segments
numSegments = 1 #100

lamCentre, theCentre = numpy.pi/4., numpy.pi/4.
radius = 0.2

def lamFunc(ts):
    return lamCentre + radius*numpy.cos(5*numpy.pi/3.*ts - 0.2)
def theFunc(ts):
    return theCentre + radius*numpy.sin(5*numpy.pi/3.*ts - 0.2)

ts = numpy.linspace(0., 1., numSegments + 1)
lams = lamFunc(ts)
thes = theFunc(ts)

# expression for the contour in lon-lat coordinates
line = numpy.array([(lams[i], thes[i]) for i in range(numSegments + 1)]).reshape(numSegments + 1, 2)
fc.setLine(line)

# for plotting
pline = PiecewiseLinearLine(lamFunc, theFunc)
pline.save('lineCubedSphere4.vtk')

# compute the flux
totFlux = fc.computeFlux()

print('Total flux: {}'.format(totFlux))
exactFlux  = psi(line[-1][0], line[-1][1]) - psi(line[0][0], line[0][1])
print('Exact flux: {} error = {}'.format(exactFlux, totFlux - exactFlux))

# write stream function
dataWriter = NodalFunctionWriter(grid, streamFuncExact, name='psi')
dataWriter.save('psiCubedSphere.vtk')


