from igFluxCalculator import FluxCalculator
from igCubedSphere import CubedSphere
import numpy

"""
Example showing how to compute the flux across stream function grad lambda, 
ie the 2-form is d lambda ^ dr.  The result across the north pole should be 2*pi. 
"""

def psi(x, y):
    """
    Stream function
    """
    return x

def integralFunction(xa, ya, xb, yb):
    """
    Compute flux attached to edge 
    x ~ longitude
    y ~ latitude
    a: starting point
    b: ending point
    """
    return psi(xb, yb) - psi(xa, ya)

n = 20
cs = CubedSphere(n)
grid = cs.getUnstructuredGrid()

fc = FluxCalculator(grid, integralFunction)

# number of contour segments
numSegments = 10

# increment in lambda
dx = 2 * numpy.pi / numSegments

# fixed latitude
the = 0.9* numpy.pi/2.

# expression for the contour in lon-lat coordinates
line = numpy.array([(0. + i*dx, the) for i in range(numSegments + 1)]).reshape(numSegments + 1, 2)
fc.setLine(line)

# compute the flux
totFlux = fc.computeFlux()

print('Total flux: {}'.format(totFlux))
print('line 2 flux map: {}'.format(fc.polyline2Flux))

