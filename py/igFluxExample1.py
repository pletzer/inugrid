from igFluxCalculator import FluxCalculator
from igCubedSphere import CubedSphere
import numpy

n = 5
cs = CubedSphere(n)
grid = cs.getUnstructuredGrid()

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

fc = FluxCalculator(grid, integralFunction)
numSegments = 10
dx = 2 * numpy.pi / numSegments
the = 0.9* numpy.pi/2.
line = numpy.array([(i*dx, the) for i in range(numSegments + 1)]).reshape(numSegments + 1, 2)
fc.setLine(line)
totFlux = fc.computeFlux()

print('Total flux: {}'.format(totFlux))

