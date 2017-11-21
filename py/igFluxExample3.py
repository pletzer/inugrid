from igFluxCalculator import FluxCalculator
from igCubedSphere import CubedSphere
import numpy

"""
Example showing how to compute the flux across stream function grad lambda, 
ie the 2-form is d lambda ^ dr.
"""

angle = 0.0 # 0.3 #numpy.pi/6.
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

n = 20
cs = CubedSphere(n)
cs.save('cs.vtk')
grid = cs.getUnstructuredGrid()

fc = FluxCalculator(grid, integralFunction)

# number of contour segments
numSegments = 100

# start longitude 	
lamMin, lamMax = 0.0, 0.3*2*numpy.pi

# increment in lambda
dx = (lamMax - lamMin) / numSegments

# fixed latitude, close to the north pole
the = 1.0 # 0.9* numpy.pi/2.

# expression for the contour in lon-lat coordinates
line = numpy.array([(lamMin + i*dx, the) for i in range(numSegments + 1)]).reshape(numSegments + 1, 2)
fc.setLine(line)

# compute the flux
totFlux = fc.computeFlux()

print('Total flux: {}'.format(totFlux))
exactFlux  = psi(line[-1][0], line[-1][1]) - psi(line[0][0], line[0][1])
print('Exact flux: {}'.format(exactFlux))
#if hasattr(fc, 'polyline2Flux'): print('line 2 flux map: {}'.format(fc.polyline2Flux))

