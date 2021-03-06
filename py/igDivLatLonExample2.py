import igLatLon
import vtk
from math import sqrt, cos, sin, pi, atan2
import numpy
from igDivFilter import DivFilter

"""
Compute closed line integral of d * d phi with 
phi = (1 - alpha sin(lambda)) * cos(theta)

* d lambda = d theta ^ dr / cos(theta)
* d theta = - cos(theta) d lambda ^ dr
"""

EPS = 1.e-14

alpha = 1.0

def getMinDLambda(lam0, lam1):
    # handle day line issue
    dlam = lam1 - lam0
    a = [abs(dlam + i*2*pi) for i in (-1, 0, 1)]
    index = numpy.argmin(a)
    return dlam + (index - 1)*2*pi


def getIntegral(xa, xb, ya, yb):
    """
    Compute the value attached to an edge
    x is longitude
    y is latitude
    """
    dy = yb - ya
    dx = getMinDLambda(xa, xb)
    if abs(dy) > EPS:
        if abs(dx) > EPS:
            res = -((alpha*dy*(-sin(xa) + sin(xb)))/dx) - \
                  (dx*((cos(2*ya) - cos(2*yb))/(-dy) + \
                  (alpha*(sin(xa - 2*ya) - sin(xb - 2*yb)))/(-dx + 2*dy) + \
                  (alpha*(-sin(xa + 2*ya) + sin(xb + 2*yb)))/(-dx - 2*dy)))/4.
        else:
            # xa == xb
            res = -(alpha*dy*cos(xa))
    else:
        # ya == yb
        res = -(((dx)*(dx - alpha*cos(xa) + alpha*cos(xb))*cos(ya)*sin(ya))/(-dx))
    return res

n = 20
# equivalent number of cells for lat-lon grid
ntot = 20**2 * 6
nlat = int(sqrt(ntot/2.))
nlon = ntot // nlat
cs = igLatLon.LatLon(nlat, nlon)
grid = cs.getUnstructuredGrid()

fltr = DivFilter(grid)
fltr.applyIntegrals(getIntegral)

# save/show
cs.save('divLatLon2.vtk')
cs.show()

