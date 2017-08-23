import igCubedSphere
import vtk
import numpy
from igDivFilter import DivFilter

"""
Compute cosed line integral of d * d phi with 
phi = (1 - alpha sin(lambda)) * cos(theta)
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
cs = igCubedSphere.CubedSphere(n)
grid = cs.getUnstructuredGrid()

fltr = DivFilter(grid)
fltr.applyIntegrals(getIntegral)

# save/show
cs.save('divCubedSphere2.vtk')
cs.show()
