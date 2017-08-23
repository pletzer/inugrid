import igLatLon
import vtk
from math import sqrt, cos, sin, pi, atan2
import numpy

"""
Compute cosed line integral of d * d phi with 
phi = (1 - alpha sin(lambda)) * cos(theta)
"""


EPS = 1.e-14

alpha = 1.0


def getLambdaTheta(xx):
    # get lon/lat in radiants
    x, y, z = xx
    lam = atan2(y, x)
    rho = sqrt(x*x + y*y)
    the = atan2(z, rho)
    return lam, the


def getMinDLambda(lam0, lam1):
    # handle day line issue
    dlam = lam1 - lam0
    a = [abs(dlam + i*2*pi) for i in (-1, 0, 1)]
    index = numpy.argmin(a)
    return dlam + (index - 1)*2*pi


def getIntegral(xa, xb, ya, yb):
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

numCells = grid.GetNumberOfCells()
divData = numpy.zeros((numCells,), numpy.float64)
vecData = numpy.zeros((numCells, 3), numpy.float64)

grid.GetCellData().SetActiveScalars('cell_areas')
cellAreas = grid.GetCellData().GetScalars()

# iterate over the cells
points = grid.GetPoints()
for cellId in range(numCells):
    cell = grid.GetCell(cellId)
    ptIds = cell.GetPointIds()
    numPoints = ptIds.GetNumberOfIds()
    #compute the closed loop integral
    divVal = 0.0
    vecVal = numpy.zeros((3,), numpy.float64)
    # iterate over the edges
    for i0 in range(numPoints):

        i1 = (i0 + 1) % numPoints

        ptId0, ptId1 = ptIds.GetId(i0), ptIds.GetId(i1)

        xx0 = numpy.array(points.GetPoint(ptId0))
        xx1 = numpy.array(points.GetPoint(ptId1))
        rr = 0.5*(xx0 + xx1)

        # retreat by a tiny bit in order to capture multivalued jumps 
        #x1 = x0 + (x1 - x0)*(1. - EPS)
        #y1 = y0 + (y1 - y0)*(1. - EPS)
        #z1 = z0 + (z1 - z0)*(1. - EPS)

        lam0, the0 = getLambdaTheta(xx0)
        lam1, the1 = getLambdaTheta(xx1)

        fluxIntegral = getIntegral(lam0, lam1, the0, the1)
        divVal += fluxIntegral

        # compute the vector at cell centre
        length = numpy.sqrt(numpy.dot(xx1 - xx0, xx1 - xx0))
        vecVal += 0.5 * fluxIntegral * numpy.cross(xx1 - xx0, rr)/length

    cellArea = cellAreas.GetComponent(cellId, 0)
    divData[cellId] = divVal / cellArea
    vecData[cellId, :] = vecVal

# attach cell centred values to the grid
divArray = vtk.vtkDoubleArray()
divArray.SetName('integral_star_d_phi_over_area')
divArray.SetNumberOfComponents(1)
divArray.SetNumberOfTuples(numCells)
save = 1
divArray.SetVoidArray(divData, numCells, save)

# attach cell centred vector field
vecArray = vtk.vtkDoubleArray()
vecArray.SetName('2-form')
vecArray.SetNumberOfComponents(3)
vecArray.SetNumberOfTuples(numCells)
vecArray.SetVoidArray(vecData, numCells*3, save)

grid.GetCellData().AddArray(divArray)
grid.GetCellData().AddArray(vecArray)

# save/show
cs.save('divLatLon2.vtk')
cs.show()
