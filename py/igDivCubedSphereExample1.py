import igCubedSphere
import vtk
from math import sqrt, cos, sin, pi, atan2
import numpy

"""
Compute cosed line integral of d * d phi with 
phi = (1 - alpha lambda) * cos(theta)
"""


EPS = 1.e-14

alpha = 0.0


def getLambdaTheta(x, y, z):
    # get lon/lat in radiants
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
        return -(alpha*(-ya + yb)) - \
          (dx*(-2*(-1 + alpha*xa)*(ya - yb)*cos(2*ya) + \
            2*(-1 + alpha*xb)*(ya - yb)*cos(2*yb) + \
            alpha*(-dx)*(sin(2*ya) - sin(2*yb))))/(8.*(ya - yb)**2)
    else:
        return -(dx*(-(cos(ya)*sin(ya)) + \
                (alpha*xa*cos(ya)*sin(ya))/2. + \
                (alpha*xb*cos(ya)*sin(ya))/2.))

n = 10
cs = igCubedSphere.CubedSphere(n)
grid = cs.getUnstructuredGrid()

numCells = grid.GetNumberOfCells()
divData = numpy.zeros((numCells,), numpy.float64)

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
    for i0 in range(numPoints):

        i1 = (i0 + 1) % numPoints

        ptId0, ptId1 = ptIds.GetId(i0), ptIds.GetId(i1)

        x0, y0, z0 = points.GetPoint(ptId0)
        x1, y1, z1 = points.GetPoint(ptId1)

        # retreat by a tiny bit in order to capture multivalued jumps 
        #x1 = x0 + (x1 - x0)*(1. - EPS)
        #y1 = y0 + (y1 - y0)*(1. - EPS)
        #z1 = z0 + (z1 - z0)*(1. - EPS)

        lam0, the0 = getLambdaTheta(x0, y0, z0)
        lam1, the1 = getLambdaTheta(x1, y1, z1)


        divVal += getIntegral(lam0, lam1, the0, the1)

    cellArea = cellAreas.GetComponent(cellId, 0)
    divData[cellId] = divVal / cellArea

# attach cell centred values to the grid
dataArray = vtk.vtkDoubleArray()
dataArray.SetName('integral_star_d_phi')
dataArray.SetNumberOfComponents(1)
dataArray.SetNumberOfTuples(numCells)
save = 1
dataArray.SetVoidArray(divData, numCells, save)

grid.GetCellData().SetScalars(dataArray)

# save/show
cs.save('divCubedSphere1.vtk')
cs.show()
