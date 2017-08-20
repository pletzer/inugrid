import igLatLon
import vtk
import math
import numpy

EPS = 1.e-14

alpha = 0.0

def getLambdaTheta(x, y, z):
    # get lon/lat in radiants
    lam = math.atan2(y, x)
    rho = math.sqrt(x*x + y*y)
    the = math.atan2(z, rho)
    return lam, the


def getMinDLambda(lam0, lam1):
    # handle day line issue
    dlam = lam1 - lam0
    a = [abs(dlam + i*2*math.pi) for i in (-1, 0, 1)]
    index = numpy.argmin(a)
    return dlam + (index - 1)*2*math.pi

def getCosThetaDLambda(lam0, the0, lam1, the1):
    dLam = getMinDLambda(lam0, lam1)
    dThe = the1 - the0
    if abs(dThe < EPS):
        # same theta
        return dLam * math.cos(0.5*(the0 + the1))
    else:
        # normal case
        return dLam * (math.sin(the1) - math.sin(the0)) / dThe


def getSinTwoThetaDLambda(lam0, the0, lam1, the1):
    dLam = getMinDLambda(lam0, lam1)
    dThe = the1 - the0
    if abs(dThe < EPS):
        # same theta
        return dLam * math.sin(the0 + the1)
    else:
        # normal case
        return - 0.5 * dLam * (math.cos(2*the1) - math.cos(2*the0)) / dThe



nlats, nlons = 11, 21
cs = igLatLon.LatLon(nlats, nlons)
grid = cs.getUnstructuredGrid()

numCells = grid.GetNumberOfCells()
divData = numpy.zeros((numCells,), numpy.float64)

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

        divVal += getCosThetaDLambda(lam0, the0, lam1, the1)
        divVal -= 0.5 * alpha * getSinTwoThetaDLambda(lam0, the0, lam1, the1)


    divData[cellId] = divVal

# attach cell centred values to the grid
dataArray = vtk.vtkDoubleArray()
dataArray.SetNumberOfComponents(1)
dataArray.SetNumberOfTuples(numCells)
save = 1
dataArray.SetVoidArray(divData, numCells, save)

grid.GetCellData().SetScalars(dataArray)

# save/show
cs.save('divLatLon1.vtk')
cs.show()
