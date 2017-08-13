import igCubedSphere
import vtk
import math
import numpy

EPS = 1.e-14

def getMinDLambda(lam0, lam1):
    # handle day line issue
    dlam = lam1 - lam0
    a = [abs(dlam + i*2*math.pi) for i in (-1, 0, 1)]
    index = numpy.argmin(a)
    return dlam + (index[0] - 1)*2*math.pi

def getLambdaTheta(x, y, z):
    # get lon/lat in radiants
    lam = math.atan2(y, x)
    rho = math.sqrt(x*x + y*y)
    the = math.atan2(z, rho)
    return lam, the    

def psi(x, y, z):
    lam, the = getLambdaTheta(x, y, z)
    # ~ d theta/ 2 pi
    return -lam/(2.*math.pi)

def getRLambdaIntegral(x0, y0, z0, x1, y1, z1):
    # integral cos(the) dr ^ dlam
    lam0, the0 = getLambdaTheta(x0, y0, z0)
    lam1, the1 = getLambdaTheta(x1, y1, z1)
    dlam = getMinDLambda(lam0, lam1)
    dthe = the1 - the0
    if abs(dthe) < EPS:
        the = 0.5*(the0 + the0)
        return dlam * math.cos(the)
    else:
        return dlam * (math.sin(the1) - math.sin(the0)) / dthe

def getThetaRIntegral(x0, y0, z0, x1, y1, z1):
    # integral sin(lam) dtheta ^ dr
    lam0, the0 = getLambdaTheta(x0, y0, z0)
    lam1, the1 = getLambdaTheta(x1, y1, z1)
    dlam = getMinDLambda(lam0, lam1)
    dthe = the1 - the0
    if abs(dlam) < EPS:
        lam = 0.5*(lam0 + lam1)
        return dthe * math.sin(lam)
    else:
        return - dthe * (math.cos(lam1) - math.cos(lam0)) / dlam

n = 21
cs = igCubedSphere.CubedSphere(n)
grid = cs.getUnstructuredGrid()

numCells = grid.GetNumberOfCells()
divData = numpy.zeros((numCells,), numpy.float64)

# iterate over the cells
points = grid.GetPoints()
for cellId in range(numCells):
    cell = grid.GetCell(cellId)
    ptIds = cell.GetPointIds()
    numPoints = ptIds.GetNumberOfIds()
    # compute the exterior derivative (div)
    divVal = 0.0
    for i0 in range(numPoints):
        i1 = (i0 + 1) % numPoints
        ptId0, ptId1 = ptIds.GetId(i0), ptIds.GetId(i1)
        x0, y0, z0 = points.GetPoint(ptId0)
        x1, y1, z1 = points.GetPoint(ptId1)
        # retreat by a tiny bit in order to capture multivalued jumps 
        x1 = x0 + (x1 - x0)*(1. - EPS)
        y1 = y0 + (y1 - y0)*(1. - EPS)
        z1 = z0 + (z1 - z0)*(1. - EPS)
        # compute the stream function for the start/end points
        psi0 = psi(x0, y0, z0)
        psi1 = psi(x1, y1, z1)
        divVal += psi1 - psi0
    divData[cellId] = divVal

# attach cell centred values to the grid
dataArray = vtk.vtkDoubleArray()
dataArray.SetNumberOfComponents(1)
dataArray.SetNumberOfTuples(numCells)
save = 1
dataArray.SetVoidArray(divData, numCells, save)

grid.GetCellData().SetScalars(dataArray)

# save/show
cs.save('div.vtk')
cs.show()
