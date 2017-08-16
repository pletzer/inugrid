import igCubedSphere
import vtk
import math
import numpy

EPS = 1.e-14

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


def integralDPhiDLambdaOverCosTheta(lama, lamb, thea, theb):
    dthe = theb - thea
    dlam = getMinDLambda(lama, lamb)
    if abs(dlam) > EPS:
        return dthe*(math.cos(lamb) - math.cos(lama))/dlam
    else:
        # lama == lamb:
        return -dthe*math.sin(lama)


def integralDPhiDThetaCosTheta(lama, lamb, thea, theb):
    dthe = theb - thea
    dlam = getMinDLambda(lama, lamb)
    la, lb = lama, lamb
    ta, tb = thea, theb
    return 0.5*dlam* (-2*dthe*math.cos(la)*math.cos(2*ta) + 2*dthe*math.cos(lb)*math.cos(2*tb) - \
                      2*dlam*math.cos(ta)*math.sin(la)*math.sin(ta) + dlam*math.sin(lb)*math.sin(2*tb))/ \
                        (dlam**2 - 4*dthe**2)

def phi(lam, the):
    return math.cos(lam) * math.sin(the)


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

        lama, thea = getLambdaTheta(x0, y0, z0)
        lamb, theb = getLambdaTheta(x1, y1, z1)

        divVal += integralDPhiDLambdaOverCosTheta(lama, lamb, thea, theb)
        divVal += integralDPhiDThetaCosTheta(lama, lamb, thea, theb)

    divData[cellId] = divVal

# attach cell centred values to the grid
dataArray = vtk.vtkDoubleArray()
dataArray.SetNumberOfComponents(1)
dataArray.SetNumberOfTuples(numCells)
save = 1
dataArray.SetVoidArray(divData, numCells, save)

grid.GetCellData().SetScalars(dataArray)

# save/show
cs.save('div2.vtk')
cs.show()
