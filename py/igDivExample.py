import igCubedSphere
import vtk
import math
import numpy

def psi(x, y, z):
    lam = math.atan2(y, x)
    rho = math.sqrt(x*x + y*y)
    the = math.atan2(z, rho)
    # ~ d theta/ 2 pi
    return -lam/(2.*math.pi)

n = 21
cs = igCubedSphere.CubedSphere(n)
grid = cs.getUnstructuredGrid()

numCells = grid.GetNumberOfCells()
divData = numpy.zeros((numCells,), numpy.float64)

# iterate over the cells
eps = 1.e-14
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
        x1 = x0 + (x1 - x0)*(1. - eps)
        y1 = y0 + (y1 - y0)*(1. - eps)
        z1 = z0 + (z1 - z0)*(1. - eps)
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
