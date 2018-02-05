import igCubedSphere
import vtk
import numpy
from igDivFilter import DivFilter
from igLandOcean import LandOcean
from igPiecewiseLinearLine import PiecewiseLinearLine
import math

"""
Compute closed line integral of d * d phi with 
phi = (1 - alpha sin(lambda)) * cos(theta)
"""


EPS = 1.e-14
M = 3
N = 1


def psi(x, y):
    """
	Stream function 
	@param x longitude 
	@param y latitude
	"""
    return numpy.cos(y)*(numpy.cos(2*M*x) * numpy.sin(N*(3*y-x))**2)


def getIntegral(xa, xb, ya, yb):
    """
    @param x is longitude
    @param y is latitude
    """
    return psi(xb, yb) - psi(xa, ya)

def getLambdaTheta(point):
    x, y, z = point
    rho = math.sqrt(x**2 + y**2)
    the = math.atan2(z, rho)
    lam = math.atan2(y, x)
    return lam, the


n = 50

actors = []

cs = igCubedSphere.CubedSphere(n, radius=1.01)
grid = cs.getUnstructuredGrid()
numCells = grid.GetNumberOfCells()
points = grid.GetPoints()

ptIdEdges2Val = {}
ptIds = vtk.vtkIdList()
p0 = numpy.array([0., 0., 0.])
p1 = numpy.array([0., 0., 0.])
for cellId in range(numCells):
    grid.GetCellPoints(cellId, ptIds)
    numPoints = ptIds.GetNumberOfIds()
    for i0 in range(numPoints):
        i1 = (i0 + 1) % numPoints
        ptId0 = ptIds.GetId(i0)
        ptId1 = ptIds.GetId(i1)
        points.GetPoint(ptId0, p0)
        points.GetPoint(ptId1, p1)
        # determine orientation of the edge (i0, i1)
        # north is positive
        # if not then east
        lam0, the0 = getLambdaTheta(p0)
        lam1, the1 = getLambdaTheta(p1)
        e01 = (ptId0, ptId1)
        orientation = 1
        dTheta, dLambda = the1 - the0, lam1 - lam0
        if dTheta < 0 or (Theta == 0. and dLambda < 0):
            e01 = (ptId1, ptId0)
            orientation = -1
        dPsi = psi(lam1, the1) - psi(lam0, the0)
        # store
        ptIdEdges2Val[e01] = orientation * dPsi

# create the edges
edges = vtk.vtkCellArray()
numEdges = len(ptIdEdges2Val)
values = vtk.vtkDoubleArray()
values.SetNumberOfComponents(1)
values.SetNumberOfTuples(numEdges)
edges.Allocate(numEdges)
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(2)
index = 0
for e01, val in ptIdEdges2Val.items():
    ptIds.SetId(0, e01[0])
    ptIds.SetId(1, e01[1])
    edges.InsertNextCell(ptIds)
    values.SetTuple(index, val)
    index += 1

poly = vtk.vtkPolyData() # may need to allocate?
poly.SetPoints(points)
poly.SetLines(edges)
poly.Build()

# write to file
writer = vtk.vtkPolyDataWriter()
writer.SetFileName('edges.vtk')
writer.SetInputData(poly)
writer.Update()


