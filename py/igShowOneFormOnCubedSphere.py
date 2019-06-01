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
M = 2
N = 1


def psi(x, y):
    """
	Stream function 
	@param x longitude 
	@param y latitude
	"""
    return numpy.cos(y)*(numpy.sin(2*x - y))
    #return numpy.cos(y)*(numpy.cos(2*M*x) * numpy.sin(N*(3*y-x))**2)


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


# resolution of cubed-sphere 
n = 20

cs = igCubedSphere.CubedSphere(n, radius=1.01)
grid = cs.getUnstructuredGrid()
numCells = grid.GetNumberOfCells()
points = grid.GetPoints()

# create the edges to value map. Use a dict to avoid duplication 
# of edges when iterating over cells (not the most efficient)
ptIdEdges2Val = {}
ptIds = vtk.vtkIdList()
p0 = numpy.array([0., 0., 0.])
p1 = numpy.array([0., 0., 0.])
for cellId in range(numCells):
    grid.GetCellPoints(cellId, ptIds)
    numPoints = ptIds.GetNumberOfIds()
    for i0 in range(numPoints):
        i1 = (i0 + 1) % numPoints # close the loop
        # get the points of the edge
        ptId0 = ptIds.GetId(i0)
        ptId1 = ptIds.GetId(i1)
        points.GetPoint(ptId0, p0)
        points.GetPoint(ptId1, p1)
        # determine orientation of the edge (i0, i1)
        # north is positive
        # if lat diff is zero then east is positive
        # might have to be careful about dateline (let's not worry at this point)
        lam0, the0 = getLambdaTheta(p0)
        lam1, the1 = getLambdaTheta(p1)
        e01 = (ptId0, ptId1)
        orientation = 1
        dTheta, dLambda = the1 - the0, lam1 - lam0
        if dTheta < 0 or (abs(dTheta) < EPS and dLambda < 0):
            # swap the point indices
            e01 = (ptId1, ptId0)
            orientation = -1
        dPsi = psi(lam1, the1) - psi(lam0, the0)
        # store
        ptIdEdges2Val[e01] = orientation * dPsi

numEdges = len(ptIdEdges2Val)

actors = []

# create the edges
edges = vtk.vtkCellArray()
edges.Allocate(numEdges)
values = vtk.vtkDoubleArray()
values.SetNumberOfComponents(1)
values.SetNumberOfTuples(numEdges)
values.SetName('1-form')
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(2) # probably not needed
index = 0
for e01, val in ptIdEdges2Val.items():
    ptIds.SetId(0, e01[0])
    ptIds.SetId(1, e01[1])
    edges.InsertNextCell(ptIds)
    values.SetTuple(index, [val,])
    index += 1

poly = vtk.vtkPolyData() # may need to allocate?
poly.SetPoints(points)
poly.SetLines(edges)
#poly.GetCellData().AddArray(values)

# write to file
writer = vtk.vtkPolyDataWriter()
writer.SetFileName('edges.vtk')
writer.SetInputData(poly)
writer.Update()

lut = vtk.vtkLookupTable()
dmin, dmax = values.GetRange()
lut.SetHueRange(0.667, 0.)
lut.SetRange(dmin, dmax)
lut.Build()

# show
tube = vtk.vtkTubeFilter()
tube.SetInputData(poly)
tube.SetRadius(0.01)
tube.SetNumberOfSides(5)
tube.Update()
tubeMapper = vtk.vtkPolyDataMapper()
tubeMapper.SetInputConnection(tube.GetOutputPort())
tubeMapper.SetLookupTable(lut)
tubeMapper.SetUseLookupTableScalarRange(1)
tubeMapper.ScalarVisibilityOn()
tubeActor = vtk.vtkActor()
tubeActor.SetMapper(tubeMapper)
actors.append(tubeActor)

mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(tube.GetOutputPort())
mapper.SetInputData(poly)
mapper.SetLookupTable(lut)
poly.GetCellData().SetActiveScalars('1-form')
mapper.SetUseLookupTableScalarRange(1)
actor = vtk.vtkActor()
actor.SetMapper(mapper)

actors.append(actor)

# add sphere
sphereSource = vtk.vtkSphereSource()
sphereSource.SetRadius(1.0)
sphereSource.SetThetaResolution(128)
sphereSource.SetPhiResolution(64)
sphereMapper = vtk.vtkPolyDataMapper()
sphereMapper.SetInputConnection(sphereSource.GetOutputPort())
sphereActor = vtk.vtkActor()
sphereActor.SetMapper(sphereMapper)
sphereActor.GetProperty().SetColor(0., 0., 0.)
actors.append(sphereActor)

# add color bar
colorbar = vtk.vtkScalarBarActor()
colorbar.SetLookupTable(lut)
colorbar.SetTitle(values.GetName())
colorbar.GetLabelTextProperty().SetFontSize(14)
#actors.append(colorbar)

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
iren = vtk.vtkRenderWindowInteractor()
renWin.AddRenderer(ren)
iren.SetRenderWindow(renWin)
for a in actors:
    ren.AddActor(a)

camera = vtk.vtkCamera()
camera.SetPosition(0., -4., 0.)
camera.SetFocalPoint(0., 0., 0.)
#ren.SetActiveCamera(camera)

ren.SetBackground(1., 1., 1.)
renWin.SetSize(1200, 1200)
iren.Initialize()
renWin.Render()
iren.Start()



