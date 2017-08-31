import igCubedSphere
import vtk
import numpy
from igDivFilter import DivFilter
from igLandOcean import LandOcean
from igPiecewiseLinearLine import PiecewiseLinearLine

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

n = 10

actors = []

cs = igCubedSphere.CubedSphere(n)
grid = cs.getUnstructuredGrid()

fltr = DivFilter(grid)
fltr.applyIntegrals(getIntegral)

# save/show
grid.GetCellData().RemoveArray('cell_areas')
grid.GetCellData().RemoveArray('integral_star_d_phi_over_area')
grid.GetCellData().RemoveArray('1-form')
grid.GetPointData().RemoveArray('2-form-approx')
##cs.save('t.vtk')
#cs.show()

# add a nive background
pnt = LandOcean(textureFile="2k_earth_daymap.jpeg", radius=0.99)
actors = pnt.actors

# show vector field 
#grid.GetCellData()
grid.GetCellData().SetActiveVectors("2-form")
#vecData = grid.GetCellData().GetVectors("2-form")
#print vecData
#for i in range(vecData.GetNumberOfTuples()):
#    print i, vecData.GetTuple(i)
arrowSource = vtk.vtkArrowSource()

cellCenters = vtk.vtkCellCenters()
cellCenters.SetVertexCells(1)
cellCenters.SetInputData(grid)

glyph = vtk.vtkGlyph3D()
glyph.SetVectorModeToUseVector()
glyph.SetScaleModeToScaleByVector()
glyph.SetSourceConnection(arrowSource.GetOutputPort())
glyph.SetInputConnection(cellCenters.GetOutputPort())
glyph.SetScaleFactor(0.1)
glyph.Update()

glyphMapper = vtk.vtkPolyDataMapper()
glyphMapper.SetInputConnection(glyph.GetOutputPort())

glyphActor = vtk.vtkActor()
glyphActor.SetMapper(glyphMapper)
actors.append(glyphActor)

# line
nt = 6
def lamFunc(ts):
	return -numpy.pi + ts*2*numpy.pi

def latFunc(ts):
	return -60.*(numpy.pi/180.)*numpy.ones((nt,), numpy.float64)


line = PiecewiseLinearLine(lamFunc, latFunc, nt=nt, radius=1.05)

tubes = vtk.vtkTubeFilter()
tubes.SetRadius(0.05)
tubes.SetInputData(line.poly)

lineMapper = vtk.vtkPolyDataMapper()
lineMapper.SetInputConnection(tubes.GetOutputPort())
lineMapper.Update()
lineActor = vtk.vtkActor()
lineActor.SetMapper(lineMapper)
lineActor.GetProperty().SetColor(1, 0, 0)
actors.append(lineActor)

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
iren = vtk.vtkRenderWindowInteractor()
renWin.AddRenderer(ren)
iren.SetRenderWindow(renWin)
for a in actors:
    ren.AddActor(a)

ren.SetBackground(0.5, 0.5, 0.5)
renWin.SetSize(900, 600)
iren.Initialize()
renWin.Render()
iren.Start()

