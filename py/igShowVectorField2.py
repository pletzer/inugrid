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
M = 2
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

n = 20

actors = []

cs = igCubedSphere.CubedSphere(n, radius=1.01)
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


# show vector field 
grid.GetCellData().SetActiveVectors("2-form")
arrowSource = vtk.vtkArrowSource()
arrowSource.SetShaftRadius(0.03)
#arrowSource.SetHeight(0.3)
#arrowSource.SetResolution(8)

cellCenters = vtk.vtkCellCenters()
cellCenters.SetVertexCells(1)
cellCenters.SetInputData(grid)

glyph = vtk.vtkGlyph3D()
glyph.SetVectorModeToUseVector()
glyph.SetScaleModeToScaleByVector()
glyph.SetColorModeToColorByVector()
glyph.SetSourceConnection(arrowSource.GetOutputPort())
glyph.SetInputConnection(cellCenters.GetOutputPort())
glyph.SetScaleFactor(0.5)
glyph.Update()

glyphMapper = vtk.vtkPolyDataMapper()
glyphMapper.SetInputConnection(glyph.GetOutputPort())
glyphMapper.SetUseLookupTableScalarRange(1)


lut = glyphMapper.GetLookupTable()
lut.SetHueRange(0.667, 0.)
#lut.SetValueRange(0., 0.5)
lut.Build()

glyphActor = vtk.vtkActor()
glyphActor.SetMapper(glyphMapper)
#glyphActor.GetProperty().SetColor(1.0, 0.5, 0.5) #(218./255., 165./255., 32./255.)
actors.append(glyphActor)

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
iren = vtk.vtkRenderWindowInteractor()
renWin.AddRenderer(ren)
iren.SetRenderWindow(renWin)
for a in actors:
    ren.AddActor(a)

camera = vtk.vtkCamera()
camera.SetPosition(0., -4., 1.)
camera.SetFocalPoint(0., 0., 0.)
ren.SetActiveCamera(camera)

ren.SetBackground(1, 1, 1)
renWin.SetSize(600, 600)
iren.Initialize()
renWin.Render()
iren.Start()

