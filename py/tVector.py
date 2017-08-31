import numpy
import vtk

nx, ny, nz = 4, 3, 2
nx1, ny1, nz1 = nx + 1, ny + 1, nz + 1
numPts = nx1 * ny1 * nz1

vecData = vtk.vtkDoubleArray()
vecData.SetNumberOfComponents(3)
vecData.SetNumberOfTuples(numPts)

ptData = vtk.vtkDoubleArray()
ptData.SetNumberOfComponents(3)
ptData.SetNumberOfTuples(numPts)

xmin, xmax = 0.0, 1.0
ymin, ymax = 0.0, 1.5
zmin, zmax = 0.0, 2.0
dx = (xmax - xmin)/float(nx)
dy = (ymax - ymin)/float(ny)
dz = (zmax - zmin)/float(nz)

index = 0
for k in range(nz1):
	for j in range(ny1):
		for i in range(nx1):
			x = xmin + i*dx
			y = ymin + j*dy
			z = zmin + k*dz
			ptData.SetTuple(index, (x, y, z))
			vx = -y
			vy = +x
			vz = 0.0
			print vx, vy, vz
			vecData.SetTuple(index, (vx, vy, vz))
			index += 1

pts = vtk.vtkPoints()
pts.SetData(ptData)

grid = vtk.vtkStructuredGrid()
grid.SetDimensions(nx1, ny1, nz1)
grid.SetPoints(pts)
grid.GetPointData().SetVectors(vecData)

arrow = vtk.vtkArrowSource()

glyph = vtk.vtkGlyph3D()
glyph.SetVectorModeToUseVector()
glyph.SetScaleModeToScaleByVector()
glyph.SetSourceConnection(arrow.GetOutputPort())
glyph.SetInputData(grid)

"""
glyph.DebugOn()
glyph.SetVectorMode(1)
glyph.OrientOn()
glyph.SetColorModeToColorByVector()
glyph.SetScaleFactor(0.2)
glyph.Update()
print glyph
"""


mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(glyph.GetOutputPort())
mapper.Update()
print mapper

actor = vtk.vtkActor()
actor.SetMapper(mapper)

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
iren = vtk.vtkRenderWindowInteractor()
renWin.AddRenderer(ren)
iren.SetRenderWindow(renWin)
ren.AddActor(actor)
ren.SetBackground(0.5, 0.5, 0.5)
renWin.SetSize(900, 600)
iren.Initialize()
renWin.Render()
iren.Start()
