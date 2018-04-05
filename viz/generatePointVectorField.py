import vtk
import numpy


points = vtk.vtkPoints()
data = vtk.vtkDoubleArray()
grid = vtk.vtkUnstructuredGrid()
writer = vtk.vtkUnstructuredGridWriter()

nx, ny = 2, 1
xmin, xmax = 0., 2.
ymin, ymax = 0., 1.

def vecFun(x, y, z):
	return numpy.array([x, 0., 0.])

nx1, ny1 = nx + 1, ny + 1

# build the points
points.SetNumberOfPoints(nx1 * ny1)
k = 0
z = 0.
for i in range(nx1):
	x = xmin + i*(xmax - xmin)/float(nx)
	for j in range(ny1):
		y = ymin + j*(ymax - ymin)/float(ny)
		points.InsertPoint(k, (x, y, z))
		k += 1

# build the data
data.SetNumberOfComponents(3)
data.SetNumberOfTuples(nx1 * ny1)
k = 0
for i in range(nx1):
	for j in range(ny1):
		x, y, z = points.GetPoint(k)
		data.SetTuple(k, vecFun(x, y, z))
		k += 1

# build the grid
numCells = nx * ny
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(4)

grid.Allocate(numCells, 1)

for i in range(nx):
	for j in range(ny):
		ptIds.SetId(0, (i+0)*ny1 + (j+0))
		ptIds.SetId(1, (i+1)*ny1 + (j+0))
		ptIds.SetId(2, (i+1)*ny1 + (j+1))
		ptIds.SetId(3, (i+0)*ny1 + (j+1))
		grid.InsertNextCell(vtk.VTK_QUAD, ptIds)

grid.SetPoints(points)
grid.GetPointData().SetScalars(data)


# save to file
writer.SetInputData(grid)
writer.SetFileName('pointVector.vtk')
writer.Update()