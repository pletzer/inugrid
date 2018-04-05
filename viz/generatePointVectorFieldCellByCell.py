import vtk
import numpy


points = vtk.vtkPoints()
data = vtk.vtkDoubleArray()
grid = vtk.vtkUnstructuredGrid()
writer = vtk.vtkUnstructuredGridWriter()

nx, ny = 3, 4
xmin, xmax = 0., 2.
ymin, ymax = 0., 1.

dx, dy = (xmax - xmin)/float(nx), (ymax - ymin)/float(ny)

def vecFun(x, y, z):
	return numpy.array([x, 0., 0.])

nx1, ny1 = nx + 1, ny + 1

# build the points
# 4 points per cell
points.SetNumberOfPoints(nx * ny * 4)
k = 0
z = 0.

for i in range(nx):
	x = xmin + i*dx
	for j in range(ny):
		y = ymin + j*dy
		points.InsertPoint(k, (x, y, z))
		k += 1
		points.InsertPoint(k, (x + dx, y, z))
		k += 1
		points.InsertPoint(k, (x + dx, y + dy, z))
		k += 1
		points.InsertPoint(k, (x, y + dy, z))
		k += 1

# build the data
data.SetNumberOfComponents(3)
data.SetNumberOfTuples(nx * ny * 4)
k = 0
for i in range(nx):
	for j in range(ny):
		for el in range(4):
			x, y, z = points.GetPoint(k)
			data.SetTuple(k, vecFun(x, y, z))
			k += 1

# build the grid
numCells = nx * ny
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(4)

grid.Allocate(numCells, 1)

k = 0
for i in range(nx):
	for j in range(ny):
		for el in range(4):
			ptIds.SetId(el, k)
			k += 1
		grid.InsertNextCell(vtk.VTK_QUAD, ptIds)

grid.SetPoints(points)
grid.GetPointData().SetScalars(data)


# save to file
writer.SetInputData(grid)
writer.SetFileName('pointVectorCellByCell.vtk')
writer.Update()