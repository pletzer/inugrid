import vtk
import numpy
import math


points = vtk.vtkPoints()
data = vtk.vtkDoubleArray()
grid = vtk.vtkUnstructuredGrid()
writer = vtk.vtkUnstructuredGridWriter()

# number of cells (radial, poloidal)
nr, nt = 5, 12
rmin, rmax = 1.e-6, 2.
tmin, tmax = 0., 2.*numpy.pi

dr, dt = (rmax - rmin)/float(nr), (tmax - tmin)/float(nt)

def psiFunc(x, y, z):
	theta = math.atan2(y, x)
	return theta

def gradXi1(x, y, z):
	r = math.sqrt(x*x + y*y)
	rHat = numpy.array([x/r, y/r, 0.])
	return rHat/dr

def gradXi2(x, y, z):
	r = math.sqrt(x*x + y*y)
	tHat = numpy.array([-y/r, x/r, 0.])
	return tHat/(r * dt)

nr1, nt1 = nr + 1, nt + 1

# build the points and the field
# 4 points per cell
points.SetNumberOfPoints(nr * nt * 4)
data.SetNumberOfComponents(3)
data.SetNumberOfTuples(nr * nt * 4)

z = 0. # on z = 0 plane (2D)
k = 0
for i in range(nr):
	r0 = rmin + i*dr
	r1 = r0 + dr
	for j in range(nt):
		t0 = tmin + j*dt
		t1 = t0 + dt

		x00, y00 = r0*numpy.cos(t0), r0*numpy.sin(t0)
		x10, y10 = r1*numpy.cos(t0), r1*numpy.sin(t0)
		x11, y11 = r1*numpy.cos(t1), r1*numpy.sin(t1)
		x01, y01 = r0*numpy.cos(t1), r0*numpy.sin(t1)

		# field values on the edges
		vs0 = psiFunc(x10, y10, z) - psiFunc(x00, y00, z)
		vs1 = psiFunc(x11, y11, z) - psiFunc(x01, y01, z)
		v0s = psiFunc(x01, y01, z) - psiFunc(x00, y00, z)
		v1s = psiFunc(x11, y11, z) - psiFunc(x10, y10, z)

		points.InsertPoint(k, (x00, y00, z))
		data.SetTuple(k, vs0*gradXi1(x00, y00, z) + v0s*gradXi2(x00, y00, z))
		k += 1

		points.InsertPoint(k, (x10, y10, z))
		data.SetTuple(k, vs0*gradXi1(x10, y10, z) + v1s*gradXi2(x10, y10, z))
		k += 1

		points.InsertPoint(k, (x11, y11, z))
		data.SetTuple(k, vs1*gradXi1(x11, y11, z) + v1s*gradXi2(x11, y11, z))
		k += 1

		points.InsertPoint(k, (x01, y01, z))
		data.SetTuple(k, vs1*gradXi1(x01, y01, z) + v0s*gradXi2(x01, y01, z))
		k += 1

# build the grid
numCells = nr * nt
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(4)

grid.Allocate(numCells, 1)

k = 0
for i in range(nr):
	for j in range(nt):
		for el in range(4):
			ptIds.SetId(el, k)
			k += 1
		grid.InsertNextCell(vtk.VTK_QUAD, ptIds)

grid.SetPoints(points)
grid.GetPointData().SetScalars(data)


# save to file
writer.SetInputData(grid)
writer.SetFileName('edgeVectorCellByCell.vtk')
writer.Update()