import vtk
import numpy
import math


points = vtk.vtkPoints()
vecDataOffset = vtk.vtkDoubleArray()
vecDataExact = vtk.vtkDoubleArray()
psiData = vtk.vtkDoubleArray()
grid = vtk.vtkUnstructuredGrid()
writer = vtk.vtkUnstructuredGridWriter()

# number of cells (radial, poloidal)
nr, nt = 10, 32
rmin, rmax = 0.01, 2.
tmin, tmax = 0., 2.*numpy.pi
amplitude = 0.3

dr, dt = (rmax - rmin)/float(nr), (tmax - tmin)/float(nt)

def psiFunc(r, t):
	return t + amplitude*r*math.cos(t)

def dPsiDr(r, t):
	return amplitude*math.cos(t)

def dPsiDtOverR(r, t):
	return 1.0/r - amplitude*math.sin(t)

def gradXi1(x, y, z):
	return rHat(x, y, z)/dr

def gradXi2(x, y, z):
	r = math.sqrt(x*x + y*y)
	return tHat(x, y, z)/(r * dt)

def rHat(x, y, z):
	r = math.sqrt(x*x + y*y)
	return numpy.array([x/r, y/r, 0.])

def tHat(x, y, z):
	r = math.sqrt(x*x + y*y)
	return numpy.array([-y/r, x/r, 0.])


nr1, nt1 = nr + 1, nt + 1

points.SetNumberOfPoints(nr1 * nt1)

vecDataOffset.SetNumberOfComponents(3)
vecDataOffset.SetNumberOfTuples(nr1 * nt1)
vecDataOffset.SetName('velocityOffset')

vecDataExact.SetNumberOfComponents(3)
vecDataExact.SetNumberOfTuples(nr1 * nt1)
vecDataExact.SetName('velocityExact')

psiData.SetNumberOfComponents(1)
psiData.SetNumberOfTuples(nr1 * nt1)
psiData.SetName('potential')

z = 0. # on z = 0 plane (2D)
k = 0
for i in range(nr1):
	r0 = rmin + i*dr
	r1 = r0 + dr
	rMid = r0 + 0.5*dr

	for j in range(nt1):
		t0 = tmin + j*dt
		t1 = t0 + dt
		tMid = t0 + 0.5*dt

		x00, y00 = r0*numpy.cos(t0), r0*numpy.sin(t0)
		x10, y10 = r1*numpy.cos(t0), r1*numpy.sin(t0)
		x11, y11 = r1*numpy.cos(t1), r1*numpy.sin(t1)
		x01, y01 = r0*numpy.cos(t1), r0*numpy.sin(t1)

		psi = psiFunc(r0, t0)

		# exact velocity at the cell corners
		vExact = dPsiDr(r0, t0)*rHat(x00, y00, z) + dPsiDtOverR(r0, t0)*tHat(x00, y00, z)

		xs0 = 0.5*(x00 + x10)
		ys0 = 0.5*(y00 + y10)
		x0s = 0.5*(x00 + x01)
		y0s = 0.5*(y00 + y01)
		vOffset = dPsiDr(rMid, t0)*rHat(xs0, ys0, z) + \
		          dPsiDtOverR(r0, tMid)*tHat(x0s, y0s, z)

		points.InsertPoint(k, (x00, y00, z))
		vecDataOffset.SetTuple(k, vOffset)
		vecDataExact.SetTuple(k, vExact)
		psiData.SetTuple(k, [psi])
		k += 1


# build the grid
numCells = nr * nt
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(4)

grid.Allocate(numCells, 1)

k = 0
for i in range(nr):
	for j in range(nt):

		k = nt1*i + j
		ptIds.SetId(0, k)
		k += nt1
		ptIds.SetId(1, k)
		k += 1
		ptIds.SetId(2, k)
		k -= nt1
		ptIds.SetId(3, k)

		grid.InsertNextCell(vtk.VTK_QUAD, ptIds)

grid.SetPoints(points)
grid.GetPointData().AddArray(vecDataOffset)
grid.GetPointData().AddArray(vecDataExact)
grid.GetPointData().AddArray(psiData)


# save to file
writer.SetInputData(grid)
writer.SetFileName('offsetVector.vtk')
writer.Update()