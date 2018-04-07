import vtk
import numpy
import math


points = vtk.vtkPoints()
vecData = vtk.vtkDoubleArray()
vecDataExact = vtk.vtkDoubleArray()
vecDataAvgCell = vtk.vtkDoubleArray()
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

# build the points and the field
# 4 points per cell
points.SetNumberOfPoints(nr * nt * 4)

vecData.SetNumberOfComponents(3)
vecData.SetNumberOfTuples(nr * nt * 4)
vecData.SetName('velocityEdge')

vecDataExact.SetNumberOfComponents(3)
vecDataExact.SetNumberOfTuples(nr * nt * 4)
vecDataExact.SetName('velocityExact')

vecDataAvgCell.SetNumberOfComponents(3)
vecDataAvgCell.SetNumberOfTuples(nr * nt)
vecDataAvgCell.SetName('velocityAvgCell')

psiData.SetNumberOfComponents(1)
psiData.SetNumberOfTuples(nr * nt * 4)
psiData.SetName('potential')

z = 0. # on z = 0 plane (2D)
k = 0
kCell = 0
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

		psi00 = psiFunc(r0, t0)
		psi10 = psiFunc(r1, t0)
		psi11 = psiFunc(r1, t1)
		psi01 = psiFunc(r0, t1)


		# integrated field values on the edges
		vs0 = psi10 - psi00 # bottom
		vs1 = psi11 - psi01 # top
		v0s = psi01 - psi00 # left
		v1s = psi11 - psi10 # right

		# edge vector velocities evaluated at the cell corners
		v00 = vs0*gradXi1(x00, y00, z) + v0s*gradXi2(x00, y00, z)
		v10 = vs0*gradXi1(x10, y10, z) + v1s*gradXi2(x10, y10, z)
		v11 = vs1*gradXi1(x11, y11, z) + v1s*gradXi2(x11, y11, z)
		v01 = vs1*gradXi1(x01, y01, z) + v0s*gradXi2(x01, y01, z)

		# exact velocity at the cell corners
		v00Exact = dPsiDr(r0, t0)*rHat(x00, y00, z) + dPsiDtOverR(r0, t0)*tHat(x00, y00, z)
		v10Exact = dPsiDr(r1, t0)*rHat(x10, y10, z) + dPsiDtOverR(r1, t0)*tHat(x10, y10, z)
		v11Exact = dPsiDr(r1, t1)*rHat(x11, y11, z) + dPsiDtOverR(r1, t1)*tHat(x11, y11, z)
		v01Exact = dPsiDr(r0, t1)*rHat(x01, y01, z) + dPsiDtOverR(r0, t1)*tHat(x01, y01, z)

		points.InsertPoint(k, (x00, y00, z))
		vecData.SetTuple(k, v00)
		vecDataExact.SetTuple(k, v00Exact)
		psiData.SetTuple(k, [psi00])
		k += 1

		points.InsertPoint(k, (x10, y10, z))
		vecData.SetTuple(k, v10)
		vecDataExact.SetTuple(k, v10Exact)
		psiData.SetTuple(k, [psi10])
		k += 1

		points.InsertPoint(k, (x11, y11, z))
		vecData.SetTuple(k, v11)
		vecDataExact.SetTuple(k, v11Exact)
		psiData.SetTuple(k, [psi11])
		k += 1

		points.InsertPoint(k, (x01, y01, z))
		vecData.SetTuple(k, v01)
		vecDataExact.SetTuple(k, v01Exact)
		psiData.SetTuple(k, [psi01])
		k += 1

		# cell average
		xCell = 0.25*(x00 + x10 + x11 + x01)
		yCell = 0.25*(y00 + y10 + y11 + y01)

		# average of the north/south edge fields. Note: vs0 and vs1 are the 
		# line integrals of the field over the edge so need to divide
		# by the edge length (dr)
		vAvgCell = 0.5*(vs0/dr + vs1/dr)*rHat(xCell, yCell, z)
		# average of the west/east edge fields. Edge length is r*dt
		vAvgCell += 0.5*(v0s/(r0*dt) + v1s/(r1*dt))*tHat(xCell, yCell, z)

		vecDataAvgCell.SetTuple(kCell, vAvgCell)
		kCell += 1

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
grid.GetPointData().AddArray(vecData)
grid.GetPointData().AddArray(vecDataExact)
grid.GetCellData().AddArray(vecDataAvgCell)
grid.GetPointData().AddArray(psiData)


# save to file
writer.SetInputData(grid)
writer.SetFileName('edgeVectorCellByCell.vtk')
writer.Update()