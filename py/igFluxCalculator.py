import vtk
import numpy

class BasisFunctionIntegral:

    def __init__(self, xisA, xisB):
        # difference
        dXis = xisB - xisA
        # average
        aXis = 0.5*(xisA + xisB)

        self.value = {0: dXis[0]*(1.0 - aXis[2]),
                      1: dXis[1]*aXis[0],
                      2: -dXis[0]*aXis[1],
                      3: -dXis[1]*(1.0 - aXis[0]),}

    def __call__(self, index):
        return self.value[index]

class CellLineIntersector:

    def __init__(self,):
        self.pts = vtk.vtkPoints()
        self.pts.SetNumberOfPoints(4)
        self.quad = vtk.vtkQuad()
        self.quad.SetPoints(self.pts)
        self.pA = numpy.zeros((3,), numpy.float64)
        self.pB = numpy.zeros((3,), numpy.float64)
        self.x = numpy.zeros((3,), numpy.flaot64)
        self.tol = 1.e-10

    def setLine(self, lamA, thetA, lamB, theB):
        self.pA[:] = lamA, theA, 0.0
        self.pB[:] = lamB, theB, 0.0

    def setCell(self, lam0, the0, lam1, the1, lam2, the2, lam3, the3):
        self.pts.InsertPoint(0, lam0, the0, 0.0)
        self.pts.InsertPoint(1, lam1, the1, 0.0)
        self.pts.InsertPoint(2, lam2, the2, 0.0)
        self.pts.InsertPoint(3, lam3, the3, 0.0)

    def findIntersection(self, xiBeg, xiEnd, xis):
        t = vtk.mutable(-1.0)
        subId = vtk.mutable(0)
        res = self.quad.IntersectWithLine(self.pA, self.pB, self.tol, self.x, subId)
        return (res != 0)


class FluxCalculator:

    def __init__(self, grid, integralFunction):

    	self.grid = grid
        self.fieldName = ''
        self.xyzLine = []
        self.totalFlux = 0.0

        self.cellLoc = vtk.vtkCellLocator()
        self.cellLoc.SetDataSet(self.grid.GetOutput())
        self.cellLoc.BuildLocator()

    def setLine(self, lamThes):
    	n = len(lamThes)
        self.xyzLine = numpy.zeros((n, 3), numpy.float64)
        for i in range(n):
            lam, the = lamThes[i, :]
            self.xyzLine[i, :] = self._getXYZFromLambdaTheta(lam, the)

    def computeFlux(self):

        self.totalFlux = 0.0
        xiBeg = numpy.zeros((2,), numpy.float64) # 2D
        xiEnd = numpy.zeros((2,), numpy.float64) # 2D
        intersector = CellLineIntersector()
        nSegs = self.xyzLine.shape[0] - 1
        # iterate over the segments of the line
    	for iSeg in range(nSegs):
            xyzA = self.xyzLine[iSeg, :]
            xyzB = self.xyzLine[iSeg + 1, :]
            lamA, theA = self._getLambdaThetaFromXYZ(xyzA)
            lamB, theB = self._getLambdaThetaFromXYZ(xyzB)
            intersector.setLine(lamA, thetA, lamB, theB)


            # iterate over the grid cells that are (likely) intersected by the line segment
            # xyzA -> xyzB
            for cellId in self._findCells(xyzA, xyzB):
                cell = self.grid.GetCell(cellId)
                ptIds = cell.GetPointIds()

                ptId0 = ptIds.GetId(0)
                ptId1 = ptIds.GetId(1)
                ptId2 = ptIds.GetId(2)
                ptId3 = ptIds.GetId(3)

                xyz0 = self.grid.GetPoint(ptId0)
                xyz1 = self.grid.GetPoint(ptId1)
                xyz2 = self.grid.GetPoint(ptId2)
                xyz3 = self.grid.GetPoint(ptId3)

                lam0, the0 = self._getLambdaThetaFromXYZ(xyz0)
                lam1, the1 = self._getLambdaThetaFromXYZ(xyz1)
                lam2, the2 = self._getLambdaThetaFromXYZ(xyz2)
                lam3, the3 = self._getLambdaThetaFromXYZ(xyz3)

                intersector.setCell(lam0, the0, lam1, the1, lam2, the2, lam3, the3)
                isIntersecting = intersector.findIntersection(xiBeg, xiEnd)
                if isIntersecting:

                    basisIntegrator = BasisFunctionIntegral(xiBeg, xiEnd)

                    # iterate over the edges
                    for i0 in range(numPts):
                        i1 = (i0 + 1) % numPts
                        ptId0, ptId1 = ptIds.GetId(i0), ptIds.GetId(i1)
                        xyz0, xyz1 = self.grid.GetPoint(ptId0), self.grid.GetPoint(ptId1)
                        lam0, the0 = self._getLambdaThetaFromXYZ(xyz0)
                        lam1, the1 = self._getLambdaThetaFromXYZ(xyz1)
                        faceFlux = integralFunction(lam0, the0, lam1, the1)
                        self.totalFlux += faceFlux * basisIntegrator(i0)

        return self.totalFlux


    def _findCells(self, xyz0, xyz1):
    	cellIds = vtk.vtkIdList()
        tol = 1.e-3
        self.cellLoc.FindCellsAlongLine(xyz0, xyz1, tol, cellIds)
        numCells = cellids.GetNumberOfIds()
        return [cellIds.GetId(i) for i in range(numCells)]

    def _getXYZFromLambdaTheta(self, lam, the):
        rho = cos(the)
        x = rho * cos(lam)
        y = rho * sin(lam)
        z = sin(the)
        return x, y, z

    def _getLambdaThetaFromXYZ(self, xyz):
        x, y, z = xyz
        rho = sqrt(x**2 + y**2)
        the = atan2(z, rho)
        lam = atan2(y, x)
        return lam, the

