import vtk
import numpy

class BasisFunctionIntegral:
    """
    Class to interpolate fluxes between faces
    """

    def __init__(self, xisA, xisB):
        """
        Constructor:
        @param xisA starting parametric position (2-vector)
        @param xisB ending parametric position (2-vector)
        """
        # difference
        dXis = xisB - xisA
        # average
        aXis = 0.5*(xisA + xisB)

        self.value = {0: dXis[0]*(1.0 - aXis[1]),
                      1: dXis[1]*aXis[0],
                      2: -dXis[0]*aXis[1],           # negative sign 
                      3: -dXis[1]*(1.0 - aXis[0]),}  # negative sign

    def __call__(self, index):
        """
        Evaluate the integral for basis function "index"
        @param index in the range(0, 4)
        @return integral
        """
        return self.value[index]


class CellLineIntersector:
    """
    Class that finds all the cells that are intersected by a line
    """

    def __init__(self,):
        """
        Constructor
        """
        self.quad = vtk.vtkQuad()
        self.pts = self.quad.GetPoints()

        # starting point of the line
        self.pA = numpy.zeros((3,), numpy.float64)
        # ending point of the line
        self.pB = numpy.zeros((3,), numpy.float64)

        # intersection point
        self.x = numpy.zeros((3,), numpy.float64)

        # parametric coordinates
        self.xi = numpy.zeros((3,), numpy.float64)

        # point closest to the intersection
        self.closestPoint = numpy.zeros((3,), numpy.float64)

        # interpolation weights, 4 values for a quad
        self.weights = numpy.zeros((4,), numpy.float64)

        # tolerance for decising if a line intersects with a cell
        self.tol = 1.e-14

    def setLine(self, lamA, theA, lamB, theB):
        """
        Set the line
        @param lamA longitude of starting point (in radiant)
        @param theA latitude of starting point (in radiant)
        @param lamB longitude of ending point
        @param theB latitude of ending point
        """
        self.pA[:] = lamA, theA, 0.0
        self.pB[:] = lamB, theB, 0.0

    def setCell(self, lam0, the0, lam1, the1, lam2, the2, lam3, the3):
        """
        Set the cell vertices
        @param lam0 longitude vertex
        @param the0 latitude vertex
        @param lam1 longitude vertex
        @param the1 latitude vertex
        @param lam2 longitude vertex
        @param the2 latitude vertex
        @param lam3 longitude vertex
        @param the3 latitude vertex
        """
        self.pts.InsertPoint(0, lam0, the0, 0.0)
        self.pts.InsertPoint(1, lam1, the1, 0.0)
        self.pts.InsertPoint(2, lam2, the2, 0.0)
        self.pts.InsertPoint(3, lam3, the3, 0.0)

    def findIntersection(self, xiBeg, xiEnd):
        """
        Find intersection of quad with line
        @param xiBeg the starting parametric coordinates (output)
        @param xiEnd the ending parametric coordinates (output)
        @return True if an intersection was found
        """
        # parametric coordinate along the line (0 <= t <= 1)
        t = vtk.mutable(-1.0)
        subId = vtk.mutable(-1)

        dist = vtk.mutable(0.0)
        res1 = 0
        res2 = 0
        pA = self.pA.copy()
        pB = self.pB.copy()

        # find if starting point is inside
        insideA = self.quad.EvaluatePosition(pA, self.closestPoint, subId, self.xi, dist, self.weights)
        if insideA == 1:
            # self.pA is inside cell
            xiBeg[:] = self.xi[:2]
        else:
            # self.pA is outside the cell
            res1 = self.quad.IntersectWithLine(pA, pB, self.tol, t, self.x, self.xi, subId)
            if res1:
                xiBeg[:] = self.xi[:2]
                # move the starting point up for the next search
                pA = self.x + 2*self.tol*(self.pB - self.pA)

        # find if ending position is inside
        insideB = self.quad.EvaluatePosition(pB, self.closestPoint, subId, self.xi, dist, self.weights)
        if insideB == 1:
            # self.pB is inside cell
            xiEnd[:] = self.xi[:2]
        else:
            res2 = self.quad.IntersectWithLine(pA, pB, self.tol, t, self.x, self.xi, subId)
            if res2:
                xiEnd[:] = self.xi[:2]

        return (insideA == 1 or res1) and (insideB == 1 or res2)


class FluxCalculator:
    """
    Class to compute flux across a segmented line
    """

    def __init__(self, grid, integralFunction):
        """
        Constructor
        @param grid instance of vtkUnstructuredGrid
        @param integralFunction function of (lam0, lam1, the0, the1)
        """
    	self.grid = grid
        # array of x, y, z positions along the line
        self.xyzLine = []

        self.totalFlux = 0.0

        # to find the cells of self.grid that are intersected by the line
        self.cellLoc = vtk.vtkCellLocator()
        self.cellLoc.SetDataSet(self.grid.GetOutput())
        self.cellLoc.BuildLocator()

    def setLine(self, lamThes):
        """
        Set the line
        @param lamThes [(lam, the), ...]
        """
    	n = len(lamThes)
        # convert to x, y, z coordinates
        self.xyzLine = numpy.zeros((n, 3), numpy.float64)
        for i in range(n):
            lam, the = lamThes[i, :]
            self.xyzLine[i, :] = self._getXYZFromLambdaTheta(lam, the)

    def computeFlux(self):
        """
        Compute the flux
        @return flux
        """
        self.totalFlux = 0.0

        # parametric strat end points in the cell
        xiBeg = numpy.zeros((2,), numpy.float64) # 2D
        xiEnd = numpy.zeros((2,), numpy.float64) # 2D

        # to find the intersection between a quad and a line
        intersector = CellLineIntersector()

        # number of segments
        nSegs = self.xyzLine.shape[0] - 1

        # iterate over the segments of the line
    	for iSeg in range(nSegs):
            xyzA = self.xyzLine[iSeg, :]
            xyzB = self.xyzLine[iSeg + 1, :]

            # get the lon/lat
            lamA, theA = self._getLambdaThetaFromXYZ(xyzA)
            lamB, theB = self._getLambdaThetaFromXYZ(xyzB)
            intersector.setLine(lamA, thetA, lamB, theB)

            # iterate over the grid cells that are (likely) intersected by the line segment
            # xyzA -> xyzB
            for cellId in self._findCells(xyzA, xyzB):

                cell = self.grid.GetCell(cellId)
                ptIds = cell.GetPointIds()

                # get the vertices of the quad
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

                # compute xiBeg and xiEnd
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

                        # assumes counterclockwise orientation
                        faceFlux = integralFunction(lam0, lam1, the0, the1)

                        # update the fluxes
                        self.totalFlux += faceFlux * basisIntegrator(i0)

        return self.totalFlux


    def _findCells(self, xyz0, xyz1):
        """
        Find all the cells intersected by the line that goes through xyz0 and xyz1
        @param xyz0 starting point
        @param xyz1 ending point
        @return list of cells
        """
    	cellIds = vtk.vtkIdList()
        tol = 1.e-3
        self.cellLoc.FindCellsAlongLine(xyz0, xyz1, tol, cellIds)
        numCells = cellids.GetNumberOfIds()
        return [cellIds.GetId(i) for i in range(numCells)]

    def _getXYZFromLambdaTheta(self, lam, the):
        """
        Convert from lon/lat to x, y, z
        @param lam longitude in radiant
        @param the latitude in radiant
        @return Cartesian coordinates
        """
        rho = cos(the)
        x = rho * cos(lam)
        y = rho * sin(lam)
        z = sin(the)
        return x, y, z

    def _getLambdaThetaFromXYZ(self, xyz):
        """
        Convert from Cartesian to lon/lat coordinates
        @param xyz point 
        @return longitude, latitude
        """
        x, y, z = xyz
        rho = sqrt(x**2 + y**2)
        the = atan2(z, rho)
        lam = atan2(y, x)
        return lam, the

###############################################################################
def testBasisFunctionIntegral0():
    xiA, xiB = numpy.array((0., 0.)), numpy.array((1., 0.))
    bi = BasisFunctionIntegral(xiA, xiB)
    for i in range(4):
        print('testBasisFunctionIntegral0: i = {} weight = {}'.format(i, bi(i)))

def testBasisFunctionIntegral1():
    xiA, xiB = numpy.array((1., 0.)), numpy.array((1., 1.))
    bi = BasisFunctionIntegral(xiA, xiB)
    for i in range(4):
        print('testBasisFunctionIntegral1: i = {} weight = {}'.format(i, bi(i)))

def testBasisFunctionIntegral2():
    xiA, xiB = numpy.array((1., 1.)), numpy.array((0., 1.))
    bi = BasisFunctionIntegral(xiA, xiB)
    for i in range(4):
        print('testBasisFunctionIntegral2: i = {} weight = {}'.format(i, bi(i)))

def testBasisFunctionIntegral3():
    xiA, xiB = numpy.array((0., 1.)), numpy.array((0., 0.))
    bi = BasisFunctionIntegral(xiA, xiB)
    for i in range(4):
        print('testBasisFunctionIntegral3: i = {} weight = {}'.format(i, bi(i)))

def testCellLineIntersector0():

    lam0, the0 = 0., 0.
    lam1, the1 = 1., 0.
    lam2, the2 = 1., 1.
    lam3, the3 = 0., 1.

    lamA, theA = -0.5, 0.5
    lamB, theB = 0.5, 0.5

    cli = CellLineIntersector()
    cli.setCell(lam0, the0,
                lam1, the1, 
                lam2, the2, 
                lam3, the3)
    cli.setLine(lamA, theA,
                lamB, theB)

    xiBeg = numpy.zeros((2,), numpy.float64)
    xiEnd = numpy.zeros((2,), numpy.float64)

    found = cli.findIntersection(xiBeg, xiEnd)

    print('testCellLineIntersector0: found = {} xiBeg = {} xiEnd = {}'.format(found, xiBeg, xiEnd))

if __name__ == '__main__':
    testBasisFunctionIntegral0()
    testBasisFunctionIntegral1()
    testBasisFunctionIntegral2()
    testBasisFunctionIntegral3()

    testCellLineIntersector0()

