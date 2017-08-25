import vtk
import numpy
import math
from igBasisFunctionIntegral import BasisFunctionIntegral
from igCellLineIntersector import CellLineIntersector


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
        self.integralFunction = integralFunction

        # array of x, y, z positions along the line
        self.xyzLine = []

        self.totalFlux = 0.0

        # to find the cells of self.grid that are intersected by the line
        self.cellLoc = vtk.vtkCellLocator()
        self.cellLoc.SetDataSet(self.grid)
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
            intersector.setLine(lamA, theA, lamB, theB)

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
                    numPts = 4
                    for i0 in range(numPts):
                        i1 = (i0 + 1) % numPts

                        ptId0, ptId1 = ptIds.GetId(i0), ptIds.GetId(i1)
                        xyz0, xyz1 = self.grid.GetPoint(ptId0), self.grid.GetPoint(ptId1)

                        lam0, the0 = self._getLambdaThetaFromXYZ(xyz0)
                        lam1, the1 = self._getLambdaThetaFromXYZ(xyz1)

                        # assumes counterclockwise orientation
                        faceFlux = self.integralFunction(lam0, lam1, the0, the1)

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
        numCells = cellIds.GetNumberOfIds()
        return [cellIds.GetId(i) for i in range(numCells)]

    def _getXYZFromLambdaTheta(self, lam, the):
        """
        Convert from lon/lat to x, y, z
        @param lam longitude in radiant
        @param the latitude in radiant
        @return Cartesian coordinates
        """
        rho = math.cos(the)
        x = rho * math.cos(lam)
        y = rho * math.sin(lam)
        z = numpy.sin(the)
        return x, y, z

    def _getLambdaThetaFromXYZ(self, xyz):
        """
        Convert from Cartesian to lon/lat coordinates
        @param xyz point 
        @return longitude, latitude
        """
        x, y, z = xyz
        rho = math.sqrt(x**2 + y**2)
        the = math.atan2(z, rho)
        lam = math.atan2(y, x)
        return lam, the

###############################################################################

def psi(x, y):
    # stream function
    # x: longitude
    # y: latitude
    return math.sin(x)*math.cos(y)

# define form
def edgeIntegral(xa, xb, ya, yb):
    """
    Compute the value attached to an edge
    x is longitude
    y is latitude
    """
    return psi(xb, yb) - psi(xa, ya)

def testClosed():
    from igLatLon import LatLon

    # create grid
    nlat, nlon = 10, 20
    coord = LatLon(numLats=nlat, numLons=nlon)
    grd = coord.getUnstructuredGrid()

    # compute flux
    fc = FluxCalculator(grd, edgeIntegral)

    lamA, theA = -math.pi, math.pi/5.
    lamB, theB = math.pi, math.pi/5.
    line = numpy.array([(lamA, theA), (lamB, theB)], numpy.float64).reshape(2, 2)
    fc.setLine(line)

    totFlux = fc.computeFlux()
    exact = psi(lamB, theB) - psi(lamA, theA)
    print('testClosed: total flux = {} exact = {}'.format(totFlux, exact))

    # check
    assert abs(totFlux - exact) < 1.e-10

def testOpen():
    from igLatLon import LatLon

    # create grid
    nlat, nlon = 10, 20
    coord = LatLon(numLats=nlat, numLons=nlon)
    grd = coord.getUnstructuredGrid()

    # compute flux
    fc = FluxCalculator(grd, edgeIntegral)

    lamA, theA = -math.pi, math.pi/5.324325
    lamB, theB = 0.0, math.pi/5.324325
    line = numpy.array([(lamA, theA), (lamB, theB)], numpy.float64).reshape(2, 2)
    fc.setLine(line)


    totFlux = fc.computeFlux()
    exact = psi(lamB, theB) - psi(lamA, theA)
    print('testOpen: total flux = {} exact = {}'.format(totFlux, exact))

    # check
    assert abs(totFlux - exact) < 1.e-10

def testOpenSmall():
    from igLatLon import LatLon

    # create grid
    nlat, nlon = 2, 4
    coord = LatLon(numLats=nlat, numLons=nlon)
    grd = coord.getUnstructuredGrid()

    # compute flux
    fc = FluxCalculator(grd, edgeIntegral)

    lamA, theA =        0.0, 0.0
    lamB, theB = math.pi/2., 0.0
    line = numpy.array([(lamA, theA), (lamB, theB)], numpy.float64).reshape(2, 2)
    fc.setLine(line)

    totFlux = fc.computeFlux()
    exact = psi(lamB, theB) - psi(lamA, theA)
    print('testOpenSmall: total flux = {} exact = {}'.format(totFlux, exact))

    # check
    assert abs(totFlux - exact) < 1.e-10

def testOpenSmall2():
    from igLatLon import LatLon

    # create grid
    nlat, nlon = 2, 4
    coord = LatLon(numLats=nlat, numLons=nlon)
    grd = coord.getUnstructuredGrid()

    # compute flux
    fc = FluxCalculator(grd, edgeIntegral)

    lamA, theA =        0.0, 1.e-6
    lamB, theB = math.pi/2., 1.e-6
    line = numpy.array([(lamA, theA), (lamB, theB)], numpy.float64).reshape(2, 2)
    fc.setLine(line)

    totFlux = fc.computeFlux()
    exact = psi(lamB, theB) - psi(lamA, theA)
    print('testOpenSmall: total flux = {} exact = {}'.format(totFlux, exact))

    # check
    assert abs(totFlux - exact) < 1.e-10

def testOpenSmall3():
    from igLatLon import LatLon

    # create grid
    nlat, nlon = 2, 4
    coord = LatLon(numLats=nlat, numLons=nlon)
    grd = coord.getUnstructuredGrid()

    # compute flux
    fc = FluxCalculator(grd, edgeIntegral)

    piHalf = 0.5*math.pi
    lamA, theA = piHalf*0.2, piHalf*0.5
    lamB, theB = piHalf*0.8, piHalf*0.5
    line = numpy.array([(lamA, theA), (lamB, theB)], numpy.float64).reshape(2, 2)
    fc.setLine(line)

    totFlux = fc.computeFlux()
    exact = psi(lamB, theB) - psi(lamA, theA)
    print('testOpenSmall: total flux = {} exact = {}'.format(totFlux, exact))

    # check
    assert abs(totFlux - exact) < 1.e-10

if __name__ == '__main__':
    testOpenSmall()
    testOpenSmall2()
    testOpenSmall3()
    #testClosed()
