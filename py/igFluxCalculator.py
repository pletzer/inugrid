import vtk
import numpy
import math
from igBasisFunctionIntegral import BasisFunctionIntegral
from igCellLineIntersector import CellLineIntersector


class FluxCalculator:
    """
    Class to compute flux across a segmented line
    """

    # to handle floating point comparisons
    EPS = 1.2323435e-14


    def __init__(self, grid, integralFunction):
        """
        Constructor
        @param grid instance of vtkUnstructuredGrid
        @param integralFunction function of (lam0, the0, lam1, the1)
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

        # compute the total distance
        self.lineDistance = 0.0
        for i in range(n - 1):
        	p0 = self.xyzLine[i, :]
        	p1 = self.xyzLine[i + 1, :]
        	dp = p1 - p0
        	self.lineDistance += numpy.sqrt(numpy.dot(dp, dp))

    def computeFlux(self):
        """
        Compute the flux
        @return flux
        """

        self.totalFlux = 0.0
        if self.lineDistance < self.EPS:
        	return 0.0

        # line parametric coordinates, will be filled in by 
        # the line-cell intersector
        tBeg = vtk.mutable(-1.)
        tEnd = vtk.mutable(-1.)

        # parametric start/end points in the cell
        # these will be filled in by the intersector
        xiBeg = numpy.zeros((2,), numpy.float64) # 2D
        xiEnd = numpy.zeros((2,), numpy.float64) # 2D

        # to find the intersection between a quad and a line
        intersector = CellLineIntersector()

        # number of segments
        nSegs = self.xyzLine.shape[0] - 1
        self.polyline2Flux = []

        # iterate over the segments of the line
    	for iSeg in range(nSegs):
            xyzA = self.xyzLine[iSeg, :]
            xyzB = self.xyzLine[iSeg + 1, :]

            # get the lon/lat
            lamA, theA = self._getLambdaThetaFromXYZ(xyzA)
            lamB, theB = self._getLambdaThetaFromXYZ(xyzB)
            intersector.setCartesianLine(lamA, theA, lamB, theB)
            #print '*** iSeg={} (lamA, theA)={} (lamB, theB)={}'.format(iSeg, (lamA, theA), (lamB, theB))

            # (tBeg, tEnd): flux
            segment2Flux = {}

            # iterate over the grid cells that are (likely) intersected by the line segment
            # A -> B
            for cellId in self._findCells(xyzA, xyzB):

                cell = self.grid.GetCell(cellId)
                ptIds = cell.GetPointIds()

                # get the vertices of the quad
                xyz0, xyz1, xyz2, xyz3 = (self.grid.GetPoint(ptIds.GetId(i)) for i in range(4))
                lam0, the0 = self._getLambdaThetaFromXYZ(xyz0)
                lam1, the1 = self._getLambdaThetaFromXYZ(xyz1)
                lam2, the2 = self._getLambdaThetaFromXYZ(xyz2)
                lam3, the3 = self._getLambdaThetaFromXYZ(xyz3)
                #print '\t--- cellId={} verts 0={} 1={} 2={} 3={}'.format(cellId, (lam0, the0), (lam1, the1), (lam2, the2), (lam3, the3),)

                intersector.setCartesianCell(lam0, the0, lam1, the1, lam2, the2, lam3, the3)

                # compute tBeg, tEnd, xiBeg and xiEnd
                isIntersecting = intersector.findIntersection(tBeg, tEnd, xiBeg, xiEnd)
                if isIntersecting:
                    #print '\t\t... found intersection tBeg={} tEnd={} xiBeg={} xiEnd={}'.format(tBeg.get(), tEnd.get(), xiBeg, xiEnd)

                    basisIntegrator = BasisFunctionIntegral(xiBeg, xiEnd)

                    flux = 0.0

                    # iterate over the edges
                    numPts = 4
                    for i0 in range(numPts):
                        i1 = (i0 + 1) % numPts

                        ptId0, ptId1 = ptIds.GetId(i0), ptIds.GetId(i1)
                        xyz0, xyz1 = self.grid.GetPoint(ptId0), self.grid.GetPoint(ptId1)

                        lam0, the0 = self._getLambdaThetaFromXYZ(xyz0)
                        lam1, the1 = self._getLambdaThetaFromXYZ(xyz1)

                        # assumes counterclockwise orientation
                        faceFlux = self.integralFunction(lam0, the0, lam1, the1)

                        # update the fluxes
                        flux += faceFlux * basisIntegrator(i0)
                        #print '\t\t\tedge = ', i0, i1, ' face flux = ', faceFlux, ' integral of basis = ', basisIntegrator(i0), ' adding flux = ', faceFlux * basisIntegrator(i0)

                    # this prevents sub segments from assigning duplicate fluxes to different 
                    # cells when both tBeg and tEnd are the same
                    segment2Flux[tBeg.get(), tEnd.get()] = flux

            self.polyline2Flux.append(segment2Flux)

        self.totalFlux = self._addFluxes()

        return self.totalFlux


    def _addFluxes(self):
        """
        Add the fluxes, taking care of duplicates
        @return total flux
        """
        totalFlux = 0.0
        iSeg = 0
        for segment2Flux in self.polyline2Flux:
            # the starting and ending parametric coordinates of the 
            # segment always start from 0 and should finish at 1.0
            totalT = 0.0
            currentTEnd = -float('inf')
            # sort the keys
            tBegEnd = sorted(segment2Flux)
            for tBeg, tEnd in tBegEnd:
                if tBeg > currentTEnd - self.EPS:
                    totalFlux += segment2Flux[tBeg, tEnd]
                    totalT += tEnd - tBeg
                    currentTEnd = tBeg

            try:
        	   assert abs(1.0 - totalT) < 1000*self.EPS
            except:
        	   print('ERROR the integrated parametric coordinate "t" for segment {} amounts to {} != 1 (diff {})'.format(\
                     iSeg, totalT, 1. - totalT))
        	   print('      This indicates that some segments are not properly accounted for.')
        	   print('      segment2Flux = {}'.format(segment2Flux))
        	   raise RuntimeError, 'FATAL'
            iSeg += 1
        
        return totalFlux


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
def edgeIntegral(xa, ya, xb, yb):
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
    print('testOpenSmall2: total flux = {} exact = {}'.format(totFlux, exact))

    # check
    assert abs(totFlux - exact) < 1.e-6 # 1.e-10

def testOpenSmall3():
    from igLatLon import LatLon

    # create grid
    nlat, nlon = 8, 16
    coord = LatLon(numLats=nlat, numLons=nlon)
    grd = coord.getUnstructuredGrid()

    # compute flux
    fc = FluxCalculator(grd, edgeIntegral)

    piHalf = 0.5*math.pi
    lamA, theA = piHalf*0.2, piHalf*0.5 #piHalf*0.2, piHalf*0.5
    lamB, theB = piHalf*0.7, piHalf*0.5 # piHalf*0.8, piHalf*0.5
    line = numpy.array([(lamA, theA), (lamB, theB)], numpy.float64).reshape(2, 2)
    fc.setLine(line)

    totFlux = fc.computeFlux()
    exact = psi(lamB, theB) - psi(lamA, theA)
    print('testOpenSmall3: total flux = {} exact = {}'.format(totFlux, exact))

    # check
    assert abs(totFlux - exact) < 0.02

def testOpenSmall3TwoSegments():
    from igLatLon import LatLon

    # create grid
    nlat, nlon = 8, 16
    coord = LatLon(numLats=nlat, numLons=nlon)
    grd = coord.getUnstructuredGrid()

    # compute flux
    fc = FluxCalculator(grd, edgeIntegral)

    piHalf = 0.5*math.pi
    lamA, theA = piHalf*0.2, piHalf*0.5
    lamB, theB = piHalf*0.21, piHalf*0.5
    lamC, theC = piHalf*0.7, piHalf*0.5
    line = numpy.array([(lamA, theA), (lamB, theB), (lamC, theC)], numpy.float64).reshape(3, 2)
    fc.setLine(line)

    totFlux = fc.computeFlux()
    exact = psi(lamC, theC) - psi(lamA, theA)
    print('testOpenSmall3TwoSegments: total flux = {} exact = {}'.format(totFlux, exact))

    # check
    assert abs(totFlux - exact) < 0.02

def testOpenSmall3ThreeSegments():
    from igLatLon import LatLon

    # create grid
    nlat, nlon = 8, 16
    coord = LatLon(numLats=nlat, numLons=nlon)
    grd = coord.getUnstructuredGrid()

    # compute flux
    fc = FluxCalculator(grd, edgeIntegral)

    piHalf = 0.5*math.pi
    lamA, theA = piHalf*0.2, piHalf*0.5
    lamB, theB = piHalf*0.21, piHalf*0.5
    lamC, theC = piHalf*0.6, piHalf*0.5
    lamD, theD = piHalf*0.7, piHalf*0.5
    line = numpy.array([(lamA, theA), (lamB, theB), (lamC, theC), (lamD, theD)], numpy.float64).reshape(4, 2)
    fc.setLine(line)

    totFlux = fc.computeFlux()
    exact = psi(lamD, theD) - psi(lamA, theA)
    print('testOpenSmall3ThreeSegments: total flux = {} exact = {}'.format(totFlux, exact))

    # check
    assert abs(totFlux - exact) < 0.02


def testZero():
    from igLatLon import LatLon

    # create grid
    nlat, nlon = 8, 16
    coord = LatLon(numLats=nlat, numLons=nlon)
    grd = coord.getUnstructuredGrid()

    # compute flux
    fc = FluxCalculator(grd, edgeIntegral)

    piHalf = 0.5*math.pi
    lamA, theA = piHalf*0.2, piHalf*0.5
    lamB, theB = piHalf*0.2, piHalf*0.5 # piHalf*0.8, piHalf*0.5
    line = numpy.array([(lamA, theA), (lamB, theB)], numpy.float64).reshape(2, 2)
    fc.setLine(line)

    totFlux = fc.computeFlux()
    exact = psi(lamB, theB) - psi(lamA, theA)
    print('testZero: total flux = {} exact = {}'.format(totFlux, exact))

    # check
    assert abs(totFlux - exact) < 1.e-10

def testCubedSphere1():
    from igCubedSphere import CubedSphere

    n = 10
    cs = CubedSphere(n)
    grid = cs.getUnstructuredGrid()

    fc = FluxCalculator(grid, edgeIntegral)

    numSegments = 1

    piHalf = 0.5*math.pi
    lamA, theA = 0.2 * piHalf, 0.5*piHalf
    lamB, theB = 0.7 * piHalf, 0.5*piHalf

    dLam, dThe = (lamB - lamA)/float(numSegments), (theB - theA)/float(numSegments)

    # expression for the contour in lon-lat coordinates
    line = numpy.array([(lamA + i*dLam, theA + i*dThe) for i in range(numSegments + 1)]).reshape(numSegments + 1, 2)
    fc.setLine(line)

    totFlux = fc.computeFlux()

    print('testCubedSphere1: Total flux: {}'.format(totFlux))
    exact = psi(lamB, theB) - psi(lamA, theA)
    assert abs(totFlux - exact) < 0.02

def testCubedSphere2():
    from igCubedSphere import CubedSphere

    n = 10
    cs = CubedSphere(n)
    grid = cs.getUnstructuredGrid()

    fc = FluxCalculator(grid, edgeIntegral)

    numSegments = 2

    piHalf = 0.5*math.pi
    lamA, theA = 0.2 * piHalf, 0.5*piHalf
    lamB, theB = 0.7 * piHalf, 0.5*piHalf

    dLam, dThe = (lamB - lamA)/float(numSegments), (theB - theA)/float(numSegments)

    # expression for the contour in lon-lat coordinates
    line = numpy.array([(lamA + i*dLam, theA + i*dThe) for i in range(numSegments + 1)]).reshape(numSegments + 1, 2)
    fc.setLine(line)

    totFlux = fc.computeFlux()

    print('testCubedSphere2: Total flux: {}'.format(totFlux))
    exact = psi(lamB, theB) - psi(lamA, theA)
    assert abs(totFlux - exact) < 0.02

def testCubedSphere3():
    from igCubedSphere import CubedSphere

    n = 10
    cs = CubedSphere(n)
    grid = cs.getUnstructuredGrid()

    fc = FluxCalculator(grid, edgeIntegral)

    numSegments = 20

    piHalf = 0.5*math.pi
    lamA, lamB = -numpy.pi, numpy.pi
    theA, theB = 0.9*piHalf, 0.9*piHalf

    dLam, dThe = (lamB - lamA)/float(numSegments), (theB - theA)/float(numSegments)

    # expression for the contour in lon-lat coordinates
    line = numpy.array([(lamA + i*dLam, theA + i*dThe) for i in range(numSegments + 1)]).reshape(numSegments + 1, 2)
    fc.setLine(line)

    totFlux = fc.computeFlux()
    print('fc.polyline2Flux = {}'.format(fc.polyline2Flux))

    exact = psi(lamB, theB) - psi(lamA, theA)
    print('testCubedSphere3: Total flux: {} exact: {}'.format(totFlux, exact))
    assert abs(totFlux - exact) < 0.02


if __name__ == '__main__':
    testCubedSphere3()

    testClosed()
    testOpenSmall()
    testOpenSmall2()
    testZero()
    testOpenSmall3()
    testOpenSmall3TwoSegments()
    testOpenSmall3ThreeSegments()
    testCubedSphere1()
    testCubedSphere2()
