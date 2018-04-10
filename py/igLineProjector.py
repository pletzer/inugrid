import vtk
import numpy
import math
from igBasisFunctionIntegral import BasisFunctionIntegral
from igCellLineIntersector import CellLineIntersector


class LineProjector:
    """
    Class to compute the projection of a 1-form onto a target line
    """

    # to handle floating point comparisons
    EPS = 1.2323435e-14


    def __init__(self, grid, integralFunction):
        """
        Constructor
        @param grid instance of vtkUnstructuredGrid
        @param integralFunction function of (lam0, the0, lam1, the1) 
                                used to set the edge values
        """
        self.grid = grid
        self.integralFunction = integralFunction

        # array of positions
        self.lineSegments = []

        self.totalProjection = 0.0

        # to find the cells of self.grid that are intersected by the line
        self.cellLoc = vtk.vtkCellLocator()
        self.cellLoc.SetDataSet(self.grid)
        self.cellLoc.BuildLocator()


    def setLine(self, lamThes):
        """
        Set the closed line
        @param lamThes [(lam, the), ...]
        """
        self.lineSegments = numpy.array([(lt[0], lt[1], 0.0) for lt in lamThes])


    def project(self):
        """
        Compute the projection
        @return value
        """

        self.totalProjection = 0.0

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
        nSegs = self.lineSegments.shape[0] - 1

        # contribution from each line segment
        self.contributions = []

        # iterate over the segments of the line
    	for iSeg in range(nSegs):

            # starting point of segment
            a = self.lineSegments[iSeg    , :]
            # end point of segement
            b = self.lineSegments[iSeg + 1, :]

            # get the lon/lat
            lamA, theA = a[:2]
            lamB, theB = b[:2]
            intersector.setLine(lamA, theA, lamB, theB)

            # (tBeg, tEnd): integral
            segment2Integral = {}

            # iterate over the grid cells that are (likely) intersected 
            # by the line segment a -> b
            for cellId in self._findCells(a, b):

                cell = self.grid.GetCell(cellId)
                ptIds = cell.GetPointIds()

                # get the vertices of the quad
                lt0, lt1, lt2, lt3 = (self.grid.GetPoint(ptIds.GetId(i)) for i in range(4))
                lam0, the0 = lt0[:2]
                lam1, the1 = lt1[:2]
                lam2, the2 = lt2[:2]
                lam3, the3 = lt3[:2]

                intersector.setSphericalCell(lam0, the0, lam1, the1, lam2, the2, lam3, the3)

                # compute tBeg, tEnd, xiBeg and xiEnd
                isIntersecting = intersector.findIntersection(tBeg, tEnd, xiBeg, xiEnd)
                if isIntersecting:

                    basisIntegrator = BasisFunctionIntegral(xiBeg, xiEnd)

                    integral = 0.0

                    # iterate over the edges
                    numPts = 4
                    for i0 in range(numPts):

                        i1 = (i0 + 1) % numPts

                        ptId0, ptId1 = ptIds.GetId(i0), ptIds.GetId(i1)
                        lt0, lt1 = self.grid.GetPoint(ptId0), self.grid.GetPoint(ptId1)

                        lam0, the0 = lt0[:2]
                        lam1, the1 = lt1[:2]

                        # assume counterclockwise orientation
                        lineIntegral= self.integralFunction(lam0, the0, lam1, the1)

                        # update the line integral
                        integral += lineIntegral * basisIntegrator(i0)

                    # this prevents sub segments from assigning duplicate line integrals to different 
                    # cells when both tBeg and tEnd are the same
                    segment2Integral[tBeg.get(), tEnd.get()] = integral

            self.contributions.append(segment2Integral)

        self.totalProjection = self._addContributions()

        return self.totalProjection


    def _findCells(self, a, b):
        """
        Find all the cells intersected by the line that goes through two points
        @param a starting point
        @param b end point
        @return list of cells
        """
        cellIds = vtk.vtkIdList()
        tol = 1.e-3
        self.cellLoc.FindCellsAlongLine(a, b, tol, cellIds)
        numCells = cellIds.GetNumberOfIds()
        return [cellIds.GetId(i) for i in range(numCells)]


    def _addContributions(self):
        """
        Add the contributions, taking care of duplicates
        @return total line integral
        """

        totalIntegral = 0.0

        iSeg = 0
        for segment2Integral in self.contributions:
            # the starting and ending parametric coordinates of the 
            # segment always start from 0 and should finish at 1.0
            totalT = 0.0
            currentTEnd = -float('inf')
            # sort the keys
            tBegEnd = sorted(segment2Integral)
            for tBeg, tEnd in tBegEnd:
                if tBeg > currentTEnd - self.EPS:
                    totalIntegral += segment2Integral[tBeg, tEnd]
                    totalT += tEnd - tBeg
                    currentTEnd = tBeg

            try:
        	   assert abs(1.0 - totalT) < 1000*self.EPS
            except:
        	   print('ERROR the integrated parametric coordinate "t" for segment {} amounts to {} != 1 (diff {})'.format(\
                     iSeg, totalT, 1. - totalT))
        	   print('      This indicates that some segments are not properly accounted for.')
        	   print('      segment2Integral = {}'.format(segment2Integral))
        	   raise RuntimeError, 'FATAL'
            iSeg += 1
        
        return totalIntegral


###############################################################################

def psi(lon, lat):
    # potential function
    # lon: longitude in rad
    # lat: latitude in rad
    return lon/(2.0 * numpy.pi)

# define form
def edgeIntegral(lona, lata, lonb, latb):
    """
    Compute the value attached to an edge
    """
    return psi(lonb, latb) - psi(lona, lata)

def testPole():

    from igLatLon import LatLon

    # create grid
    nlat, nlon = 2, 4
    sphere = LatLon(numLats=nlat, numLons=nlon, radius=1.0, 
                  coords='spherical')
    sphere.save('testPole.vtk')
    grd = sphere.getUnstructuredGrid()

    # compute vorticity
    vort = LineProjector(grd, edgeIntegral)

    lamA, theA =      0.0, math.pi/5.
    lamB, theB =  math.pi, math.pi/5.
    vort.setLine([(lamA, theA), (lamB, theB)])

    totVort = vort.project()
    exact = psi(lamB, theB) - psi(lamA, theA)
    print('testPole: total integral = {} exact = {}'.format(totVort, exact))

    # check
    assert abs(totVort - exact) < 1.e-10


if __name__ == '__main__':
    testPole()
