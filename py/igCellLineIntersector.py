import vtk
import numpy

class CellLineIntersector:
    """
    Class that finds all the intersection points between a line and a cell
    """

    def __init__(self, radius= 1.0):
        """
        Constructor
        @param radius radius of the earth 
        """

        self.radius = radius

        # intersections are computed in 3d between lines and faces, so need 
        # to create a 3d hexahedron
        self.cell = vtk.vtkHexahedron()
        self.pts = self.cell.GetPoints()

        # starting/ending points of the line
        self.pA = numpy.zeros((3,), numpy.float64)
        self.pB = numpy.zeros((3,), numpy.float64)

        # intersection point
        self.intersectPt = numpy.zeros((3,), numpy.float64)

        # parametric coordinates (3d)
        self.xi = numpy.zeros((3,), numpy.float64)

        # point closest to the intersection
        self.closestPoint = numpy.zeros((3,), numpy.float64)

        # interpolation weights, 8 values for hex
        self.weights = numpy.zeros((8,), numpy.float64)

        # tolerance for deciding if a line intersects with a cell
        self.tol = 1.e-14

    def setSphericalLine(self, lamA, theA, lamB, theB):
        """
        Set the line
        @param lamA longitude of starting point (in radiant)
        @param theA latitude of starting point (in radiant)
        @param lamB longitude of ending point
        @param theB latitude of ending point
        """
        self.pA = numpy.array([lamA, theA, 0.0])
        self.pB = numpy.array([lamB, theB, 0.0])

    def setCartesianLine(self, lamA, theA, lamB, theB):
        """
        Set the line
        @param lamA longitude of starting point (in radiant)
        @param theA latitude of starting point (in radiant)
        @param lamB longitude of ending point
        @param theB latitude of ending point
        """
        self.pA = self._getXYZ(lamA, theA)
        self.pB = self._getXYZ(lamB, theB)

    def setSphericalCell(self, lam0, the0, lam1, the1, lam2, the2, lam3, the3):
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

        # lower level
        self.pts.InsertPoint(0, (lam0, the0, -1.))
        self.pts.InsertPoint(1, (lam1, the1, -1.))
        self.pts.InsertPoint(2, (lam2, the2, -1.))
        self.pts.InsertPoint(3, (lam3, the3, -1.))

        # upper level
        self.pts.InsertPoint(4, (lam0, the0, 1.))
        self.pts.InsertPoint(5, (lam1, the1, 1.))
        self.pts.InsertPoint(6, (lam2, the2, 1.))
        self.pts.InsertPoint(7, (lam3, the3, 1.))


    def setCartesianCell(self, lam0, the0, lam1, the1, lam2, the2, lam3, the3):
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

        # lower level
        self.pts.InsertPoint(0, self._getXYZ(lam0, the0, -0.5))
        self.pts.InsertPoint(1, self._getXYZ(lam1, the1, -0.5))
        self.pts.InsertPoint(2, self._getXYZ(lam2, the2, -0.5))
        self.pts.InsertPoint(3, self._getXYZ(lam3, the3, -0.5))

        # upper level
        self.pts.InsertPoint(4, self._getXYZ(lam0, the0, 0.5))
        self.pts.InsertPoint(5, self._getXYZ(lam1, the1, 0.5))
        self.pts.InsertPoint(6, self._getXYZ(lam2, the2, 0.5))
        self.pts.InsertPoint(7, self._getXYZ(lam3, the3, 0.5))


    def isPointInsideCell(self, pt, xi):
        """
        Check if point is inside cell
        @param pt point
        @param xi parametric coordinates
        @return True is inside, false otherwise
        """
        subId = vtk.mutable(0)
        dist = vtk.mutable(0.0)
        res = self.cell.EvaluatePosition(pt, self.closestPoint, subId, xi, dist, self.weights)
        return res


    def findParametric(self, tStart, t, xi):
        """
        Find intersection of hex with line
        @paran tStart starting value for the ray
        @param t starting parametric line coordinate 0 <= t <= 1 (output)
        @param xi parametric coordinates (output)
        @return True if an intersection was found
        """
        subId = vtk.mutable(0)
        tPrime = vtk.mutable(-1.)
        pStart = self.pA + tStart*(self.pB - self.pA)
        res = self.cell.IntersectWithLine(pStart, self.pB, self.tol, 
                                           tPrime, self.intersectPt, xi, subId)
        # modify the parametric coord to start from self.pA
        t.set(tStart + tPrime.get()*(1.0 - tStart))
        return res

    def findParametricIntersection(self, tBeg, xiBeg, tEnd, xiEnd):
        """
        """
        t = vtk.mutable(-1)
        xi = numpy.zeros((3,), numpy.float64)

        # compute the parametric coordinates xiBeg and xiEnd
        if self.isPointInsideCell(self.pA, xi):
            tBeg.set(0.0)
            xiBeg[:] = xi[:2]
            #print '....... a is inside cell tBeg = {} xiBeg = {}'.format(tBeg.get(), xiBeg)
        else:
            found = self.findParametric(0.0, t, xi)
            if found:
                tBeg.set(t.get())
                xiBeg[:] = xi[:2]
                #print '....... 1st intersection tBeg = {} xiBeg = {}'.format(tBeg.get(), xiBeg)
            else:
                print 'No intersection and a is not inside cell!!!'

        if self.isPointInsideCell(self.pB, xi):
            tEnd.set(1.0)
            xiEnd[:] = xi[:2]
            #print '....... b is inside cell tEnd = {} xiEnd = {}'.format(tEnd.get(), xiEnd)
        else:
            # reset the starting point
            tStart  = tBeg.get() + self.tol
            found = self.findParametric(tStart, t, xi)
            if found:
                tEnd.set(t.get())
                xiEnd[:] = xi[:2]
                #print '....... 2nd intersection tEnd = {} xiEnd = {}'.format(tEnd.get(), xiEnd)
            else:
                print 'No intersection and b is not inside cell!!!'


    def findIntersection(self, tStart, tEnd, xiBeg, xiEnd):
        """
        Find intersection of hex with line
        @param tStart starting parametric line coordinate 0 <= tStart <= 1 (output)
        @param tEnd ending parametric line coordinate 0 <= tEnd <= 1 (output)
        @param xiBeg the starting parametric coordinates (output)
        @param xiEnd the ending parametric coordinates (output)
        @return True if an intersection was found
        """

        tStart.set(0.0)
        tEnd.set(1.0)

        # parametric position
        tprime = vtk.mutable(-1.0)
        # not used
        subId = vtk.mutable(0)
        # not used
        dist = vtk.mutable(0.0)

        hasIntersection = False

        res1 = 0
        res2 = 0
        pA = self.pA.copy()
        pB = self.pB.copy()

        # find if starting point is inside
        insideA = self.cell.EvaluatePosition(pA, self.closestPoint, subId, self.xi, dist, self.weights)
        if insideA == 1:
            # self.pA is inside cell
            tStart.set(0.0)
            xiBeg[:] = self.xi[:2]
        else:
            # self.pA is outside the cell
            self.xi *= 0 # initialize
            self.xi -= 1.0
            res1 = self.cell.IntersectWithLine(pA, pB, self.tol, tStart, self.intersectPt, self.xi, subId)
            if res1:
                # the parametric coords return by IntersectWithLine don't look right
                self.cell.EvaluatePosition(self.intersectPt, self.closestPoint, subId, self.xi, dist, self.weights)
                xiBeg[:] = self.xi[:2]
                # move the starting point up for the next search
                pA = self.intersectPt + 2*self.tol*(self.pB - self.pA)

        # find if ending position is inside
        insideB = self.cell.EvaluatePosition(pB, self.closestPoint, subId, self.xi, dist, self.weights)
        if insideB == 1:
            # self.pB is inside cell
            tEnd.set(1.0)
            xiEnd[:] = self.xi[:2]
        else:
            res2 = self.cell.IntersectWithLine(pA, pB, self.tol, tprime, self.intersectPt, self.xi, subId)
            # correct the parametric coordinate to account for moving pA to the previous intersection
            tEnd.set(tStart.get() + (1. - tStart.get())*tprime.get())
            if res2:
                self.cell.EvaluatePosition(self.intersectPt, self.closestPoint, subId, self.xi, dist, self.weights)
                xiEnd[:] = self.xi[:2]

        hasIntersection = (insideA == 1 or res1) and (insideB == 1 or res2)
        return hasIntersection


    def _getXYZ(self, lam, the, elev=0.0):
        """
        Get the Cartesian coordinates from lon-lat
        @param lam longitude in radian
        @param the latitude in radian
        @param elev elevation in normalized units
        """
        a = self.radius*(1.0 + elev)
        cos_the = numpy.cos(the)
        x = a*cos_the*numpy.cos(lam)
        y = a*cos_the*numpy.sin(lam)
        z = a*numpy.sin(the)
        return numpy.array([x, y, z])

###############################################################################

def testEasy():

    # standard, easy intersection

    tBeg, tEnd = vtk.mutable(-1.0), vtk.mutable(-1.0)

    lam0, the0 = 0., 0.
    lam1, the1 = 0.01, 0.
    lam2, the2 = 0.01, 0.01
    lam3, the3 = 0., 0.01

    lamA, theA = -0.005, 0.005
    lamB, theB = 0.005, 0.005

    cli = CellLineIntersector()
    cli.setCartesianCell(lam0, the0,
                         lam1, the1, 
                         lam2, the2, 
                         lam3, the3)
    cli.setCartesianLine(lamA, theA,
                         lamB, theB)

    xiBeg = numpy.zeros((2,), numpy.float64)
    xiEnd = numpy.zeros((2,), numpy.float64)

    found = cli.findIntersection(tBeg, tEnd, xiBeg, xiEnd)
    assert abs(xiBeg[0] - 0.0) < 1.e-5 and abs(xiBeg[1] - 0.5) < 1.e-5
    assert abs(xiEnd[0] - 0.5) < 1.e-5 and abs(xiEnd[1] - 0.5) < 1.e-5

def test2Intersections():

    # 2 intersections

    tBeg, tEnd = vtk.mutable(-1.0), vtk.mutable(-1.0)

    lam0, the0 = 0., 0.
    lam1, the1 = 0.01, 0.
    lam2, the2 = 0.01, 0.01
    lam3, the3 = 0., 0.01

    lamA, theA = -0.005, 0.005
    lamB, theB = 0.015, 0.005

    cli = CellLineIntersector()
    cli.setCartesianCell(lam0, the0,
                         lam1, the1, 
                         lam2, the2, 
                         lam3, the3)
    cli.setCartesianLine(lamA, theA,
                         lamB, theB)

    xiBeg = numpy.zeros((2,), numpy.float64)
    xiEnd = numpy.zeros((2,), numpy.float64)

    found = cli.findIntersection(tBeg, tEnd, xiBeg, xiEnd)
    assert abs(xiBeg[0] - 0.0) < 1.e-4 and abs(xiBeg[1] - 0.5) < 1.e-4
    assert abs(xiEnd[0] - 1.0) < 1.e-4 and abs(xiEnd[1] - 0.5) < 1.e-4

def testInside():

    # segment is fully inside the cell

    tBeg, tEnd = vtk.mutable(-1.0), vtk.mutable(-1.0)

    lam0, the0 = 0., 0.
    lam1, the1 = 0.01, 0.
    lam2, the2 = 0.01, 0.01
    lam3, the3 = 0., 0.01

    lamA, theA = 0.001, 0.005
    lamB, theB = 0.009, 0.005

    cli = CellLineIntersector()
    cli.setCartesianCell(lam0, the0,
                         lam1, the1, 
                         lam2, the2, 
                         lam3, the3)
    cli.setCartesianLine(lamA, theA,
                         lamB, theB)

    xiBeg = numpy.zeros((2,), numpy.float64)
    xiEnd = numpy.zeros((2,), numpy.float64)

    found = cli.findIntersection(tBeg, tEnd, xiBeg, xiEnd)
    assert abs(xiBeg[0] - 0.1) < 1.e-5 and abs(xiBeg[1] - 0.5) < 1.e-5
    assert abs(xiEnd[0] - 0.9) < 1.e-5 and abs(xiEnd[1] - 0.5) < 1.e-5

def testOutside():

    # segment is fully outside of the cell

    tBeg, tEnd = vtk.mutable(-1.0), vtk.mutable(-1.0)

    lam0, the0 = 0., 0.
    lam1, the1 = 0.01, 0.
    lam2, the2 = 0.01, 0.01
    lam3, the3 = 0., 0.01

    lamA, theA = -0.02, -0.005
    lamB, theB = +0.03, -0.002

    cli = CellLineIntersector()
    cli.setCartesianCell(lam0, the0,
                         lam1, the1, 
                         lam2, the2, 
                         lam3, the3)
    cli.setCartesianLine(lamA, theA,
                         lamB, theB)

    xiBeg = numpy.zeros((2,), numpy.float64)
    xiEnd = numpy.zeros((2,), numpy.float64)

    found = cli.findIntersection(tBeg, tEnd, xiBeg, xiEnd)
    assert not found 

def testTangent():

    # segment is tangent to cell

    tBeg, tEnd = vtk.mutable(-1.0), vtk.mutable(-1.0)

    lam0, the0 = 0., 0.
    lam1, the1 = 0.01, 0.
    lam2, the2 = 0.01, 0.01
    lam3, the3 = 0., 0.01

    lamA, theA = -0.02, -0.0
    lamB, theB = +0.03, -0.0

    cli = CellLineIntersector()
    cli.setCartesianCell(lam0, the0,
                         lam1, the1, 
                         lam2, the2, 
                         lam3, the3)
    cli.setCartesianLine(lamA, theA,
                         lamB, theB)

    xiBeg = numpy.zeros((2,), numpy.float64)
    xiEnd = numpy.zeros((2,), numpy.float64)

    found = cli.findIntersection(tBeg, tEnd, xiBeg, xiEnd)
    print('testTangent: xiBeg = {} xiEnd = {}'.format(xiBeg, xiEnd))

def testTangent2():

    # segment is tangent to cell

    tBeg, tEnd = vtk.mutable(-1.0), vtk.mutable(-1.0)

    piHalf = numpy.pi/2.

    # cell
    lam0, the0 = 0.025*piHalf, 0.05*piHalf
    lam1, the1 = 0.050*piHalf, 0.05*piHalf
    lam2, the2 = 0.050*piHalf, 0.1*piHalf
    lam3, the3 = 0.025*piHalf, 0.1*piHalf

    lamA, theA = 0.02*piHalf, 0.05*piHalf
    lamB, theB = 0.07*piHalf, 0.05*piHalf

    cli = CellLineIntersector()
    cli.setCartesianCell(lam0, the0,
                         lam1, the1, 
                         lam2, the2, 
                         lam3, the3)
    cli.setCartesianLine(lamA, theA,
                         lamB, theB)

    xiBeg = numpy.zeros((2,), numpy.float64)
    xiEnd = numpy.zeros((2,), numpy.float64)

    found = cli.findIntersection(tBeg, tEnd, xiBeg, xiEnd)
    print('testTangent2: tBeg = {} tEnd = {} xiBeg = {} xiEnd = {}'.format(tBeg.get(), tEnd.get(), xiBeg, xiEnd))
    assert abs(tBeg.get() - 0.1) < 1.e-3
    assert abs(tEnd.get() - 0.6) < 1.e-3


if __name__ == '__main__':
    testEasy()
    test2Intersections()
    testInside()
    testOutside()
    testTangent()
    testTangent2()


