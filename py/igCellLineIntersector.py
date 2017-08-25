import vtk
import numpy

class CellLineIntersector:
    """
    Class that finds all the cells that are intersected by a line
    """

    def __init__(self,):
        """
        Constructor
        """
        # intersecvtions are computed in 3d between lines and faces, so need 
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

        # lower level
        self.pts.InsertPoint(0, lam0, the0, -0.5)
        self.pts.InsertPoint(1, lam1, the1, -0.5)
        self.pts.InsertPoint(2, lam2, the2, -0.5)
        self.pts.InsertPoint(3, lam3, the3, -0.5)

        # upper level
        self.pts.InsertPoint(4, lam0, the0, 0.5)
        self.pts.InsertPoint(5, lam1, the1, 0.5)
        self.pts.InsertPoint(6, lam2, the2, 0.5)
        self.pts.InsertPoint(7, lam3, the3, 0.5)


    def findIntersection(self, xiBeg, xiEnd):
        """
        Find intersection of hex with line
        @param xiBeg the starting parametric coordinates (output)
        @param xiEnd the ending parametric coordinates (output)
        @return True if an intersection was found
        """
        # parametric coordinate along the line (0 <= t <= 1)
        t = vtk.mutable(0.0)
        subId = vtk.mutable(0)

        dist = vtk.mutable(0.0)
        res1 = 0
        res2 = 0
        pA = self.pA.copy()
        pB = self.pB.copy()

        # find if starting point is inside
        insideA = self.cell.EvaluatePosition(pA, self.closestPoint, subId, self.xi, dist, self.weights)
        if insideA == 1:
            # self.pA is inside cell
            xiBeg[:] = self.xi[:2]
        else:
            # self.pA is outside the cell
            self.xi *= 0 # initialize
            self.xi -= 1.0
            res1 = self.cell.IntersectWithLine(pA, pB, self.tol, t, self.intersectPt, self.xi, subId)
            if res1:
                # it seems the parametric coords return by IntersectWithLine don't look right
                self.cell.EvaluatePosition(self.intersectPt, self.closestPoint, subId, self.xi, dist, self.weights)
                xiBeg[:] = self.xi[:2]
                # move the starting point up for the next search
                pA = self.intersectPt + 2*self.tol*(self.pB - self.pA)

        # find if ending position is inside
        insideB = self.cell.EvaluatePosition(pB, self.closestPoint, subId, self.xi, dist, self.weights)
        if insideB == 1:
            # self.pB is inside cell
            xiEnd[:] = self.xi[:2]
        else:
            res2 = self.cell.IntersectWithLine(pA, pB, self.tol, t, self.intersectPt, self.xi, subId)
            if res2:
                self.cell.EvaluatePosition(self.intersectPt, self.closestPoint, subId, self.xi, dist, self.weights)
                xiEnd[:] = self.xi[:2]

        return (insideA == 1 or res1) and (insideB == 1 or res2)

###############################################################################

def testCellLineIntersector0():
    # standard, easy intersection

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
    assert abs(xiBeg[0] - 0.0) < 1.e-10 and abs(xiBeg[1] - 0.5) < 1.e-10
    assert abs(xiEnd[0] - 0.5) < 1.e-10 and abs(xiEnd[1] - 0.5) < 1.e-10

def testCellLineIntersector2Intersections():
    # 2 intersection

    lam0, the0 = 0., 0.
    lam1, the1 = 1., 0.
    lam2, the2 = 1., 1.
    lam3, the3 = 0., 1.

    lamA, theA = -0.5, 0.5
    lamB, theB = 1.5, 0.5

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
    assert abs(xiBeg[0] - 0.0) < 1.e-10 and abs(xiBeg[1] - 0.5) < 1.e-10
    assert abs(xiEnd[0] - 1.0) < 1.e-10 and abs(xiEnd[1] - 0.5) < 1.e-10

def testCellLineIntersectorInside():
    # segment is fully inside the cell

    lam0, the0 = 0., 0.
    lam1, the1 = 1., 0.
    lam2, the2 = 1., 1.
    lam3, the3 = 0., 1.

    lamA, theA = 0.1, 0.5
    lamB, theB = 0.9, 0.5

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
    assert abs(xiBeg[0] - 0.1) < 1.e-10 and abs(xiBeg[1] - 0.5) < 1.e-10
    assert abs(xiEnd[0] - 0.9) < 1.e-10 and abs(xiEnd[1] - 0.5) < 1.e-10

def testCellLineIntersectorOutside():
    # segment is fully outside of the cell

    lam0, the0 = 0., 0.
    lam1, the1 = 1., 0.
    lam2, the2 = 1., 1.
    lam3, the3 = 0., 1.

    lamA, theA = -2.0, -0.5
    lamB, theB = +3.0, -0.2

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
    assert not found 


if __name__ == '__main__':
    testCellLineIntersector0()
    testCellLineIntersector2Intersections()
    testCellLineIntersectorInside()
    testCellLineIntersectorOutside()


