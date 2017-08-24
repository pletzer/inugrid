import vtk
import numpy

class FluxCalculator:

    def __init__(self, grid, integralFunction):
    	self.grid = grid
        self.fieldName = ''
        self.xyzLine = []
        self.totalFlux = 0.0

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
                    # iterate over the edges
                    for i0 in range(numPts):
                        i1 = (i0 + 1) % numPts
                        ptId0, ptId1 = ptIds.GetId(i0), ptIds.GetId(i1)
                        xyz0, xyz1 = self.grid.GetPoint(ptId0), self.grid.GetPoint(ptId1)
                        lam0, the0 = self._getLambdaThetaFromXYZ(xyz0)
                        lam1, the1 = self._getLambdaThetaFromXYZ(xyz1)
                        faceFlux = integralFunction(lam0, the0, lam1, the1)
                        self.totalFlux += faceFlux * self.basisFluxIntegral(xiBeg, xiEnd, i0)

        return self.totalFlux


    def _findCells(self, xyz0, xyz1):
    	pass

    def _addCellContributionToFlux(self, cellId):
    	pass

    def _addEdgeContributionToFlux(self, cellId, edgeId):
        pass

    def _getXYZFromLambdaTheta(self, lam, the):
        pass

    def _getLambdaThetaFromXYZ(self, xyz):
        pass