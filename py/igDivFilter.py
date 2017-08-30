import vtk
import numpy
from math import sqrt, cos, sin, pi, atan2


def getLambdaTheta(xx):
    """
    Get lon/lat in radiants from the Cartesian positions
    @param xx array of x, y and z
    @return lambda, theta
    """
    x, y, z = xx
    lam = atan2(y, x)
    rho = sqrt(x*x + y*y)
    the = atan2(z, rho)
    return lam, the


class DivFilter:

    def __init__(self, grid):
        self.grid = grid

    def applyIntegrals(self, integralFunction):

        numCells = self.grid.GetNumberOfCells()
        numPoints = self.grid.GetNumberOfPoints()
        self.divData = numpy.zeros((numCells,), numpy.float64)
        self.vec1Data = numpy.zeros((numCells, 3), numpy.float64)
        self.vec2Data = numpy.zeros((numCells, 3), numpy.float64)
        self.vec2DataApprox = numpy.zeros((numPoints, 3), numpy.float64)

        self.grid.GetCellData().SetActiveScalars('cell_areas')
        cellAreas = self.grid.GetCellData().GetScalars()

        # iterate over the cells
        points = self.grid.GetPoints()
        for cellId in range(numCells):
            cell = self.grid.GetCell(cellId)
            ptIds = cell.GetPointIds()
            npts = ptIds.GetNumberOfIds()
            #compute the closed loop integral
            divVal = 0.0
            vec1Val = numpy.zeros((3,), numpy.float64)
            vec2Val = numpy.zeros((3,), numpy.float64)

            # edge vector field attached to nodes. There are three ways of doing this, hence three
            # averages
            #
            #  |      |  __o  o__
            #  o__  __o    |  |   
            #   0    1    2     3
            vec2ValApprox0 = numpy.zeros((3,), numpy.float64)
            vec2ValApprox1 = numpy.zeros((3,), numpy.float64)
            vec2ValApprox2 = numpy.zeros((3,), numpy.float64)
            vec2ValApprox3 = numpy.zeros((3,), numpy.float64)

            for i0 in range(npts):
                i1 = (i0 + 1) % npts
                ptId0, ptId1 = ptIds.GetId(i0), ptIds.GetId(i1)
                xx0 = numpy.array(points.GetPoint(ptId0))
                xx1 = numpy.array(points.GetPoint(ptId1))
                rr = 0.5*(xx0 + xx1)

                # retreat by a tiny bit in order to capture multivalued jumps 
                #x1 = x0 + (x1 - x0)*(1. - EPS)
                #y1 = y0 + (y1 - y0)*(1. - EPS)
                #z1 = z0 + (z1 - z0)*(1. - EPS)

                lam0, the0 = getLambdaTheta(xx0)
                lam1, the1 = getLambdaTheta(xx1)

                fluxIntegral = integralFunction(lam0, lam1, the0, the1)
                divVal += fluxIntegral

                # compute the vector (1 and 2 forms) at cell centre
                length = numpy.sqrt(numpy.dot(xx1 - xx0, xx1 - xx0))
                vec1Val += 0.5 * fluxIntegral * (xx1 - xx0) / length
                vecCross = numpy.cross(xx1 - xx0, rr)/length
                vec2Val += 0.5 * fluxIntegral * vecCross

                vec2ValApprox0 += (i0 == 0 or i0 == 3) * fluxIntegral * vecCross
                vec2ValApprox1 += (i0 == 0 or i0 == 1) * fluxIntegral * vecCross
                vec2ValApprox2 += (i0 == 1 or i0 == 2) * fluxIntegral * vecCross
                vec2ValApprox2 += (i0 == 2 or i0 == 3) * fluxIntegral * vecCross

            cellArea = cellAreas.GetComponent(cellId, 0)
            # put on cells
            self.divData[cellId] = divVal / cellArea
            self.vec1Data[cellId, :] = vec1Val
            self.vec2Data[cellId, :] = vec2Val

            # put vector on vertices. We go through every cell and assign a vector for each the points
            # of the cell. Note that points will get multiple values from all the cells that sahre the 
            # node. The final value assigned to the node will depend on the order we loop through cells. 
            ptId0, ptId1, ptId2, ptId3 = ptIds.GetId(0), ptIds.GetId(1), ptIds.GetId(2), ptIds.GetId(3)
            self.vec2DataApprox[ptId0, :] = vec2ValApprox0
            self.vec2DataApprox[ptId1, :] = vec2ValApprox1
            self.vec2DataApprox[ptId2, :] = vec2ValApprox2
            self.vec2DataApprox[ptId3, :] = vec2ValApprox3

        # attach cell centred values to the grid
        self.divArray = vtk.vtkDoubleArray()
        self.divArray.SetName('integral_star_d_phi_over_area')
        self.divArray.SetNumberOfComponents(1)
        self.divArray.SetNumberOfTuples(numCells)
        save = 1
        self.divArray.SetVoidArray(self.divData, numCells, save)

        # attach cell centred vector field
        self.vec1Array = vtk.vtkDoubleArray()
        self.vec1Array.SetName('1-form')
        self.vec1Array.SetNumberOfComponents(3)
        self.vec1Array.SetNumberOfTuples(numCells)
        self.vec1Array.SetVoidArray(self.vec1Data, numCells*3, save)

        self.vec2Array = vtk.vtkDoubleArray()
        self.vec2Array.SetName('2-form')
        self.vec2Array.SetNumberOfComponents(3)
        self.vec2Array.SetNumberOfTuples(numCells)
        self.vec2Array.SetVoidArray(self.vec2Data, numCells*3, save)

        self.vec2ArrayApprox = vtk.vtkDoubleArray()
        self.vec2ArrayApprox.SetName('2-form-approx')
        self.vec2ArrayApprox.SetNumberOfComponents(3)
        self.vec2ArrayApprox.SetNumberOfTuples(numPoints)
        self.vec2ArrayApprox.SetVoidArray(self.vec2DataApprox, numPoints*3, save)

        self.grid.GetCellData().AddArray(self.divArray)
        self.grid.GetCellData().AddArray(self.vec1Array)
        self.grid.GetCellData().AddArray(self.vec2Array)
        self.grid.GetPointData().AddArray(self.vec2ArrayApprox)
