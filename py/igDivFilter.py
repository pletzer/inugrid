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
        self.divData = numpy.zeros((numCells,), numpy.float64)
        self.vec1Data = numpy.zeros((numCells, 3), numpy.float64)
        self.vec2Data = numpy.zeros((numCells, 3), numpy.float64)

        self.grid.GetCellData().SetActiveScalars('cell_areas')
        cellAreas = self.grid.GetCellData().GetScalars()

        # iterate over the cells
        points = self.grid.GetPoints()
        for cellId in range(numCells):
            cell = self.grid.GetCell(cellId)
            ptIds = cell.GetPointIds()
            numPoints = ptIds.GetNumberOfIds()
            #compute the closed loop integral
            divVal = 0.0
            vec1Val = numpy.zeros((3,), numpy.float64)
            vec2Val = numpy.zeros((3,), numpy.float64)
            for i0 in range(numPoints):
                i1 = (i0 + 1) % numPoints
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
                vec2Val += 0.5 * fluxIntegral * numpy.cross(xx1 - xx0, rr)/length

            cellArea = cellAreas.GetComponent(cellId, 0)
            self.divData[cellId] = divVal / cellArea
            self.vec1Data[cellId, :] = vec1Val
            self.vec2Data[cellId, :] = vec2Val

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

        self.grid.GetCellData().AddArray(self.divArray)
        self.grid.GetCellData().AddArray(self.vec1Array)
        self.grid.GetCellData().AddArray(self.vec2Array)
