import vtk
import numpy

class GridGeometry:
    """
    A class that computes and stores all the metric coefficients in 
    an unstructured grid made of hexagons
    """

    def __init__(self, grid):
        """
        Constructor 
        @param grid vtkUnstructuredGrid instance
        """
        self.grid = grid
        self.cellLocator = vtk.vtkCellLocator()
        self.cellLocator.SetDataSet(self.grid)
        self.cellLocator.BuildLocator()
        self.cell = vtk.vtkGenericCell()
        self.subId = vtk.mutable(0)
        self.x = numpy.zeros((3,), numpy.float64)
        # hexagon
        self.weights = numpy.zeros((8,), numpy.float64)

    def findCell(self, x):
        """
        Find the cell and the parametric coordinates for a given position
        @param x target position
        @return cellId, pcoords
        """
        pcoords = numpy.zeros((3,), numpy.float64)
        cellId = self.cellLocator.FindCell(x, self.tol2, self.cell, pcoords, self.weights)
        return cellId, pcoords

    def getX(self, cellId, pcoords):
        """
        Interpolate the cell vertices to get the position
        @param cellId cell index
        @param pcoords parametric coordinates
        """
        self.grid.EvaluateLocation(self.subId, pcoords, self.x, self.weights)
        return self.x

    def getDXDXi(self, cellId, pcoords, index):
        """
        Compute the partial derivative of position with respect to the parametric coordinates
        @param cellId cell index
        @param pcoords parametric coordinates
        @param index 0, 1, or 2
        @return vector
        """
        dXi = numpy.zeros((3,), numpy.float64)
        dPlus = 1.0 - pcoords[index]
        dMnus = pcoords[index] - 0.0

        dXi[index] = dPlus
        pcoordsPlus = pcoords + dXi
        dXi[index] = dMnus
        pcoordsMnus = pcoords - dXi
        xPlus = self.getX(cellId, pcoordsPlus)
        xMnus = self.getX(cellId, pcoordsMnus)
        return xPlus - xMnus

    def getJac(self, cellId, pcoords):
        """
        Compute the Jacobian at the given pcoords position
        @param cellId cell index
        @param pcoords parametric coordinates
        @return number
        """
        dx0 = self.getDXDXi(cellId, 0)
        dx1 = self.getDXDXi(cellId, 1)
        dx2 = self.getDXDXi(cellId, 2)
        return numpy.dot(dx0, numpy.cross(dx1, dx2))

    def getGradX(self, cellId, pcoords, index):
        """
        Compute the gradient of position with respect to the parametric coordinates
        @param cellId cell index
        @param pcoords parametric coordinates
        @return vector
        """
        i1 = (index + 1) % 3
        i2 = (index + 2) % 3
        dx1 = self.getDXDXi(cellId, i1)
        dx2 = self.getDXDXi(cellId, i2)
        jac = self.getJac(cellid, pcoords)
        return numpy.cross(dx1, dx2) / jac

    def getGradX0CrossGradX1(self, cellId, pcoords, indexPerp):
        """
        Compute the cross product of the gradients of position with respect to the parametric coordinates
        @param cellId cell index
        @param pcoords parametric coordinates
        @param indexPerp the complementatary index, eg. (0, 1) -> 2, (1, 2) - > 0, (2, 0) -> 1
        @return vector
        """
        jac = self.getJac(cellid, pcoords)
        return self.getDXDXi(cellId, pcoords, indexPerp) / jac

###############################################################################

def testSphere():
    

if __name__ == '__main__': test()




