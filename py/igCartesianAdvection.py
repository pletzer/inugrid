import numpy
import scipy
import vtk

class CartesianAdvection:

    def __init__(self):
        """
        Constructor
        """
        self.grid = None
        self.psiFunc = None
        self.x0 = []
        self.cellLocator = vtk.vtkCellLocator()
        self.tol2 = 1.e-10
        self.cell = vtk.vtkGenericCell()
        self.pcoords = numpy.zeros((3,), numpy.float64)
        self.weights = numpy.zeros((8,), numpy.float64)

    def setGrid(self, grid):
        """
        Set the grid
        @param grid vtkUnstructuredGrid instance
        """
        self.grid = grid
        self.cellLocator.SetDataSet(self.grid)
        self.cellLocator.BuildLocator()

    def setStreamFunction(self, function):
        """
        Set the stream function
        @param function scalar function of position
        """
        self.psiFunc = function

    def setPosition(self, pos):
        """
        Set the target position
        @param pos numpy array of size 3
        """
        self.x0 = pos

    def advance(self, dt):
        """
        Advance the point by dt
        @param dt time interval
        """
        
        # find the cell that contains the point

        # advect the point along the velocity field until the point intersect 
        # cell face

    def tendency(self, x):
        """
        Compute the tendency
        @param x target point
        """
        cellId = self.cellLocator.FindCell(x, self.tol2, self.cell, self.pcoords, self.weights)
        if cellId == -1:
            print('ERROR: could not find cell')
            return []

        # compute the stream function at the vertices 
        points = self.grid.GetPoints()
        pointIds = self.cell.GetPointIds()
        numPoints = pointIds.GetNumberOfIds()
        psis = numpy.zeros((numPoints,), numpy.float64)
        for i in pointIds.GetNumberOfIds():
            pointId = pointIds.GetId(i)
            xyz = points.GetPoint(pointId)
            # evalaute the stream function
            psis.append(self.psiFunc(xyz))

        # compute the velocity by interpolating the edge values using the face basis functions
        xi0, xi1, xi2 = self.pcoords
        vx = (1.0 - xi1)*(psis[1] - psis[0]) + xi1*(psi[2] - psis[3])
        vy = (1.0 - xi0)*(psis[3] - psis[0]) + xi0*(psi[2] - psis[1])
        vz = 0.0
        return numpy.array([vx, vy, vz])


###############################################################################

def main():

    def psiFunction(pos):
        x, y = pos
        return numpy.cos(x)**2 + y**2

    nx, ny = 10, 30
    lx, ly = 2*numpy.pi, 1.0
    cart = CartesianGrid(nx, ny, lx, ly)
    grid = cart.getUnstructuredGrid()

    adv = CartesianGrid()
    adv.setGrid(grid)
    pos = numpy.array([numpy.pi, 0.5, 0.0]) # vtk wants 3d
    adv.setPosition(pos)
    adv.setStreamFunction(psiFunction)

    points = []
    dt = 0.1
    for i in range(nt):
        points.append(adv.advance(dt))

    print points

if __name__ == '__main__': main()

