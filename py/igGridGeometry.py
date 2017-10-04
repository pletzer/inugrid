import vtk
import numpy
import math

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
        self.cellId = -1
        self.subId = vtk.mutable(0)
        self.x = numpy.zeros((3,), numpy.float64)
        # hexagon
        self.weights = numpy.zeros((8,), numpy.float64)
        self.tol2 = 1.e-12
        self.pcoords = numpy.zeros((3,), numpy.float64)

    def findCell(self, x):
        """
        Find the cell, the parametric coordinates and the weights for a given position
        @param x target position
        @return True if the cell was found
        """
        self.x = x
        self.cellId = self.cellLocator.FindCell(self.x, self.tol2, self.cell, self.pcoords, self.weights)
        return (self.cellId >= 0)

    def getX(self):
        """
        Interpolate the cell vertices to get the position
        @return vector
        @note call this after findCell(x)
        """
        self.cell.EvaluateLocation(self.subId, self.pcoords, self.x, self.weights)
        return self.x.copy()

    def getDXDXi(self, index):
        """
        Compute the partial derivative of position with respect to the parametric coordinates
        @param index 0, 1, or 2
        @return vector
        @note call this after findCell(x)
        """

        # note that VTK stores the data in inverse order to C
        # so we define a complementary index
        compIndex = 2 - index

        # save the paramatric position
        pcoordsBase = self.pcoords.copy()

        dXi = numpy.zeros((3,), numpy.float64)
        dPlus = 1.0 - self.pcoords[compIndex]
        dMnus = self.pcoords[compIndex] - 0.0

        dXi[compIndex] = dPlus
        self.pcoords = pcoordsBase + dXi
        xPlus = self.getX()

        dXi[compIndex] = dMnus
        self.pcoords = pcoordsBase - dXi
        xMnus = self.getX()

        # reset
        self.pcoords[:] = pcoordsBase

        return xPlus - xMnus

    def getJac(self):
        """
        Compute the Jacobian at the given pcoords position
        @return number
        @note call this after findCell(x)
        """
        dx0 = self.getDXDXi(0)
        dx1 = self.getDXDXi(1)
        dx2 = self.getDXDXi(2)
        return numpy.dot(dx0, numpy.cross(dx1, dx2))

    def getGradXi(self, index):
        """
        Compute the gradient of parametric coordinates with respect to the physical position
        @return vector
        @note call this after findCell(x)
        """
        i1 = (index + 1) % 3
        i2 = (index + 2) % 3
        dx1 = self.getDXDXi(i1)
        dx2 = self.getDXDXi(i2)
        jac = self.getJac()
        return numpy.cross(dx1, dx2) / jac

    def getGradX0CrossGradX1(self, indexPerp):
        """
        Compute the cross product of the gradients of position with respect to the parametric coordinates
        @param indexPerp the complementatary index, eg. (0, 1) -> 2, (1, 2) - > 0, (2, 0) -> 1
        @return vector
        @note call this after findCell(x)
        """
        jac = self.getJac()
        return self.getDXDXi(indexPerp) / jac

###############################################################################

def testCartesian():
    from igCartesianGrid import CartesianGrid

    ns = (10, 11, 12)
    ls = (1.0, 1.1, 1.2)
    cart =CartesianGrid(ns, ls)
    grid = cart.getUnstructuredGrid()

    geom = GridGeometry(grid)

    # can we find the cell?
    target = numpy.array([0.3, 0.4, 0.5])
    geom.findCell(target)
    print('cellId for target = {} is {} (pcoords = {})'.format(target, geom.cellId, geom.pcoords))

    # can we get back the position?
    target2 = geom.getX()
    print('position obtained from interpolation: {}'.format(target2))

    jac = geom.getJac()
    print('Jacobian: {} should be 0.1**3 = 1.e-3'.format(jac))

    grad0 = geom.getGradXi(0)
    grad1 = geom.getGradXi(1)
    grad2 = geom.getGradXi(2)
    print('grad0 = {} (should be [10, 0, 0])'.format(grad0))
    print('grad1 = {} (should be [0, 10, 0])'.format(grad1))
    print('grad2 = {} (should be [0, 0, 10])'.format(grad2))


def testSphere():
    from igLatLonElv import LatLonElv

    numLons, numLats, numElvs = 16, 8, 1
    radius = 1.0
    maxRelElv = 0.5
    sph = LatLonElv(numLons, numLats, numElvs, radius=radius, maxRelElv=maxRelElv)
    grid = sph.getUnstructuredGrid()

    geom = GridGeometry(grid)

    # can we find the cell?
    x, y, z = 1., 0., 0.
    r = math.sqrt(x**2 + y**2 + z**2)
    rho = math.sqrt(x**2 + y**2)
    the = math.acos(rho/r)
    lam = math.atan2(y, x)
    dLam = 2*numpy.pi/float(numLons)
    dThe = numpy.pi/float(numLats)
    dR = radius*maxRelElv/float(numElvs)

    lamHat = numpy.array([-y/r, x/r, 0.0])
    theHat = numpy.array([-z*math.cos(lam)/r, -z*math.sin(lam)/r, rho/r])
    rHat = numpy.array([x/r, y/r, z/r])

    target = numpy.array([x, y, z])
    geom.findCell(target)
    print('cellId for target = {} is {} (pcoords = {})'.format(target, geom.cellId, geom.pcoords))

    # can we get back the position?
    target2 = geom.getX()
    print('position obtained from interpolation: {}'.format(target2))

    dx0 = geom.getDXDXi(0)
    dx1 = geom.getDXDXi(1)
    dx2 = geom.getDXDXi(2)
    print('dx0 = {} (should be about {})'.format(dx0, r*lamHat * dLam))
    print('dx1 = {} (should be about {})'.format(dx1, r*math.cos(the)*theHat * dThe))
    print('dx2 = {} (should be about {})'.format(dx2, rHat * dR))

    jac = geom.getJac()
    print('Jacobian: {}'.format(jac))

    grad0 = geom.getGradXi(0)
    grad1 = geom.getGradXi(1)
    grad2 = geom.getGradXi(2)
    print('grad0 = {} grad1 = {} grad2 = {}'.format(grad0, grad1, grad2))
    

if __name__ == '__main__':
    #testCartesian()
    testSphere()
