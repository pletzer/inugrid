import numpy
from scipy.integrate import odeint
import vtk

class Advection:
    """
    A class to advect/integrate trajectories
    """

    def __init__(self, grid):
        """
        Constructor
        """
        self.grid = grid
        self.tendency = None
        self.x0 = []


    def setVectorFieldFunction(self, vectorFunction):
        """
        Set the vector field function
        @param x position
        """
        self.tendency = vectorFunction


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
        sol = odeint(self.tendency, self.x0, dt, args=(self.grid,))
        print sol
        self.x0[:] = newX


###############################################################################

def testSimple():

    from igCartesianGrid import CartesianGrid

    def vectorFieldFunction(x, *args):
        print 'calling vectorFieldFunction with x = ', x, ' and args = ', args
        return numpy.array([1., 2., 3.])

    ns = (10, 11, 12)
    ls = (1.0, 1.1, 1.2)
    cart = CartesianGrid(ns, ls)
    grid = cart.getUnstructuredGrid()

    adv = Advection(grid)
    adv.setPosition(numpy.array([0., 0., 0.]))
    adv.setVectorFieldFunction(vectorFieldFunction)

    points = []
    dt = 0.1
    nt = 10
    for i in range(nt):
        points.append(adv.advance(dt))

    print points

if __name__ == '__main__': 
    testSimple()

