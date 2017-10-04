import numpy
from scipy.integrate import odeint
from matplotlib import pylab
from igStreamVectorField import StreamVectorField
from igCartesianGrid import CartesianGrid
from igGridGeometry import GridGeometry

"""
Stream function y^2 + sin(x)^2
v = d psi ^ dz

v is interpolated using bilinear basis functions
"""

#  create the grid
ns = (60, 30, 1)
ls = (2*numpy.pi, 2.0, 1.0)
origin = (0., -1.0, -0.5)
hs = (ls[0]/float(ns[0]), ls[1]/float(ns[1]), ls[2]/float(ns[2]))
cart = CartesianGrid(ns, ls, origin)
grid = cart.getUnstructuredGrid()

geom = GridGeometry(grid)


def velocityNodal(x):
    # -d psi/dx
    vy = -2. * numpy.sin(x[0]) * numpy.cos(x[0])
    # d psi/dy
    vx = 2. * x[1]
    return numpy.array([vx, vy, 0.])

def velocity(x, t):

    res = numpy.zeros((3,), numpy.float64)

    # find the cell the contains the point
    found = geom.findCell(x)
    if not found:
        # return zero veolocity if outside the domain
        return res

    # parametric coordinates are stored in inverse order
    xis = (geom.pcoords[2], geom.pcoords[1])

    # vertices of the cell
    pts = geom.cell.GetPoints()

    # depends on the indexing of the hex in VTK
    verts = [pts.GetPoint(i) for i in (0, 4, 7, 3)]
    velNodes = [velocityNodal(v) for v in verts]

    # bilinear interpolation of velocity at vertices
    w0 = (1. - xis[0]) * (1. - xis[1])
    w1 = (xis[0] - 0.) * (1. - xis[1])
    w2 = (xis[0] - 0.) * (xis[1] - 0.)
    w3 = (1. - xis[0]) * (xis[1] - 0.)

    res = w0*velNodes[0] + w1*velNodes[1] + w2*velNodes[2] + w3*velNodes[3]

    return res



x = numpy.array([0.1, 0.2, 0.0]) #[numpy.pi/2. + 0.01, 0., 0.])
ts = numpy.linspace(0., 1.0, 11) #numpy.linspace(0., 1000., 10001)
sol = odeint(velocity, x, ts)

pylab.plot(sol[:, 0], sol[:, 1], 'b-')
pylab.title('igStreamlineNodalExample')

pylab.show()
