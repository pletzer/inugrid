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
ls = (2*numpy.pi, 3.0, 1.0)
origin = (0., -1.5, -0.5)
hs = (ls[0]/float(ns[0]), ls[1]/float(ns[1]), ls[2]/float(ns[2]))
cart = CartesianGrid(ns, ls, origin)
grid = cart.getUnstructuredGrid()

geom = GridGeometry(grid)


def velocityExact(x, *args):
    """
    Exact velocity field
    """
    # -d psi/dx
    vy = -2. * numpy.sin(x[0]) * numpy.cos(x[0])
    # d psi/dy
    vx = 2. * x[1]
    return numpy.array([vx, vy, 0.])

def velocityNodal(x, *args):
    """
    Velocity field linear linterpolated from the exact values at the nodes
    """

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
    velNodes = [velocityExact(v) for v in verts]

    # bilinear interpolation of velocity at vertices
    w0 = (1. - xis[0]) * (1. - xis[1])
    w1 = (xis[0] - 0.) * (1. - xis[1])
    w2 = (xis[0] - 0.) * (xis[1] - 0.)
    w3 = (1. - xis[0]) * (xis[1] - 0.)

    res = w0*velNodes[0] + w1*velNodes[1] + w2*velNodes[2] + w3*velNodes[3]

    return res

# stream function
def streamFuncExact(x):
    return numpy.sin(x[0])**2 + x[1]**2

def velocityFace(x, *args):
    """
    Velocity field interpolated using the face centred basis functions and exact 
    fluxes on cell faces
    """
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
    psis = [streamFuncExact(v) for v in verts]

    vx = (1. - xis[0])*(psis[3] - psis[0]) + xis[0]*(psis[2] - psis[1])
    vx /= hs[1]

    vy = (1. - xis[1])*(psis[0] - psis[1]) + xis[1]*(psis[3] - psis[2])
    vy /= hs[0]

    res[0:2] = vx, vy
    return res

def velocityFaceAsNodal(x, *args):
    """
    Use the face basis functions to approximate the nodal vector field
    """
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
    # compute the velocity at the nodes from the face centred basis functions
    velNodes = [velocityFace(v) for v in verts]

    # bilinear interpolation of velocity at vertices
    w0 = (1. - xis[0]) * (1. - xis[1])
    w1 = (xis[0] - 0.) * (1. - xis[1])
    w2 = (xis[0] - 0.) * (xis[1] - 0.)
    w3 = (1. - xis[0]) * (xis[1] - 0.)

    res = w0*velNodes[0] + w1*velNodes[1] + w2*velNodes[2] + w3*velNodes[3]

    return res

def velocityFaceAsNodal2(x, *args):
    """
    Use the face basis functions to approximate the nodal vector field
    """
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


    dx = numpy.array([hs[0], 0., 0.])
    dy = numpy.array([0., hs[1], 0.])
    velNodes = []
    for v in verts:
        psi0 = streamFuncExact(v)
        dpsiy = (streamFuncExact(v + dy) - psi0)/hs[1]
        dpsix = (streamFuncExact(v + dx) - psi0)/hs[0]
        vx = + dpsiy
        vy = - dpsix
        velNodes.append(numpy.array([vx, vy, 0.]))

    # bilinear interpolation of velocity at vertices
    w0 = (1. - xis[0]) * (1. - xis[1])
    w1 = (xis[0] - 0.) * (1. - xis[1])
    w2 = (xis[0] - 0.) * (xis[1] - 0.)
    w3 = (1. - xis[0]) * (xis[1] - 0.)

    res = w0*velNodes[0] + w1*velNodes[1] + w2*velNodes[2] + w3*velNodes[3]

    return res


# initial condition
x = numpy.array([numpy.pi/2. + 0.1, 0.015, 0.0])

# times
ts = numpy.linspace(0., 100., 1001)

# solve
solExact = odeint(velocityExact, x, ts)
solNodal = odeint(velocityFaceAsNodal, x, ts)
solFace = odeint(velocityFace, x, ts)

# plot
fig = pylab.figure()
ax = fig.add_subplot(111)

pylab.plot(solExact[:, 0], solExact[:, 1], 'g-')
pylab.plot(solNodal[:, 0], solNodal[:, 1], 'r-')
pylab.plot(solFace[:, 0], solFace[:, 1], 'b-')
pylab.legend(['exact', 'nodal', 'face'])
pylab.plot([solExact[0, 0]], [solExact[0, 1]], 'k*')

x = numpy.linspace(origin[0], origin[0] + ls[0], 101)
y = numpy.linspace(origin[1], origin[1] + ls[1], 101)
xx, yy = numpy.meshgrid(x, y)
pylab.contour(x, y, numpy.sin(xx)**2 + yy**2, colors='k')
ax.set_aspect('equal')

#pylab.title('igStreamlineNodalExample')

pylab.show()
