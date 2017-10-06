import numpy
from scipy.integrate import odeint
from matplotlib import pylab
from igStreamVectorField import StreamVectorField
from igCartesianGrid import CartesianGrid
from igGridGeometry import GridGeometry

"""
Stream function y^2 + sin(x)^2
v = d psi ^ dz
"""

#  create the grid
ns = (2*12, 2*6, 1) #(2*30, 2*15, 1)
ls = (6., 3., 1.) #(2*numpy.pi, 3.0, 1.0)
origin = (0, -1.5, -0.5) # (0., -2.0, -0.5)
hs = (ls[0]/float(ns[0]), ls[1]/float(ns[1]), ls[2]/float(ns[2]))
cart = CartesianGrid(ns, ls, origin)
grid = cart.getUnstructuredGrid()

geom = GridGeometry(grid)

angle = 0.1 # 0.3 #numpy.pi/6.
cos_angle = numpy.cos(angle)
sin_angle = numpy.sin(angle)

# stream function
def streamFuncExact(x):
    # apply rotation of the coordinates to make it more interesting
    xp = cos_angle * x[0] - sin_angle * x[1]
    yp = sin_angle * x[0] + cos_angle * x[1]
    return numpy.sin(xp)**2 + yp**2


def velocityExact(x, *args):
    """
    Exact velocity field
    """
    #zero = numpy.zeros((3,), numpy.float64)
    #for i in range(3):
    #    if x[i] < origin[i] or x[i] > origin[i] + ls[i]:
    #        return zero

    xp = cos_angle * x[0] - sin_angle * x[1]
    yp = sin_angle * x[0] + cos_angle * x[1]

    # -d psi/dx'
    vyp = -2. * numpy.sin(xp) * numpy.cos(xp)

    # d psi/dy'
    vxp = 2. * yp

    vy = vyp * cos_angle - vxp * sin_angle
    vx = vxp * cos_angle + vyp * sin_angle

    return numpy.array([vx, vy, 0.])

def velocityNodal(x, *args):
    """
    Velocity field linearly linterpolated using the exact values at the nodes
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

    # stream function at the vertices
    psis = [streamFuncExact(v) for v in verts]

    # integrate the flux and divide by cell width to get an approximation of the derivative
    vx = (1. - xis[0])*(psis[3] - psis[0]) + xis[0]*(psis[2] - psis[1])
    vx /= hs[1]

    vy = (1. - xis[1])*(psis[0] - psis[1]) + xis[1]*(psis[3] - psis[2])
    vy /= hs[0]

    res[0:2] = vx, vy
    #print 'vx, hs, xi0, dpsi_west, dpsi_east = ', vx, hs, xis[0], psis[3] - psis[0], psis[2] - psis[1]
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

    # vertices on the lower z plane
    verts = [pts.GetPoint(i) for i in (0, 4, 7, 3)]

    # displacement vectors in the x and y directions
    dx = numpy.array([hs[0], 0., 0.])
    dy = numpy.array([0., hs[1], 0.])

    # compute the velocity by computing finite differences in the 
    # downwind direction. Assumes a uniform mesh!
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

# time steps
ts = numpy.linspace(0., 10.0, 101)

# plot
fig = pylab.figure()
ax = fig.add_subplot(111)

# vary initial conditions
ystart = 0.1
for xstart in numpy.linspace(1., 1.4, 3):

    # initial condition
    x = numpy.array([xstart, ystart, 0.]) 


    # solve
    solExact = odeint(velocityExact, x, ts)
    solNodal = odeint(velocityFaceAsNodal2, x, ts)
    solFace = odeint(velocityFace, x, ts)


    pylab.plot(solExact[:, 0], solExact[:, 1], 'g-')
    pylab.plot(solNodal[:, 0], solNodal[:, 1], 'r-')
    pylab.plot(solFace[:, 0], solFace[:, 1], 'b-')
    pylab.plot([solExact[0, 0]], [solExact[0, 1]], 'k*', markersize=10)


pylab.legend(['exact', 'nodal', 'face'])

# plot the stream function
x = numpy.linspace(origin[0], origin[0] + ls[0], 101)
y = numpy.linspace(origin[1], origin[1] + ls[1], 101)
xx, yy = numpy.meshgrid(x, y)
xxp = cos_angle*xx - sin_angle*yy
yyp = sin_angle*xx + cos_angle*yy
pylab.contour(x, y, numpy.sin(xxp)**2 + yyp**2, colors='k')
ax.set_aspect('equal')

#pylab.title('igStreamlineNodalExample')

pylab.show()
