import numpy
from scipy.integrate import odeint
from matplotlib import pylab


"""
Stream function y^2 + sin(x)^2
v = d psi ^ dz
"""

def velocity(x, t):

	# -d psi/dx
    vy = -2. * numpy.sin(x[0]) * numpy.cos(x[0])

    # d psi/dy
    vx = 2. * x[1]

    #print x, vx, vy
    return numpy.array([vx, vy, 0.])



x = numpy.array([0.1, 0.2, 0.0]) # [numpy.pi/2. + 0.01, 0., 0.])
ts = numpy.linspace(0., 1., 11) #numpy.linspace(0., 1000., 10001)
sol = odeint(velocity, x, ts)


pylab.plot(sol[:, 0], sol[:, 1], 'b-')
pylab.title('igStreamlineExample1')

pylab.show()
