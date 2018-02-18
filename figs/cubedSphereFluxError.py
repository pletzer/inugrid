from matplotlib import pylab
pylab.rcParams['font.size'] = 20

import numpy

numCells = numpy.array([24, 96, 600, 2400, 15000, 60000])
error = numpy.array([-0.35, -0.11, -6.9e-4, -4.6e-4, -3.3e-4, -1.2e-4])

pylab.loglog(numCells, abs(error), 'o', markeredgecolor='k', markerfacecolor='w')
pylab.loglog(numCells, abs(error), 'b-')
# h ~ sqrt(1/numCells)
# error is ~ h^2
pylab.loglog(numCells, 10./numCells, 'k--')

pylab.title('Flux error on cubed sphere')
pylab.xlabel('number of horizontal grid cells')
pylab.ylabel('error')

pylab.show()