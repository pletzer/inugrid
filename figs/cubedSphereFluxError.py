from matplotlib import pylab
pylab.rcParams['font.size'] = 20

import numpy

# 20 segements

numCells = numpy.array([24, 150, 600, 2400, 15000, 60000, 240000,])
error = numpy.array([0.11, 0.014, 0.0043, 4.4e-5, 1.8e-4, 1.2e-5, 5.5e-6])

pylab.loglog(numCells, abs(error), 'o', markeredgecolor='k', markerfacecolor='w')
pylab.loglog(numCells, abs(error), 'b-')
# h ~ sqrt(1/numCells)
# error is ~ h^2
pylab.loglog(numCells, 10./numCells, 'k--')

pylab.title('Flux error on cubed sphere')
pylab.xlabel('number of horizontal grid cells')
pylab.ylabel('error')

pylab.show()
