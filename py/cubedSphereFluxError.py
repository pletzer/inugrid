from matplotlib import pylab
import matplotlib
import numpy

font = {'family' : 'normal',
        'size'   : 20}

matplotlib.rc('font', **font)

ns = [2, 5, 10, 20, 50,]
ntots = [6*n**2 for n in ns]

# 100 segments
error100 = [0.110302180825, 0.0136772740688, 0.00426548330072, 4.35606850407e-05, 0.000180361841866]
quadratic = 0.4*numpy.array([1./float(nt) for nt in ntots])

# total number of cells is 6*n**2
pylab.loglog(ns, error100, 'm-',)
pylab.loglog(ns, quadratic, 'k--')
pylab.legend(['100 segments', '1/N^2'])
pylab.plot(ns, error100, 'kx')
pylab.ylabel('error')
pylab.xlabel('tot num cells')
pylab.title('Cubed sphere flux error')
pylab.show()
