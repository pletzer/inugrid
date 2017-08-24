import vtk
import numpy

class BasisFunctionIntegral:
    """
    Class to interpolate fluxes between faces
    """

    def __init__(self, xisA, xisB):
        """
        Constructor:
        @param xisA starting parametric position (2-vector)
        @param xisB ending parametric position (2-vector)
        """
        # difference
        dXis = xisB - xisA
        # average
        aXis = 0.5*(xisA + xisB)

        self.value = {0: dXis[0]*(1.0 - aXis[1]),
                      1: dXis[1]*aXis[0],
                      2: -dXis[0]*aXis[1],           # negative sign 
                      3: -dXis[1]*(1.0 - aXis[0]),}  # negative sign

    def __call__(self, index):
        """
        Evaluate the integral for basis function "index"
        @param index in the range(0, 4)
        @return integral
        """
        return self.value[index]

###############################################################################

def testBasisFunctionIntegral0():
    xiA, xiB = numpy.array((0., 0.)), numpy.array((1., 0.))
    bi = BasisFunctionIntegral(xiA, xiB)
    for i in range(4):
        print('testBasisFunctionIntegral0: i = {} weight = {}'.format(i, bi(i)))

def testBasisFunctionIntegral1():
    xiA, xiB = numpy.array((1., 0.)), numpy.array((1., 1.))
    bi = BasisFunctionIntegral(xiA, xiB)
    for i in range(4):
        print('testBasisFunctionIntegral1: i = {} weight = {}'.format(i, bi(i)))

def testBasisFunctionIntegral2():
    xiA, xiB = numpy.array((1., 1.)), numpy.array((0., 1.))
    bi = BasisFunctionIntegral(xiA, xiB)
    for i in range(4):
        print('testBasisFunctionIntegral2: i = {} weight = {}'.format(i, bi(i)))

def testBasisFunctionIntegral3():
    xiA, xiB = numpy.array((0., 1.)), numpy.array((0., 0.))
    bi = BasisFunctionIntegral(xiA, xiB)
    for i in range(4):
        print('testBasisFunctionIntegral3: i = {} weight = {}'.format(i, bi(i)))

if __name__ == '__main__':
    testBasisFunctionIntegral0()
    testBasisFunctionIntegral1()
    testBasisFunctionIntegral2()
    testBasisFunctionIntegral3()
