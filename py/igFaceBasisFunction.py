import numpy
import vtk

class FaceBasisFunction:

    def __init__(self, theta0, dLambda, dTheta, dRadius, n1, n2, n3,):

        ntot = (n1 + 1) * (n2 + 1) * (n3 + 1)
        dXi1 = 1.0/float(n1)
        dXi2 = 1.0/float(n2)
        dXi3 = 1.0/float(n3)

        self.ptArray = vtk.vtkDoubleArray()
        self.ptArray.SetNumberOfComponents(3)
        self.ptArray.SetNumberOfTuples(ntot)

        self.data = vtk.vtkDoubleArray()
        self.data.SetNumberOfComponents(3) # vector
        self.data.SetNumberOfTuples(ntot)

        index = 0
        for i3 in range(n3 + 1):
            for i2 in range(n2 + 1):
                for i1 in range(n1 + 1):

                    lam = 0.0 + i3*dLambda
                    the = theta0 + i2*dTheta
                    rad = 1.0 + i1*dRadius
                    x = rad * numpy.cos(the) * numpy.cos(lam)
                    y = rad * numpy.cos(the) * numpy.cos(lam)
                    z = rad * numpy.sin(the)
                    self.ptArray.SetTuple(index, (x, y, z))

                    xi1 = 0.0 + i1*dXi1
                    xi2 = 0.0 + i2*dXi2
                    xi3 = 0.0 + i3*dXi3
                    basisFunc = - 8.*(xi2 - 0.75)*(xi3 - 0.25)*(xi1 - 0.5)
                    rho = numpy.sqrt(x**2 + y**2)
                    thetaHat = numpy.array([-z*numpy.cos(lam), -z*numpy.sin(lam), rho]) / rad
                    rHat = numpy.array([x, y, z]) / rad
                    basisVec = numpy.cross(thetaHat, rHat) / rad
                    vx, vy, vz = basisFunc*basisVec
                    self.data.SetTuple(index, (vx, vy, vz))

                    index += 1

        self.pts = vtk.vtkPoints()
        self.pts.SetData(self.ptArray)

        self.grid = vtk.vtkStructuredGrid()
        self.SetPoints(self.pts)
        self.grid.GetPointData().SetScalars(self.data)

    def save(self, filename):
        writer = vtk.vtkStructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.grid)
        writer.Update()

###############################################################################
def main():
    fb = FaceBasisFunction(theta0=numpy.pi/5., dLambda=numpy.pi/8., dTheta=numpy.pi/8., dRadius=0.2, n1=10, n2=10, n3=10)
    fb.save('faceBasis.vtk')
