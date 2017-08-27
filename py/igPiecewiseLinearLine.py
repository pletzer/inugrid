import vtk
import numpy

class PiecewiseLinearLine:

	def __init__(self, lambdaFunction, thetaFunction, nt=11, radius=1.0):
		"""
		Create a polyline 
		@param lambdaFunction lambda (longitude) function of t (0 <= t <= 1)
		@param thetaFunction theta (latitude) functuon of t (0 <= t <= 1)
		"""
        self.ptData = vtk.DoubleArray()
        self.pts = vtk.vtkPoints()
		self.line = vtk.vtkPolyLine()
        self.grid = vtk.vtkUnstructuredGrid()

		ts = numpy.linspace(0., 1., nt)
		lams = lambdaFunction(ts)
		thes = thetaFunction(ts)
		xs = radius * numpy.cos(thes) * numpy.cos(lams)
		ys = radius * numpy.cos(thes) * numpy.sin(lams)
		zs = radius * numpy.sin(thes)
        self.xyz = numpy.zeros((nt, 3), numpy.float64)
        self.xyz[:, 0] = xs
        self.xyz[:, 1] = ys
        self.xyz[:, 2] = zs

        self.ptData.SetNumberOfComponents(3)
        self.ptData.SetNumberOfTuples(nt)
        self.ptData.SetVoidArray(self.xyz, nt*3, 1)

        self.pts.SetNumberOfPoints(nt)
        self.pts.SetData(self.ptData)

        ptIds = self.line.GetPointIds()
        ptIds.SetNumberOfIds(nt)
        for i in range(ptIds.GetNumberOfIds()):
            segs.SetPointId(i, i)

        self.grid.SetPoints(self.pts)
        self.grid.Allocate(1, 1)
        # one cell
        self.grid.InsertNextCell(self.line)


	def save(self, filename):
		"""
        Write unstructured grid to file
        @param filename file name
        """
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.grid)
        writer.Update()

###############################################################################
def main():
    lam0, the0 = 165.0, -40.0
    a, b = 20.0, 15.0
    def lamFunc(ts):
        return lam0 + a*numpy.cos(ts)
    def theFunc(ts):
        return the0 + b*numpy.sin(ts)

    pl = PiecewiseLinearLine(lamFunc, theFunc, nt=21)
    pl.save('line.vtk')

    
if __name__ == '__main__':
    main()


