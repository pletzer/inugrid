import numpy
import vtk

class CartesianGrid:

    def __init__(self, nx, ny, lx, ly):

        # this is tp convert from structured to unstructured grid
        self.appendGrids = vtk.vtkAppendFilter()

        # create grid
        nx1, ny1 = nx + 1, ny + 1
        xs = numpy.linspace(0., lx, nx1)
        ys = numpy.linspace(0., ly, ny1)
        xxs, yys = numpy.meshgrid(xs, ys)

        xxs = xxs.flat
        yys = yys.flat

        # coordinates
        self.xyz = numpy.zeros((nx1*ny1, 3), numpy.float64)
        self.xyz[:, 0] = xxs
        self.xyz[:, 1] = yys
        self.xyz[:, 2] = 0.0

        # create the VTK unstructured grid
        self.vxyz = vtk.vtkDoubleArray()
        self.vxyz.SetNumberOfComponents(3)
        ntot = nx1 * ny1
        self.vxyz.SetNumberOfTuples(ntot)
        self.vxyz.SetVoidArray(self.xyz, 3*ntot, 1)

        self.pts = vtk.vtkPoints()
        self.pts.SetNumberOfPoints(ntot)
        self.pts.SetData(self.vxyz)

        self.sgrid = vtk.vtkStructuredGrid()
        self.sgrid.SetDimensions(nx1, ny1, 1)
        self.sgrid.SetPoints(self.pts)

        self.appendGrids.AddInputData(self.sgrid)
        self.appendGrids.Update()
        self.grid = self.appendGrids.GetOutput()

    def getUnstructuredGrid(self):
        return self.grid


    def save(self, filename):
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.grid)
        writer.Update()


    def show(self):
        gridMapper = vtk.vtkDataSetMapper()
        gridMapper.SetInputData(self.grid)

        gridActor = vtk.vtkActor()
        gridActor.SetMapper(gridMapper)
        #gridActor.GetProperty().SetColor(93./255., 173./255., 226./255.)

        light = vtk.vtkLight()
        light.SetFocalPoint(0., 0., 0)
        light.SetPosition(0.5, 0.3, 1.) #(1., 1., 0.)
        #light.SetSpecularColor(0.8, 0.2, 0.0)
        #light.SetDiffuseColor(0., 0.2, 0.8)
        light.SetConeAngle(0.2)
        light.SetIntensity(1.0)

        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        # add the actors to the renderer, set the background and size
        ren.AddActor(gridActor)
        ren.AddLight(light)
        ren.SetBackground(1, 1, 1)
        renWin.SetSize(640, 640)
        iren.Initialize()
        renWin.Render()
        iren.Start()


#############################################################################
def test():
    nx, ny = 10, 11
    lx, ly = 2*numpy.pi, 1.0
    cart = CartesianGrid(nx, ny, lx, ly)
    grid = cart.getUnstructuredGrid()
    cart.save('cart.vtk')
    cart.show()

if __name__ == '__main__':
    test()

