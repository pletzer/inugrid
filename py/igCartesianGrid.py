import numpy
import vtk

class CartesianGrid:

    def __init__(self, ns, ls, origin=(0., 0., 0.)):
        """
        Create Cartesian grid
        @param ns number of cells in x, y, and z (3-tuple)
        @param ls domain sizes in x, y, and z
        """

        # to convert from structured to unstructured grid
        self.appendGrids = vtk.vtkAppendFilter()

        # create grid
        lx, ly, lz = ls
        nx, ny, nz = ns
        nx1, ny1, nz1 = nx + 1, ny + 1, nz + 1
        ntot = nx1 * ny1 * nz1
        xs = numpy.linspace(origin[0], origin[0] + lx, nx1)
        ys = numpy.linspace(origin[1], origin[1] + ly, ny1)
        zs = numpy.linspace(origin[2], origin[2] + lz, nz1)
        xxs, yys, zzs = numpy.meshgrid(xs, ys, zs, sparse=False, indexing='ij')

        # coordinates
        self.xyz = numpy.zeros((ntot, 3), numpy.float64)
        self.xyz[:, 0] = xxs.flat
        self.xyz[:, 1] = yys.flat
        self.xyz[:, 2] = zzs.flat

        # create the VTK unstructured grid
        self.vxyz = vtk.vtkDoubleArray()
        self.vxyz.SetNumberOfComponents(3)
        self.vxyz.SetNumberOfTuples(ntot)
        self.vxyz.SetVoidArray(self.xyz, 3*ntot, 1)

        self.pts = vtk.vtkPoints()
        self.pts.SetNumberOfPoints(ntot)
        self.pts.SetData(self.vxyz)

        self.sgrid = vtk.vtkStructuredGrid()
        # VTK requires indexing to be in reverse order
        self.sgrid.SetDimensions(nz1, ny1, nx1)
        self.sgrid.SetPoints(self.pts)

        # convert the structured grid into an unstructured grid
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
    ns = (10, 11, 12)
    ls = (1., 1.1, 1.2)
    cart = CartesianGrid(ns, ls)
    grid = cart.getUnstructuredGrid()    
    #cart.save('cart.vtk')
    cart.show()

if __name__ == '__main__':
    test()

