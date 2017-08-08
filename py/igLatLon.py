import numpy
import vtk

class LatLon:

    def __init__(self, numLats, numLons, radius=1.0):

        # this is tp convert from structured to unstructured grid
        self.appendGrids = vtk.vtkAppendFilter()

        # create grid
        numLats1, numLons1 = numLats + 1, numLons + 1
        lats = numpy.linspace(-numpy.pi/2., numpy.pi/2., numLats1)
        lons = numpy.linspace(0., 2*numpy.pi, numLons1)
        llats, llons = numpy.meshgrid(lats, lons)
        llats = llats.flat
        llons = llons.flat

        # coordinates
        self.xyz = numpy.zeros((numLats1*numLons1, 3), numpy.float64)
        self.xyz[:, 0] = radius*numpy.cos(llats)*numpy.cos(llons)
        self.xyz[:, 1] = radius*numpy.cos(llats)*numpy.sin(llons)
        self.xyz[:, 2] = radius*numpy.sin(llats)

        # create the VTK unstructired grid
        self.vxyz = vtk.vtkDoubleArray()
        self.vxyz.SetNumberOfComponents(3)
        ntot = numLats1 * numLons1
        self.vxyz.SetNumberOfTuples(ntot)
        self.vxyz.SetVoidArray(self.xyz, 3*ntot, 1)

        self.pts = vtk.vtkPoints()
        self.pts.SetNumberOfPoints(ntot)
        self.pts.SetData(self.vxyz)

        self.sgrid = vtk.vtkStructuredGrid()
        self.sgrid.SetDimensions(numLats1, numLons1, 1)
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
        gridActor.GetProperty().SetColor(93./255., 173./255., 226./255.)

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
    numLat, numLon = 16, 32
    ll = LatLon(numLat, numLon)
    grid = ll.getUnstructuredGrid()
    ll.save('ll.vtk')
    ll.show()

if __name__ == '__main__':
    test()

