import numpy
import vtk
import igAreas

class LatLonElv:

    def __init__(self, numLats, numLons, numElvs, radius=1.0, maxRelElv = 0.1):

        # this is tp convert from structured to unstructured grid
        self.appendGrids = vtk.vtkAppendFilter()

        # create grid
        numLats1, numLons1, numElvs1 = numLats + 1, numLons + 1, numElvs + 1
        lats = numpy.linspace(-numpy.pi/2., numpy.pi/2., numLats1)
        lons = numpy.linspace(0., 2*numpy.pi, numLons1)
        elvs = numpy.linspace(0., maxRelElv, numElvs1)
        llons, llats, eelvs = numpy.meshgrid(lons, lats, elvs, sparse=False, indexing='ij')

        llats = llats.flat
        llons = llons.flat
        eelvs = eelvs.flat

        # coordinates
        ntot = numLats1 * numLons1 * numElvs1
        self.xyz = numpy.zeros((ntot, 3), numpy.float64)
        rr = eelvs.copy()
        rr += 1.0
        rr *= radius
        rrho = rr * numpy.cos(llats)
        self.xyz[:, 0] = rrho * numpy.cos(llons)
        self.xyz[:, 1] = rrho * numpy.sin(llons)
        self.xyz[:, 2] = rr * numpy.sin(llats)

        # create the VTK unstructured grid
        self.vxyz = vtk.vtkDoubleArray()
        self.vxyz.SetNumberOfComponents(3)
        self.vxyz.SetNumberOfTuples(ntot)
        self.vxyz.SetVoidArray(self.xyz, 3*ntot, 1)

        self.pts = vtk.vtkPoints()
        self.pts.SetNumberOfPoints(ntot)
        self.pts.SetData(self.vxyz)

        self.sgrid = vtk.vtkStructuredGrid()
        self.sgrid.SetDimensions(numLons1, numLats1, numElvs1)
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
    numLats, numLons, numElvs = 8, 16, 1
    ll3 = LatLonElv(numLats, numLons, numElvs,)
    grid = ll3.getUnstructuredGrid()
    ll3.save('ll3.vtk')
    ll3.show()

if __name__ == '__main__':
    test()

