import numpy
import vtk
import igAreas

class LatLon:

    def __init__(self, numLats, numLons, radius=1.0):

        # this is tp convert from structured to unstructured grid
        self.appendGrids = vtk.vtkAppendFilter()

        # create grid
        numLats1, numLons1 = numLats + 1, numLons + 1
        lats = numpy.linspace(-numpy.pi/2., numpy.pi/2., numLats1)
        lons = numpy.linspace(0., 2*numpy.pi, numLons1)
        llons, llats = numpy.meshgrid(lons, lats)

        # add a tiny perturbation  to avoid issues with computing
        # fluxes
        #eps = 1.e-3
        #llons += eps*numpy.cos(llats)*numpy.sin(llons)
        #llats += eps*numpy.cos(llats)*numpy.sin(llons)

        llats = llats.flat
        llons = llons.flat

        # coordinates
        self.xyz = numpy.zeros((numLats1*numLons1, 3), numpy.float64)
        rrho = radius*numpy.cos(llats)
        self.xyz[:, 0] = rrho*numpy.cos(llons)
        self.xyz[:, 1] = rrho*numpy.sin(llons)
        self.xyz[:, 2] = radius*numpy.sin(llats)

        # compute the cell areas, enforce positiveness (not sure why all the areas are negative)
        self.areas = -igAreas.getCellAreas(self.xyz, n0=numLats1, n1=numLons1)
        self.vareas = vtk.vtkDoubleArray()
        self.vareas.SetName('cell_areas')
        self.vareas.SetNumberOfComponents(1)
        self.vareas.SetNumberOfTuples(numLats * numLons)
        self.vareas.SetVoidArray(self.areas, numLats * numLons, 1)

        # create the VTK unstructured grid
        self.vxyz = vtk.vtkDoubleArray()
        self.vxyz.SetNumberOfComponents(3)
        ntot = numLats1 * numLons1
        self.vxyz.SetNumberOfTuples(ntot)
        self.vxyz.SetVoidArray(self.xyz, 3*ntot, 1)

        self.pts = vtk.vtkPoints()
        self.pts.SetNumberOfPoints(ntot)
        self.pts.SetData(self.vxyz)

        self.sgrid = vtk.vtkStructuredGrid()
        self.sgrid.SetDimensions(numLons1, numLats1, 1)
        self.sgrid.SetPoints(self.pts)

        self.sgrid.GetCellData().SetScalars(self.vareas)

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

        actors = []

        gridMapper = vtk.vtkDataSetMapper()
        gridMapper.SetInputData(self.grid)

        data = self.grid.GetCellData().GetScalars()
        if data:
            lut = vtk.vtkLookupTable()
            lut.SetHueRange(0., 0.666)
            dmin, dmax = data.GetRange()
            dmin = 0.0
            lut.SetTableRange(dmin, dmax)
            lut.Build()

            cbar = vtk.vtkScalarBarActor()
            cbar.SetLookupTable(lut)
            actors.append(cbar)

            gridMapper.SetLookupTable(lut)
            gridMapper.SetUseLookupTableScalarRange(1)

        gridActor = vtk.vtkActor()
        gridActor.SetMapper(gridMapper)
        #gridActor.GetProperty().SetColor(93./255., 173./255., 226./255.)
        actors.append(gridActor)

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
        # add the actors to the renderer, set the background and size
        for a in actors:
            ren.AddActor(a)
        #ren.AddLight(light)
        ren.SetBackground(1, 1, 1)
        renWin.SetSize(640, 640)
        iren.Initialize()
        renWin.Render()
        iren.Start()


#############################################################################
def test():
    numLat, numLon = 10, 20
    ll = LatLon(numLat, numLon)
    grid = ll.getUnstructuredGrid()
    ll.save('ll.vtk')
    ll.show()

if __name__ == '__main__':
    test()

