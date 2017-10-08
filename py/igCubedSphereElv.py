import numpy
import vtk
import igAreas

class CubedSphereElv:

    def __init__(self, numCellsPerTile, numElvs, radius=1.0, maxRelElv=0.1):


        self.appendGrids = vtk.vtkAppendFilter()
        #appendGrids.MergePointsOn()
        self.appendGrids.SetOutputPointsPrecision(1) # double

        # create tile grids
        numPointsPerTile = numCellsPerTile + 1
        numElvs1 = numElvs + 1
        ntot = numPointsPerTile**2 * numElvs1

        us = numpy.linspace(0., 1., numPointsPerTile)
        vs = numpy.linspace(0., 1., numPointsPerTile)
        ws = numpy.linspace(0., 1., numElvs1)
        uu, vv, ww = numpy.meshgrid(us, vs, ws, sparse=False, indexing='ij')

        ww += 1.0

        # box is [0, 1] x [0, 1], let's fit a sphere inside the box
        centre = numpy.array([0.5, 0.5, 0.5])

        self.xyzList = []
        self.tileXyzList = []
        self.tilePtsList = []
        self.tileGridList = []

        # iterate over the space dimensions
        for dim0 in range(3):

            # low or high side
            for pm in range(-1, 2, 2):

                # normal vector, pointing out
                normal = numpy.zeros((3,), numpy.float64)
                normal[dim0] = pm

                # indices of the dimensions on the tile
                dim1 = (dim0 + 1) % 3
                dim2 = (dim0 + 2) % 3

                # coordinates
                xyz = numpy.zeros((ntot, 3), numpy.float64)

                # grid on the box's side/tile 
                xyz[:, dim0] = (pm + 1.0)/2.
                xyz[:, dim0] *= ww.flat
                xyz[:, dim1] = uu.flat
                xyz[:, dim2] = vv.flat
                # fix the vertex ordering so the area points outwards
                if pm > 0:
                    xyz[:, dim1] *= -1.0
                    xyz[:, dim1] += 1.0

                # project the vertices onto sphere
                for i in range(3):
                    xyz[:, i] -= centre[i]
                dist = numpy.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2)
                for i in range(3):
                    # normalize
                    xyz[:, i] /= dist
                    # extend to the sphere's surface
                    xyz[:, i] *= radius

                # create the VTK unstructired grid
                tileXyz = vtk.vtkDoubleArray()
                tileXyz.SetNumberOfComponents(3)
                tileXyz.SetNumberOfTuples(ntot)
                tileXyz.SetVoidArray(xyz, 3*ntot, 1)

                tilePts = vtk.vtkPoints()
                tilePts.SetNumberOfPoints(ntot)
                tilePts.SetData(tileXyz)

                tileGrid = vtk.vtkStructuredGrid()
                # inverse order, first index varies fastest
                tileGrid.SetDimensions(numElvs + 1, numPointsPerTile, numPointsPerTile)
                tileGrid.SetPoints(tilePts)

                self.appendGrids.AddInputData(tileGrid)

                self.xyzList.append(xyz)
                self.tileXyzList.append(tileXyz)
                self.tilePtsList.append(tilePts)
                self.tileGridList.append(tileGrid)

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
            lut.SetHueRange(0.666, 0.)
            dmin, dmax = data.GetRange()
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
        light.SetPosition(1., 1., 0.)
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
    numCells = 5
    cs3d = CubedSphereElv(numCells, numElvs=2)
    grid = cs3d.getUnstructuredGrid()
    cs3d.save('cs.vtk')
    cs3d.show()

if __name__ == '__main__':
    test()

