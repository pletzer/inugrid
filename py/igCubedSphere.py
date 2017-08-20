import numpy
import vtk

class CubedSphere:

    def __init__(self, numCellsPerTile, radius=1.0):


        self.appendGrids = vtk.vtkAppendFilter()
        #appendGrids.MergePointsOn()
        self.appendGrids.SetOutputPointsPrecision(1) # double

        # create tile grids
        numPointsPerTile = numCellsPerTile + 1
        us = numpy.linspace(0., 1., numPointsPerTile)
        vs = numpy.linspace(0., 1., numPointsPerTile)
        uu, vv = numpy.meshgrid(us, vs)
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
                xyz = numpy.zeros((numPointsPerTile*numPointsPerTile, 3), 
                                  numpy.float64)
                # grid on the box's side/tile 
                xyz[:, dim0] = (pm + 1.0)/2.
                xyz[:, dim1] = uu.flat
                xyz[:, dim2] = vv.flat
                # project the vertices onto sphere
                for i in range(3):
                    xyz[:, i] -= centre[i]
                dist = numpy.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2)
                for i in range(3):
                    # normalize
                    xyz[:, i] /= dist
                    # extend to the sphere's surface
                    xyz[:, i] *= radius

                # compute the cell areas
                areas = self.getCellAreas(xyz)
                tileAreas = vtk.vtkDoubleArray()
                tileAreas.SetNumberOfComponents(1)
                tileAreas.SetNumberOfTuples(numCellsPerTile * numCellsPerTile)
                tileAreas.SetVoidArray(areas, numCellsPerTile * numCellsPerTile, 1)

                ntot = numPointsPerTile**2

                # create the VTK unstructired grid
                tileXyz = vtk.vtkDoubleArray()
                tileXyz.SetNumberOfComponents(3)
                tileXyz.SetNumberOfTuples(ntot)
                tileXyz.SetVoidArray(xyz, 3*ntot, 1)

                tilePts = vtk.vtkPoints()
                tilePts.SetNumberOfPoints(ntot)
                tilePts.SetData(tileXyz)

                tileGrid = vtk.vtkStructuredGrid()
                tileGrid.SetDimensions(numPointsPerTile, numPointsPerTile, 1)
                tileGrid.SetPoints(tilePts)
                tileGrid.GetCellData().SetScalars(tileAreas)

                self.appendGrids.AddInputData(tileGrid)

                self.xyzList.append(xyz)
                self.tileXyzList.append(tileXyz)
                self.tilePtsList.append(tilePts)
                self.tileGridList.append(tileGrid)

        self.appendGrids.Update()
        self.grid = self.appendGrids.GetOutput()

    def getCornerCoords(self, xyz, coordIndex):
        numPointsPerTile = int(numpy.sqrt(xyz.shape[0]))
        xx = xyz[:, coordIndex].reshape(numPointsPerTile, numPointsPerTile)
        xx0 = xx[:-1, :-1]
        xx1 = xx[1:, :-1]
        xx2 = xx[1:, 1:]
        xx3 = xx[:-1, 1:]
        return xx0, xx1, xx2, xx3


    def getTriangleAreas(self, xx0, xx1, xx2, yy0, yy1, yy2, zz0, zz1, zz2):
        dxx10, dyy10, dzz10 = xx1 - xx0, yy1 - yy0, zz1 - zz0
        dxx20, dyy20, dzz20 = xx2 - xx0, yy2 - yy0, zz2 - zz0
        areasxx = dyy10*dzz20 - dyy20*dzz10
        areasyy = dzz10*dxx20 - dzz20*dxx10
        areaszz = dxx10*dyy20 - dxx20*dyy10
        xx = (xx0 + xx1 + xx2)/3.
        yy = (yy0 + yy1 + yy2)/3.
        zz = (zz0 + zz1 + zz2)/3.
        rr = numpy.sqrt(xx**2 + yy**2 + zz**2)
        return (areasxx*xx + areasyy*yy + areaszz*zz)/rr


    def getCellAreas(self, xyz):
        xx0, xx1, xx2, xx3 = self.getCornerCoords(xyz, 0)
        yy0, yy1, yy2, yy3 = self.getCornerCoords(xyz, 1)
        zz0, zz1, zz2, zz3 = self.getCornerCoords(xyz, 2)
        areas = self.getTriangleAreas(xx0, xx1, xx3, yy0, yy1, yy3, zz0, zz1, zz3)
        areas += self.getTriangleAreas(xx2, xx3, xx1, yy2, yy3, yy1, zz2, zz3, zz1)
        return areas


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
    numCells = 10
    cs = CubedSphere(numCells)
    grid = cs.getUnstructuredGrid()
    cs.save('cs.vtk')
    cs.show()

if __name__ == '__main__':
    test()

