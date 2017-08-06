import vtk
import numpy

class Intersector:
    """
    Compute the intersection between two grids
    """

    TOL = 1.e-8

    def __init__(self, grid1, grid2):
        self.cell2ToCell1 = self.getIntersections(grid1, grid2)
        self.cell1ToCell2 = self.getIntersections(grid2, grid1)

    def getCells1(self):
        return self.cell2ToCell1

    def getCells2(self):
        return self.cell1ToCell2

    def getIntersections(self, grida, gridb):

        cellBToCellAList = {}
        cellLocator = vtk.vtkCellLocator()
        cellLocator.SetDataSet(grida)
        cellLocator.BuildLocator()
        ptIds = vtk.vtkIdList()
        p0 = numpy.zeros((3,), numpy.float64)
        p1 = numpy.zeros((3,), numpy.float64)
        cellAList = vtk.vtkIdList()
        # iterate over the gridb cells
        for iCellB in range(gridb.GetNumberOfCells()):
            cellb = gridb.GetCell(iCellB)
            # iterate over the edges of this cell
            for iEdgeb in range(cellb.GetNumberOfEdges()):
                edge = cellb.GetEdge(iEdgeb)
                # find the intersection of this edge with grida
                pts = edge.GetPoints()
                p0[:] = pts.GetPoint(0)
                p1[:] = pts.GetPoint(1)
                #print('start/end points: {} {}'.format(p0, p1))
                cellLocator.FindCellsAlongLine(p0, p1, self.TOL, cellAList)
                numCells = cellAList.GetNumberOfIds()
                if numCells > 0 and iCellB not in cellBToCellAList:
                    cellBToCellAList[iCellB] = set()
                for i in range(numCells):
                    #print('numCells = {} i = {} adding cell {}'.format(numCells, i, cellAList.GetId(i)))
                    cellBToCellAList[iCellB].add(cellAList.GetId(i))
        return cellBToCellAList

#####################################################################################

def createTriangulation(pts):
    vPts = vtk.vtkPoints()
    numPoints = len(pts)
    vPts.SetNumberOfPoints(numPoints)
    for i in range(numPoints):
        vPts.SetPoint(i, pts[i])
    poly = vtk.vtkPolyData()
    poly.SetPoints(vPts)
    delny = vtk.vtkDelaunay2D()
    delny.SetInputData(poly)
    delny.Update()
    # convert the output of delny to a vtkUnstructuredGrid
    appendFilter = vtk.vtkAppendFilter()
    appendFilter.AddInputData(delny.GetOutput())
    appendFilter.Update()
    grid = vtk.vtkUnstructuredGrid()
    grid.ShallowCopy(appendFilter.GetOutput())
    return grid

def test0():
    grid = createTriangulation([(0., 0., 0.), (1., 0., 0.), (0., 1., 0.)])

def test1():
    # 2 triangles, they overlap
    grid1 = createTriangulation([(0., 0., 0.),
                                 (1., 0., 0.),
                                 (0., 1., 0.)])
    grid2 = createTriangulation([(0., 0., 0.),
                                 (1., 0., 0.),
                                 (1., 1., 0.)])

    insectr = Intersector(grid1, grid2)
    print('test1: Cells of grid 2 intersecting cells of grid 1: {}'.format(insectr.getCells1()))
    print('test1: Cells of grid 1 intersecting cells of grid 2: {}'.format(insectr.getCells2()))

def test2():
    # 2 triangles, no overlap but close by
    grid1 = createTriangulation([(0., 0., 0.),
                                 (1., 0., 0.),
                                 (0., 1., 0.)])
    grid2 = createTriangulation([(0., 1.2, 0.),
                                 (1., 0.2, 0.),
                                 (1., 1.2, 0.)])

    insectr = Intersector(grid1, grid2)
    print('test2: Cells of grid 2 intersecting cells of grid 1: {}'.format(insectr.getCells1()))
    print('test2: Cells of grid 1 intersecting cells of grid 2: {}'.format(insectr.getCells2()))

def test3():
    # 2 triangles, no overlap and far away
    grid1 = createTriangulation([(0., 0., 0.),
                                 (1., 0., 0.),
                                 (0., 1., 0.)])
    grid2 = createTriangulation([(0., 2.2, 0.),
                                 (1., 1.2, 0.),
                                 (1., 2.2, 0.)])

    insectr = Intersector(grid1, grid2)
    print('test3: Cells of grid 2 intersecting cells of grid 1: {}'.format(insectr.getCells1()))
    print('test3: Cells of grid 1 intersecting cells of grid 2: {}'.format(insectr.getCells2()))


if __name__ == '__main__':
    test0()
    test1()
    test2()
    test3()

