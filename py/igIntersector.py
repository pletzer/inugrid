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
				cellLocator.FindCellsAlongLine(p0, p1, TOL, cellAList)
				numCells = cellAList.GetNumberOfIds()
				for i in range(numCells):
					cellBToCellAList = cellBToCellAList.get(iCellB, []) + \
						[cellAList.GetId(i)]

		return cellBToCellAList

#####################################################################################


def test1():
	# create grid 1, one triangle
	pts1 = vtk.vtkPoints()
	pts1.SetNumberOfPoints(3)
	pts1.SetPoint(0, [0., 0., 0.])
	pts1.SetPoint(1, [1., 0., 0.])
	pts1.SetPoint(2, [0., 1., 0.])
	grid1 = vtk.vtkUnstructuredGrid()
	grid1.Allocate()
	pt1Ids = vtk.vtkIdList()
	pt1Ids.SetNumberOfIds(3)
	pt1Ids.SetId(0, 0)
	pt1Ids.SetId(1, 1)
	pt1Ids.SetId(2, 2)
	grid1.InsertNextCell(vtk.VTK_TRIANGLE, pt1Ids)

	# create grid 2, another triangle
	pts2 = vtk.vtkPoints()
	pts2.SetNumberOfPoints(3)
	pts2.SetPoint(0, [0., 0., 0.])
	pts2.SetPoint(1, [1., 0., 0.])
	pts2.SetPoint(2, [1., 1., 0.])
	grid2 = vtk.vtkUnstructuredGrid()
	grid2.Allocate()
	pt2Ids = vtk.vtkIdList()
	pt2Ids.SetNumberOfIds(3)
	pt2Ids.SetId(0, 0)
	pt2Ids.SetId(1, 1)
	pt2Ids.SetId(2, 2)
	grid2.InsertNextCell(vtk.VTK_TRIANGLE, pt2Ids)

	insectr = Intersector(grid1, grid2)

if __name__ == '__main__': 
	test1()

