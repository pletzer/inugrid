import numpy
import vtk

class CubedSphere:

	def __init__(self, numCellsPerTile, radius=1.0):


		appendGrids = vtk.vtkAppendFilter()
		#appendGrids.MergePointsOn()
		appendGrids.SetOutputPointsPrecision(1) # double

		# create tile grids
		numPointsPerTile = numCellsPerTile + 1
		us = numpy.linspace(0., 1., numPointsPerTile)
		vs = numpy.linspace(0., 1., numPointsPerTile)
		uu, vv = numpy.meshgrid(us, vs)
		xyz = numpy.zeros((numPointsPerTile*numPointsPerTile, 3), numpy.float64)
		# box is [0, 1] x [0, 1], let's fit a sphere inside the box
		centre = numpy.array([0.5, 0.5, 0.5])

		self.tilesXyz = []
		self.tilesPts = []
		self.tilesGrid = []
		
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

				# grid on the box's side/tile 
				xyz[:, dim0] = (pm + 1.0)/2.
				xyz[:, dim1] = uu.flat
				xyz[:, dim2] = vv.flat
				#print xyz

				"""
				# project the vertices onto sphere
				for i in range(3):
					xyz[:, i] -= centre[i]
				dist = numpy.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2)
				for i in range(3):
					# normalize
					xyz[:, i] /= dist
					# extend to the sphere's surface
					xyz[:, i] *= radius
				"""

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
				#print tileGrid

				self.tilesXyz.append(tileXyz)
				self.tilesPts.append(tilePts)
				self.tilesGrid.append(tileGrid)

				appendGrids.AddInputData(tileGrid)

		appendGrids.Update()
		self.grid = appendGrids.GetOutput()
		print self.grid


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

		ren = vtk.vtkRenderer()
		renWin = vtk.vtkRenderWindow()
		renWin.AddRenderer(ren)
		iren = vtk.vtkRenderWindowInteractor()
		iren.SetRenderWindow(renWin)
		# add the actors to the renderer, set the background and size
		ren.AddActor(gridActor)
		ren.SetBackground(0.5, 0.5, 0.5)
		renWin.SetSize(400, 300)
		iren.Initialize()
		renWin.Render()
		iren.Start()


###################################################################################
def test():
	numCells = 2
	cs = CubedSphere(numCells)
	cs.save('cs.vtk')
	cs.show()

if __name__ == '__main__':
	test()

