import vtk
import numpy

class NodalFunctionWriter:

    def __init__(self, ugrid, nodalFunc, name='nodalFunction'):
        """
        Constructor
        @param ugrid vtkUnstructuredGrid instance
        @param nodalFunc
        """
        self.grid = ugrid
        pts = self.grid.GetPoints()
        numPoints = pts.GetNumberOfPoints()
        data = vtk.vtkDoubleArray()
        data.SetName(name)
        data.SetNumberOfComponents(1)
        data.SetNumberOfTuples(numPoints)
        for i in range(numPoints):
            xyz = pts.GetPoint(i)
            f = nodalFunc(xyz)
            data.SetTuple(i, (f,))
        # attach to the grid
        self.grid.GetPointData().AddArray(data)

    def save(self, filename):
        """
        Save data in file
        @param filename file name
        """
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.grid)
        writer.Update()

