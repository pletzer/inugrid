import vtk
import numpy

class PolyLineWriter:

    def __init__(self, xyzArray):
        """
        Constructor
        @param xyzArray [(x,y,z), ...]
        """
        self.appendFilter = vtk.vtkAppendFilter()

        self.vpts = vtk.vtkPoints()
        self.line = vtk.vtkPolyLine()
        self.ug = vtk.vtkUnstructuredGrid()

        numPoints = xyzArray.shape[0]
        self.vpts.SetNumberOfPoints(numPoints)

        self.ptIds = self.line.GetPointIds()
        self.ptIds.SetNumberOfIds(numPoints)
        index = 0
        for x, y, z in xyzArray:
            self.vpts.InsertPoint(index, x, y, z)
            self.ptIds.SetId(index, index)
            index += 1

        # one cell
        self.ug.InsertNextCell(self.line.GetCellType(), self.ptIds)
        self.ug.SetPoints(self.vpts)

        self.appendFilter.AddInputData(self.ug)
        self.appendFilter.Update()
        self.ugrid = self.appendFilter.GetOutput()

    def save(self, filename):
        """
        Save data in file
        @param filename file name
        """
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.ugrid)
        writer.Update()

