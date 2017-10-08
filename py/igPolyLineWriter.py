import vtk
import numpy

class PolyLineWriter:

    def __init__(self):
        """
        Constructor
        """
        self.appendFilter = vtk.vtkAppendFilter()

        self.vpts = vtk.vtkPoints()
        self.ugrid = vtk.vtkUnstructuredGrid()

        self.numPoints = 0
        self.lines = []

    def addLine(self, xyzArray):
        """
        Add a line
        @param xyzArray
        """
        self.lines.append(xyzArray)
        self.numPoints += xyzArray.shape[0]


    def build(self):
        """
        Create unstructured grid
        """
        index = 0
        self.vpts.SetNumberOfPoints(self.numPoints)
        self.vlineList = []
        self.ptIdsList = []
        for line in self.lines:
            npts = line.shape[0]
            vline = vtk.vtkPolyLine()
            ptIds = vline.GetPointIds()
            ptIds.SetNumberOfIds(npts)
            for i in range(npts):
                self.vpts.SetPoint(index, line[i, :])
                ptIds.SetId(i, index)
                index += 1

            self.ugrid.InsertNextCell(vline.GetCellType(), ptIds)
            self.vlineList.append(vline)
            self.ptIdsList.append(ptIds)
        
        self.ugrid.SetPoints(self.vpts)


    def save(self, filename):
        """
        Save data in file
        @param filename file name
        """
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.ugrid)
        writer.Update()

