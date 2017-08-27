import shapefile # must have pyshp installed
import numpy as np
import vtk
        

class Coastlines:

    def __init__(self, radius=1.0):
        """
        Constructor
        @param radius earth's radius
        """
    
        self.sf = shapefile.Reader('ne_10m_coastline')

        self.pts = []
        self.lines = []
        self.ugrids = []
        self.appendFilter = vtk.vtkAppendFilter()

        for s in self.sf.shapes():

            numPoints = len(s.points)

            # skip some points (not sure what they represent)
            if numPoints < 100:
                # skip
                continue

            vpts = vtk.vtkPoints()
            line = vtk.vtkPolyLine()
            ug = vtk.vtkUnstructuredGrid()

            vpts.SetNumberOfPoints(numPoints)

            segs = line.GetPointIds()
            segs.SetNumberOfIds(numPoints + 1)
            index = 0
            for p in s.points:
                lam, the = p[0]*np.pi/180., p[1]*np.pi/180.
                x = radius*np.cos(the)*np.cos(lam)
                y = radius*np.cos(the)*np.sin(lam)
                z = radius*np.sin(the)
                vpts.InsertPoint(index, x, y, z)
                segs.SetId(index, index)
                index += 1
            segs.SetId(index, 0) # close the loop
            ug.InsertNextCell(line.GetCellType(), segs)
            ug.SetPoints(vpts)

            self.pts.append(vpts)
            self.lines.append(line)
            self.ugrids.append(ug)

            self.appendFilter.AddInputData(ug)

        self.appendFilter.Update()
        self.ugrid = self.appendFilter.GetOutput()

    def show(self):
        """
        Show the cloastline
        """
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(self.ugrid)
        actor = vtk.vtkActor()
        actor.GetProperty().SetColor(0, 0, 0)
        actor.SetMapper(mapper)

        renderer = vtk.vtkRenderer()
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)
        renderer.AddActor(actor)
        renderer.SetBackground(.8, 0.8, 0.8)
        renderWindow.Render()
        renderWindowInteractor.Start()

    def save(self, filename):
        """
        Save data in file
        @param filename file name
        """
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.ugrid)
        writer.Update()


###############################################################################
def test():
    cl = Coastlines()
    cl.show()
    cl.save('coastlines.vtk')

if __name__ == '__main__':
    test()