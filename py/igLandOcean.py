#!/usr/bin/env python

import re
import operator
import vtk
import numpy


# r, g, b, alpha colormaps
COLORMAPS = {
    'bone': [(x/20., x/20., x/20., x/20.) for x in range(21)],
    'hot':  [(x/4.,0., 0., x/8.0) for x in range(5)] + \
            [(1.,x/4., 0., 0.5) for x in range(5)] + \
            [(1.,1., x/4., 0.5 + x/8.) for x in range(5)], 
    'rainbow': [(1-x/4.,0.,1.,x/8.) for x in range(5)] + \
               [(0.,x/4.,1.,0.5) for x in range(5)] + \
               [(0.,1.,1.-x/4.,0.5) for x in range(5)] + \
               [(x/4.,1.,0.,0.5) for x in range(5)] + \
               [(1.,1.-x/4.,0.,0.5 + x/8.) for x in range(5)],
    'PuBu': [(1-x/4., 1-x/4., 1., x/8.) for x in range(5)] + \
            [(0., 0., 1-x/4., 0.5 + x/8.) for x in range(5)]
}

class LandOcean:

    def __init__(self, textureFile='texture_land_ocean_ice_960_2.png'):
        """
        Constructor
        @param textureFile PNG file with texture data
        """
        self.globe = vtk.vtkTexturedSphereSource()
        self.texture = vtk.vtkTexture()
        self.reader = None
        if re.search(r'jpe?g', textureFile.split('.')[-1]):
            self.reader = vtk.vtkJPEGReader()
        elif textureFile.split('.')[-1] == 'png':
            self.reader = vtk.vtkPNGReader()
        self.mapper = vtk.vtkPolyDataMapper()
        a = vtk.vtkActor()
        self.actors = [a]

        # connect
        a.SetMapper(self.mapper)
        a.SetTexture(self.texture)

        self.mapper.SetInputConnection(self.globe.GetOutputPort())
        self.texture.SetInputConnection(self.reader.GetOutputPort())
        self.reader.SetFileName(textureFile)

        self.globe.SetThetaResolution(128)
        self.globe.SetPhiResolution(64)
        self.globe.SetRadius(1.0)

    def setData(self, lats, lons, data, colormap='rainbow', text=''):
        """
        Set data 
        @param lats 1-D array of latitudes (deg)
        @param lons 1-D array of longitudes (deg)
        @param data 2-D array 
        """

        self.text = text

        # mid position for the camera
        self.midLat = reduce(operator.add, lats, 0.0)
        self.midLon = reduce(operator.add, lons, 0.0)
        if len(lats) > 0:
            self.midLat /= float(len(lats))
        if len(lons) > 0:
            self.midLon /= float(len(lons))

        nlat1, nlon1 = len(lats), len(lons)
        self.pointArray = vtk.vtkDoubleArray()
        self.dataArray = vtk.vtkDoubleArray()
        self.points = vtk.vtkPoints()
        self.grid = vtk.vtkStructuredGrid()
        self.fltr = vtk.vtkGeometryFilter()
        self.dataMapper = vtk.vtkPolyDataMapper()

        npoints = nlat1*nlon1
        self.pointArray.SetNumberOfTuples(npoints)
        self.pointArray.SetNumberOfComponents(3)
        llats = numpy.outer(lats, numpy.ones(nlon1))
        llons = numpy.outer(numpy.ones(nlat1), lons)
        self.xyz = numpy.zeros( (npoints, 3), numpy.float64 )
        radius = 1.05
        xx = radius*numpy.cos(llats*numpy.pi/180.)*numpy.cos(llons*numpy.pi/180.)
        yy = radius*numpy.cos(llats*numpy.pi/180.)*numpy.sin(llons*numpy.pi/180.)
        zz = radius*numpy.sin(llats*numpy.pi/180.)
        self.xyz[:, 0] = xx.flat
        self.xyz[:, 1] = yy.flat
        self.xyz[:, 2] = zz.flat

        self.pointArray.SetVoidArray(self.xyz, npoints*3, 1)

        self.dataArray.SetNumberOfTuples(npoints)
        self.dataArray.SetNumberOfComponents(1)
        self.dataArray.SetVoidArray(data, npoints, 1)

        self.points.SetNumberOfPoints(npoints)
        self.points.SetDataTypeToDouble()
        self.points.SetData(self.pointArray)

        # fastest varying index goes first
        self.grid.SetDimensions(nlon1, nlat1, 1)
        self.grid.SetPoints(self.points)
        self.grid.GetPointData().SetScalars(self.dataArray)

        if vtk.VTK_MAJOR_VERSION >= 6:
            self.fltr.SetInputData(self.grid)
        else:
            self.fltr.SetInput(self.grid)
        self.dataMapper.SetInputConnection(self.fltr.GetOutputPort())
        a = vtk.vtkActor()
        a.SetMapper(self.dataMapper)
        self.actors.append(a)

        dataMin = min(data.flat)
        dataMax = max(data.flat)
        self.dataMapper.SetScalarRange(dataMin, dataMax)

        # build the color map
        lut = self.dataMapper.GetLookupTable()
        rgbaValues = COLORMAPS[colormap]
        lut.SetNumberOfTableValues(len(rgbaValues))
        for i in range(len(rgbaValues)):
            lut.SetTableValue(i, rgbaValues[i][0], rgbaValues[i][1], rgbaValues[i][2], rgbaValues[i][3])
        lut.SetRange(dataMin, dataMax)

    def addText(self, text, xpos=20, ypos=30, size=20):
        a = vtk.vtkTextActor()
        a.SetInput(text)
        a.SetPosition(xpos, ypos)
        a.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        txtProp = a.GetTextProperty()
        txtProp.SetFontSize(size)
        txtProp.SetFontFamilyToArial()
        self.actors.append(a)

    def render(self, width=960, height=640):

        self.__addImage(solarfsCfg.LOGOFILE, xpos=0.05, ypos=0.05, xpixSize=width//10)

        self.ren = vtk.vtkRenderer()
        self.renWin = vtk.vtkRenderWindow()
        self.iren = vtk.vtkRenderWindowInteractor()

        # add color bar
        if hasattr(self, 'dataMapper'):
            colorbar = vtk.vtkScalarBarActor()
            colorbar.SetLookupTable(self.dataMapper.GetLookupTable())
            colorbar.SetTitle(self.text)
            colorbar.GetLabelTextProperty().SetFontSize(solarfsCfg.COLORBAR_TEXT_SIZE)
            self.actors.append(colorbar)

        self.camera = vtk.vtkCamera()
        self.camera.OrthogonalizeViewUp()
        self.camera.SetFocalPoint(0., 0., 0.)
        self.camera.Roll(270. - self.midLon)
        radius = 5.0
        x = radius*numpy.cos(self.midLat*numpy.pi/180.)*numpy.cos(self.midLon*numpy.pi/180.)
        y = radius*numpy.cos(self.midLat*numpy.pi/180.)*numpy.sin(self.midLon*numpy.pi/180.)
        z = radius*numpy.sin(self.midLat*numpy.pi/180.)
        self.camera.SetPosition(x, y, z)

        self.ren.SetActiveCamera(self.camera)

        self.renWin.AddRenderer(self.ren)
        self.iren.SetRenderWindow(self.renWin)
        for a in self.actors:
            self.ren.AddActor(a)

        self.ren.SetBackground(solarfsCfg.WINDOW_BACKGROUND)

        self.renWin.SetSize(width, height)
        self.iren.Initialize()
        self.renWin.Render()
        self.iren.Start()

    def __addImage(self, imageFile, xpos, ypos, xpixSize):
        self.imageReader = None
        suffix = imageFile.split('.')[-1].lower()
        if suffix == 'png':
            self.imageReader = vtk.vtkPNGReader()
        elif re.search(r'jp?g', suffix):
            self.imageReader = vtk.vtkJPEGReader()
        self.imageReader.SetFileName(imageFile)

        #self.imageResizer = vtk.vtkImageResize()
        self.imageResizer = vtk.vtkImageResample()
        self.imageResizer.SetInputConnection(self.imageReader.GetOutputPort())
        self.imageResizer.SetAxisMagnificationFactor(0, 0.25)
        self.imageResizer.SetAxisMagnificationFactor(1, 0.25)
        #xp, yp, zp = self.imageResizer.GetOutputDimensions()
        #self.imageResizer.SetOutputDimensions(xpixSize, 0.8*xpixSize*yp//xp, xpixSize*zp//xp)

        self.imageMapper = vtk.vtkImageMapper()
        self.imageMapper.SetInputConnection(self.imageResizer.GetOutputPort())

        self.imageMapper.SetColorWindow(255)
        self.imageMapper.SetColorLevel(127.5)

        a = vtk.vtkActor2D()
        a.SetMapper(self.imageMapper)
        a.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        a.GetPositionCoordinate().SetValue(xpos, ypos)
        self.actors.append(a)

#########################################################
def test():
    pnt = Painter()
    lats = numpy.arange(-20, 60.1, 1.)
    lons = numpy.arange(80,  130, 1.)
    llats = numpy.outer(lats, numpy.ones((len(lons),), numpy.float64))
    llons = numpy.outer(numpy.ones((len(lats),), numpy.float64), lons)
    data = llons
    pnt.setData(lats, lons, data, colormap='bone')
    pnt.addText("hello", xpos=0.1, ypos=0.9, size=50)
    pnt.render()

if __name__ == '__main__': test()

