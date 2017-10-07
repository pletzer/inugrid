import vtk
import numpy
from igGridGeometry import GridGeometry

class StreamVectorField:
    """
    Functor to compute a vector field from a stream function, ie d psi ^ d z 
    where z is the perpendicular direction. It is assumed that psi = psi(x, y) does
    not depend on z.
    """

    def __init__(self, grid):
        """
        Constructor
        @param grid vtkUnstructuredGrid instance
        """
        self.grid = grid

        # to compute the metric quantities
        self.geom = GridGeometry(grid)

        self.psiIndices = {
            (0 ,-1, -1): (0, 3),
            (1, -1, -1): (4, 7),
            (-1, 0, -1): (4, 0),
            (-1, 1, -1): (7, 3),
        }


    def setStreamFunction(self, streamFunction):
        """
        Set the stream function
        @param streamFunction function
        """
        self.streamFunc = streamFunction

    def __call__(self, x, *args):
        """
        Overload of () operator
        @param x position 
        @return vector field at above position
        """

        res = numpy.zeros((3,), numpy.float64)

        # call findcell and check if point is in domain, if not just return 0s
        if not self.geom.findCell(x):
            return res

        # get the stream function at the cell vertices
        verts = self.geom.getVertices()

        # evaluate the stream function at the vertices
        psis = [self.streamFunc(vert) for vert in verts]

        # get the parametric coordinates of the point in the cell
        xis = self.geom.pcoords

        # grad xi x grad xi vectors
        dS0 = self.geom.getGradXiCrossGradXi(0)
        dS1 = self.geom.getGradXiCrossGradXi(1)
        dS2 = self.geom.getGradXiCrossGradXi(2)
        dSs = [dS0, dS1, dS2]
        zHat = numpy.array([0., 0., 1.])

        # iterate over the faces
        for dim in (0, 1):

            # pick the vector that is perpendicular to the face
            dS = dSs[dim]

            # the parametric coordinate in the direction perpendicular to the face
            xi = xis[dim]

            # iterate over low/high sides
            for lh in (0, 1):

                # set the face ID
                face = [-1, -1, -1]
                face[dim] = lh

                ptId0, ptId1 = self.psiIndices[tuple(face)]
                flux = psis[ptId1] - psis[ptId0]

                # interpolate the flux at this location
                weight = (1 - lh)*(1.0 - xi) + lh*xi

                res += dS * weight * flux

        return res

###############################################################################

def test():
    pass

if __name__ == '__main__':
    test()