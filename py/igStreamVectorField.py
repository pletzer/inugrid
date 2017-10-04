import vtk
import numpy
from igGridGeometry import GridGeometry

class StreamVectorField:
    """
    Functor to compute a vector field from a stream function, ie d psi ^ d z 
    where z is the perpendicular direction.
    """

    def __init__(self, grid):
        """
        Constructor
        @param grid vtkUnstructuredGrid instance
        """
        self.grid = grid

        # to compute the metric quantities
        self.geom = GridGeometry(grid)

        # maps a face to a list of vertex Ids (hexagon)
        self.faces = {
            (0 ,-1, -1): [0, 3, 7, 4],
            (1, -1, -1): [1, 2, 6, 5],
            (-1, 0, -1): [1, 0, 4, 5],
            (-1, 1, -1): [2, 3, 7, 6],
            (-1, -1, 0): [0, 1, 2, 3],
            (-1, -1, 1): [4, 5, 6, 7],
        }


    def setStreamFunction(self, streamFunction):
        """
        Set the stream function
        @param streamFunction function
        """
        self.streamFunc = streamFunction

    def __call__(self, x):
        """
        Overload of () operator
        @param x position 
        @return vector field at above position
        """

        res = numy.zeros((3,), numpy.float64)

        # call findcell and check if point is in domain, if not just return 0s
        if not self.geom.findCell(x):
            return self.res

        # get the stream function at the cell vertices
        verts = self.geom.getVertices()

        # evaluate the stream function at the vertices
        psis = [self.streamFunc(vert) in verts]

        # get the parametric coordinates of the point in the cell
        xis = self.geom.pcoords

        # grad xi x grad xi vectors
        dS0 = self.geom.getGradX0CrossGradX1(self, 0)
        dS1 = self.geom.getGradX0CrossGradX1(self, 1)
        dS2 = self.geom.getGradX0CrossGradX1(self, 2)
        dSs = [dS0, dS1, dS2]

        # iterate over the faces
        for dim in (0, 1, 2):

            # pick the vector that is perpendicular to the face
            dS = dSs[dim]

            # the parametric coordinate in the direction perpendicular to the face
            xi = xis[dim]

            # iterate over low/high sides
            for lh in (0, 1):

                # set the face ID
                face = [-1, -1, -1]
                face[dim] = lh

                # get the vertex indices for that face
                ptIds = self.faces[tuple(face)]

                # iterate over the edges to compute the flux associated with this face
                # by performing a loop integral (ie using Stokes's theorem). Note that
                # the stream function is expected to be constant in the last coordinate
                # so the line integral is expected to be zero along the last, vertical 
                # coordinate. Along the other directions, it is just the difference
                # between the psi's
                flux = 0.0
                for i0 in range(len(ptIds)):
                    i1 = (i0 + 1) % 4
                    ptId0 = ptIds[i0]
                    ptId1 = ptIds[i1]
                    # the vertical coordinate, expect z0 == z1 along edges where 
                    # psi changes
                    z0 = verts[ptId0][2]
                    z1 = verts[ptId1][2]
                    zmid = 0.5*(z0 + z1)

                    psi0 = psi[ptId0]
                    psi1 = psi[ptId1]

                    # now add the contribution
                    flux += zmid*(psi1 - psi0)
                
                # compute the vector basis at this location by interpolating 
                # the low/high side values
                weight = (1 - lh)*(1.0 - xi) + lh*xi

                self.res += dS * weight * flux

        return res

###############################################################################

def test():
    pass

if __name__ == '__main__':
    test()