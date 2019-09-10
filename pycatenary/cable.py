from . import catenary
import numpy as np


class MooringLine:
    """
    Class to create a mooring line between anchor and fairlead points
    Parameters
    ----------
    w: double
        submerged weight [N/m]
    EA: double
        axial stiffness
    L: double
        unstretched line length [m]
    anchor: np.ndarray
        anchor coordinates
    fairlead: np.ndarray
        fairlead coordinates
    """
    count = 0
    def __init__(
            self,
            L,
            w,
            EA=0.,
            anchor=None,
            fairlead=None,
            nd=3,
            tol=1e-6,
            maxiter=10000,
            floor=True,
    ):
        self.__class__.count += 1
        self.nd = nd
        self.name = 'cable_'+str(self.count)
        self.w = get_array(w)
        #self.E = E
        #self.A0 = A0
        self.EA = get_array(EA)
        self.L = get_array(L)
        self.tol = tol
        self.maxiter = maxiter
        if anchor is None:
            self.anchor = np.zeros(nd)  # coordinates of anchor
        else:
            self.anchor = np.array(anchor)
        if fairlead is None:
            self.fairlead = np.zeros(nd)  # coordinates of fairlead
        else:
            self.fairlead = np.array(fairlead)
        self.anchor_coords_system = np.eye(3)
        self.setDirectionDistance()
        self.X = np.zeros(nd)  # local coords with anchor at origin
        self.x = np.zeros(nd)  # local coords (Ls, y)
        self.Tf = np.zeros(nd)  # horizontal pretension
        self.Ta = np.zeros(nd)
        self.Ls = 0.  # lifted line length
        self.e = 0.
        self.a = 1.
        self.x0 = None
        self.broken = False
        self.floor = floor
        if EA == 0.:
            self.catenary = catenary.CatenaryRigid(
                L=self.L,
                w=self.w,
            )
        else:
            self.catenary = catenary.CatenaryElastic(
                L=self.L,
                w=self.w,
                EA=self.EA,
            )

        #print(self.fairlead)
        self.setDirectionDistance()

    def setDirectionDistance(self):
        if self.nd == 3:
            self.distance = np.sqrt(np.sum((self.fairlead[:2]-self.anchor[:2])**2))
            self.direction = (self.fairlead[:2]-self.anchor[:2])/self.distance
        elif self.nd == 2:
            self.direction = np.array([1., 0.])
            self.distance = np.abs(self.fairlead[0]-self.anchor[0])

    def setVariables(self):
        self.setDirectionDistance()
        if self.nd == 2:
            h = self.fairlead[1]-self.anchor[1]
        elif self.nd == 3:
            h = self.fairlead[2]-self.anchor[2]
        d = self.distance
        cat = self.catenary
        cat.getState(d, h, self.L)

        self.Ls = self.catenary.Ls
        if self.nd == 2:
            # self.Tf = Tf
            # self.Ta = Ta
            self.s2xyz = lambda s: self.anchor+self.catenary.s2xy(s)
            self.ds2xyz = self.catenary.ds2xy
        if self.nd == 3:
            Tf3D = np.zeros(3)
            Ta3D = np.zeros(3)
            # horizontal component from 2D to 3D on x-y
            Ta3D[:2] = Tf[0]*self.direction
            Tf3D[:2] = Tf[0]*self.direction
            self.y = lambda x: y_coords(x)*self.direction
            self.s2xyz = lambda s: self.anchor+[s_coords(s)[0]*self.direction[0], s_coords(s)[0]*self.direction[1], s_coords(s)[1]]
            self.ds2xyz = lambda s: np.array([
                cat.ds2xyz(s)[0]*self.direction[0],
                cat.ds2xyz(s)[0]*self.direction[1],
                cat.ds2xyz(s)[1],
            ])
            Tf3D[2] = Tf[1]
            Ta3D[2] = Ta[1]
            self.Tf = Tf3D
            self.Ta = Ta3D

    def setAnchorCoords(self, coords):
        """
        Sets coordinates of anchor
        :param coords: coordinates of anchor
        """
        self.anchor[:] = np.array(coords)
        self.setDirectionDistance()

    def setFairleadCoords(self, coords):
        """
        Sets coordinates of fairlead
        :param coords: coordinates of fairlead
        """
        self.fairlead[:] = np.array(coords)
        self.setDirectionDistance()

    def setCoords(self):
        transf = np.linalg.inv(self.anchor_coords_system)
        anchor = np.dot(self.anchor, transf)
        fairlead = np.dot(self.anchor, transf)

    def getTension(self, s):
        """
        Gives the tension using line length and the coordinates of the anchor and fairlead
        :param tol: tolerance for calculation of tension
        :param maxit: maximum of iterations to calculate tension
        (!) works for seabed flat and perpendicular to gravity
        """
        if 0 <= s <= self.L:
            return self.catenary.getTension(s)


def get_array(x):
    if np.isscalar(x):
        x = np.array([x])
    else:
        x = np.asarray(x)
    return x
