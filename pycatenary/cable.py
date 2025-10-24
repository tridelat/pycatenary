import numpy as np

from . import catenary


class MooringLine:
    """Class to create a mooring line

    Parameters
    ----------
    L: double
        unstretched line length [m]
    w: double
        submerged weight [N/m]
    EA: double
        axial stiffness
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
        EA=None,
        anchor=None,
        fairlead=None,
        nd=3,
        floor=True,
    ):
        self.__class__.count += 1
        self.L = get_array(L)
        self.w = get_array(w)
        self.EA = EA
        if self.EA is not None:
            self.EA = get_array(self.EA)
        self.nd = nd
        self.name = "cable_" + str(self.count)
        if anchor is None:
            self.anchor = np.zeros(nd)  # coordinates of anchor
        else:
            self.anchor = np.array(anchor)
        if fairlead is None:
            self.fairlead = np.zeros(nd)  # coordinates of fairlead
        else:
            self.fairlead = np.array(fairlead)
        self.anchor_coords_system = np.eye(3)
        self.floor = floor
        if self.EA is None:
            self.catenary = catenary.CatenaryRigid(self)
        else:
            self.catenary = catenary.CatenaryElastic(self)
        self._setDirectionDistance()

    def computeSolution(self):
        """Computes solution of the catenary"""
        self.catenary.getState(
            d=self.distance_h,
            h=self.distance_v,
            floor=self.floor,
        )

    def s2xyz(self, s):
        """Gives xyz coordinates along line

        Parameters
        ----------
        s: double
            distance along line (from anchor)
        """
        Lt = np.sum(self.L)
        assert (
            0.0 <= s <= Lt
        ), f"Cannot get position for s = {s} (should be 0.0 <= s <= L = {Lt})."
        if not self.fairlead_above_anchor:
            return self.fairlead + self._transformVector(
                self.catenary.s2xy(Lt - s)
            )
        else:
            return self.anchor + self._transformVector(self.catenary.s2xy(s))

    def getTension(self, s):
        """Gives tension along line

        Parameters
        ----------
        s: double
            distance along line (from anchor)
        """
        Lt = np.sum(self.L)
        assert (
            0.0 <= s <= Lt
        ), f"Cannot get tension for s = {s} (should be 0.0 <= s <= L = {Lt})."
        if not self.fairlead_above_anchor:
            return self._transformVector(self.catenary.getTension(Lt - s))
        else:
            return self._transformVector(self.catenary.getTension(s))

    def getTensionFairlead(self):
        """Returns tension at fairlead."""
        return self.getTension(np.sum(self.L))

    def getTensionAnchor(self):
        """Returns tension at anchor."""
        return self.getTension(0.0)

    def plot(self, npoints=100):
        """Plots line from anchor to fairlead"""
        self.plot3D(npoints=npoints)

    def plot2D(self, npoints=100):
        """Plots line from anchor to fairlead in 2D"""
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        xyzs = []
        dd = []
        hh = []
        ss = np.linspace(0.0, np.sum(self.L), npoints)
        for s in ss:
            xyz = self.s2xyz(s)
            xyzs.append(xyz)
            dd.append(np.sqrt(xyz[0] ** 2 + xyz[1] ** 2))
            hh.append(xyz[2])
        ax.plot(dd, hh)
        ax.grid("both")
        ax.set_xlabel("d")
        ax.set_ylabel("h")
        plt.show()

    def plot3D(self, npoints=100):
        """Plots line from anchor to fairlead in 3D"""
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        xyzs = []
        xx = []
        yy = []
        zz = []
        ss = np.linspace(0.0, np.sum(self.L), npoints)
        for s in ss:
            xyz = self.s2xyz(s)
            xyzs.append(xyz)
            xx.append(xyz[0])
            yy.append(xyz[1])
            zz.append(xyz[2])
        ax.plot(xx, yy, zz)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        plt.show()

    def _setDirectionDistance(self):
        if self.nd == 3:
            self.distance_h = np.sqrt(
                np.sum((self.fairlead[:2] - self.anchor[:2]) ** 2)
            )
            self.distance_v = np.abs(self.fairlead[2] - self.anchor[2])
            self.fairlead_above_anchor = self.fairlead[2] - self.anchor[2] > 0
            self.direction = (
                self.fairlead[:2] - self.anchor[:2]
            ) / self.distance_h
        elif self.nd == 2:
            if self.fairlead[0] - self.anchor[0] > 0:
                self.direction = np.array([1.0, 0.0])
            else:
                self.direction = np.array([-1.0, 0.0])
            self.distance_h = np.abs(self.fairlead[0] - self.anchor[0])
            self.distance_v = np.abs(self.fairlead[1] - self.anchor[1])
            self.fairlead_above_anchor = self.fairlead[1] - self.anchor[1] > 0
        if not self.fairlead_above_anchor:
            self.direction = -self.direction

    def _transformVector(self, vector):
        """Transforms a 2D vector back in 3D (or 2D) according to direction

        Note that it is assumed that gravity acts in the Y direction in 2D,
        and Z direction in 3D
        """
        assert len(vector) == 2, "must be 2D vector"
        if self.nd == 2:
            vector[0] *= self.direction[0]
            return np.array([vector[0], vector[1], 0.0])
        else:
            vector3D = np.zeros(3)
            vector3D[0] = vector[0] * self.direction[0]
            vector3D[1] = vector[0] * self.direction[1]
            vector3D[2] = vector[1]
            return vector3D

    def setAnchorCoords(self, coords):
        """Sets coordinates of anchor

        Parameters
        ----------
        coords: array
            coordinates of anchor
        """
        self.anchor[:] = np.array(coords)
        self._setDirectionDistance()

    def setFairleadCoords(self, coords):
        """Sets coordinates of fairlead

        Parameters
        ----------
        coords: array
            coordinates of fairlead
        """
        self.fairlead[:] = np.array(coords)
        self._setDirectionDistance()


def get_array(x):
    if np.isscalar(x):
        x = np.array([x])
    else:
        x = np.asarray(x)
    return x
