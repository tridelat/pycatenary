from typing import Optional, Sequence, Union

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
        L: Union[float, Sequence[float]],
        w: Union[float, Sequence[float]],
        EA: Optional[Union[float, Sequence[float]]] = None,
        anchor: Optional[Sequence[float]] = None,
        fairlead: Optional[Sequence[float]] = None,
        nd: int = 3,
        floor: bool = True,
    ) -> None:
        self.__class__.count += 1
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
        if EA is None:
            self.catenary = catenary.CatenaryRigid(L=L, w=w, floor=floor)
        else:
            self.catenary = catenary.CatenaryElastic(
                L=L, w=w, EA=EA, floor=floor
            )
        self._setDirectionDistance()

    def updateAxialStiffness(self, EA: Union[float, Sequence[float]]) -> None:
        if isinstance(self.catenary, catenary.CatenaryElastic):
            EA = get_array(EA)
            old_len = len(self.catenary.EA)
            if len(EA) != old_len:
                raise ValueError(
                    f"Length of new EA is {len(EA)} (should be {old_len})."
                )
            self.catenary.EA = EA
        else:
            raise ValueError(
                "Catenary is not elastic, cannot update axial stiffness."
            )

    def computeSolution(self) -> None:
        """Computes solution of the catenary"""
        self.catenary.getState(
            d=self.distance_h,
            h=self.distance_v,
        )

    def s2xyz(self, s: float) -> np.ndarray:
        """Gives xyz coordinates along line

        Parameters
        ----------
        s: double
            distance along line (from anchor)
        """
        Lt = np.sum(self.catenary.L)
        assert (
            0.0 <= s <= Lt
        ), f"Cannot get position for s = {s} (should be 0.0 <= s <= L = {Lt})."
        if not self.fairlead_above_anchor:
            return self.fairlead + self._transformVector2D(
                self.catenary.s2xy(Lt - s)
            )
        else:
            return self.anchor + self._transformVector2D(self.catenary.s2xy(s))

    def getTension(self, s: float) -> np.ndarray:
        """Gives tension along line

        Parameters
        ----------
        s: double
            distance along line (from anchor)
        """
        Lt = np.sum(self.catenary.L)

        if not self.fairlead_above_anchor:
            return self._transformVector2D(self.catenary.getTension(Lt - s))
        else:
            return self._transformVector2D(self.catenary.getTension(s))

    def getTensionFairlead(self) -> np.ndarray:
        """Returns tension at fairlead."""
        return self.getTension(np.sum(self.catenary.L))

    def getTensionAnchor(self) -> np.ndarray:
        """Returns tension at anchor."""
        return self.getTension(0.0)

    def plot(
        self,
        npoints: int = 100,
        show_tension: bool = True,
        colormap: str = "viridis",
    ) -> None:
        """Plots line from anchor to fairlead.

        Parameters
        ----------
        npoints: int, optional
            Number of points along the line, by default 100.
        show_tension: bool, optional
            If True, color the line by tension magnitude, by default True.
        colormap: str, optional
            Matplotlib colormap name, by default "viridis".
        """
        if self.nd == 2:
            self.plot2D(
                npoints=npoints, show_tension=show_tension, colormap=colormap
            )
        else:
            self.plot3D(
                npoints=npoints, show_tension=show_tension, colormap=colormap
            )

    def plot2D(
        self,
        npoints: int = 100,
        show_tension: bool = True,
        colormap: str = "viridis",
    ) -> None:
        """Plots line from anchor to fairlead in 2D.

        Parameters
        ----------
        npoints: int, optional
            Number of points along the line, by default 100.
        show_tension: bool, optional
            If True, color the line by tension magnitude, by default True.
        colormap: str, optional
            Matplotlib colormap name, by default "viridis".
        """
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        xyzs = list()
        dd = list()
        hh = list()
        tensions = list()
        ss = np.linspace(0.0, np.sum(self.catenary.L), npoints)

        for s in ss:
            xyz = self.s2xyz(s)
            tension = self.getTension(s)
            xyzs.append(xyz)
            if self.nd == 2:
                dd.append(xyz[0])
                hh.append(xyz[1])
            else:
                dd.append(np.linalg.norm(xyz[:2] - self.anchor[:2]))
                hh.append(xyz[2])
            tensions.append(tension)

        if show_tension:
            from matplotlib.collections import LineCollection

            tension_magnitudes = np.linalg.norm(np.array(tensions), axis=1)
            # create segments
            points = np.array([dd, hh]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            # make a line with tension-based colors
            lc = LineCollection(segments, cmap=colormap, linewidths=2)
            lc.set_array(tension_magnitudes)
            line = ax.add_collection(lc)
            # add colorbar
            cbar = plt.colorbar(line, ax=ax)
            cbar.set_label("Tension Magnitude [N]")
        else:
            ax.plot(dd, hh)

        ax.grid("both")
        if self.nd == 2:
            ax.set_xlabel("x [m]")
            ax.set_ylabel("y [m]")
            ax.plot(self.anchor[0], self.anchor[1], "ko")
            ax.plot(self.fairlead[0], self.fairlead[1], "ko")
        else:
            ax.set_xlabel("distance from anchor [m]")
            ax.set_ylabel("z [m]")

            ax.plot(0.0, self.anchor[2], "ko")
            ax.plot(
                np.linalg.norm(self.fairlead[:2] - self.anchor[:2]),
                self.fairlead[2],
                "ko",
            )
        plt.show()

    def plot3D(
        self,
        npoints: int = 100,
        show_tension: bool = True,
        colormap: str = "viridis",
    ) -> None:
        """Plots line from anchor to fairlead in 3D.

        Parameters
        ----------
        npoints: int, optional
            Number of points along the line, by default 100.
        show_tension: bool, optional
            If True, color the line by tension magnitude, by default True.
        colormap: str, optional
            Matplotlib colormap name, by default "viridis".
        """
        import matplotlib.pyplot as plt

        if self.nd == 2:
            raise ValueError("3D plot not available for 2D cables.")

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        xyzs = list()
        xx = list()
        yy = list()
        zz = list()
        tensions = list()
        ss = np.linspace(0.0, np.sum(self.catenary.L), npoints)

        for s in ss:
            xyz = self.s2xyz(s)
            tension = self.getTension(s)
            xyzs.append(xyz)
            xx.append(xyz[0])
            yy.append(xyz[1])
            zz.append(xyz[2])
            tensions.append(tension)

        if show_tension:
            from mpl_toolkits.mplot3d.art3d import Line3DCollection

            tension_magnitudes = np.linalg.norm(np.array(tensions), axis=1)
            # create segments
            points = np.array([xx, yy, zz]).T.reshape(-1, 1, 3)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            # make a line with tension-based colors
            lc = Line3DCollection(segments, cmap=colormap, linewidths=2)
            lc.set_array(tension_magnitudes)
            line = ax.add_collection(lc)
            # add colorbar
            cbar = plt.colorbar(line, ax=ax, shrink=0.5, aspect=5)
            cbar.set_label("Tension Magnitude [N]")
        else:
            ax.plot(xx, yy, zz)

        ax.plot(self.anchor[0], self.anchor[1], self.anchor[2], "ko")
        ax.plot(self.fairlead[0], self.fairlead[1], self.fairlead[2], "ko")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        plt.show()

    def _setDirectionDistance(self) -> None:
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

    def _transformVector2D(self, vector: Sequence[float]) -> np.ndarray:
        """Transforms a 2D vector back in 3D (or 2D) according to direction

        Note that it is assumed that gravity acts in the Y direction in 2D,
        and Z direction in 3D
        """
        assert (
            len(vector) == 2
        ), f"Length of input vector is {len(vector)} (should be 2)."
        if self.nd == 2:
            return np.array([vector[0] * self.direction[0], vector[1]])
        elif self.nd == 3:
            vector3D = np.zeros(3)
            vector3D[0] = vector[0] * self.direction[0]
            vector3D[1] = vector[0] * self.direction[1]
            vector3D[2] = vector[1]
            return vector3D
        else:
            raise RuntimeError(f"Dimension nd = {self.nd} (should be 2 or 3).")

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


def get_array(x: Union[float, Sequence[float]]) -> np.ndarray:
    if np.isscalar(x):
        x = np.array([x])
    else:
        x = np.asarray(x)
    return x
