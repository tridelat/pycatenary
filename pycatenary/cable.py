from typing import Optional, Sequence, Union

import numpy as np

from . import catenary
from .utils import deprecated


class MooringLine:
    """Class to create a mooring line.

    Mooring lines can be elastic or rigid, and multisegmented or not.
    If the line is multisegmented, the properties must be given in a list
    going from the anchor to the fairlead.

    Parameters
    ----------
    fairlead: sequence of floats
        Fairlead coordinates [x, y, z] (3D) or [x, y] (2D).
    anchor: sequence of floats
        Anchor coordinates [x, y, z] (3D) or [x, y] (2D).
    L: float, or sequence of floats
        Unstretched line length [m].
        If a list is provided, it is assumed to be a multisegmented cable.
    w: float, or sequence of floats
        Submerged weight [N/m].
        If a list is provided, it must match the length of the L list.
    EA: float, or sequence of floats
        Axial stiffness [N].
        Must be provided if the cable is elastic.
        If EA is None, the cable is assumed to be rigid.
        If a list is provided, it must match the length of the L list.
    floor: bool
        If True, the floor is assumed to be at the anchor level.
        If fairlead is below anchor, the floor will be at fairlead level.
    """

    def __init__(
        self,
        fairlead: Sequence[float],
        anchor: Sequence[float],
        L: Union[float, Sequence[float]],
        w: Union[float, Sequence[float]],
        EA: Optional[Union[float, Sequence[float]]] = None,
        floor: bool = True,
    ) -> None:
        self._nd = len(fairlead)
        assert (
            len(anchor) == self._nd
        ), "Anchor and fairlead vectors must have the same length."
        if not 2 <= self._nd <= 3:
            raise ValueError(
                "Invalid fairlead or anchor vector length"
                "(should be 2 or 3)."
            )
        self._anchor = np.array(anchor)
        self._fairlead = np.array(fairlead)
        if EA is None:
            self.catenary = catenary.CatenaryRigid(L=L, w=w, floor=floor)
        else:
            self.catenary = catenary.CatenaryElastic(
                L=L, w=w, EA=EA, floor=floor
            )
        self._set_direction_distance()

    def update_axial_stiffness(
        self, EA: Union[float, Sequence[float]]
    ) -> None:
        """Updates the axial stiffness of the cable.

        Parameters
        ----------
        EA: float, or sequence of floats
            Axial stiffness [N].
            Must be of the same length as the original EA list.
        """
        if isinstance(self.catenary, catenary.CatenaryElastic):
            EA = catenary.get_array(EA)
            old_len = len(self.catenary.EA)
            if len(EA) != old_len:
                raise ValueError(
                    f"Length of new EA is {len(EA)} (should be {old_len})."
                )
            if self.catenary._has_reversed_properties:
                self.catenary.EA = EA[::-1]
            else:
                self.catenary.EA[:] = EA
        else:
            raise ValueError(
                "Catenary is rigid, cannot update axial stiffness."
            )

    def compute_solution(self) -> None:
        """Computes solution of the catenary.

        It is computed according to current anchor and fairlead positions."""
        self.catenary.compute_solution(
            d=self.distance_h,
            h=self.distance_v,
        )

    def get_position(
        self, s: float, from_fairlead: bool = False
    ) -> np.ndarray:
        """Returns position at a given distance along line from 0 to L.

        Parameters
        ----------
        s: float
            Distance along line [m].
        from_fairlead: bool, optional
            False (by default): distance is measured from the anchor,
            True: distance is measured from the fairlead.

        Returns
        -------
        position: np.ndarray
            Position [x, y, z] (3D) or [x, y] (2D).
        """
        if from_fairlead:
            s = np.sum(self.catenary.L) - s
        return self._s2xyz(s)

    def _s2xyz(self, s: float) -> np.ndarray:
        """Returns xyz coordinates at a given distance line from anchor.

        Parameters
        ----------
        s: float
            Distance along line (from anchor) [m].

        Returns
        -------
        xyz: np.ndarray
            xyz coordinates [x, y, z] (3D) or [x, y] (2D).
        """
        if not self._fairlead_above_anchor:
            return self._fairlead + self._transform_vector_2d(
                self.catenary.s2xy(np.sum(self.catenary.L) - s)
            )
        else:
            return self._anchor + self._transform_vector_2d(
                self.catenary.s2xy(s)
            )

    def get_tension(self, s: float, from_fairlead: bool = False) -> np.ndarray:
        """Returns tension at a given distance along line from 0 to L.

        Parameters
        ----------
        s: float
            Distance along line [m].
        from_fairlead: bool, optional
            False (by default): distance is measured from the anchor,
            True: distance is measured from the fairlead.

        Returns
        -------
        tension: np.ndarray
            Tension vector [N].
        """
        if from_fairlead:
            s = np.sum(self.catenary.L) - s

        if not self._fairlead_above_anchor:
            tension = self._transform_vector_2d(
                self.catenary.get_tension(np.sum(self.catenary.L) - s)
            )
        else:
            tension = self._transform_vector_2d(self.catenary.get_tension(s))

        if self._nd == 2:
            tension[0] = abs(tension[0])
        elif self._nd == 3:
            tension[0] = abs(tension[0])
            tension[1] = abs(tension[1])
        return tension

    def get_fairlead_force(self) -> np.ndarray:
        """Returns force at fairlead.

        Returns
        -------
        force: np.ndarray
            Tension vector [N].
        """
        return self._transform_vector_2d(self.catenary.get_force_end_of_line())

    def get_anchor_force(self) -> np.ndarray:
        """Returns force at anchor.

        Returns
        -------
        force: np.ndarray
            Tension vector [N].
        """
        return self._transform_vector_2d(
            self.catenary.get_force_beginning_of_line()
        )

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
        if self._nd == 2:
            self.plot_2d(
                npoints=npoints, show_tension=show_tension, colormap=colormap
            )
        else:
            self.plot_3d(
                npoints=npoints, show_tension=show_tension, colormap=colormap
            )

    def plot_2d(
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

        Raises
        ------
        ImportError
            If matplotlib is not installed.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "matplotlib is required for plotting. Install with: "
                "pip install matplotlib or pip install pycatenary[plotting]"
            )

        fig = plt.figure()
        ax = fig.add_subplot(111)
        xyzs = list()
        dd = list()
        hh = list()
        tensions = list()
        ss = np.linspace(0.0, np.sum(self.catenary.L), npoints)

        for s in ss:
            xyz = self._s2xyz(s)
            tension = self.get_tension(s)
            xyzs.append(xyz)
            if self._nd == 2:
                dd.append(xyz[0])
                hh.append(xyz[1])
            else:
                dd.append(np.linalg.norm(xyz[:2] - self._anchor[:2]))
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
            cbar.set_label("Tension Magnitude")
        else:
            ax.plot(dd, hh)

        ax.grid("both")
        if self._nd == 2:
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.plot(self._anchor[0], self._anchor[1], "ko")
            ax.plot(self._fairlead[0], self._fairlead[1], "ko")
        else:
            ax.set_xlabel("distance from anchor")
            ax.set_ylabel("z")
            ax.plot(0.0, self._anchor[2], "ko")
            ax.plot(
                np.linalg.norm(self._fairlead[:2] - self._anchor[:2]),
                self._fairlead[2],
                "ko",
            )
        # add tension information
        anchor_tension = np.linalg.norm(self.get_anchor_force())
        fairlead_tension = np.linalg.norm(self.get_fairlead_force())
        ax.set_title(
            f"Tensions: Fairlead {fairlead_tension:.3e} | "
            f"Anchor {anchor_tension:.3e}"
        )
        plt.show()

    def plot_3d(
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

        Raises
        ------
        ImportError
            If matplotlib is not installed.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "matplotlib is required for plotting. Install with: "
                "pip install matplotlib or pip install pycatenary[plotting]"
            )

        if self._nd == 2:
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
            xyz = self._s2xyz(s)
            tension = self.get_tension(s)
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
            cbar.set_label("Tension Magnitude")
        else:
            ax.plot(xx, yy, zz)

        ax.plot(self._anchor[0], self._anchor[1], self._anchor[2], "ko")
        ax.plot(self._fairlead[0], self._fairlead[1], self._fairlead[2], "ko")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_zlim(bottom=min(zz), top=max(zz))
        # add tension information
        anchor_tension = np.linalg.norm(self.get_anchor_force())
        fairlead_tension = np.linalg.norm(self.get_fairlead_force())
        ax.set_title(
            f"Tensions: Fairlead {fairlead_tension:.3e} | "
            f"Anchor {anchor_tension:.3e}"
        )
        plt.show()

    def _set_direction_distance(self) -> None:
        """Sets the direction and distance between the anchor and the fairlead

        For internal use only, do not call this method directly."""
        if (
            self._fairlead[0] - self._anchor[0] == 0.0
            and self._fairlead[1] - self._anchor[1] == 0.0
            and self._fairlead[2] - self._anchor[2] == 0.0
        ):
            raise ValueError("Anchor and fairlead are at the same position.")
        if self._nd == 3:
            self.distance_h = np.sqrt(
                np.sum((self._fairlead[:2] - self._anchor[:2]) ** 2)
            )
            self.distance_v = np.abs(self._fairlead[2] - self._anchor[2])
            self._fairlead_above_anchor = (
                self._fairlead[2] - self._anchor[2] > 0
            )
            self.direction = (
                self._fairlead[:2] - self._anchor[:2]
            ) / self.distance_h
        elif self._nd == 2:
            if self._fairlead[0] - self._anchor[0] > 0:
                self.direction = np.array([1.0, 0.0])
            else:
                self.direction = np.array([-1.0, 0.0])
            self.distance_h = np.abs(self._fairlead[0] - self._anchor[0])
            self.distance_v = np.abs(self._fairlead[1] - self._anchor[1])
            self._fairlead_above_anchor = (
                self._fairlead[1] - self._anchor[1] > 0
            )
        if not self._fairlead_above_anchor:
            self.direction = -self.direction

        # reverse properties for catenary if necessary
        if (
            self._fairlead_above_anchor
            and self.catenary._has_reversed_properties
        ):
            self.catenary._reverse_properties()
        elif (
            not self._fairlead_above_anchor
            and not self.catenary._has_reversed_properties
        ):
            self.catenary._reverse_properties()

    def _transform_vector_2d(self, vector: Sequence[float]) -> np.ndarray:
        """Transforms a 2D vector back in 3D (or 2D) according to direction

        Note that it is assumed that gravity acts in the Y direction in 2D,
        and Z direction in 3D
        """
        assert (
            len(vector) == 2
        ), f"Length of input vector is {len(vector)} (should be 2)."
        if self._nd == 2:
            return np.array([vector[0] * self.direction[0], vector[1]])
        elif self._nd == 3:
            vector3D = np.zeros(3)
            vector3D[0] = vector[0] * self.direction[0]
            vector3D[1] = vector[0] * self.direction[1]
            vector3D[2] = vector[1]
            return vector3D
        else:
            raise RuntimeError(
                f"Dimension nd = {self._nd} (should be 2 or 3)."
            )

    def set_anchor_position(self, position: Sequence[float]) -> None:
        """Sets coordinates of anchor.

        Parameters
        ----------
        position: sequence of floats
            Anchor position [x, y, z] (3D) or [x, y] (2D).
        """
        self._anchor[:] = np.array(position)
        self._set_direction_distance()

    def get_anchor_position(self) -> np.ndarray:
        """Returns coordinates of anchor.

        Returns
        -------
        position: np.ndarray
            Anchor position [x, y, z] (3D) or [x, y] (2D).
        """
        return np.array(self._anchor)

    def set_fairlead_position(self, position: Sequence[float]) -> None:
        """Sets coordinates of fairlead.

        Parameters
        ----------
        coords: sequence of floats
            Fairlead position [x, y, z] (3D) or [x, y] (2D).
        """
        self._fairlead[:] = np.array(position)
        self._set_direction_distance()

    def get_fairlead_position(self) -> np.ndarray:
        """Returns coordinates of fairlead.

        Returns
        -------
        position: np.ndarray
            Fairlead position [x, y, z] (3D) or [x, y] (2D).
        """
        return np.array(self._fairlead)

    @deprecated("computeSolution is deprecated, use compute_solution instead.")
    def computeSolution(self) -> None:
        return self.compute_solution()

    @deprecated("getTension is deprecated, use get_tension instead.")
    def getTension(self, s: float) -> np.ndarray:
        return self.get_tension(s)

    @deprecated(
        "setAnchorCoords is deprecated, use set_anchor_position instead."
    )
    def setAnchorCoords(self, coords: Sequence[float]) -> None:
        return self.set_anchor_position(coords)

    @deprecated(
        "setFairleadCoords is deprecated, use set_fairlead_position instead."
    )
    def setFairleadCoords(self, coords: Sequence[float]) -> None:
        return self.set_fairlead_position(coords)

    @deprecated("plot2D is deprecated, use plot_2d instead.")
    def plot2D(
        self,
        npoints: int = 100,
        show_tension: bool = True,
        colormap: str = "viridis",
    ) -> None:
        return self.plot_2d(npoints, show_tension, colormap)

    @deprecated("plot3D is deprecated, use plot_3d instead.")
    def plot3D(
        self,
        npoints: int = 100,
        show_tension: bool = True,
        colormap: str = "viridis",
    ) -> None:
        return self.plot_3d(npoints, show_tension, colormap)

    @deprecated("s2xyz is deprecated, use get_position instead.")
    def s2xyz(self, s: float) -> np.ndarray:
        return self.get_position(s)
