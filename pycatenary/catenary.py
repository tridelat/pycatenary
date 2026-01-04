import warnings
from abc import ABC, abstractmethod
from typing import Sequence, Union

import numpy as np

from . import root_finding


def get_array(x: Union[float, Sequence[float]]) -> np.ndarray:
    """Converts input to numpy array.

    Parameters
    ----------
    x: Union[float, Sequence[float]]
        Input to convert to numpy array.

    Returns
    -------
    x: np.ndarray
        Numpy array of the input.
    """
    if np.isscalar(x):
        x = np.array([x])
    else:
        x = np.asarray(x)
    return x


class CatenaryBase(ABC):
    """Base class for catenaries.

    Parameters
    ----------
    L: Union[float, Sequence[float]]
        Unstretched line length [m].
        If a list is provided, it is assumed to be a multisegmented cable.
    w: Union[float, Sequence[float]]
        Submerged weight [N/m].
        If a list is provided, it must match the length of the L list.
    floor: bool
        If True, the floor is assumed to be at the anchor level.
    """

    def __init__(
        self,
        L: Union[float, Sequence[float]],
        w: Union[float, Sequence[float]],
        floor: bool = True,
    ) -> None:
        # unstretched line length
        self.L = get_array(L)
        # submerged weight
        self.w = get_array(w)
        # check if lengths are the same
        if len(self.L) != len(self.w):
            raise ValueError("Length of L and w vectors must be the same.")
        # floor
        self.floor = floor
        # elongation
        self.e = np.zeros_like(self.L)
        # lifted line length
        self.Ls = np.zeros_like(self.L)
        # horizontal distance
        self.d = 0.0
        # vertical distance
        self.h = 0.0
        # catenary a
        self.a = 0.0
        # horizontal span
        self.x0 = 0.0
        # maximum number of iterations
        self.maxit = 1000
        # tolerance
        self.tol = 1e-10
        # first guess for bisection (int1)
        self.bisection_int1 = 1e-6
        # first guess for bisection (int2)
        self.bisection_int2 = 1e6
        # offset for x
        self._x_offset = 0.0
        # offset for y
        self._y_offset = 0.0
        # offset for s
        self._s_offset = 0.0
        # properties not reversed initially
        self._has_reversed_properties = False

    def get_tension(self, s: float) -> np.ndarray:
        """Returns tension at a given distance along the line from the anchor.

        Parameters
        ----------
        s: float
            Distance along line [m].

        Returns
        -------
        tension: np.ndarray
            Tension vector [N].
        """
        s0 = self.d - self.x0
        # total line lengths
        Lt = np.sum(self.L)  # unstretched
        Lst = np.sum(self.Ls)  # unstretched
        Lset = Lst + np.sum(self.e)  # stretched
        if self.x0 == 0.0:  # line straight to seabed
            # length of line on floor
            L_floor = Lt - np.sum(self.Ls)
            if s < L_floor:
                Ts = np.array([0.0, 0.0])
            else:
                # average w
                w_av = np.sum(self.w * self.Ls) / Lst
                Tv = (s - L_floor) * w_av
                Ts = np.array([0.0, Tv])
        elif Lt >= s >= s0:  # s in lifted line part
            # average w
            w_av = np.sum(self.w * self.Ls) / Lst
            # horizontal tension
            Th = self.a * w_av * (Lst / Lset)
            # vertical tension
            dydx = np.sinh((self.s2xy(s)[0] - self._x_offset - s0) / self.a)
            angle = np.arctan(dydx)
            Tv = Th * np.tan(angle)
            # Tv assumed always negative
            Tv = np.abs(Tv)
            # tension at point
            Ts = np.array([Th, Tv])
        elif 0 <= s < s0:  # s on floor
            # average w
            w_av = np.sum(self.w * self.Ls) / Lst
            # horizontal tension
            Th = self.a * w_av * (Lst / Lset)
            # tension at point
            Ts = np.array([Th, 0.0])
        else:
            raise RuntimeError(
                f"Cannot get tension for s = {s} (should be 0.0 <= s <= {Lt})."
            )
        return Ts

    def s2xy(self, s: float) -> np.ndarray:
        """Returns [x,y] coords at a given distance along line.

        Parameters
        ----------
        s: float
            Distance along line [m].

        Returns
        -------
        xy: np.ndarray
            [x,y] coordinates.
        """
        s0 = self.d - self.x0
        Lt = np.sum(self.L)
        if self.x0 == 0.0:  # line straight to seabed
            # length of line on floor
            L_floor = Lt - np.sum(self.Ls)
            if s < L_floor:
                x = s * self.d / L_floor
                y = 0.0
            else:
                x = self.d
                y = s - L_floor + self._get_elongation_at_s(s)
        elif (
            0.0 <= s < s0 and self.floor
        ):  # line partly lifted, with s on floor
            x = s
            y = 0.0 - self._y_offset
        elif 0.0 <= s <= Lt:  # s in lifted line part
            s += self._get_elongation_at_s(s)
            # add offset from catenary
            s = s + self._s_offset
            # calculate x and y coordinates
            a = self.a
            x = s0 + a * np.arcsinh((s - s0) / a)
            y = a * np.cosh((x - s0) / a)
        else:
            raise RuntimeError(
                f"Cannot get coords for s = {s} (should be 0.0 <= s <= {Lt})."
            )
        xy = np.array([x + self._x_offset, y + self._y_offset])
        return xy

    def plot(
        self,
        npoints: int = 100,
        show_tension: bool = True,
        colormap: str = "viridis",
    ) -> None:
        """Plots catenary in 2D from (0, 0) to (d, h).

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

        fig, ax = plt.subplots()
        xys = list()
        xx = list()
        yy = list()
        tensions = list()
        ss = np.linspace(0.0, np.sum(self.L), npoints)

        for s in ss:
            xy = self.s2xy(s)
            tension = self.get_tension(s)
            xys.append(xy)
            xx.append(xy[0])
            yy.append(xy[1])
            tensions.append(tension)

        if show_tension:
            from matplotlib.collections import LineCollection

            tension_magnitudes = np.linalg.norm(np.array(tensions), axis=1)
            # create segments
            points = np.array([xx, yy]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            # make a line with tension-based colors
            lc = LineCollection(segments, cmap=colormap, linewidths=2)
            lc.set_array(tension_magnitudes)
            line = ax.add_collection(lc)
            # add colorbar
            cbar = plt.colorbar(line, ax=ax)
            cbar.set_label("Tension Magnitude")
            # Set axis limits to show the line
            ax.set_xlim(min(xx), max(xx))
            ax.set_ylim(min(yy), max(yy))
        else:
            ax.plot(xx, yy)

        ax.grid("both")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        plt.show()

    def _get_elongation_at_s(self, s: float) -> float:
        """Returns total elongation at a given distance along line.

        Parameters
        ----------
        s: float
            Distance along line [m].

        Returns
        -------
        elongation: float
            Total elongation [m].
        """
        for ii in range(len(self.L)):
            if s <= np.sum(self.L[: ii + 1]):
                if self.Ls[ii] > 0:
                    # distance along segment
                    s_segment = self.L[ii] - (np.sum(self.L[: ii + 1]) - s)
                    if s_segment > self.L[ii] - self.Ls[ii]:
                        # distance along lifted part of segment
                        s_segment = s_segment - (self.L[ii] - self.Ls[ii])
                        if s_segment > 0.0:  # in lifted part
                            # line lifted at s --> elongation
                            s_frac = s_segment / self.Ls[ii]
                            return np.sum(self.e[:ii]) + s_frac * self.e[ii]
                        else:  # line not lifted at s --> no elongation
                            return 0.0
                    else:
                        # line not lifted at s --> no elongation
                        return 0.0
                else:
                    # line not lifted yet at s --> no elongation
                    return 0.0

    @abstractmethod
    def compute_solution(self, d: float, h: float) -> None:
        """Abstract method to calculate the catenary solution.

        This method must be implemented by subclasses to define the specific
        catenary calculation algorithm (rigid or elastic).

        Parameters
        ----------
        d: float
            Horizontal distance between anchor and fairlead [m].
        h: float
            Vertical distance between anchor and fairlead [m].
        """
        pass

    def _reverse_properties(self) -> None:
        """Reverses the properties of the catenary (for internal use)."""
        self.L = self.L[::-1]
        self.w = self.w[::-1]
        self.e = self.e[::-1]
        self.Ls = self.Ls[::-1]
        self._has_reversed_properties = not self._has_reversed_properties

    def get_force_beginning_of_line(self):
        Lt = np.sum(self.L)
        if not self._has_reversed_properties:
            s = 0.0
        else:
            s = Lt

        force = self.get_tension(s)

        s0 = self.d - self.x0
        if Lt >= s >= s0:  # s in lifted line part
            if s == Lt:
                force[0] = -force[0]
            if not (s == 0.0 and self._s_offset > 0.0):  # fully lifted line
                # inverse sign of vertical component in all cases
                # apart from case where line is fully lifted and not hanging
                force[1] = -force[1]
        return force

    def get_force_end_of_line(self):
        Lt = np.sum(self.L)
        if not self._has_reversed_properties:
            s = Lt
        else:
            s = 0.0

        force = self.get_tension(s)

        s0 = self.d - self.x0
        if Lt >= s >= s0:  # s in lifted line part
            # inverse sign of horizontal component if on right side of catenary
            if s == Lt:
                force[0] = -force[0]
            if not (s == 0.0 and self._s_offset > 0.0):  # fully lifted line
                # inverse sign of vertical component in all cases
                # apart from case where line is fully lifted and not hanging
                force[1] = -force[1]
        return force

    # Deprecated camelCase methods with warnings
    def getTension(self, s: float) -> np.ndarray:
        import warnings

        warnings.warn(
            "getTension is deprecated, use get_tension.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.get_tension(s)

    def getState(self, d: float, h: float) -> None:
        warnings.warn(
            "getState is deprecated, use compute_solution.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.compute_solution(d, h)


class CatenaryRigid(CatenaryBase):
    """A class for rigid catenary.

    Parameters
    ----------
    L: Union[float, Sequence[float]]
        Unstretched line length [m].
        If a list is provided, it is assumed to be a multisegmented cable.
    w: Union[float, Sequence[float]]
        Submerged weight [N/m].
        If a list is provided, it must match the length of the L list.
    floor: bool
        If True, the floor is assumed to be at the anchor level.
    """

    def __init__(
        self,
        L: Union[float, Sequence[float]],
        w: Union[float, Sequence[float]],
        floor: bool = True,
    ) -> None:
        super(CatenaryRigid, self).__init__(L=L, w=w, floor=floor)

    def compute_solution(self, d: float, h: float) -> None:
        """Calculates the solution for rigid catenary.

        Parameters
        ----------
        d: float
            Horizontal distance between anchor and fairlead [m].
        h: float
            Vertical distance between anchor and fairlead [m].
        """
        self.d = d
        self.h = h
        tol = self.tol
        maxit = self.maxit
        L = self.L
        self.e = np.zeros(len(L))
        floor = self.floor
        Ls = np.zeros(len(L))
        Lt = np.sum(L)
        x_offset = 0.0
        y_offset = 0.0
        s_offset = 0.0
        a = 1.0
        x0 = 0.0
        a2f = np.sqrt(h**2 + d**2)  # distance between anchor and fairlead
        if Lt + tol < a2f:
            raise RuntimeError(
                f"Cannot find solution for rigid line of length ({Lt})"
                f" inferior to distance between anchor and fairlead "
                f"({a2f}), delta={Lt - a2f} < tol={tol}."
            )
        if floor is False:
            a = root_finding.nofloor_rigid(d=d, h=h, L=L, tol=tol, maxit=maxit)
            x0 = d
            Ls[:] = L
            Lst = np.sum(Ls + self.e)
            xx = 0.5 * (a * np.log((Lst + h) / (Lst - h)) - d)
            xy = 0.5 * (a * np.log((Lst + h) / (Lst - h)) + d)
            x_offset = -xx
            y_offset = h - a * np.cosh(xy / a)
            s_offset = a * np.sinh(xx / a)
        else:
            # cable straight to seabed:
            if np.sum(L) + tol >= h + d:
                # no horizontal tension
                a = 0.0
                x0 = 0.0
                Lst = 0.0
                for ii in reversed(range(len(L))):
                    if Lst >= h:
                        break
                    Lst += L[ii]
                    if Lst < h:
                        Ls[ii] = L[ii]
                    if Lst >= h:
                        Ls[ii] = L[ii] - (Lst - h)
            else:
                # check if line is partly or fully lifted
                f = lambda a: a * (np.cosh(d / a) - 1) - h
                a = root_finding.bisection(
                    f,
                    int1=self.bisection_int1,
                    int2=self.bisection_int2,
                    tol=tol,
                    maxit=maxit,
                )
                Ls0 = a * np.sinh(
                    d / a
                )  # maximum line length to be fully lifted
                # get actual line length assuming fully lifted (from a)
                Ls1 = np.sum(L)
                if Ls1 > Ls0:  # partly lifted
                    a, Ls = root_finding.partly_lifted_rigid(
                        d=d, h=h, L=L, maxit=maxit, tol=tol
                    )
                    x0 = a * np.arccosh(1 + h / a)
                    y_offset = -a
                elif Ls1 <= Ls0:  # fully lifted
                    a = root_finding.fully_lifted_rigid(
                        d=d,
                        h=h,
                        L=L,
                        maxit=maxit,
                        tol=tol,
                        must_converge=False,
                    )
                    if a is np.nan:
                        raise RuntimeError(
                            "Line is too stretched, cannot solve catenary."
                        )
                    Ls[:] = L
                    Lst = Lt
                    x0 = d
                    xx = 0.5 * (a * np.log((Lst + h) / (Lst - h)) - d)
                    xy = 0.5 * (a * np.log((Lst + h) / (Lst - h)) + d)
                    x_offset = -xx
                    y_offset = h - a * np.cosh(xy / a)
                    s_offset = a * np.sinh(xx / a)
        self.Ls = Ls
        self._x_offset = x_offset
        self._y_offset = y_offset
        self._s_offset = s_offset
        self.a = a
        self.x0 = x0


class CatenaryElastic(CatenaryBase):
    """A class for elastic catenary

    Parameters
    ----------
    L: Union[float, Sequence[float]]
        Unstretched line length [m].
        If a list is provided, it is assumed to be a multisegmented cable.
    w: Union[float, Sequence[float]]
        Submerged weight [N/m].
        If a list is provided, it must match the length of the L list.
    EA: Union[float, Sequence[float]]
        Axial stiffness [N].
        If a list is provided, it must match the length of the L list.
    floor: bool
        If True, the floor is assumed to be at the anchor level.
    """

    def __init__(
        self,
        L: Union[float, Sequence[float]],
        w: Union[float, Sequence[float]],
        EA: Union[float, Sequence[float]] = None,
        floor: bool = True,
    ) -> None:
        super(CatenaryElastic, self).__init__(L=L, w=w, floor=floor)
        # axial stiffness
        self.EA = get_array(EA)
        # check if lengths are the same
        if len(self.L) != len(self.EA):
            raise ValueError("Length of L and EA vectors must be the same.")

    def _reverse_properties(self) -> None:
        """Reverses the properties of the catenary (for internal use)."""
        super(CatenaryElastic, self)._reverse_properties()
        self.EA = self.EA[::-1]

    def compute_solution(self, d: float, h: float) -> None:
        """Calculates the solution for elastic catenary.

        Parameters
        ----------
        d: float
            Horizontal distance between anchor and fairlead [m].
        h: float
            Vertical distance between anchor and fairlead [m].
        """
        self.d = d
        self.h = h
        tol = self.tol
        maxit = self.maxit
        L = self.L
        w = self.w
        EA = self.EA
        floor = self.floor

        Lt = np.sum(L)  # total unstretched line length
        Ls = np.zeros(len(L))  # unstretched lifted line length
        e = np.zeros(len(L))  # stretching
        x_offset = 0.0
        y_offset = 0.0
        s_offset = 0.0

        diff = tol + 1

        if floor is False:
            a, e = root_finding.nofloor_elastic(
                d=d, h=h, L=L, w=w, EA=EA, tol=tol, maxit=maxit
            )
            x0 = d
            Ls[:] = L
            Lst = np.sum(Ls + e)
            xx = 0.5 * (a * np.log((Lst + h) / (Lst - h)) - d)
            xy = 0.5 * (a * np.log((Lst + h) / (Lst - h)) + d)
            x_offset = -xx
            y_offset = h - a * np.cosh(xy / a)
            s_offset = a * np.sinh(xx / a)
        else:
            # cable straight to seabed:
            # find tension and stretching
            for i in reversed(range(len(L))):
                if np.sum(Ls + e) + tol >= h:
                    break
                else:
                    Ls[i] = L[i]
                    for j in range(i, len(L)):
                        e[j] = (
                            (w[j] * Ls[j] / 2.0 + np.sum(w[i:j] * Ls[i:j]))
                            * Ls[j]
                            / EA[j]
                        )
                    if np.sum(Ls + e) >= h:
                        Lhi_low = 0
                        Lhi_high = L[i]
                        while diff > tol:
                            Ls[i] = (Lhi_low + Lhi_high) / 2.0
                            for j in range(i, len(L)):
                                e[j] = (
                                    (
                                        w[j] * Ls[j] / 2.0
                                        + np.sum(w[i:j] * Ls[i:j])
                                    )
                                    * Ls[j]
                                    / EA[j]
                                )
                            if np.sum(Ls + e) > h:
                                Lhi_high = Ls[i]
                            elif np.sum(Ls + e) <= h:
                                Lhi_low = Ls[i]
                            diff = np.abs(np.sum(Ls + e) - h)
            # check if cable straight to seabed is solution
            if np.sum(L + e) + tol >= h + d:
                # no horizontal tension
                a = 0.0
                x0 = 0.0
            else:
                # check if line is partly or fully lifted
                f = lambda a: a * (np.cosh(d / a) - 1) - h
                a = root_finding.bisection(
                    f,
                    self.bisection_int1,
                    self.bisection_int2,
                    tol=tol,
                    maxit=maxit,
                )
                Ls0 = a * np.sinh(
                    d / a
                )  # maximum line length to be fully lifted
                # get actual line length assuming fully lifted (from a)
                H = a * np.sum(w * L) / Lt
                Va = 0
                for i in range(len(e)):
                    e[i] = (
                        np.sqrt(
                            H**2
                            + (Va + np.sum(w[:i] * L[:i]) + w[i] * L[i] / 2.0)
                            ** 2
                        )
                        * L[i]
                        / EA[i]
                    )
                Ls1 = Lt + np.sum(e)
                if Ls1 > Ls0:  # partly lifted
                    a, e, Lsu = root_finding.partly_lifted_elastic(
                        d=d, h=h, L=L, w=w, EA=EA, maxit=maxit, tol=tol
                    )
                    Ls[:] = Lsu
                    x0 = a * np.arccosh(1 + h / a)
                    y_offset = -a
                elif Ls1 <= Ls0:  # fully lifted
                    x0 = d
                    Ls[:] = L
                    a, e = root_finding.fully_lifted_elastic(
                        d=d,
                        h=h,
                        L=L,
                        w=w,
                        EA=EA,
                        int1=a,
                        maxit=maxit,
                        tol=tol,
                        must_converge=False,
                    )
                    if a is np.nan:
                        raise RuntimeError(
                            "Line is too stretched, cannot solve catenary."
                        )
                    Lst = np.sum(Ls + e)
                    xx = 0.5 * (a * np.log((Lst + h) / (Lst - h)) - d)
                    xy = 0.5 * (a * np.log((Lst + h) / (Lst - h)) + d)
                    x_offset = -xx
                    y_offset = h - a * np.cosh(xy / a)
                    s_offset = a * np.sinh(xx / a)
        self.Ls = Ls
        self._x_offset = x_offset
        self._y_offset = y_offset
        self._s_offset = s_offset
        self.a = a
        self.e = e
        self.x0 = x0
