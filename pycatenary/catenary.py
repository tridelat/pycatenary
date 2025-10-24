from typing import Sequence, Union

import numpy as np

from . import utils


def get_array(x: Union[float, Sequence[float]]) -> np.ndarray:
    if np.isscalar(x):
        x = np.array([x])
    else:
        x = np.asarray(x)
    return x


class CatenaryBase(object):
    """Base class for catenaries

    Parameters
    ----------
    L: Union[float, Sequence[float]]
        unstretched line length [m]
    w: Union[float, Sequence[float]]
        submerged weight [N/m]
    floor: bool
        if True, the floor is assumed to be at the anchor level
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

    def getTension(self, s: float) -> np.ndarray:
        s0 = self.d - self.x0
        # total line lengths
        Lt = np.sum(self.L)  # unstretched
        Lst = np.sum(self.Ls)  # unstretched
        Lset = Lst + np.sum(self.e)  # stretched
        if Lt >= s >= s0:
            # average w
            w_av = np.sum(self.w * self.Ls) / Lst
            # horizontal tension
            Th = self.a * w_av * (Lst / Lset)
            # reverse sign for Th if s > 0 for catenary
            if s + self._s_offset > 0.0:
                Th = -Th
            # vertical tension
            dydx = np.sinh((self.s2xy(s)[0] - self._x_offset - s0) / self.a)
            angle = np.arctan(dydx)
            Tv = Th * np.tan(angle)
            # Tv assumed always negative
            Tv = -np.abs(Tv)
            # tension at point
            Ts = np.array([Th, Tv])
        elif 0 <= s < s0:
            Ts = np.array([0.0, 0.0])
        else:
            raise RuntimeError(
                f"Cannot get tension for s = {s} (should be 0.0 <= s <= {Lt})."
            )
        return Ts

    def s2xy(self, s: float) -> np.ndarray:
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
        """
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        xys = list()
        xx = list()
        yy = list()
        tensions = list()
        ss = np.linspace(0.0, np.sum(self.L), npoints)

        for s in ss:
            xy = self.s2xy(s)
            tension = self.getTension(s)
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
        for ii in range(len(self.L)):
            if s <= np.sum(self.L[: ii + 1]):
                s_frac = 1 - (np.sum(self.L[: ii + 1]) - s) / self.L[ii]
                return np.sum(self.e[:ii]) + s_frac * self.e[ii]
        raise RuntimeError("Could not calculate elongation along line.")


class CatenaryRigid(CatenaryBase):
    """A class for rigid catenary

    Parameters
    ----------
    L: Union[float, Sequence[float]]
        unstretched line length [m]
    w: Union[float, Sequence[float]]
        submerged weight [N/m]
    floor: bool
        if True, the floor is assumed to be at the anchor level
    """

    def __init__(
        self,
        L: Union[float, Sequence[float]],
        w: Union[float, Sequence[float]],
        floor: bool = True,
    ) -> None:
        super(CatenaryRigid, self).__init__(L=L, w=w, floor=floor)

    def getState(self, d: float, h: float) -> None:
        """Calculates the solution for rigid catenary

        Parameters
        ----------
        d: double
            horizontal distance between anchor and fairlead
        h: double
            vertical distance between anchor and fairlead
        floor: bool
            if True, the floor is assumed to be at the anchor level
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
            a = utils.nofloor_rigid(d=d, h=h, L=L, tol=tol, maxit=maxit)
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
                        Ls[ii] = Lst - h
            else:
                # check if line is partly or fully lifted
                f = lambda a: a * (np.cosh(d / a) - 1) - h
                a = utils.bisection(
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
                    a, Ls = utils.partly_lifted_rigid(
                        d=d, h=h, L=L, maxit=maxit, tol=tol
                    )
                    x0 = a * np.arccosh(1 + h / a)
                    y_offset = -a
                elif Ls1 <= Ls0:  # fully lifted
                    a = utils.fully_lifted_rigid(
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
        unstretched line length [m]
    w: Union[float, Sequence[float]]
        submerged weight [N/m]
    EA: Union[float, Sequence[float]]
        axial stiffness
    floor: bool
        if True, the floor is assumed to be at the anchor level
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

    def getState(self, d: float, h: float) -> None:
        """Calculates the solution for elastic catenary

        Parameters
        ----------
        d: double
            horizontal distance between anchor and fairlead
        h: double
            vertical distance between anchor and fairlead
        floor: bool
            if True, the floor is assumed to be at the anchor level
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
            a, e = utils.nofloor_elastic(
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
                a = utils.bisection(
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
                    a, e, Lsu = utils.partly_lifted_elastic(
                        d=d, h=h, L=L, w=w, EA=EA, maxit=maxit, tol=tol
                    )
                    Ls[:] = Lsu
                    x0 = a * np.arccosh(1 + h / a)
                    y_offset = -a
                elif Ls1 <= Ls0:  # fully lifted
                    x0 = d
                    Ls[:] = L
                    a, e = utils.fully_lifted_elastic(
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
