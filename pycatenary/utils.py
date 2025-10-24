import warnings
from typing import Callable, Sequence

import numpy as np

# ignore overflow warnings
# cosh can overflow (e.g. inf expected for cosh(d/a) when d/a is very large)
np.seterr(over="ignore")

int1_default = 1e-2
int2_default = 1e10
maxit_default = 1000
tol_default = 1e-6


def get_root_a(
    L: float,
    d: float,
    h: float,
    a0: float = 1.0,
    tol: float = tol_default,
    maxit: int = maxit_default,
    int1: float = int1_default,
    int2: float = int2_default,
) -> float:
    """Returns the initial guess for the catenary a parameter.

    Parameters
    ----------
    L: float
        Unstretched line length [m].
    d: float
        Horizontal distance between anchor and fairlead [m].
    h: float
        Vertical distance between anchor and fairlead [m].
    a0: float
        Initial guess for the catenary a parameter.
    tol: float
        Tolerance for the root finding algorithm.
    maxit: int
        Maximum number of iterations for the root finding algorithm.
    int1: float
        Lower bound for the bisection algorithm.
    int2: float
        Upper bound for the bisection algorithm.

    Returns
    -------
    a: float
        Initial guess for the catenary a parameter.
    """
    g = lambda a: 2.0 * a * np.sinh(d / (2.0 * a)) - np.sqrt(
        L**2.0 - h**2.0
    )
    dg = (
        lambda a: 2.0 * np.sinh(d / (2.0 * a)) - d * np.cosh(d / (2.0 * a)) / a
    )
    a = newton_raphson(f=g, df=dg, x0=a0, tol=tol, maxit=maxit)
    if np.isnan(a) or a < 0:
        a = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
    return a


def newton_raphson(
    f: Callable[[float], float],
    df: Callable[[float], float],
    x0: float,
    tol: float = tol_default,
    maxit: int = maxit_default,
    must_converge: bool = True,
) -> float:
    """Newton-Raphson root finding algorithm (for transcendental equations).

    Parameters
    ----------
    f: function
        Function to find the root of.
    df: function
        Derivative of the function f (df/dx).
    x0: float
        Initial guess of x.
    tol: float
        Tolerance for the root finding algorithm.
    maxit: int
        Maximum number of iterations for the root finding algorithm.
    must_converge: bool
        If True, an error will be raised if the algorithm does not converge.
        If False, a warning will be issued if the algorithm does not converge.

    Returns
    -------
    x: float
        Root of the function.
    """
    x_prev = x0
    x = x0 - f(x0) / df(x0)
    err = np.abs(x - x_prev)
    niter = 0
    while err > tol and niter < maxit:
        niter += 1
        x_prev = x
        x = x - f(x) / df(x)
        err = np.abs(x - x_prev)
    if maxit <= niter:
        if must_converge:
            raise RuntimeError("Newton-Raphson did not converge!")
        else:
            warnings.warn("Newton-Raphson did not converge!")
        x = np.nan
    return x


def bisection(
    f: Callable[[float], float],
    int1: float,
    int2: float,
    tol: float = tol_default,
    maxit: int = maxit_default,
    must_converge: bool = True,
) -> float:
    """Bisection root finding algorithm (for transcendental equations).

    Parameters
    ----------
    f: function
        Function to find the root of.
    int1: float
        Lower bound for the bisection algorithm.
    int2: float
        Upper bound for the bisection algorithm.
    tol: float
        Tolerance for the root finding algorithm.
    maxit: int
        Maximum number of iterations for the root finding algorithm.

    Returns
    -------
    x: float
        Root of the function.
    """
    err = np.abs(int2 - int1) / 2.0
    niter = 0
    while err > tol and niter < maxit:
        niter += 1
        x = (int1 + int2) / 2.0
        if np.sign(f(x)) == np.sign(f(int1)):
            int1 = x
        else:
            int2 = x
        err = np.abs(int2 - int1) / 2.0
    if maxit <= niter:
        if must_converge:
            raise RuntimeError("Bisection did not converge!")
        else:
            warnings.warn("Bisection did not converge!")
        x = np.nan
    return x


def integrate_tension(
    s1: float, s2: float, w: float, Ha: float, Va: float
) -> float:
    """Helper function to integrate tension along a cable segment.

    Parameters
    ----------
    s1: float
        Start of the cable segment [m].
    s2: float
        End of the cable segment [m].
    w: float
        Submerged weight [N/m].
    Ha: float
        Horizontal tension at s1 [N].
    Va: float
        Vertical tension at s1 [N].

    Returns
    -------
    Ts: float
        Integrated tension.
    """
    Ts1 = (
        1
        / (2 * w)
        * (
            (Va + w * s1) * np.sqrt(Ha**2 + (Va + w * s1) ** 2)
            + Ha**2 * np.arcsinh((Va + w * s1) / Ha)
        )
    )
    Ts2 = (
        1
        / (2 * w)
        * (
            (Va + w * s2) * np.sqrt(Ha**2 + (Va + w * s2) ** 2)
            + Ha**2 * np.arcsinh((Va + w * s2) / Ha)
        )
    )
    return Ts2 - Ts1


def nofloor_rigid(
    d: float,
    h: float,
    L: Sequence[float],
    tol: float = tol_default,
    maxit: int = maxit_default,
    int1: float = int1_default,
    int2: float = int2_default,
    must_converge: bool = True,
) -> float:
    """Returns catenary shape for rigid cable with no floor.

    Parameters
    ----------
    d: float
        Horizontal distance between anchor and fairlead [m].
    h: float
        Vertical distance between anchor and fairlead [m].
    L: Sequence[float]
        Unstretched line length [m].
    tol: float
        Tolerance for the root finding algorithm.
    maxit: int
        Maximum number of iterations for the root finding algorithm.
    int1: float
        Lower bound for the bisection algorithm.
    int2: float
        Upper bound for the bisection algorithm.
    must_converge: bool
        If True, an error will be raised if the algorithm does not converge.
        If False, a warning will be issued if the algorithm does not converge.

    Returns
    -------
    a: float
        Catenary shape parameter.
    """
    Lt = np.sum(L)
    g = lambda a: 2 * a * np.sinh(d / (2 * a)) - np.sqrt(Lt**2 - h**2)
    a0 = bisection(
        f=g,
        int1=int1,
        int2=int2,
        tol=tol,
        maxit=maxit,
        must_converge=must_converge,
    )
    return a0


def nofloor_elastic(
    d: float,
    h: float,
    L: Sequence[float],
    w: Sequence[float],
    EA: Sequence[float],
    tol=tol_default,
    maxit: int = maxit_default,
    int1: float = int1_default,
    int2: float = int2_default,
    must_converge: bool = True,
) -> tuple[float, np.ndarray]:
    """Returns catenary solution for elastic cable with no floor.

    Parameters
    ----------
    d: float
        Horizontal distance between anchor and fairlead [m].
    h: float
        Vertical distance between anchor and fairlead [m].
    L: Sequence[float]
        Unstretched line length [m].
    w: Sequence[float]
        Submerged weight [N/m].
    EA: Sequence[float]
        Axial stiffness [N].
    tol: float
        Tolerance for the root finding algorithm.
    maxit: int
        Maximum number of iterations for the root finding algorithm.
    int1: float
        Lower bound for the bisection algorithm.
    int2: float
        Upper bound for the bisection algorithm.
    must_converge: bool
        If True, an error will be raised if the algorithm does not converge.
        If False, a warning will be issued if the algorithm does not converge.

    Returns
    -------
    a: float
        Catenary shape parameter.
    e: np.ndarray
        Elongation of the cable segments [m].
    """
    Lt = np.sum(L)  # total length of cable
    w_av = np.sum(w * L / Lt)  # average weight of cable
    e = np.zeros(len(L))  # stretching of cable segments
    diff = tol + 1
    niter = 0
    while diff > tol and niter < maxit:
        niter += 1
        Lte = np.sum(L + e)
        g = lambda a: 2 * a * np.sinh(d / (2 * a)) - np.sqrt(Lte**2 - h**2)
        a = bisection(
            f=g,
            int1=int1,
            int2=int2,
            tol=tol,
            maxit=maxit,
            must_converge=must_converge,
        )
        # HACK: not real elongation if multi-segmented line here
        T = np.sqrt((a * w_av) ** 2 + np.sum(w * L))
        e = T * L / EA
        et = np.sum(e)
        Lte_check = Lt + et  # store new Ls value as calculated with stretching
        diff = np.abs(Lte - Lte_check)
    return a, e


def fully_lifted_elastic(
    d: float,
    h: float,
    L: Sequence[float],
    w: Sequence[float],
    EA: Sequence[float],
    tol: float = tol_default,
    maxit: int = maxit_default,
    int1: float = int1_default,
    int2: float = int2_default,
    must_converge: bool = True,
) -> tuple[float, np.ndarray]:
    """Returns catenary solution for fully lifted elastic cable.

    Parameters
    ----------
    d: float
        Horizontal distance between anchor and fairlead [m].
    h: float
        Vertical distance between anchor and fairlead [m].
    L: Sequence[float]
        Unstretched line length [m].
    w: Sequence[float]
        Submerged weight [N/m].
    EA: Sequence[float]
        Axial stiffness [N].
    tol: float
        Tolerance for the root finding algorithm.
    maxit: int
        Maximum number of iterations for the root finding algorithm.
    int1: float
        Lower bound for the bisection algorithm.
    int2: float
        Upper bound for the bisection algorithm.
    must_converge: bool
        If True, an error will be raised if the algorithm does not converge.
        If False, a warning will be issued if the algorithm does not converge.

    Returns
    -------
    a: float
        Catenary shape parameter.
    e: np.ndarray
        Elongation of the cable segments [m].
    """
    Ls_tot = Le_tot = 0
    Lt = np.sum(L)  # total length of cable
    w_av = np.sum(w * L / Lt)  # average weight of cable
    e = np.zeros(len(L))  # stretching of cable segments

    t_high = h / d
    t_low = 0.0

    diff = 1.0
    niter = 0
    a = 1.0
    while diff > tol and niter < maxit:
        niter += 1
        t = (t_low + t_high) / 2.0
        angle = np.arctan(t)
        # transcendental equation
        g = (
            lambda a: a
            * (np.cosh(d / a + np.arcsinh(t)) - np.cosh(np.arcsinh(t)))
            - h
        )
        a = bisection(
            f=g,
            int1=int1,
            int2=int2,
            tol=tol,
            maxit=maxit,
            must_converge=must_converge,
        )
        # dg = lambda a: np.cosh(d / a + np.arcsinh(t)) - d / a * np.sinh(
        #     d / a + np.arcsinh(t)
        # )
        # a = newton_raphson(f=g, df=dg, x0=a, tol=tol, maxit=maxit)
        # if a is np.nan:
        #    a = bisection(f=g, int1=1., int2=100000, tol=tol, maxit=maxit)
        # get new total Ls from solution a
        Ls_tot = np.sqrt((2 * a * np.sinh(d / (2 * a))) ** 2 + h**2)
        # get new stretching from solution a
        Ta = a * w_av / np.cos(angle) * Lt / Ls_tot
        Ha = Ta * np.cos(angle)
        Va = Ta * np.sin(angle)

        # compute elongation
        for i in range(len(e)):
            # integrate tension
            T_int = integrate_tension(
                0, L[i], w[i], Ha, Va + np.sum(w[:i] * L[:i])
            )
            # get elongation
            e[i] = T_int / EA[i]

        et = np.sum(e)
        Le_tot = Lt + et
        diff = np.abs(Le_tot - Ls_tot)
        if Le_tot > Ls_tot:
            t_high = t
        elif Le_tot < Ls_tot:
            t_low = t
    return a, e


def fully_lifted_rigid(
    d: float,
    h: float,
    L: Sequence[float],
    tol=tol_default,
    maxit: int = maxit_default,
    int1: float = int1_default,
    int2: float = int2_default,
    must_converge: bool = True,
) -> float:
    """Returns catenary solution for fully lifted rigid cable.

    Parameters
    ----------
    d: float
        Horizontal distance between anchor and fairlead [m].
    h: float
        Vertical distance between anchor and fairlead [m].
    L: Sequence[float]
        Unstretched line length [m].
    tol: float
        Tolerance for the root finding algorithm.
    maxit: int
        Maximum number of iterations for the root finding algorithm.
    int1: float
        Lower bound for the bisection algorithm.
    int2: float
        Upper bound for the bisection algorithm.
    must_converge: bool
        If True, an error will be raised if the algorithm does not converge.
        If False, a warning will be issued if the algorithm does not converge.

    Returns
    -------
    a: float
        Catenary shape parameter.
    """
    g = lambda a: 2.0 * a * np.sinh(d / (2.0 * a)) - np.sqrt(
        np.sum(L) ** 2.0 - h**2.0
    )
    dg = (
        lambda a: 2.0 * np.sinh(d / (2.0 * a)) - d * np.cosh(d / (2.0 * a)) / a
    )
    a0 = bisection(
        f=g,
        int1=int1,
        int2=int2,
        tol=tol,
        maxit=maxit,
        must_converge=must_converge,
    )
    a1 = newton_raphson(
        f=g, df=dg, x0=a0, tol=tol, maxit=maxit, must_converge=False
    )
    if np.isnan(a1) or a1 < 0:
        a = a0
    else:
        a = a1
    return a


def partly_lifted_elastic(
    d: float,
    h: float,
    L: Sequence[float],
    w: Sequence[float],
    EA: Sequence[float],
    tol: float = tol_default,
    maxit: int = maxit_default,
    int1: float = int1_default,
    int2: float = int2_default,
    must_converge: bool = True,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Returns catenary solution for partly lifted elastic cable.

    Parameters
    ----------
    d: float
        Horizontal distance between anchor and fairlead [m].
    h: float
        Vertical distance between anchor and fairlead [m].
    L: Sequence[float]
        Unstretched line length [m].
    w: Sequence[float]
        Submerged weight [N/m].
    EA: Sequence[float]
        Axial stiffness [N].
    tol: float
        Tolerance for the root finding algorithm.
    maxit: int
        Maximum number of iterations for the root finding algorithm.
    int1: float
        Lower bound for the bisection algorithm.
    int2: float
        Upper bound for the bisection algorithm.
    must_converge: bool
        If True, an error will be raised if the algorithm does not converge.
        If False, a warning will be issued if the algorithm does not converge.

    Returns
    -------
    a: float
        Catenary shape parameter.
    e: np.ndarray
        Elongation of the cable segments [m].
    Lsu: np.ndarray
        Lifted line lengths of the cable segments [m].
    """
    diff = 1.0
    niter = 0
    a = 1.0
    e = np.zeros(len(L))
    Ls = np.zeros(len(L))
    Lt = np.sum(L)
    x0_high = d
    x0_low = 0
    et = 0
    Ls = 0
    Lsu = np.zeros(len(L))
    while np.abs(diff) > tol and niter < maxit:
        niter += 1
        x0 = (x0_low + x0_high) / 2.0
        g = lambda a: a * (np.cosh(x0 / a) - 1.0) - h
        a = bisection(
            f=g,
            int1=int1,
            int2=int2,
            tol=tol,
            maxit=maxit,
            must_converge=must_converge,
        )
        # dg = lambda a: np.cosh(d / a + np.arcsinh(t)) - d / a * np.sinh(
        #     d / a + np.arcsinh(t)
        # )
        Ls = h * np.sqrt(1 + 2 * a / h)
        Lns_tot_check = 0
        ground = d - x0
        lifted = False
        for i in range(len(L)):
            if lifted is False:
                Lsu[i] = 0
                Lns_tot_check += L[i]
                if Lns_tot_check > ground:
                    Lsu[i] = Lns_tot_check - ground
                    lifted = True
            else:
                Lsu[i] = L[i]
        w_av = np.sum(w * Lsu) / Ls
        H = a * w_av
        # compute elongation
        for i in range(len(L)):
            # integrate tension
            T_int = integrate_tension(
                0, Lsu[i], w[i], H, np.sum(w[:i] * Lsu[:i])
            )
            # get elongation
            e[i] = T_int / EA[i]
        et = np.sum(e)
        X0 = Lt + et - Ls
        diff = X0 + x0 - d
        if diff > 0:
            x0_high = x0
        elif diff < 0:
            x0_low = x0
    return a, e, Lsu


def partly_lifted_rigid(
    d: float,
    h: float,
    L: Sequence[float],
    tol: float = tol_default,
    maxit: int = maxit_default,
    int1: float = int1_default,
    int2: float = int2_default,
    must_converge: bool = True,
) -> tuple[float, np.ndarray]:
    """Returns catenary solution for partly lifted rigid cable.

    Parameters
    ----------
    d: float
        Horizontal distance between anchor and fairlead [m].
    h: float
        Vertical distance between anchor and fairlead [m].
    L: Sequence[float]
        Unstretched line length [m].
    tol: float
        Tolerance for the root finding algorithm.
    maxit: int
        Maximum number of iterations for the root finding algorithm.
    int1: float
        Lower bound for the bisection algorithm.
    int2: float
        Upper bound for the bisection algorithm.
    must_converge: bool
        If True, an error will be raised if the algorithm does not converge.
        If False, a warning will be issued if the algorithm does not converge.

    Returns
    -------
    a: float
        Catenary shape parameter.
    Lsu: np.ndarray
        Lifted line lengths of the cable segments [m].
    """
    diff = 1.0
    niter = 0
    a = 1.0
    Lsu = np.zeros(len(L))
    Lt = np.sum(L)
    x0_high = d
    x0_low = 0
    while np.abs(diff) > tol and niter < maxit:
        niter += 1
        x0 = (x0_low + x0_high) / 2.0
        g = lambda a: a * (np.cosh(x0 / a) - 1.0) - h
        # dg = lambda a: np.cosh(d / a + np.arcsinh(t)) - d / a * np.sinh(
        #     d / a + np.arcsinh(t)
        # )
        a = bisection(
            f=g,
            int1=int1,
            int2=int2,
            tol=tol,
            maxit=maxit,
            must_converge=must_converge,
        )
        Ls = h * np.sqrt(1 + 2 * a / h)
        Lns_tot_check = 0
        ground = d - x0
        lifted = False
        for i in range(len(L)):
            if lifted is False:
                Lsu[i] = 0
                Lns_tot_check += L[i]
                if Lns_tot_check > ground:
                    Lsu[i] = Lns_tot_check - ground
                    lifted = True
            else:
                Lsu[i] = L[i]
        x0 = a * np.arccosh(1 + h / a)
        X0 = Lt - Ls
        diff = X0 + x0 - d
        if diff > 0:
            x0_high = x0
        elif diff < 0:
            x0_low = x0
    return a, Lsu


def straight_elastic(
    d: float,
    h: float,
    L: Sequence[float],
    w: Sequence[float],
    EA: Sequence[float],
    H_low: float = 0,
    H_high: float = 1e10,
    tol: float = tol_default,
    maxit: int = maxit_default,
    must_converge: bool = True,
) -> tuple[float, np.ndarray]:
    """Returns horizontal tension and elongation for straight elastic cable.

    Parameters
    ----------
    d: float
        Horizontal distance between anchor and fairlead [m].
    h: float
        Vertical distance between anchor and fairlead [m].
    L: Sequence[float]
        Unstretched line length [m].
    w: Sequence[float]
        Submerged weight [N/m].
    EA: Sequence[float]
        Axial stiffness [N].
    H_low: float
        Lower bound for the horizontal tension at anchor [N].
    H_high: float
        Upper bound for the horizontal tension at anchor [N].
    tol: float
        Tolerance for the root finding algorithm.
    maxit: int
        Maximum number of iterations for the root finding algorithm.
    must_converge: bool
        If True, an error will be raised if the algorithm does not converge.
        If False, a warning will be issued if the algorithm does not converge.

    Returns
    -------
    H: float
        Horizontal tension at anchor [N].
    e: np.ndarray
        Elongation of the cable segments [m].
    """
    Lt = np.sum(L)  # total length of cable
    assert Lt <= np.sqrt(d**2 + h**2)
    e = np.zeros(len(L))  # stretching of cable segments
    et = (d**2 + h**2) ** 0.5 - Lt  # stretching to reach minimum length
    angle = np.arctan(h / d)  # angle

    et_check = et
    H = (H_low + H_high) / 2.0  # guess of horizontal tension at anchor
    diff = 1.0
    niter = 0
    while diff > tol and niter < maxit:
        niter += 1
        if et_check > et:
            H_high = H
        elif et_check < et:
            H_low = H
        H = (H_low + H_high) / 2.0
        Va = H * np.tan(angle)  # tension at anchor
        for i in range(len(e)):
            e[i] = (
                (
                    (H**2 + (Va + np.sum(w[i:] * L[i:]) + w[i] * L[i]) ** 2)
                    ** 0.5
                )
                * L[i]
                / EA[i]
            )
        et_check = np.sum(e)
        diff = np.abs(et_check - et)
    return H, e
