import numpy as np

# ignore overflow warnings
# cosh can overflow (e.g. inf expected for cosh(d/a) when d/a is very large)
np.seterr(over='ignore')

int1_default = 1e-2
int2_default = 1e10
maxit_default = 1000
tol_default = 1e-6

def get_root_a(L, d, h, a0=1., tol=tol_default, maxit=maxit_default, int1=int1_default, int2=int2_default):
    g = lambda a: 2.*a*np.sinh(d/(2.*a))-np.sqrt((L)**2.-h**2.)
    dg = lambda a: 2.*np.sinh(d/(2.*a))-d*np.cosh(d/(2.*a))/a
    a = newton_raphson(f=g, df=dg, x0=a0, tol=tol, maxit=maxit)
    if np.isnan(a) or a < 0:
        a = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
    return a

def newton_raphson(f, df, x0, tol=tol_default, maxit=maxit_default):
    """Root finding algorithm (for transcendental equations)

    Parameters
    ----------
    f: function
        must be a function (so f(x) = 0 returns the required x)
    df: function
        derivative of the function f (df/dx)
    x0: double
        initial guess of x
    tol: double
        tolerance
    maxit: int
        maximum number of iterations

    Returns
    -------
    x: double
        root
    """
    x_prev = x0
    x = x0-f(x0)/df(x0)
    err = np.abs(x-x_prev)
    niter = 0
    while err > tol and niter < maxit:
        niter += 1
        x_prev = x
        x = x-f(x)/df(x)
        err = np.abs(x-x_prev)
    #print('Newton-Raphson: iterations', niter, ', solution', x, ', err', err)
    if maxit <= niter:
        print('did not converge!')
        x = np.nan
    return x

def bisection(f, int1, int2, tol=tol_default, maxit=maxit_default):
    """Root finding algorithm (for transcendental equations)

    Parameters
    ----------
    f: function
        must be a function (so f(x) = 0 returns the required x)
    int1: double
        lower end value
    int2: double
        lower end value
    tol: double
        tolerance
    maxit: int
        maximum number of iterations

    Returns
    -------
    x: double
        root
    """
    err = np.abs(int2-int1)/2.
    niter = 0
    while err > tol and niter < maxit:
        niter += 1
        x = (int1+int2)/2.
        if np.sign(f(x)) == np.sign(f(int1)):
            int1 = x
        else:
            int2 = x
        err = np.abs(int2-int1)/2.
    #print('Bisection: iterations', niter, ', solution', x, ', err', err)
    if maxit <= niter:
        print('did not converge!')
        x = np.nan
    return x

def nofloor_rigid(d, h, L, tol=tol_default, maxit=maxit_default, int1=int1_default, int2=int2_default):
    Lt = np.sum(L)
    g = lambda a: 2*a*np.sinh(d/(2*a))-np.sqrt(Lt**2-h**2)
    a0 = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
    return a0

def nofloor_elastic(d, h, L, w, EA, tol=tol_default, maxit=maxit_default, int1=int1_default, int2=int2_default):
    Lt = np.sum(L)  # total length of cable
    w_av = np.sum(w*L/Lt) # average weight of cable
    e = np.zeros(len(L))  # stretching of cable segments
    diff = tol+1
    niter = 0
    while diff > tol and niter < maxit:
        niter += 1
        Lte = np.sum(L+e)
        g = lambda a: 2*a*np.sinh(d/(2*a))-np.sqrt(Lte**2-h**2)
        a = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
        #
        T = np.sqrt((a*w_av)**2+np.sum(w*L))
        et = T*L/EA
        Lte_check = Lt+et  # store new Ls value as calculated with stretching
        diff = np.abs(Lte-Lte_check)
    # HACK: not real elongation if multi-segmented line here
    e[:] = et*L/Lt
    return a, e

def fully_lifted_elastic(d, h, L, w, EA, tol=tol_default, maxit=maxit_default, int1=int1_default, int2=int2_default):
    Ls_tot = Le = 0
    Lt = np.sum(L)  # total length of cable
    w_av = np.sum(w*L/Lt) # average weight of cable
    e = np.zeros(len(L))  # stretching of cable segments

    t_high = h/d
    t_low = 0.

    diff = 1.
    niter = 0
    a = 1.
    while diff > tol and niter < maxit:
        niter += 1
        t = (t_low+t_high)/2.
        angle = np.arctan(t)
        # transcendental equation
        g = lambda a: a*(np.cosh(d/a+np.arcsinh(t))-np.cosh(np.arcsinh(t)))-h
        dg = lambda a: np.cosh(d/a+np.arcsinh(t))-d/a*np.sinh(d/a+np.arcsinh(t))
        a = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
        #a = newton_raphson(f=g, df=dg, x0=a, tol=tol, maxit=maxit)
        #if a is np.nan:
        #    a = bisection(f=g, int1=1., int2=100000, tol=tol, maxit=maxit)
        # get new total Ls from solution a
        Ls_tot = np.sqrt((2*a*np.sinh(d/(2*a)))**2+h**2)
        # get new stretching from solution a
        Ta = a*w_av/np.cos(angle)*Lt/Ls_tot
        Ha = Ta*np.cos(angle)
        Va = Ta*np.sin(angle)
        for i in range(len(e)):
            e[i] = np.sqrt(Ha**2+(Va+(np.sum(w[:i]*L[:i])+w[i]*L[i]))**2)*L[i]/EA[i]
        et = np.sum(e)
        Le = Lt+et  # store new Ls value as calculated with stretching
        x0 = d
        diff = np.abs(Le-Ls_tot)
        if Le > Ls_tot:
            t_high = t
        elif Le < Ls_tot:
            t_low = t
    return a, e

def fully_lifted_rigid(d, h, L, tol=tol_default, maxit=maxit_default, int1=int1_default, int2=int2_default):
    g = lambda a: 2.*a*np.sinh(d/(2.*a))-np.sqrt(np.sum(L)**2.-h**2.)
    dg = lambda a: 2.*np.sinh(d/(2.*a))-d*np.cosh(d/(2.*a))/a
    a0 = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
    a1 = newton_raphson(f=g, df=dg, x0=a0, tol=tol, maxit=maxit)
    if np.isnan(a1) or a1 < 0:
        a = a0
    else:
        a = a1
    return a

def partly_lifted_elastic(d, h, L, w, EA, tol=tol_default, maxit=maxit_default, int1=int1_default, int2=int2_default):
    diff = 1.
    niter = 0
    a = 1.
    e = np.zeros(len(L))
    Ls = np.zeros(len(L))
    Lt = np.sum(L)
    x0_high = d
    x0_low = 0
    Ls_tot = 0
    Lsu_tot = 0
    et = 0
    Lse = 0
    Ls = 0
    Lsu = np.zeros(len(L))
    while np.abs(diff) > tol and niter < maxit:
        niter += 1
        x0 = (x0_low+x0_high)/2.
        g = lambda a: a*(np.cosh(x0/a)-1.)-h
        dg = lambda a: np.cosh(d/a+np.arcsinh(t))-d/a*np.sinh(d/a+np.arcsinh(t))
        a = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
        Ls = h*np.sqrt(1+2*a/h)
        Lns_tot_check = 0
        ground = d-x0
        lifted = False
        for i in (range(len(L))):
            if lifted is False:
                Lsu[i] = 0
                Lns_tot_check += L[i]
                if Lns_tot_check > ground:
                    Lsu[i] = Lns_tot_check-ground
                    lifted = True
            else:
                Lsu[i] = L[i]
        Lsu_tot = np.sum(Lsu)
        w_av = np.sum(w[i]*Lsu[i])/Ls
        H = a*w_av
        for i in range(len(L)):
            e[i] = np.sqrt(H**2+(np.sum(w[i:]*Lsu[i:])+w[i]*Lsu[i])**2)*Lsu[i]/EA[i]
        et = np.sum(e)
        Lse = Lsu_tot+et
        X0 = Lt+et-Ls
        diff = X0+x0-d
        if diff > 0:
            x0_high = x0
        elif diff < 0:
            x0_low = x0
    return a, e, Lsu

def partly_lifted_rigid(d, h, L, tol=tol_default, maxit=maxit_default, int1=int1_default, int2=int2_default):
    diff = 1.
    niter = 0
    a = 1.
    Lsu = np.zeros(len(L))
    Lt = np.sum(L)
    x0_high = d
    x0_low = 0
    while np.abs(diff) > tol and niter < maxit:
        niter += 1
        x0 = (x0_low+x0_high)/2.
        g = lambda a: a*(np.cosh(x0/a)-1.)-h
        dg = lambda a: np.cosh(d/a+np.arcsinh(t))-d/a*np.sinh(d/a+np.arcsinh(t))
        a = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
        Ls = h*np.sqrt(1+2*a/h)
        Lns_tot_check = 0
        ground = d-x0
        lifted = False
        for i in (range(len(L))):
            if lifted is False:
                Lsu[i] = 0
                Lns_tot_check += L[i]
                if Lns_tot_check > ground:
                    Lsu[i] = Lns_tot_check-ground
                    lifted = True
            else:
                Lsu[i] = L[i]
        x0 = a*np.arccosh(1+h/a)
        X0 = (Lt-Ls)
        diff = X0+x0-d
        if diff > 0:
            x0_high = x0
        elif diff < 0:
            x0_low = x0
    return a, Lsu

def straight_elastic(d, h, L, w, EA, H_low=0, H_high=1e10, tol=tol_default, maxit=maxit_default):

    Lt = np.sum(L)  # total length of cable
    assert Lt <= np.sqrt(d**2+h**2)
    w_av = np.sum(w*L/Lt) # average weight of cable
    e = np.zeros(len(L))  # stretching of cable segments
    et = (d**2+h**2)**0.5-Lt  # stretching to reach minimum length
    angle = np.arctan(h/d)  # angle

    et_check = et
    H = (H_low+H_high)/2.  # guess of horizontal tension at anchor
    diff = 1.
    niter = 0
    while diff > tol and niter < maxit:
        niter += 1
        if et_check > et:
            H_high = H
        elif et_check < et:
            H_low = H
        H = (H_low+H_high)/2.
        Va = H*np.tan(angle)  # tension at anchor
        for i in range(len(e)):
            e[i] = ((H**2+(Va+np.sum(w[i:]*L[i:])+w[i]*L[i])**2)**0.5)*L[i]/EA[i]
        et_check = np.sum(e)
        diff = np.abs(et_check-et)
    return H, e
