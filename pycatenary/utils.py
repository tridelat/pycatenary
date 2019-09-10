import numpy as np

def get_root_a(L, d, h, a0=1., tol=1e-6, maxit=1000, int1=1e-6, int2=1e6):
    g = lambda a: 2.*a*np.sinh(d/(2.*a))-np.sqrt((L)**2.-h**2.)
    dg = lambda a: 2.*np.sinh(d/(2.*a))-d*np.cosh(d/(2.*a))/a
    a = newton_raphson(f=g, df=dg, x0=a0, tol=tol, maxit=maxit)
    if np.isnan(a) or a < 0:
        a = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
    return a

def newton_raphson(f, df, x0, tol=1e-6, maxit=1000):
    """
    Root finding algorithm (for transcendental equations)
    :param f: must be a function (so f(x) = 0 returns the required x)
    :param df: derivative of the function f (df/dx)
    :param x0: initial guess of x
    :param tol: tolerance
    :param maxit: maximum number of iterations
    :return: x
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

def bisection(f, int1, int2, tol=1e-6, maxit=1000):
    """
    Root finding algorithm (for transcendental equations).
    A solution is found between a user-provided interval.
    :param f: must be a function (so f(x) = 0 returns the required x)
    :param int1: lower end value
    :param int2: higher end value
    :param tol: tolerance
    :param maxit: maximum number of iterations
    :return: x (solution)
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

def nofloor_rigid(d, h, L, tol=1e-6, maxit=1000):
    g = lambda a: 2.*a*np.arcsinh((L/2.)/a)-d
    a0 = bisection(f=g, int1=0.01, int2=1e15, tol=tol, maxit=maxit)
    return a0

def nofloor_elastic(d, h, L, w, EA, tol=1e-6, maxit=1000, int1=1e-6, int2=1e6):
    e = w*L/EA
    Le = L+e
    g = lambda a: 2.*a*np.arcsinh((Le/2.)/a)-d
    a0 = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
    return a0, e

def fully_lifted_elastic(d, h, L, w, EA, tol=1e-6, maxit=1000, int1=1e-6, int2=1e6):
    Ls_tot = Le = 0
    Lt = np.sum(L)  # total length of cable
    w_av = np.sum(w*L/Lt) # average weight of cable
    e = np.zeros(len(L))  # stretching of cable segments

    t_high = h/d
    t_low = 0

    diff = 1.
    niter = 0
    a = 1.
    while diff > tol and niter < maxit:
        niter += 1
        if Le > Ls_tot:
            t_high = t
        elif Le < Ls_tot:
            t_low = t
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
            e[i] = np.sqrt(Ha**2+(Va+(np.sum(w[:i]*L[:i])+w[i]*L[i]/2.))**2)*L[i]/EA[i]
        et = np.sum(e)
        Le = Lt+et  # store new Ls value as calculated with stretching
        x0 = d
        diff = np.abs(Le-Ls_tot)
    return a, e

def fully_lifted_rigid(d, h, L, tol=1e-6, maxit=1000, int1=1e-6, int2=1e6):
    g = lambda a: 2.*a*np.sinh(d/(2.*a))-np.sqrt(np.sum(L)**2.-h**2.)
    dg = lambda a: 2.*np.sinh(d/(2.*a))-d*np.cosh(d/(2.*a))/a
    a0 = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
    a1 = newton_raphson(f=g, df=dg, x0=a0, tol=tol, maxit=maxit)
    if np.isnan(a1) or a1 < 0:
        a = a0
    else:
        a = a1
    return a

def partly_lifted_elastic(d, h, L, w, EA, tol=1e-6, maxit=1000, int1=1e-6, int2=1e6):
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
    while diff > tol and niter < maxit:
        niter += 1
        if Ls < Lse:
            x0_high = x0
        elif Ls > Lse:
            x0_low = x0
        x0 = (x0_low+x0_high)/2.
        lifted = False
        g = lambda a: a*(np.cosh(x0/a)-1.)-h
        dg = lambda a: np.cosh(d/a+np.arcsinh(t))-d/a*np.sinh(d/a+np.arcsinh(t))
        a = bisection(f=g, int1=int1, int2=int2, tol=tol, maxit=maxit)
        Ls = h*np.sqrt(1+2*a/h)
        Lns_tot_check = 0
        ground = d-x0
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
            e[i] = np.sqrt(H**2+(np.sum(w[i:]*Lsu[i:])+w[i]*Lsu[i]/2.)**2)*Lsu[i]/EA[i]
        et = np.sum(e)
        Lse = Lsu_tot+et
        diff = np.abs(Ls-Lse)
    return a, e, Lsu

def partly_lifted_rigid(d, h, L, tol=1e-6, maxit=1000):
    diff = 1.
    niter = 0
    a = 1.
    Ls = np.zeros(len(L))
    Lt = np.sum(L)
    while diff > tol and niter < maxit:
        niter += 1
        Ls_tot = h*np.sqrt(1+2*a/h)  # lifted line length
        x0 = a*np.arccosh(1+h/a)
        X0 = (Lt-Ls_tot)+x0
        if Ls_tot > Lt:
            diff = 0  # cannot find solution; line must be fully lifted
            a = Ls = np.nan
        else:
            diff = np.abs(X0-d)
        if diff > tol:
            a = a*((d/X0)**(d/x0))
    return a, Ls

def straight_elastic(d, h, L, w, EA, H_low=0, H_high=1e10, tol=1e-6, maxit=1000):

    Lt = np.sum(L)  # total length of cable
    w_av = np.sum(w*L/Lt) # average weight of cable
    e = np.zeros(len(L))  # stretching of cable segments
    et = (d**2+h**2)**0.5-Lt  # stretching to reach minimum length
    angle = np.arctan(h/d)  # angle

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
            e[i] = ((H**2+(Va+w[i:]*L[w:]+w[i]*L[i]/2.)**2)**0.5)*L[i]/EA[i]
        et_check = np.sum(e)
    return H, e
