from . import utils
import numpy as np
import sys

def get_array(x):
    if np.isscalar(x):
        x = np.array([x])
    else:
        x = np.asarray(x)
    return x

class CatenaryBase(object):
    """Base class for catenaries

    Parameters
    ----------
    line: pycatenary.cable.MooringLine
        line holding properties necessary for calculation of catenary
    """
    def __init__(self, line):
        self.line = line
        # horizontal distance
        self.d = 0.
        # vertical distance
        self.h = 0.
        # catenary a
        self.a = 0.
        # elongation
        self.e = 0.
        # horizontal span
        self.x0 = 0.
        # lifted line length
        self.Ls = 0.
        # submerged weight
        self.maxit = 1000
        # tolerance
        self.tol = 1e-10
        # first guess for bisection (int1)
        self.bisection_int1 = 1e-6
        # first guess for bisection (int2)
        self.bisection_int2 = 1e6
        # offset for x
        self._x_offet = 0.
        # offset for y
        self._y_offset = 0.
        # offset for s
        self._s_offset = 0.

    def getTension(self, s):
        s0 = self.d-self.x0
        # total line lengths
        Lt = np.sum(self.line.L) # unstretched
        Lst = np.sum(self.Ls) # unstretched
        Lset = Lst+np.sum(self.e) # stretched
        if Lt >= s >= s0:
            # average w
            w_av = np.sum(self.line.w*self.Ls)/Lst
            # horizontal tension
            Th = self.a*w_av*(Lst/Lset)
            # vertical tension
            ds0 = self.ds2xy(0)
            angle0 = np.arctan(ds0[1]/ds0[0])
            # vertical tension at anchor
            Tv_a = Th*np.tan(angle0)
            # vertical tension at point
            Tv = Tv_a+w_av*(s-s0)
            # tension at point
            Ts = np.array([Th, Tv])
        elif 0 <= s < s0:
            Ts = np.array([0., 0.])
        else:
            assert False, 'wrong value for s'
        return Ts

    def s2xy(self, s):
        s0 = self.d-self.x0
        s = s+self._s_offset
        s = s*np.sum(self.Ls+self.e)/np.sum(self.Ls)
        a = self.a
        if s < s0 and self.line.floor:
            x = s
            y = 0.-self._y_offset
        else:
            x = s0+a*np.arcsinh((s-s0)/a)
            y = a*np.cosh((x-s0)/a)
        xy = np.array([x+self._x_offset, y+self._y_offset])
        return xy

    def ds2xy(self, s):
        s0 = self.d-self.x0
        s = s+self._s_offset
        a = self.a
        if s < s0 and self.line.floor:
            x = 1.
            y = 0.
        else:
            x = a/np.sqrt(a**2+s**2)
            y = s/np.sqrt(a**2+s**2)
        xy = np.array([x, y])
        return xy

    def plot(self, npoints=100):
        """Plots catenary in 2D from (0, 0) to (d, h)
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        xys = []
        xx = []
        yy = []
        ss = np.linspace(0., np.sum(self.line.L), npoints)
        for s in ss:
            xy = self.s2xy(s)
            xys.append(xy)
            xx.append(xy[0])
            yy.append(xy[1])
        ax.plot(xx, yy)
        print('anchor: {anchor}'.format(anchor=str(xys[0])))
        print('fairlead: {fairlead}'.format(fairlead=str(xys[-1])))
        fig.show()


class CatenaryRigid(CatenaryBase):
    """A class for rigid catenary

    Parameters
    ----------
    line: pycatenary.cable.MooringLine
        line holding properties necessary for calculation of catenary
    """
    def __init__(self, line):
        super(CatenaryRigid, self).__init__(line)

    def getState(self, d, h, floor=True):
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
        L = self.line.L
        Ls = np.zeros(len(L))
        Lt = np.sum(L)
        x_offset = 0.
        y_offset = 0.
        s_offset = 0.
        a = 1.
        x0 = 0.
        if floor is False:
            a = utils.nofloor_rigid(d=d, h=h, L=L, tol=tol, maxit=maxit)
            x0 = d
            Ls[:] = L
            Lst = np.sum(Ls+self.e)
            xx = 0.5*(a*np.log((Lst+h)/(Lst-h))-d)
            xy = 0.5*(a*np.log((Lst+h)/(Lst-h))+d)
            x_offset = -xx
            y_offset = h-a*np.cosh(xy/a)
            s_offset = a*np.sinh(xx/a)
        else:
            # cable straight to seabed:
            if L >= h+d:
                # no horizontal tension
                Lst = 0.
                for ii in range(len(L)):
                    index = -1-i
                    Lst += L[index]
                    if Ltot < h:
                        Ls[index] = L[index]
                    if Ltot >= h:
                        Ls[index] = Ltot-L[index+1:]
                a = 0.
                x0 = 0.
            else:
                # check if line is partly or fully lifted
                f = lambda a: a*(np.cosh(d/a)-1)-h
                a = utils.bisection(f, int1=self.bisection_int1, int2=self.bisection_int2, tol=tol, maxit=maxit)
                Ls0 = a*np.sinh(d/a)  # maximum line length to be fully lifted
                # get actual line length assuming it is fully lifted (from parameter a)
                Ls1 = np.sum(L)
                if Ls1 > Ls0:  # partly lifted
                    a, Ls = utils.partly_lifted_rigid(d=d, h=h, L=L, maxit=maxit, tol=tol)
                    x0 = a*np.arccosh(1+h/a)
                    y_offset = -a
                elif Ls1 <= Ls0:  # fully lifted
                    a = utils.fully_lifted_rigid(d=d, h=h, L=L, maxit=maxit, tol=tol)
                    Ls[:] = L
                    Lst = Lt
                    x0 = d
                    xx = 0.5*(a*np.log((Lst+h)/(Lst-h))-d)
                    xy = 0.5*(a*np.log((Lst+h)/(Lst-h))+d)
                    x_offset = -xx
                    y_offset = h-a*np.cosh(xy/a)
                    s_offset = a*np.sinh(xx/a)
                    if a is np.nan:
                        raise RuntimeError('The line is too stretched/straight, cannot find catenary shape')
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
    line: pycatenary.cable.MooringLine
        line holding properties necessary for calculation of catenary
    """
    def __init__(self, line):
        super(CatenaryElastic, self).__init__(line)

    def getState(self, d, h, floor=True):
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
        L = self.line.L
        minL = np.sqrt(h**2+d**2)+tol  # minimum line length

        L = get_array(self.line.L)  # unstretched line length
        w = get_array(self.line.w)  # submerged weight
        EA = get_array(self.line.EA)  # axial stiffness

        Lt = np.sum(L)  # total unstretched line length
        Ls = np.zeros(len(L))  # lifted line length
        Lsu = np.zeros(len(L))  # unstretched lifted line length
        e = np.zeros(len(L))  # stretching
        et = np.sum(e)  # total stretching
        Le = Lt+et  # stretched line length
        x_offset = 0.
        y_offset = 0.
        s_offset = 0.

        diff = tol+1
        niter = 0

        if floor is False:
            a, e = utils.nofloor_elastic(d=d, h=h, L=L, w=w, EA=EA, tol=tol, maxit=maxit)
            # a, e = utils.fully_lifted_elastic(d=d, h=h, L=L, w=w, EA=EA, int1=a, maxit=maxit, tol=tol)
            x0 = d
            Ls[:] = L
            Lst = np.sum(Ls+e)
            xx = 0.5*(a*np.log((Lst+h)/(Lst-h))-d)
            xy = 0.5*(a*np.log((Lst+h)/(Lst-h))+d)
            x_offset = -xx
            y_offset = h-a*np.cosh(xy/a)
            s_offset = a*np.sinh(xx/a)
        else:
            # cable straight to seabed:
            # find tension and stretching
            for i in reversed(range(len(L))):
                Lsu[i] = L[i]
                for j in range(i, len(L)):
                    e[i] = (w[i]*Lsu[i]/2.+np.sum(w[:i]*Lsu[:i]))*Lsu[i]/EA[i]
                    Ls[i] = Lsu[i]+e[i]
                if np.sum(Ls) >= h:
                    Lhi_low = 0
                    Lhi_high = L[i]
                    while diff > tol:
                        Lsu[i] = (Lhi_low+Lhi_high)/2.
                        for j in range(i, len(L)):
                            e[i] = (w[i]*Lsu[i]/2.+np.sum(w[:i]*Lsu[:i]))*Lsu[i]/EA[i]
                            Ls[i] = Lsu[i]+e[i]
                        if np.sum(Ls) > h:
                            Lhi_high = Lsu[i]
                        elif np.sum(Ls) < h:
                            Lhi_low = Lsu[i]
                        diff = np.abs(np.sum(Ls)-h)
            # check if cable straight to seabed is solution
            if np.sum(L-Lsu+Ls)+tol >= h+d:
                # no horizontal tension
                a = 0
                Ls = h
                Le = Lt+et
                x0 = 0.
            else:
                # check if line is partly or fully lifted
                f = lambda a: a*(np.cosh(d/a)-1)-h
                a = utils.bisection(f, self.bisection_int1, self.bisection_int2, tol=tol, maxit=maxit)
                Ls0 = a*np.sinh(d/a)  # maximum line length to be fully lifted
                # get actual line length assuming it is fully lifted (from parameter a)
                H = a*np.sum(w*L)/Lt
                Va = 0
                for i in range(len(e)):
                    e[i] = np.sqrt(H**2+(Va+np.sum(w[:i]*Lsu[:i])+w[i]*Lsu[i]/2.)**2)*Lsu[i]/EA[i]
                Ls1 = Lt+np.sum(e)
                if Ls1 > Ls0:  # partly lifted
                    a, e, Lsu = utils.partly_lifted_elastic(d=d, h=h, L=L, w=w, EA=EA, maxit=maxit, tol=tol)
                    x0 = a*np.arccosh(1+h/a)
                    y_offset = -a
                elif Ls1 <= Ls0:  # fully lifted
                    x0 = d
                    Ls[:] = L
                    a, e = utils.fully_lifted_elastic(d=d, h=h, L=L, w=w, EA=EA, int1=a, maxit=maxit, tol=tol)
                    if a is np.nan:  # assume line is straight
                        raise RuntimeError('The line is too stretched/straight, cannot find catenary shape')
                        # H, e = utils.straight_elastic(d=d, h=h, L=L, w=w, EA=EA, maxit=maxit, tol=tol)
                        # a = H/np.sum(w*L)
                    Lst = np.sum(Ls+e)
                    xx = 0.5*(a*np.log((Lst+h)/(Lst-h))-d)
                    xy = 0.5*(a*np.log((Lst+h)/(Lst-h))+d)
                    x_offset = -xx
                    y_offset = h-a*np.cosh(xy/a)
                    s_offset = a*np.sinh(xx/a)
        self.Ls = Ls
        self._x_offset = x_offset
        self._y_offset = y_offset
        self._s_offset = s_offset
        self.a = a
        self.e = e
        self.x0 = x0
