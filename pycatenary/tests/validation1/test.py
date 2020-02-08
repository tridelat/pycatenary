import os
import unittest
import numpy as np
import numpy.testing as npt
from pycatenary import cable


def csv2array(filename, delimiter=',', names=None):
    fname = os.path.join(os.path.dirname(__file__), filename)
    return np.genfromtxt(fname, delimiter=delimiter, names=names)


class TestCatenaryValidation(unittest.TestCase):

    def test_elastic(self):
        # results to compare to
        pos = csv2array('xpos.csv')
        T = csv2array('T.csv')
        Tx = csv2array('Tx.csv')
        Ty = csv2array('Ty.csv')

        # make mooring line
        length = 6.98
        w = 1.036
        EA = 560e3
        anchor1 = [0.,0., 0.]
        fairlead1 = [5.3, 0., 2.65]
        l1 = cable.MooringLine(L=length, w=w, EA=EA, anchor=anchor1, fairlead=fairlead1)
        Tfs = []
        Txs = []
        Tys = []
        x0s = []
        for x in pos:
            l1.setFairleadCoords(np.array([x, 0., 2.65]))
            l1.computeSolution()
            TT = l1.getTension(length)
            Tfs += [np.linalg.norm(TT)]
            Txs += [TT[0]]
            Tys += [TT[1]]
            x0s += [l1.catenary.x0]
        L2 = np.linalg.norm(np.array(Tfs)-T)/len(T)
        npt.assert_almost_equal(L2, 0.0013525612473249868)

    def test_rigid(self):
        # results to compare to
        pos = csv2array('xpos.csv')
        T = csv2array('T.csv')
        Tx = csv2array('Tx.csv')
        Ty = csv2array('Ty.csv')

        # make mooring line
        length = 6.98
        w = 1.036
        EA = None
        anchor1 = [0.,0., 0.]
        fairlead1 = [5.3, 0., 2.65]
        l1 = cable.MooringLine(L=length, w=w, EA=EA, anchor=anchor1, fairlead=fairlead1)
        Tfs = []
        Txs = []
        Tys = []
        x0s = []
        for x in pos:
            l1.setFairleadCoords(np.array([x, 0., 2.65]))
            l1.computeSolution()
            TT = l1.getTension(length)
            Tfs += [np.linalg.norm(TT)]
            Txs += [TT[0]]
            Tys += [TT[1]]
            x0s += [l1.catenary.x0]
            print(l1.catenary.Ls, TT)
        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots()
        # ax.plot(pos, Tfs, 'b-')
        # ax.plot(pos, T, 'r--')
        # plt.show()
        L2 = np.linalg.norm(np.array(Tfs)-T)/len(T)
        npt.assert_almost_equal(L2, 0.0015298813773942969)


if __name__ == "__main__":
    unittest.main()
