import unittest
import numpy as np
import numpy.testing as npt
from pycatenary import cable

class TestCatenaryValidation(unittest.TestCase):

    def test_validation1(self):
        # results to compare to
        pos = np.genfromtxt('xpos.csv', delimiter=',')
        T = np.genfromtxt('T.csv', delimiter=',')
        Tx = np.genfromtxt('Tx.csv', delimiter=',')
        Ty = np.genfromtxt('Ty.csv', delimiter=',')

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
        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots()
        # ax.plot(pos, Tfs, 'b-')
        # ax.plot(pos, T, 'r--')
        # plt.show()
        L2 = np.linalg.norm(np.array(Tfs)-T)/len(T)
        npt.assert_almost_equal(0.001352561247, L2)


if __name__ == "__main__":
    unittest.main()
