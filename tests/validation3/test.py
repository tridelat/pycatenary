import os
import unittest

import numpy as np
import numpy.testing as npt

from pycatenary import cable

# number of points along line for testing
XPOS = [-0.91, -0.7, -0.51, -0.44, -0.3, -0.18, 0.0, 0.19, 0.35, 0.5, 0.9]


def csv2array(
    filename: str, delimiter: str = ",", names: bool = None
) -> np.ndarray:
    fname = os.path.join(os.path.dirname(__file__), filename)
    return np.genfromtxt(fname, delimiter=delimiter, names=names)


def array2csv(
    filename: str, array: np.ndarray, delimiter: str = ",", names: bool = None
) -> None:
    fname = os.path.join(os.path.dirname(__file__), filename)
    np.savetxt(fname, array, delimiter=delimiter, header=names, comments="")


def rotate_vector_2d(vector: np.ndarray, angle: float) -> np.ndarray:
    return np.array(
        [
            vector[0] * np.cos(angle) - vector[1] * np.sin(angle),
            vector[0] * np.sin(angle) + vector[1] * np.cos(angle),
            vector[2],
        ]
    )


def get_mooring_line(
    elastic: bool = True, floor: bool = True, heavy_section: bool = False
) -> tuple[cable.MooringLine, cable.MooringLine]:
    """Returns a mooring line instance.

    Parameters
    ----------
    elastic: bool
        If True, the cable is elastic.
    floor: bool
        If True, the floor is assumed to be at the anchor level.
    """
    L = [5.672, 0.126, 4.0, 0.259]
    w = np.array([0.402, 1.558, 0.00425, 1.529]) * 9.81

    if elastic:
        EA = [2.050e3, 3.636e6, 10.873e3, 6.464e6]
    else:
        EA = None

    fairlead1 = [0.142, 0.0, 5.486]
    fairlead2 = rotate_vector_2d(fairlead1, 2.0 * np.pi / 3.0)
    anchor1 = [7.083, 0.0, 0.0]
    anchor2 = rotate_vector_2d(anchor1, 2.0 * np.pi / 3.0)

    mooring1 = cable.MooringLine(
        fairlead=fairlead1, anchor=anchor1, L=L, w=w, EA=EA, floor=floor
    )
    mooring2 = cable.MooringLine(
        fairlead=fairlead2, anchor=anchor2, L=L, w=w, EA=EA, floor=floor
    )

    return (mooring1, mooring2)


class TestMultisegmentedWEC(unittest.TestCase):
    def setUp(self):
        self.save_test2ref = False  # whether to save test results to ref file
        self.compare_test = True  # compare test results to ref

    def test_elastic(self):
        ref_filename = "elastic.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            T1_ref = np.column_stack((ref["T1x"], ref["T1y"], ref["T1z"]))
            T2_ref = np.column_stack((ref["T2x"], ref["T2y"], ref["T2z"]))

        # create mooring line
        mooring1, mooring2 = get_mooring_line(elastic=True, floor=True)
        fairlead1 = mooring1.get_fairlead_position()
        fairlead2 = mooring2.get_fairlead_position()

        T1s = np.zeros((len(XPOS), 3))
        T2s = np.zeros((len(XPOS), 3))
        for ii, x in enumerate(XPOS):
            mooring1.set_fairlead_position(fairlead1 + [x, 0.0, 0.0])
            mooring1.compute_solution()
            T1s[ii] = mooring1.get_fairlead_force()

            mooring2.set_fairlead_position(fairlead2 + [x, 0.0, 0.0])
            mooring2.compute_solution()
            T2s[ii] = mooring2.get_fairlead_force()

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T1s[ii, 0], T1_ref[ii, 0])
                npt.assert_almost_equal(T1s[ii, 1], T1_ref[ii, 1])
                npt.assert_almost_equal(T1s[ii, 2], T1_ref[ii, 2])
                npt.assert_almost_equal(T2s[ii, 0], T2_ref[ii, 0])
                npt.assert_almost_equal(T2s[ii, 1], T2_ref[ii, 1])
                npt.assert_almost_equal(T2s[ii, 2], T2_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((XPOS, T1s, T2s))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="xpos,T1x,T1y,T1z,T2x,T2y,T2z",
            )

    def test_rigid(self):
        ref_filename = "rigid.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            T1_ref = np.column_stack((ref["T1x"], ref["T1y"], ref["T1z"]))
            T2_ref = np.column_stack((ref["T2x"], ref["T2y"], ref["T2z"]))

        # create mooring line
        mooring1, mooring2 = get_mooring_line(elastic=False, floor=True)
        fairlead1 = mooring1.get_fairlead_position()
        fairlead2 = mooring2.get_fairlead_position()

        T1s = np.zeros((len(XPOS), 3))
        T2s = np.zeros((len(XPOS), 3))
        for ii, x in enumerate(XPOS):
            mooring1.set_fairlead_position(fairlead1 + [x, 0.0, 0.0])
            mooring1.compute_solution()
            T1s[ii] = mooring1.get_fairlead_force()

            mooring2.set_fairlead_position(fairlead2 + [x, 0.0, 0.0])
            mooring2.compute_solution()
            T2s[ii] = mooring2.get_fairlead_force()

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T1s[ii, 0], T1_ref[ii, 0])
                npt.assert_almost_equal(T1s[ii, 1], T1_ref[ii, 1])
                npt.assert_almost_equal(T1s[ii, 2], T1_ref[ii, 2])
                npt.assert_almost_equal(T2s[ii, 0], T2_ref[ii, 0])
                npt.assert_almost_equal(T2s[ii, 1], T2_ref[ii, 1])
                npt.assert_almost_equal(T2s[ii, 2], T2_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((XPOS, T1s, T2s))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="xpos,T1x,T1y,T1z,T2x,T2y,T2z",
            )

    def test_elastic_nofloor(self):
        ref_filename = "elastic_nofloor.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            T1_ref = np.column_stack((ref["T1x"], ref["T1y"], ref["T1z"]))
            T2_ref = np.column_stack((ref["T2x"], ref["T2y"], ref["T2z"]))

        # create mooring line
        mooring1, mooring2 = get_mooring_line(elastic=True, floor=False)
        fairlead1 = mooring1.get_fairlead_position()
        fairlead2 = mooring2.get_fairlead_position()

        T1s = np.zeros((len(XPOS), 3))
        T2s = np.zeros((len(XPOS), 3))
        for ii, x in enumerate(XPOS):
            mooring1.set_fairlead_position(fairlead1 + [x, 0.0, 0.0])
            mooring1.compute_solution()
            T1s[ii] = mooring1.get_fairlead_force()

            mooring2.set_fairlead_position(fairlead2 + [x, 0.0, 0.0])
            mooring2.compute_solution()
            T2s[ii] = mooring2.get_fairlead_force()

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T1s[ii, 0], T1_ref[ii, 0])
                npt.assert_almost_equal(T1s[ii, 1], T1_ref[ii, 1])
                npt.assert_almost_equal(T1s[ii, 2], T1_ref[ii, 2])
                npt.assert_almost_equal(T2s[ii, 0], T2_ref[ii, 0])
                npt.assert_almost_equal(T2s[ii, 1], T2_ref[ii, 1])
                npt.assert_almost_equal(T2s[ii, 2], T2_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((XPOS, T1s, T2s))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="xpos,T1x,T1y,T1z,T2x,T2y,T2z",
            )

    def test_rigid_nofloor(self):
        ref_filename = "rigid_nofloor.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            T1_ref = np.column_stack((ref["T1x"], ref["T1y"], ref["T1z"]))
            T2_ref = np.column_stack((ref["T2x"], ref["T2y"], ref["T2z"]))

        # create mooring line
        mooring1, mooring2 = get_mooring_line(elastic=False, floor=False)
        fairlead1 = mooring1.get_fairlead_position()
        fairlead2 = mooring2.get_fairlead_position()

        T1s = np.zeros((len(XPOS), 3))
        T2s = np.zeros((len(XPOS), 3))
        for ii, x in enumerate(XPOS):
            mooring1.set_fairlead_position(fairlead1 + [x, 0.0, 0.0])
            mooring1.compute_solution()
            T1s[ii] = mooring1.get_fairlead_force()

            mooring2.set_fairlead_position(fairlead2 + [x, 0.0, 0.0])
            mooring2.compute_solution()
            T2s[ii] = mooring2.get_fairlead_force()

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T1s[ii, 0], T1_ref[ii, 0])
                npt.assert_almost_equal(T1s[ii, 1], T1_ref[ii, 1])
                npt.assert_almost_equal(T1s[ii, 2], T1_ref[ii, 2])
                npt.assert_almost_equal(T2s[ii, 0], T2_ref[ii, 0])
                npt.assert_almost_equal(T2s[ii, 1], T2_ref[ii, 1])
                npt.assert_almost_equal(T2s[ii, 2], T2_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((XPOS, T1s, T2s))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="xpos,T1x,T1y,T1z,T2x,T2y,T2z",
            )
