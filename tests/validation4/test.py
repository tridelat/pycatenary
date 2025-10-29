import os
import unittest

import numpy as np
import numpy.testing as npt

from pycatenary import cable

# number of points along line for testing
NPOINTS = 101


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
    elastic: bool = True, multisegmented: bool = False
) -> tuple[cable.MooringLine, cable.MooringLine]:
    """Returns a mooring line instance.

    Parameters
    ----------
    elastic: bool
        If True, the cable is elastic.
    floor: bool
        If True, the floor is assumed to be at the anchor level.
    """

    length = 20.0 if not multisegmented else [5.0, 10.0, 2.5, 2.5]
    w = 1.962 if not multisegmented else [1.962] * 4
    if not elastic:
        EA = None
    else:
        # make it very extensible on purpose
        EA = 1e1 if not multisegmented else [1e1] * 4

    # define properties of cable
    mooring = cable.MooringLine(
        fairlead=[17.69, 0.0],
        anchor=[0.0, 0.0],
        L=length,
        w=w,
        EA=EA,
        floor=False,
    )

    return mooring


class TestHangingCable(unittest.TestCase):
    def setUp(self):
        self.save_test2ref = False  # whether to save test results to ref file
        self.compare_test = True  # compare test results to ref

    def test_rigid(self):
        ref_filename = "rigid.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"]))

        # create mooring line
        mooring = get_mooring_line(elastic=False, multisegmented=False)
        length = np.sum(mooring.catenary.L)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 2))
        xyz_test = np.zeros((len(ss_test), 2))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], -T_ref[-1, 1])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], -T_ref[0, 1])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,Tx,Ty",
            )

    def test_elastic(self):
        ref_filename = "elastic.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"]))

        # create mooring line
        mooring = get_mooring_line(elastic=True, multisegmented=False)
        length = np.sum(mooring.catenary.L)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 2))
        xyz_test = np.zeros((len(ss_test), 2))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], -T_ref[-1, 1])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], -T_ref[0, 1])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,Tx,Ty",
            )

    def test_rigid_multisegmented(self):
        ref_filename = "rigid.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"]))

        # create mooring line
        mooring = get_mooring_line(elastic=False, multisegmented=True)
        length = np.sum(mooring.catenary.L)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 2))
        xyz_test = np.zeros((len(ss_test), 2))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], -T_ref[-1, 1])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], -T_ref[0, 1])
