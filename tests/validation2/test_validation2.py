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
    elastic: bool = True, floor: bool = True, heavy_section: bool = False
) -> cable.MooringLine:
    """Returns a mooring line instance.

    Parameters
    ----------
    elastic: bool
        If True, the cable is elastic.
    floor: bool
        If True, the floor is assumed to be at the anchor level.
    """
    diameter = 0.185 * 1.80  # studless chain equivalent outer diameter
    area = np.pi * diameter**2 / 4.0  # area of cable
    w = (685.0 - area * 1025.0) * 9.81  # submerged weight of cable
    if elastic:
        EA = 3.27e9 if not heavy_section else [3.27e9, 3.27e9]
    else:
        EA = None

    # create cable instance
    mooring = cable.MooringLine(
        fairlead=[-58.0, 0.0, -14.0],
        anchor=[-837.6, 0.0, -200],
        L=850.0 if not heavy_section else [840.0, 10.0],
        w=w if not heavy_section else [w, w * 5],
        EA=EA,
        floor=floor,
    )

    return mooring


class TestSingleLineFOWT(unittest.TestCase):
    def setUp(self):
        self.save_test2ref = False  # whether to save test results to ref file
        self.compare_test = True  # compare test results to ref

    def test_elastic(self):
        ref_filename = "elastic.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(elastic=True, floor=True)
        length = np.sum(mooring.catenary.L)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], T_ref[-1, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[-1, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], T_ref[0, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[0, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_rigid(self):
        ref_filename = "rigid.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(elastic=False, floor=True)
        length = np.sum(mooring.catenary.L)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], T_ref[-1, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[-1, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], T_ref[0, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[0, 2])
        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_elastic_nofloor(self):
        ref_filename = "elastic_nofloor.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(elastic=True, floor=False)
        length = np.sum(mooring.catenary.L)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], T_ref[-1, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[-1, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], T_ref[0, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[0, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_rigid_nofloor(self):
        ref_filename = "rigid_nofloor.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(elastic=False, floor=False)
        length = np.sum(mooring.catenary.L)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], T_ref[-1, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[-1, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], T_ref[0, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[0, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_elastic_reversed(self):
        ref_filename = "elastic.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(elastic=True, floor=True)
        length = np.sum(mooring.catenary.L)

        # switch anchor and fairlead positions
        anchor_position = mooring.get_anchor_position()
        fairlead_position = mooring.get_fairlead_position()
        mooring.set_anchor_position(fairlead_position + 1)
        mooring.set_fairlead_position(anchor_position)
        mooring.set_anchor_position(fairlead_position)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s, from_fairlead=True)
            xyz_test[ii] = mooring.get_position(s, from_fairlead=True)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], T_ref[0, 0])
            npt.assert_almost_equal(Tf[1], T_ref[0, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[0, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Ta[1], T_ref[-1, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[-1, 2])

    def test_rigid_reversed(self):
        ref_filename = "rigid.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(elastic=False, floor=True)
        length = np.sum(mooring.catenary.L)

        # switch anchor and fairlead positions
        anchor_position = mooring.get_anchor_position()
        fairlead_position = mooring.get_fairlead_position()
        mooring.set_anchor_position(fairlead_position + 1)
        mooring.set_fairlead_position(anchor_position)
        mooring.set_anchor_position(fairlead_position)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s, from_fairlead=True)
            xyz_test[ii] = mooring.get_position(s, from_fairlead=True)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], T_ref[0, 0])
            npt.assert_almost_equal(Tf[1], T_ref[0, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[0, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Ta[1], T_ref[-1, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[-1, 2])

    def test_elastic_vertical(self):
        ref_filename = "elastic_line_too_long.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(elastic=True, floor=True)
        length = np.sum(mooring.catenary.L)

        # set fairlead position
        mooring.set_fairlead_position(
            mooring.get_fairlead_position() - np.array([150.0, 0.0, 0.0])
        )

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], T_ref[-1, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[-1, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], T_ref[0, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[0, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_rigid_vertical(self):
        ref_filename = "rigid_line_too_long.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(elastic=False, floor=True)
        length = np.sum(mooring.catenary.L)

        # set fairlead position
        mooring.set_fairlead_position(
            mooring.get_fairlead_position() - np.array([150.0, 0.0, 0.0])
        )

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], T_ref[-1, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[-1, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], T_ref[0, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[0, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_elastic_fully_lifted(self):
        ref_filename = "elastic_fully_lifted.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(elastic=True, floor=False)
        length = np.sum(mooring.catenary.L)

        # set fairlead position
        mooring.set_fairlead_position(np.array([-10.0, 0.0, 0.0]))

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], T_ref[-1, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[-1, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], T_ref[0, 1])
            npt.assert_almost_equal(Ta[2], T_ref[0, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_rigid_fully_lifted(self):
        ref_filename = "rigid_fully_lifted.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(elastic=False, floor=True)
        length = np.sum(mooring.catenary.L)

        # set fairlead position
        mooring.set_fairlead_position(np.array([-20.0, 0.0, 0.0]))

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], T_ref[-1, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[-1, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], T_ref[0, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[0, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_elastic_heavy_section(self):
        ref_filename = "elastic_heavy_section.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(
            elastic=True, floor=True, heavy_section=True
        )
        length = np.sum(mooring.catenary.L)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], T_ref[-1, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[-1, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], T_ref[0, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[0, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_elastic_heavy_section_reversed(self):
        ref_filename = "elastic.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(
            elastic=True, floor=True, heavy_section=True
        )
        length = np.sum(mooring.catenary.L)

        # switch anchor and fairlead positions
        anchor_position = mooring.get_anchor_position()
        fairlead_position = mooring.get_fairlead_position()
        mooring.set_anchor_position(fairlead_position + 1)
        mooring.set_fairlead_position(anchor_position)
        mooring.set_anchor_position(fairlead_position)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test = mooring.get_tension(s)
            xyz_test = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[0], T_ref[-(ii + 1), 0])
                npt.assert_almost_equal(T_test[1], T_ref[-(ii + 1), 1])
                npt.assert_almost_equal(T_test[2], T_ref[-(ii + 1), 2])
                npt.assert_almost_equal(xyz_test[0], xyz_ref[-(ii + 1), 0])
                npt.assert_almost_equal(xyz_test[1], xyz_ref[-(ii + 1), 1])
                npt.assert_almost_equal(xyz_test[2], xyz_ref[-(ii + 1), 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], T_ref[0, 0])
            npt.assert_almost_equal(Tf[1], T_ref[0, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[0, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Ta[1], T_ref[-1, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[-1, 2])

    def test_rigid_heavy_section(self):
        ref_filename = "rigid_heavy_section.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(
            elastic=False, floor=True, heavy_section=True
        )
        length = np.sum(mooring.catenary.L)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = mooring.get_tension(s)
            xyz_test[ii] = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], T_ref[-1, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[-1, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], T_ref[0, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[0, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_rigid_heavy_section_reversed(self):
        ref_filename = "rigid.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(
            elastic=False, floor=True, heavy_section=True
        )
        length = np.sum(mooring.catenary.L)

        # switch anchor and fairlead positions
        anchor_position = mooring.get_anchor_position()
        fairlead_position = mooring.get_fairlead_position()
        mooring.set_anchor_position(fairlead_position + 1)
        mooring.set_fairlead_position(anchor_position)
        mooring.set_anchor_position(fairlead_position)

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test = mooring.get_tension(s)
            xyz_test = mooring.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[0], T_ref[-(ii + 1), 0])
                npt.assert_almost_equal(T_test[1], T_ref[-(ii + 1), 1])
                npt.assert_almost_equal(T_test[2], T_ref[-(ii + 1), 2])
                npt.assert_almost_equal(xyz_test[0], xyz_ref[-(ii + 1), 0])
                npt.assert_almost_equal(xyz_test[1], xyz_ref[-(ii + 1), 1])
                npt.assert_almost_equal(xyz_test[2], xyz_ref[-(ii + 1), 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = mooring.get_fairlead_force()
            npt.assert_almost_equal(Tf[0], T_ref[0, 0])
            npt.assert_almost_equal(Tf[1], T_ref[0, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[0, 2])
            # compare tension at anchor
            Ta = mooring.get_anchor_force()
            npt.assert_almost_equal(Ta[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Ta[1], T_ref[-1, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[-1, 2])

    def test_elastic_120degrees(self):
        ref_filename = "elastic.csv"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        # create mooring line
        mooring = get_mooring_line(elastic=True, floor=True)
        length = np.sum(mooring.catenary.L)

        # rotate anchor and fairlead positions (120 degrees)
        fairlead_position = mooring.get_fairlead_position()
        anchor_position = mooring.get_anchor_position()
        mooring.set_fairlead_position(
            rotate_vector_2d(fairlead_position, 120 * np.pi / 180)
        )
        mooring.set_anchor_position(
            rotate_vector_2d(anchor_position, 120 * np.pi / 180)
        )

        # compute solution
        mooring.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, NPOINTS)
        for ii, s in enumerate(ss_test):
            # rotate tension and position back (-120 degrees)
            tension = mooring.get_tension(s)
            # necessary because tension is strictly positive otherwise
            tension[0] = -tension[0]
            T_test = rotate_vector_2d(tension, -120 * np.pi / 180)
            xyz_test = rotate_vector_2d(
                mooring.get_position(s), -120 * np.pi / 180
            )

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[2], xyz_ref[ii, 2])

        if self.compare_test:
            # compare tension at fairlead
            Tf = rotate_vector_2d(
                mooring.get_fairlead_force(), -120 * np.pi / 180
            )
            npt.assert_almost_equal(Tf[0], -T_ref[-1, 0])
            npt.assert_almost_equal(Tf[1], T_ref[-1, 1])
            npt.assert_almost_equal(Tf[2], -T_ref[-1, 2])
            # compare tension at anchor
            Ta = rotate_vector_2d(
                mooring.get_anchor_force(), -120 * np.pi / 180
            )
            npt.assert_almost_equal(Ta[0], T_ref[0, 0])
            npt.assert_almost_equal(Ta[1], T_ref[0, 1])
            npt.assert_almost_equal(Ta[2], -T_ref[0, 2])


if __name__ == "__main__":
    unittest.main()
