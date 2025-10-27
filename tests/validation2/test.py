import os
import unittest

import numpy as np
import numpy.testing as npt

from pycatenary import cable


def csv2array(filename, delimiter=",", names=None):
    fname = os.path.join(os.path.dirname(__file__), filename)
    return np.genfromtxt(fname, delimiter=delimiter, names=names)


def array2csv(filename, array, delimiter=",", names=None):
    fname = os.path.join(os.path.dirname(__file__), filename)
    np.savetxt(fname, array, delimiter=delimiter, header=names, comments="")


class TestCatenaryValidation(unittest.TestCase):
    def setUp(self):
        self.save_test2ref = False  # whether to save test results to ref file
        self.compare_test = True  # compare test results to ref

    def test_elastic(self):
        ref_filename = "elastic.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        length = 850  # length of line
        diameter = 0.185 * 1.80  # diameter of cable
        area = np.pi * diameter**2 / 4.0  # area of cable
        w = (685.0 - area * 1025.0) * 9.81  # submerged weight of cable

        # create cable instance
        l1 = cable.MooringLine(
            fairlead=[-58.0, 0.0, -14.0],
            anchor=[-837.6, 0.0, -200],
            L=length,
            w=w,
            EA=3.27e9,
            floor=True,
        )
        l1.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, 101)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = l1.get_tension(s)
            xyz_test[ii] = l1.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_rigid(self):
        ref_filename = "rigid.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        length = 850  # length of line
        diameter = 0.185 * 1.80  # diameter of cable
        area = np.pi * diameter**2 / 4.0  # area of cable
        w = (685.0 - area * 1025.0) * 9.81  # submerged weight of cable

        # create cable instance
        l1 = cable.MooringLine(
            fairlead=[-58.0, 0.0, -14.0],
            anchor=[-837.6, 0.0, -200],
            L=length,
            w=w,
            EA=None,
            floor=True,
        )
        l1.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, 101)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = l1.get_tension(s)
            xyz_test[ii] = l1.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_elastic_nofloor(self):
        ref_filename = "elastic_nofloor.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        length = 850  # length of line
        diameter = 0.185 * 1.80  # diameter of cable
        area = np.pi * diameter**2 / 4.0  # area of cable
        w = (685.0 - area * 1025.0) * 9.81  # submerged weight of cable

        # create cable instance
        l1 = cable.MooringLine(
            fairlead=[-58.0, 0.0, -14.0],
            anchor=[-837.6, 0.0, -200],
            L=length,
            w=w,
            EA=3.27e9,
            floor=False,
        )
        l1.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, 101)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = l1.get_tension(s)
            xyz_test[ii] = l1.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_rigid_nofloor(self):
        ref_filename = "rigid_nofloor.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        length = 850  # length of line
        diameter = 0.185 * 1.80  # diameter of cable
        area = np.pi * diameter**2 / 4.0  # area of cable
        w = (685.0 - area * 1025.0) * 9.81  # submerged weight of cable

        # create cable instance
        l1 = cable.MooringLine(
            fairlead=[-58.0, 0.0, -14.0],
            anchor=[-837.6, 0.0, -200],
            L=length,
            w=w,
            EA=None,
            floor=False,
        )
        l1.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, 101)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = l1.get_tension(s)
            xyz_test[ii] = l1.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_elastic_anchor_above(self):
        ref_filename = "elastic.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        length = 850  # length of line
        diameter = 0.185 * 1.80  # diameter of cable
        area = np.pi * diameter**2 / 4.0  # area of cable
        w = (685.0 - area * 1025.0) * 9.81  # submerged weight of cable

        # create cable instance
        l1 = cable.MooringLine(
            anchor=[-58.0, 0.0, -14.0],
            fairlead=[-837.6, 0.0, -200],
            L=length,
            w=w,
            EA=3.27e9,
            floor=True,
        )
        l1.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, 101)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = l1.get_tension(s)
            xyz_test[ii] = l1.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[-(ii + 1), 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[-(ii + 1), 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[-(ii + 1), 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[-(ii + 1), 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[-(ii + 1), 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[-(ii + 1), 2])

    def test_rigid_anchor_above(self):
        ref_filename = "rigid.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        length = 850  # length of line
        diameter = 0.185 * 1.80  # diameter of cable
        area = np.pi * diameter**2 / 4.0  # area of cable
        w = (685.0 - area * 1025.0) * 9.81  # submerged weight of cable

        # create cable instance
        l1 = cable.MooringLine(
            anchor=[-58.0, 0.0, -14.0],
            fairlead=[-837.6, 0.0, -200],
            L=length,
            w=w,
            EA=None,
            floor=True,
        )
        l1.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, 101)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = l1.get_tension(s)
            xyz_test[ii] = l1.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[-(ii + 1), 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[-(ii + 1), 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[-(ii + 1), 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[-(ii + 1), 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[-(ii + 1), 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[-(ii + 1), 2])

    def test_elastic_line_too_long(self):
        ref_filename = "elastic_line_too_long.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        length = 850.0  # length of line
        diameter = 0.185 * 1.80  # diameter of cable
        area = np.pi * diameter**2 / 4.0  # area of cable
        w = (685.0 - area * 1025.0) * 9.81  # submerged weight of cable

        # create cable instance
        l1 = cable.MooringLine(
            fairlead=[-58.0 - 150.0, 0.0, -14.0],
            anchor=[-837.6, 0.0, -200],
            L=length,
            w=w,
            EA=3.27e9,
            floor=True,
        )
        l1.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, 101)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = l1.get_tension(s)
            xyz_test[ii] = l1.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_rigid_line_too_long(self):
        ref_filename = "rigid_line_too_long.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        length = 850.0  # length of line
        diameter = 0.185 * 1.80  # diameter of cable
        area = np.pi * diameter**2 / 4.0  # area of cable
        w = (685.0 - area * 1025.0) * 9.81  # submerged weight of cable

        # create cable instance
        l1 = cable.MooringLine(
            fairlead=[-58.0 - 150.0, 0.0, -14.0],
            anchor=[-837.6, 0.0, -200],
            L=length,
            w=w,
            EA=None,
            floor=True,
        )
        l1.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, 101)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = l1.get_tension(s)
            xyz_test[ii] = l1.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_elastic_fully_lifted(self):
        ref_filename = "elastic_fully_lifted.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        length = 850  # length of line
        diameter = 0.185 * 1.80  # diameter of cable
        area = np.pi * diameter**2 / 4.0  # area of cable
        w = (685.0 - area * 1025.0) * 9.81  # submerged weight of cable

        # create cable instance
        l1 = cable.MooringLine(
            fairlead=[-0.0, 0.0, -14.0],
            anchor=[-837.6, 0.0, -200],
            L=length,
            w=w,
            EA=3.27e9,
            floor=False,
        )
        l1.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, 101)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = l1.get_tension(s)
            xyz_test[ii] = l1.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )

    def test_rigid_fully_lifted(self):
        ref_filename = "rigid_fully_lifted.txt"
        if self.compare_test:
            ref = csv2array(ref_filename, names=True, delimiter=",")
            xyz_ref = np.column_stack((ref["x"], ref["y"], ref["z"]))
            T_ref = np.column_stack((ref["Tx"], ref["Ty"], ref["Tz"]))

        length = 850  # length of line
        diameter = 0.185 * 1.80  # diameter of cable
        area = np.pi * diameter**2 / 4.0  # area of cable
        w = (685.0 - area * 1025.0) * 9.81  # submerged weight of cable

        # create cable instance
        l1 = cable.MooringLine(
            fairlead=[-58.0 - 100.0, 0.0, -14.0],
            anchor=[-837.6, 0.0, -200],
            L=length,
            w=w,
            EA=None,
            floor=True,
        )
        l1.compute_solution()

        # test for different positions of fairlead
        ss_test = np.linspace(0.0, length, 101)
        T_test = np.zeros((len(ss_test), 3))
        xyz_test = np.zeros((len(ss_test), 3))
        for ii, s in enumerate(ss_test):
            T_test[ii] = l1.get_tension(s)
            xyz_test[ii] = l1.get_position(s)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(T_test[ii, 0], T_ref[ii, 0])
                npt.assert_almost_equal(T_test[ii, 1], T_ref[ii, 1])
                npt.assert_almost_equal(T_test[ii, 2], T_ref[ii, 2])
                npt.assert_almost_equal(xyz_test[ii, 0], xyz_ref[ii, 0])
                npt.assert_almost_equal(xyz_test[ii, 1], xyz_ref[ii, 1])
                npt.assert_almost_equal(xyz_test[ii, 2], xyz_ref[ii, 2])

        if self.save_test2ref:
            stack = np.column_stack((ss_test, xyz_test, T_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="s,x,y,z,Tx,Ty,Tz",
            )


if __name__ == "__main__":
    unittest.main()
