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
        ref = csv2array(ref_filename, names=True, delimiter=",")
        xpos = ref["xpos"]
        Tf_ref = np.column_stack((ref["Tfx"], ref["Tfy"], ref["Tfz"]))
        Ls_ref = ref["Ls"]

        # make mooring line
        length = 6.98
        l1 = cable.MooringLine(
            L=length,
            w=1.036,
            EA=560e3,
            anchor=[0.0, 0.0, 0.0],
            fairlead=[5.3, 0.0, 2.65],
            floor=True,
        )

        # test for different positions of fairlead
        Tf_test = np.zeros_like(Tf_ref)
        Ls_test = np.zeros(len(Tf_ref))
        for ii, x in enumerate(xpos):
            l1.set_fairlead_position(np.array([x, 0.0, 2.65]))
            l1.compute_solution()
            # tension at fairlead
            Tf = l1.get_tension(length)
            Tf_test[ii] = Tf
            # total lifted line length
            Ls_test[ii] = np.sum(l1.catenary.Ls)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(Tf_test[ii, 0], Tf_ref[ii, 0])
                npt.assert_almost_equal(Tf_test[ii, 1], Tf_ref[ii, 1])
                npt.assert_almost_equal(Tf_test[ii, 2], Tf_ref[ii, 2])
                npt.assert_almost_equal(Ls_test[ii], Ls_ref[ii])

        if self.save_test2ref:
            stack = np.column_stack((xpos, Ls_test, Tf_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="xpos,Ls,Tfx,Tfy,Tfz",
            )

    def test_elastic_multisegmented(self):
        ref_filename = "elastic.txt"
        ref = csv2array(ref_filename, names=True, delimiter=",")
        xpos = ref["xpos"]
        Tf_ref = np.column_stack((ref["Tfx"], ref["Tfy"], ref["Tfz"]))
        Ls_ref = ref["Ls"]

        # make mooring line
        length = 6.98
        l1 = cable.MooringLine(
            L=[length / 3.0, length * 2.0 / 3.0],
            w=[1.036, 1.036],
            EA=[560e3, 560e3],
            anchor=[0.0, 0.0, 0.0],
            fairlead=[5.3, 0.0, 2.65],
            floor=True,
        )

        # test for different positions of fairlead
        for ii, x in enumerate(xpos):
            l1.set_fairlead_position(np.array([x, 0.0, 2.65]))
            l1.compute_solution()
            # tension at fairlead
            Tf = l1.get_tension(length)
            # total lifted line length
            Ls = np.sum(l1.catenary.Ls)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(Tf[0], Tf_ref[ii, 0])
                npt.assert_almost_equal(Tf[1], Tf_ref[ii, 1])
                npt.assert_almost_equal(Tf[2], Tf_ref[ii, 2])
                npt.assert_almost_equal(Ls, Ls_ref[ii])

    def test_rigid(self):
        # load reference data
        ref_filename = "rigid.txt"
        ref = csv2array(ref_filename, names=True, delimiter=",")
        xpos = ref["xpos"]
        Tf_ref = np.column_stack((ref["Tfx"], ref["Tfy"], ref["Tfz"]))
        Ls_ref = ref["Ls"]

        # make mooring line
        length = 6.98
        l1 = cable.MooringLine(
            L=length,
            w=1.036,
            EA=None,
            anchor=[0.0, 0.0, 0.0],
            fairlead=[5.3, 0.0, 2.65],
            floor=True,
        )

        # test for different positions of fairlead
        Tf_test = np.zeros_like(Tf_ref)
        Ls_test = np.zeros(len(Tf_ref))
        for ii, x in enumerate(xpos):
            l1.set_fairlead_position(np.array([x, 0.0, 2.65]))
            l1.compute_solution()
            # tension at fairlead
            Tf = l1.get_tension(length)
            Tf_test[ii] = Tf
            # total lifted line length
            Ls_test[ii] = np.sum(l1.catenary.Ls)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(Tf_test[ii, 0], Tf_ref[ii, 0])
                npt.assert_almost_equal(Tf_test[ii, 1], Tf_ref[ii, 1])
                npt.assert_almost_equal(Tf_test[ii, 2], Tf_ref[ii, 2])
                npt.assert_almost_equal(Ls_test[ii], Ls_ref[ii])

        if self.save_test2ref:
            stack = np.column_stack((xpos, Ls_test, Tf_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="xpos,Ls,Tfx,Tfy,Tfz",
            )

    def test_rigid_multisegmented(self):
        ref_filename = "rigid.txt"
        ref = csv2array(ref_filename, names=True, delimiter=",")
        xpos = ref["xpos"]
        Tf_ref = np.column_stack((ref["Tfx"], ref["Tfy"], ref["Tfz"]))
        Ls_ref = ref["Ls"]

        # make mooring line
        length = 6.98
        l1 = cable.MooringLine(
            L=[length / 3.0, length * 2.0 / 3.0],
            w=[1.036, 1.036],
            EA=None,
            anchor=[0.0, 0.0, 0.0],
            fairlead=[5.3, 0.0, 2.65],
            floor=True,
        )

        # test for different positions of fairlead
        for ii, x in enumerate(xpos):
            l1.set_fairlead_position(np.array([x, 0.0, 2.65]))
            l1.compute_solution()
            # tension at fairlead
            Tf = l1.get_tension(length)
            # total lifted line length
            Ls = np.sum(l1.catenary.Ls)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(Tf[0], Tf_ref[ii, 0])
                npt.assert_almost_equal(Tf[1], Tf_ref[ii, 1])
                npt.assert_almost_equal(Tf[2], Tf_ref[ii, 2])
                npt.assert_almost_equal(Ls, Ls_ref[ii])

    def test_rigid_nofloor(self):
        # load reference data
        ref_filename = "rigid_nofloor.txt"
        ref = csv2array(ref_filename, names=True, delimiter=",")
        xpos = ref["xpos"]
        Tf_ref = np.column_stack((ref["Tfx"], ref["Tfy"], ref["Tfz"]))
        Ta_ref = np.column_stack((ref["Tax"], ref["Tay"], ref["Taz"]))
        Ls_ref = ref["Ls"]

        # make mooring line
        length = 6.98
        l1 = cable.MooringLine(
            L=length,
            w=1.036,
            EA=None,
            anchor=[0.0, 0.0, 0.0],
            fairlead=[5.3, 0.0, 2.65],
            floor=False,
        )

        # test for different positions of fairlead
        Tf_test = np.zeros_like(Tf_ref)
        Ta_test = np.zeros_like(Ta_ref)
        Ls_test = np.zeros(len(Tf_ref))
        for ii, x in enumerate(xpos):
            l1.set_fairlead_position(np.array([x, 0.0, 2.65]))
            l1.compute_solution()
            # tension at fairlead
            Tf = l1.get_tension(length)
            Tf_test[ii] = Tf
            # tension at anchor
            Ta = l1.get_tension(0.0)
            Ta_test[ii] = Ta
            # total lifted line length
            Ls_test[ii] = np.sum(l1.catenary.Ls)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(Tf_test[ii, 0], Tf_ref[ii, 0])
                npt.assert_almost_equal(Tf_test[ii, 1], Tf_ref[ii, 1])
                npt.assert_almost_equal(Tf_test[ii, 2], Tf_ref[ii, 2])
                npt.assert_almost_equal(Ta_test[ii, 0], Ta_ref[ii, 0])
                npt.assert_almost_equal(Ta_test[ii, 1], Ta_ref[ii, 1])
                npt.assert_almost_equal(Ta_test[ii, 2], Ta_ref[ii, 2])
                npt.assert_almost_equal(Ls_test[ii], Ls_ref[ii])

        if self.save_test2ref:
            stack = np.column_stack((xpos, Ls_test, Tf_test, Ta_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="xpos,Ls,Tfx,Tfy,Tfz,Tax,Tay,Taz",
            )

    def test_rigid_nofloor_multisegmented(self):
        # load reference data
        ref_filename = "rigid_nofloor.txt"
        ref = csv2array(ref_filename, names=True, delimiter=",")
        xpos = ref["xpos"]
        Tf_ref = np.column_stack((ref["Tfx"], ref["Tfy"], ref["Tfz"]))
        Ta_ref = np.column_stack((ref["Tax"], ref["Tay"], ref["Taz"]))
        Ls_ref = ref["Ls"]

        # make mooring line
        length = 6.98
        l1 = cable.MooringLine(
            L=[length / 3.0, length * 2.0 / 3.0],
            w=[1.036, 1.036],
            EA=None,
            anchor=[0.0, 0.0, 0.0],
            fairlead=[5.3, 0.0, 2.65],
            floor=False,
        )

        # test for different positions of fairlead
        Tf_test = np.zeros_like(Tf_ref)
        Ta_test = np.zeros_like(Ta_ref)
        Ls_test = np.zeros(len(Tf_ref))
        for ii, x in enumerate(xpos):
            l1.set_fairlead_position(np.array([x, 0.0, 2.65]))
            l1.compute_solution()
            # tension at fairlead
            Tf = l1.get_tension(length)
            Tf_test[ii] = Tf
            # tension at anchor
            Ta = l1.get_tension(0.0)
            Ta_test[ii] = Ta
            # total lifted line length
            Ls_test[ii] = np.sum(l1.catenary.Ls)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(Tf_test[ii, 0], Tf_ref[ii, 0])
                npt.assert_almost_equal(Tf_test[ii, 1], Tf_ref[ii, 1])
                npt.assert_almost_equal(Tf_test[ii, 2], Tf_ref[ii, 2])
                npt.assert_almost_equal(Ta_test[ii, 0], Ta_ref[ii, 0])
                npt.assert_almost_equal(Ta_test[ii, 1], Ta_ref[ii, 1])
                npt.assert_almost_equal(Ta_test[ii, 2], Ta_ref[ii, 2])
                npt.assert_almost_equal(Ls_test[ii], Ls_ref[ii])

        if self.save_test2ref:
            stack = np.column_stack((xpos, Ls_test, Tf_test, Ta_test))
            array2csv(
                ref_filename,
                stack,
                delimiter=",",
                names="xpos,Ls,Tfx,Tfy,Tfz,Tax,Tay,Taz",
            )

    def test_rigid_nofloor_anchor_above(self):
        # load reference data
        ref_filename = "rigid_nofloor.txt"
        ref = csv2array(ref_filename, names=True, delimiter=",")
        xpos = ref["xpos"]
        Tf_ref = np.column_stack((ref["Tfx"], ref["Tfy"], ref["Tfz"]))
        Ta_ref = np.column_stack((ref["Tax"], ref["Tay"], ref["Taz"]))
        Ls_ref = ref["Ls"]

        # make mooring line
        length = 6.98
        l1 = cable.MooringLine(
            L=length,
            w=1.036,
            EA=None,
            anchor=[5.3, 0.0, 2.65],
            fairlead=[0.0, 0.0, 0.0],
            floor=False,
        )

        # test for different positions of fairlead
        for ii, x in enumerate(xpos):
            l1.set_anchor_position(np.array([x, 0.0, 2.65]))
            l1.compute_solution()
            # tension at fairlead
            Tf = l1.get_tension(0.0)
            # tension at anchor
            Ta = l1.get_tension(length)
            # total lifted line length
            Ls = np.sum(l1.catenary.Ls)

            # check solution
            if self.compare_test:
                npt.assert_almost_equal(Tf[0], Tf_ref[ii, 0])
                npt.assert_almost_equal(Tf[1], Tf_ref[ii, 1])
                npt.assert_almost_equal(Tf[2], Tf_ref[ii, 2])
                npt.assert_almost_equal(Ls, Ls_ref[ii])
                npt.assert_almost_equal(Ta[0], Ta_ref[ii, 0])
                npt.assert_almost_equal(Ta[1], Ta_ref[ii, 1])
                npt.assert_almost_equal(Ta[2], Ta_ref[ii, 2])


if __name__ == "__main__":
    unittest.main()
