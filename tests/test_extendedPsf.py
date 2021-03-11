# This file is part of pipe_tasks.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest
import numpy as np
import tempfile

from lsst.afw import image as afwImage
from lsst.pipe.tasks import extendedPsf
import lsst.utils.tests

np.random.seed(51778)


def makeExtendedPsf(nExtendedPsf=1):
    ePsfImages = [afwImage.MaskedImageF(25, 25) for _ in range(nExtendedPsf)]
    for ePsfIm in ePsfImages:
        ePsfIm.image.array += np.random.rand(25, 25)
        ePsfIm.mask.array += np.random.choice(3, size=(25, 25))
    return ePsfImages


class ExtendedPsfTestCase(lsst.utils.tests.TestCase):
    """Test ExtendedPsf.
    """
    def setUp(self):
        self.defaultEPsf = makeExtendedPsf(1)[0]
        self.constantExtendedPsf = extendedPsf.ExtendedPsf(self.defaultEPsf)
        self.regions = ["NW", "SW", "E"]
        self.regionDetectors = [list(range(10)), list(range(10, 20)), list(range(20, 40))]
        self.regionalEPsfs = makeExtendedPsf(3)

    def tearDown(self):
        del self.defaultEPsf
        del self.regions
        del self.regionDetectors
        del self.regionalEPsfs

    def testRegionalPsfAddition(self):
        # start with either an empty instance, or one containing a default
        # extended PSF
        startsEmptyExtendedPsf = extendedPsf.ExtendedPsf()
        withDefaultExtendedPsf = extendedPsf.ExtendedPsf(self.defaultEPsf)
        self.assertEqual(len(startsEmptyExtendedPsf), 0)
        self.assertEqual(len(withDefaultExtendedPsf), 1)
        # add a couple of regional PSFs
        for j in range(2):
            startsEmptyExtendedPsf.addRegionalExtendedPsf(self.regionalEPsfs[j], self.regions[j],
                                                          self.regionDetectors[j])
            withDefaultExtendedPsf.addRegionalExtendedPsf(self.regionalEPsfs[j], self.regions[j],
                                                          self.regionDetectors[j])
        self.assertEqual(len(startsEmptyExtendedPsf), 2)
        self.assertEqual(len(withDefaultExtendedPsf), 3)
        # ensure we recover the correct regional PSF
        for j in range(2):
            for det in self.regionDetectors[j]:
                # try by calling the class directly...
                regPsf0, regPsf1 = startsEmptyExtendedPsf(det), withDefaultExtendedPsf(det)
                self.assertMaskedImagesAlmostEqual(regPsf0, self.regionalEPsfs[j])
                self.assertMaskedImagesAlmostEqual(regPsf1, self.regionalEPsfs[j])
                # ... when passing on a detector number to the
                # getRegionalExtendedPsf method...
                regPsf0 = startsEmptyExtendedPsf.getRegionalExtendedPsf(detector=det)
                regPsf1 = withDefaultExtendedPsf.getRegionalExtendedPsf(detector=det)
                self.assertMaskedImagesAlmostEqual(regPsf0, self.regionalEPsfs[j])
                self.assertMaskedImagesAlmostEqual(regPsf1, self.regionalEPsfs[j])
            # ... and when passing on a region name
            regPsf0 = startsEmptyExtendedPsf.getRegionalExtendedPsf(regionName=self.regions[j])
            regPsf1 = withDefaultExtendedPsf.getRegionalExtendedPsf(regionName=self.regions[j])
            self.assertMaskedImagesAlmostEqual(regPsf0, self.regionalEPsfs[j])
            self.assertMaskedImagesAlmostEqual(regPsf1, self.regionalEPsfs[j])
        # ensure we recover the original default PSF
        self.assertMaskedImagesAlmostEqual(withDefaultExtendedPsf(), self.defaultEPsf)

    def testIO(self):
        # Test IO with a constant extended PSF
        with tempfile.NamedTemporaryFile() as f:
            self.constantExtendedPsf.writeFits(f.name)
            readEPsf = extendedPsf.ExtendedPsf.readFits(f.name)
            self.assertMaskedImagesAlmostEqual(self.constantExtendedPsf(), readEPsf())
        # Test IO with per-region extended PSFs (with default)
        perRegionEPsf0 = extendedPsf.ExtendedPsf(self.defaultEPsf)
        for j in range(3):
            perRegionEPsf0.addRegionalExtendedPsf(self.regionalEPsfs[j], self.regions[j],
                                                  self.regionDetectors[j])
        with tempfile.NamedTemporaryFile() as f:
            perRegionEPsf0.writeFits(f.name)
            readEPsf0 = extendedPsf.ExtendedPsf.readFits(f.name)
            self.assertEqual(perRegionEPsf0.detectorsFocalPlaneRegions, readEPsf0.detectorsFocalPlaneRegions)
            # check default extended PSF
            self.assertMaskedImagesAlmostEqual(perRegionEPsf0(), readEPsf0())
            # and per-region extended PSFs
            for j in range(3):
                for det in self.regionDetectors[j]:
                    regPsf0, readRegPsf0 = perRegionEPsf0(det), readEPsf0(det)
                    self.assertMaskedImagesAlmostEqual(regPsf0, readRegPsf0)
        # Test IO with a single per-region extended PSF
        perRegionEPsf1 = extendedPsf.ExtendedPsf()
        perRegionEPsf1.addRegionalExtendedPsf(self.regionalEPsfs[1], self.regions[1], self.regionDetectors[1])
        with tempfile.NamedTemporaryFile() as f:
            perRegionEPsf1.writeFits(f.name)
            readEPsf1 = extendedPsf.ExtendedPsf.readFits(f.name)
            self.assertEqual(perRegionEPsf0.detectorsFocalPlaneRegions, readEPsf0.detectorsFocalPlaneRegions)
            for det in self.regionDetectors[1]:
                regPsf1, readRegPsf1 = perRegionEPsf1(det), readEPsf1(det)
                self.assertMaskedImagesAlmostEqual(regPsf1, readRegPsf1)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
