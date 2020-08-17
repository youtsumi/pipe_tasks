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
#
"""Read preprocessed bright stars and stack them to build an extended
PSF model."""

__all__ = ["MeasureExtendedPsfTask"]

from dataclasses import dataclass
from typing import List

from lsst.pipe import base as pipeBase
from lsst.pipe.tasks.assembleCoadd import AssembleCoaddTask
import lsst.pex.config as pexConfig
from lsst.afw import math as afwMath
from lsst.afw import image as afwImage
from lsst.afw import fits as afwFits
from lsst.geom import Extent2I
from lsst.daf.base import PropertyList
from lsst.daf.butler import DeferredDatasetHandle


@dataclass
class FocalPlaneRegionExtendedPsf:
    """Single extended PSF over a focal plane region.

    The focal plane region is defined through a list
    of detectors.

    Parameters
    ----------
    extendedPsfImage : `lsst.afw.image.MaskedImageF`
        Image of the extended PSF model.
    detectorList : `list` [`int`]
        List of detector IDs that define the focal plane region over which this
        extended PSF model has been built (and can be used).
    """
    extendedPsfImage: afwImage.MaskedImageF
    detectorList: List[int]


class ExtendedPsf:
    """Extended PSF model.

    Each instance may contain a default extended PSF, a set of extended PSFs
    that correspond to different focal plane regions, or both. At this time,
    focal plane regions are always defined as a subset of detectors.

    Parameters
    ----------
    defaultExtendedPsf : `lsst.afw.image.MaskedImageF`
        Extended PSF model to be used as default (or only) extended PSF model.
    """
    def __init__(self, defaultExtendedPsf=None):
        self.defaultExtendedPsf = defaultExtendedPsf
        self.focalPlaneRegions = {}
        self.detectorsFocalPlaneRegions = {}

    def addRegionalExtendedPsf(self, extendedPsfImage, regionName, detectorList):
        """Add a new focal plane region, along wit hits extended PSF, to the
        ExtendedPsf instance.

        Parameters
        ----------
        extendedPsfImage : `lsst.afw.image.MaskedImageF`
            Extended PSF model for the region.
        regionName : `str`
            Name of the focal plane region. Will be converted to all-uppercase.
        detectorList : `list` [`int`]
            List of IDs for the detectors that define the focal plane region.
        """
        regionName = regionName.upper()
        if regionName in self.focalPlaneRegions:
            raise ValueError(f"Region name {regionName} is already used by this ExtendedPsf instance.")
        self.focalPlaneRegions[regionName] = FocalPlaneRegionExtendedPsf(extendedPsfImage=extendedPsfImage,
                                                                         detectorList=detectorList)
        for det in detectorList:
            self.detectorsFocalPlaneRegions[det] = regionName
        return None

    def __call__(self, detector=None):
        """Return the appropriate extended PSF.

        If the instance contains no extended PSF defined over focal plane
        regions, the default extended PSF will be returned regardless of
        whether a detector ID was passed as argument.

        Parameters
        ----------
        detector : `int`, optional
            Detector ID. If this instance contains extended PSFs defined over
            focal plane regions, the extended PSF model for the region that
            contains ``detector`` is returned. If not, the default extended PSF
            is returned.
        """
        if detector is None:
            if self.defaultExtendedPsf is None:
                raise ValueError("No default extended PSF available; please provide detector number.")
            return self.defaultExtendedPsf
        elif detector and not self.focalPlaneRegions:
            return self.defaultExtendedPsf
        return self.getRegionalExtendedPsf(detector=detector)

    def __len__(self):
        nRegions = len(self.focalPlaneRegions)
        if self.defaultExtendedPsf is not None:
            nRegions += 1
        return nRegions

    def getRegionalExtendedPsf(self, regionName=None, detector=None):
        """Returns the extended PSF for a focal plane region.

        The region can be identified either by name, or through a detector ID.

        Parameters
        ----------
        regionName : `str` or `None`, optional
            Name of the region for which the extended PSF should be retrieved.
            Ignored if  ``detector`` is provided. Must be provided if
            ``detector`` is None.
        detector : `int` or `None`, optional
            If provided, returns the extended PSF for the focal plane region
            that includes this detector.
        """
        if detector is None:
            if regionName is None:
                raise ValueError('One of either a regionName or a detector number must be provided.')
            return self.focalPlaneRegions[regionName].extendedPsfImage
        return self.focalPlaneRegions[self.detectorsFocalPlaneRegions[detector]].extendedPsfImage

    def writeFits(self, filename):
        """Write this object to a file.

        Parameters
        ----------
        filename : `str`
            Name of file to write.
        """
        # create primary HDU with global metadata
        metadata = PropertyList()
        if self.defaultExtendedPsf is not None:
            metadata["HAS_DEFAULT"] = True
        else:
            metadata["HAS_DEFAULT"] = False
        if self.focalPlaneRegions:
            metadata["HAS_REGIONS"] = True
            metadata["REGION_NAMES"] = list(self.focalPlaneRegions.keys())
            for region, ePsfRegion in self.focalPlaneRegions.items():
                metadata[region] = ePsfRegion.detectorList
        else:
            metadata["HAS_REGIONS"] = False
        fitsPrimary = afwFits.Fits(filename, "w")
        fitsPrimary.createEmpty()
        fitsPrimary.writeMetadata(metadata)
        fitsPrimary.closeFile()
        # write default extended PSF
        if self.defaultExtendedPsf is not None:
            defaultHduMetadata = PropertyList()
            defaultHduMetadata.update({"REGION": "DEFAULT", "EXTNAME": "IMAGE"})
            self.defaultExtendedPsf.image.writeFits(filename, metadata=defaultHduMetadata, mode='a')
            defaultHduMetadata.update({"REGION": "DEFAULT", "EXTNAME": "MASK"})
            self.defaultExtendedPsf.mask.writeFits(filename, metadata=defaultHduMetadata, mode='a')
        # write extended PSF for each focal plane region
        for j, (region, ePsfRegion) in enumerate(self.focalPlaneRegions.items()):
            metadata = PropertyList()
            metadata.update({"REGION": region, "EXTNAME": "IMAGE"})
            ePsfRegion.extendedPsfImage.image.writeFits(filename, metadata=metadata, mode='a')
            metadata.update({"REGION": region, "EXTNAME": "MASK"})
            ePsfRegion.extendedPsfImage.mask.writeFits(filename, metadata=metadata, mode='a')
        return None

    @classmethod
    def readFits(cls, filename):
        """Build an instance of this class from a file.

        Parameters
        ----------
        filename : `str`
            Name of the file to read.
        """
        # extract info from metadata
        globalMetadata = afwFits.readMetadata(filename, hdu=0)
        hasDefault = globalMetadata.getBool("HAS_DEFAULT")
        if globalMetadata.getBool("HAS_REGIONS"):
            focalPlaneRegionNames = globalMetadata.getArray("REGION_NAMES")
        else:
            focalPlaneRegionNames = []
        f = afwFits.Fits(filename, 'r')
        nExtensions = f.countHdus()
        extendedPsfParts = {}
        for j in range(1, nExtensions):
            md = afwFits.readMetadata(filename, hdu=j)
            if hasDefault and md["REGION"] == "DEFAULT":
                if md['EXTNAME'] == 'IMAGE':
                    defaultImage = afwImage.ImageF(filename, hdu=j)
                elif md['EXTNAME'] == 'MASK':
                    defaultMask = afwImage.MaskX(filename, hdu=j)
                continue
            if md['EXTNAME'] == 'IMAGE':
                extendedPsfPart = afwImage.ImageF(filename, hdu=j)
            elif md['EXTNAME'] == 'MASK':
                extendedPsfPart = afwImage.MaskX(filename, hdu=j)
            extendedPsfParts.setdefault(md['REGION'], {})[md['EXTNAME'].lower()] = extendedPsfPart
        # handle default if present
        if hasDefault:
            extendedPsf = cls(afwImage.MaskedImageF(defaultImage, defaultMask))
        else:
            extendedPsf = cls()
        # ensure we recovered an extended PSF for all focal plane regions
        if len(extendedPsfParts) != len(focalPlaneRegionNames):
            raise ValueError(f'Number of per-region extended PSFs read ({len(extendedPsfParts)}) does not '
                             'match with the number of regions recorded in the metadata '
                             f'({len(focalPlaneRegionNames)}).')
        # generate extended PSF regions mappings
        for rName in focalPlaneRegionNames:
            extendedPsfImage = afwImage.MaskedImageF(**extendedPsfParts[rName])
            detectorList = globalMetadata.getArray(rName)
            extendedPsf.addRegionalExtendedPsf(extendedPsfImage, rName, detectorList)
        # instantiate ExtendedPsf
        return extendedPsf


class StackBrightStarsConfig(pexConfig.Config):
    """Configuration parameters for MeasureExtendedPsfTask.
    """
    subregionSize = pexConfig.ListField(
        dtype=int,
        doc="Size, in pixels, of the subregions over which the stacking will be "
            "iteratively performed.",
        default=(100, 100)
    )
    stackingStatistic = pexConfig.ChoiceField(
        dtype=str,
        doc="Type of statistic to use for stacking.",
        default="MEANCLIP",
        allowed={
            "MEAN": "mean",
            "MEDIAN": "median",
            "MEANCLIP": "clipped mean",
        }
    )
    numSigmaClip = pexConfig.Field(
        dtype=float,
        doc="Sigma for outlier rejection; ignored if stackingStatistic != 'MEANCLIP'.",
        default=4
    )
    numIter = pexConfig.Field(
        dtype=int,
        doc="Number of iterations of outlier rejection; ignored if atackingStatistic != 'MEANCLIP'.",
        default=3
    )
    badMaskPlanes = pexConfig.ListField(
        dtype=str,
        doc="Mask planes that, if set, lead to associated pixels not being included in the stacking of the "
            "bright star stamps.",
        default=('BAD', 'CR', 'CROSSTALK', 'EDGE', 'NO_DATA', 'SAT', 'SUSPECT', 'UNMASKEDNAN')
    )
    doMagCut = pexConfig.Field(
        dtype=bool,
        doc="Apply magnitude cut before stacking?",
        default=False
    )
    magLimit = pexConfig.Field(
        dtype=float,
        doc="Magnitude limit, in Gaia G; all stars brighter than this value will be stacked",
        default=18
    )


class StackBrightStarsTask(pipeBase.CmdLineTask):
    """Stack bright stars together to build an extended PSF model.
    """
    ConfigClass = StackBrightStarsConfig
    _DefaultName = "stackBrightStars"

    def __init__(self, initInputs=None, *args, **kwargs):
        pipeBase.CmdLineTask.__init__(self, *args, **kwargs)

    def _setUpStacking(self, exampleStamp):
        """Configure stacking statistic and control from config fields.
        """
        statsControl = afwMath.StatisticsControl()
        statsControl.setNumSigmaClip(self.config.numSigmaClip)
        statsControl.setNumIter(self.config.numIter)
        badMasks = self.config.badMaskPlanes
        if badMasks:
            andMask = exampleStamp.mask.getPlaneBitMask(badMasks[0])
            for bm in badMasks[1:]:
                andMask = andMask | exampleStamp.mask.getPlaneBitMask(bm)
            statsControl.setAndMask(andMask)
        statsFlags = afwMath.stringToStatisticsProperty(self.config.stackingStatistic)
        return statsControl, statsFlags

    def run(self, bssRefList, regionName=None):
        """Read input bright star stamps and stack them together.

        The stacking is done iteratively over smaller areas of the final model
        image to allow for a great number of bright star stamps to be used.

        Parameters
        ----------
        bssRefList : `list` of
                `lsst.daf.butler._deferredDatasetHandle.DeferredDatasetHandle`
            List of available bright star stamps data references.
        regionName : `str`, optional
            Name of the focal plane region, if applicable. Only used for
            logging purposes, when running over multiple such regions
            (typically from `MeasureExtendedPsfTask`)
        """
        logMessage = f'Building extended PSF from stamps extracted from {len(bssRefList)} detector images '
        if regionName:
            logMessage += f'for region "{regionName}".'
        self.log.info(logMessage)
        # read in example set of full stamps
        for bssRef in bssRefList:
            if not isinstance(bssRef, DeferredDatasetHandle):
                if not bssRef.datasetExists("brightStarStamps"):
                    self.log.warn("Could not find %s %s; skipping it", "brightStarStamps", bssRef.dataId)
                    bssRefList.remove(bssRef)
                    continue
            bss = bssRef.get(datasetType="brightStarStamps", immediate=True)
            break
        exampleStamp = bss[0].stamp_im
        # create model image
        extPsf = afwImage.MaskedImageF(exampleStamp.getBBox())
        # divide model image into smaller subregions
        subregionSize = Extent2I(*self.config.subregionSize)
        subBBoxes = AssembleCoaddTask._subBBoxIter(extPsf.getBBox(), subregionSize)
        # compute approximate number of subregions
        nSubregions = int(extPsf.getDimensions()[0]/subregionSize[0] + 1)*int(
            extPsf.getDimensions()[1]/subregionSize[1] + 1)
        self.log.info(f"Stacking will performed iteratively over approximately {nSubregions} "
                      "smaller areas of the final model image.")
        # set up stacking statistic
        statsControl, statsFlags = self._setUpStacking(exampleStamp)
        # perform stacking
        for jbbox, bbox in enumerate(subBBoxes):
            allStars = None
            for bssRef in bssRefList:
                if not isinstance(bssRef, DeferredDatasetHandle):
                    if not bssRef.datasetExists("brightStarStamps"):
                        self.log.warn("Could not find %s %s; skipping it", "brightStarStamps", bssRef.dataId)
                        bssRefList.remove(bssRef)
                        continue
                readStars = bssRef.get(datasetType="brightStarStamps", parameters={'bbox': bbox})
                if self.config.doMagCut:
                    readStars = readStars.selectByMag(magMax=self.config.magLimit)
                if allStars:
                    allStars.extend(readStars)
                else:
                    allStars = readStars
            # TODO: DM-27371 add weights to bright stars for stacking
            coaddSubregion = afwMath.statisticsStack(allStars.getMaskedImages(), statsFlags, statsControl)
            extPsf.assign(coaddSubregion, bbox)
        return extPsf


class MeasureExtendedPsfConnections(pipeBase.PipelineTaskConnections,
                                    dimensions=("band", "instrument")):
    inputBrightStarStamps = pipeBase.connectionTypes.Input(
        doc="Input list of bright star collections to be stacked.",
        name="brightStarStamps",
        storageClass="BrightStarStamps",
        dimensions=("visit", "detector"),
        deferLoad=True,
        multiple=True
    )
    extendedPsf = pipeBase.connectionTypes.Output(
        doc="Extended PSF model built by stacking bright stars.",
        name="extendedPsf",
        storageClass="ExtendedPsf",
        dimensions=("band",),
    )


class MeasureExtendedPsfConfig(pipeBase.PipelineTaskConfig,
                               pipelineConnections=MeasureExtendedPsfConnections):
    """Configuration parameters for MeasureExtendedPsfTask.
    """
    stackBrightStars = pexConfig.ConfigurableField(
        target=StackBrightStarsTask,
        doc="Stack selected bright stars",
    )
    detectorsFocalPlaneRegions = pexConfig.DictField(
        keytype=int,
        itemtype=str,
        doc="Mapping from detector IDs to focal plane region names. If empty, a constant "
            "extended PSF model is built from all selected bright stars.",
        default={}
    )


class MeasureExtendedPsfTask(pipeBase.CmdLineTask):
    """Build and save extended PSF model.

    The model is built by stacking bright star stamps, extracted and
    preprocessed by
    `lsst.pipe.tasks.processBrightStars.ProcessBrightStarsTask`.
    If a mapping from detector IDs to focal plane regions is provided,
    a different extended PSF model will be built for each focal plane
    region. If not, a single, constant extended PSF model is built using
    all available data.
    """
    ConfigClass = MeasureExtendedPsfConfig
    _DefaultName = "measureExtendedPsf"

    def __init__(self, initInputs=None, *args, **kwargs):
        pipeBase.CmdLineTask.__init__(self, *args, **kwargs)
        self.makeSubtask("stackBrightStars")
        self.focalPlaneRegions = {region: [] for region in
                                  set(self.config.detectorsFocalPlaneRegions.values())}
        for det, region in self.config.detectorsFocalPlaneRegions.items():
            self.focalPlaneRegions[region].append(det)
        # make no assumption on what detector IDs should be, but if we come
        # across one where there are processed bright stars, but no
        # corresponding focal plane region, make sure we keep track of
        # it (eg to raise a warning only once)
        self.regionlessDetectors = []

    def selectDetectorRefs(self, refList):
        """Split available sets of bright star stamps according to focal plane
        regions.

        Parameters
        ----------
        refList : `list` of
                `lsst.daf.butler._deferredDatasetHandle.DeferredDatasetHandle`
            List of available bright star stamps data references.
        """
        regionRefList = {region: [] for region in self.focalPlaneRegions.keys()}
        for deferredDatasetHandle in refList:
            detectorId = deferredDatasetHandle.ref.dataId["detector"]
            if detectorId in self.regionlessDetectors:
                continue
            try:
                regionName = self.config.detectorsFocalPlaneRegions[detectorId]
            except KeyError:
                self.log.warn(f'Bright stars were available for detector {detectorId}, but it was missing '
                              'from the "detectorsFocalPlaneRegions" config field, so they will not '
                              'be used to build any of the extended PSF models')
                self.regionlessDetectors.append(detectorId)
                continue
            regionRefList[regionName].append(deferredDatasetHandle)
        return regionRefList

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputData = butlerQC.get(inputRefs)
        bssRefList = inputData['inputBrightStarStamps']
        # Handle default cause of a single region with empty detector list
        if not self.config.detectorsFocalPlaneRegions:
            self.log.info("No detector groups were provided to MeasureExtendedPsfTask; computing a single, "
                          "constant extended PSF model over all available observations.")
            outputEPsf = ExtendedPsf(self.stackBrightStars.run(bssRefList))
        else:
            outputEPsf = ExtendedPsf()
            regionRefList = self.selectDetectorRefs(bssRefList)
            for regionName, refList in regionRefList.items():
                if not refList:
                    # no valid references found
                    self.log.warn(f'No valid brightStarStamps reference found for region "{regionName}"; '
                                  'skipping it.')
                    continue
                extPsf = self.stackBrightStars.run(refList, regionName)
                outputEPsf.addRegionalExtendedPsf(extPsf, regionName, self.focalPlaneRegions[regionName])
        output = pipeBase.Struct(extendedPsf=outputEPsf)
        butlerQC.put(output, outputRefs)
