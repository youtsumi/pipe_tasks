#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import numpy as np

from lsstDebug import getDebugFrame
import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.daf.base as dafBase
import lsst.pipe.base.connectionTypes as cT
from lsst.afw.math import BackgroundList
from lsst.afw.table import SourceTable, SourceCatalog
from lsst.meas.algorithms import SubtractBackgroundTask, SourceDetectionTask, MeasureApCorrTask
from lsst.meas.algorithms.installGaussianPsf import InstallGaussianPsfTask
from lsst.meas.astrom import RefMatchTask, displayAstrometry
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask
from lsst.obs.base import ExposureIdInfo
from lsst.meas.base import SingleFrameMeasurementTask, ApplyApCorrTask, CatalogCalculationTask
from lsst.meas.deblender import SourceDeblendTask
from .measurePsf import MeasurePsfTask
from .repair import RepairTask
from .computeExposureSummaryStats import ComputeExposureSummaryStatsTask
from lsst.pex.exceptions import LengthError

__all__ = ["CharacterizeSpotsConfig", "CharacterizeSpotsTask"]


class CharacterizeSpotsConnections(pipeBase.PipelineTaskConnections,
                                   dimensions=("instrument", "exposure", "detector")):
    exposure = cT.Input(
        doc="Input exposure data",
        name="postISRCCD",
        storageClass="Exposure",
        dimensions=["instrument", "exposure", "detector"],
    )
    characterized = cT.Output(
        doc="Output characterized data.",
        name="spotExp",
        storageClass="ExposureF",
        dimensions=["instrument", "exposure", "detector"],
    )
    sourceCat = cT.Output(
        doc="Output source catalog.",
        name="spotSrc",
        storageClass="SourceCatalog",
        dimensions=["instrument", "exposure", "detector"],
    )
    backgroundModel = cT.Output(
        doc="Output background model.",
        name="spotExpBackground",
        storageClass="Background",
        dimensions=["instrument", "exposure", "detector"],
    )
    outputSchema = cT.InitOutput(
        doc="Schema of the catalog produced by CharacterizeSpots",
        name="spotSrc_schema",
        storageClass="SourceCatalog",
    )

#    def adjustQuantum(self, datasetRefMap: pipeBase.InputQuantizedConnection):
#        # Docstring inherited from PipelineTaskConnections
#        try:
#            return super().adjustQuantum(datasetRefMap)
#        except pipeBase.ScalarError as err:
#            raise pipeBase.ScalarError(
#                "CharacterizeSpotsTask can at present only be run on visits that are associated with "
#                "exactly one exposure.  Either this is not a valid exposure for this pipeline, or the "
#                "snap-combination step you probably want hasn't been configured to run between ISR and "
#                "this task (as of this writing, that would be because it hasn't been implemented yet)."
#            ) from err


class CharacterizeSpotsConfig(pipeBase.PipelineTaskConfig,
                              pipelineConnections=CharacterizeSpotsConnections):

    """!Config for CharacterizeSpotsTask"""
    doMeasurePsf = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Measure PSF? If False then for all subsequent operations use either existing PSF "
            "model when present, or install simple PSF model when not (see installSimplePsf "
            "config options)"
    )
    doWrite = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Persist results?",
    )
    doWriteExposure = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Write icExp and icExpBackground in addition to icSrc? Ignored if doWrite False.",
    )
    psfIterations = pexConfig.RangeField(
        dtype=int,
        default=2,
        min=1,
        doc="Number of iterations of detect sources, measure sources, "
            "estimate PSF. If useSimplePsf is True then 2 should be plenty; "
            "otherwise more may be wanted.",
    )
    background = pexConfig.ConfigurableField(
        target=SubtractBackgroundTask,
        doc="Configuration for initial background estimation",
    )
    detection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Detect sources"
    )
    measurement = pexConfig.ConfigurableField(
        target=SingleFrameMeasurementTask,
        doc="Measure sources"
    )
    catalogCalculation = pexConfig.ConfigurableField(
        target=CatalogCalculationTask,
        doc="Subtask to run catalogCalculation plugins on catalog"
    )
    useSimplePsf = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Replace the existing PSF model with a simplified version that has the same sigma "
        "at the start of each PSF determination iteration? Doing so makes PSF determination "
        "converge more robustly and quickly.",
    )
    installSimplePsf = pexConfig.ConfigurableField(
        target=InstallGaussianPsfTask,
        doc="Install a simple PSF model",
    )
    measurePsf = pexConfig.ConfigurableField(
        target=MeasurePsfTask,
        doc="Measure PSF",
    )
    repair = pexConfig.ConfigurableField(
        target=RepairTask,
        doc="Remove cosmic rays",
    )
    requireCrForPsf = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Require cosmic ray detection and masking to run successfully before measuring the PSF."
    )
    checkUnitsParseStrict = pexConfig.Field(
        doc="Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'",
        dtype=str,
        default="raise",
    )

    def setDefaults(self):
        super().setDefaults()
        # just detect bright stars; includeThresholdMultipler=10 seems large,
        # but these are the values we have been using
        self.detection.thresholdValue = 5.0
        self.detection.includeThresholdMultiplier = 10.0
        self.detection.doTempLocalBackground = False
        # minimal set of measurements needed to determine PSF
        self.measurement.plugins.names = [
            "base_PixelFlags",
            "base_SdssCentroid",
            "base_SdssShape",
            "base_GaussianFlux",
            "base_PsfFlux",
            "base_CircularApertureFlux",
        ]

## \addtogroup LSST_task_documentation
## \{
## \page CharacterizeSpotsTask
## \ref CharacterizeSpotsTask_ "CharacterizeSpotsTask"
## \copybrief CharacterizeSpotsTask
## \}


class CharacterizeSpotsTask(pipeBase.PipelineTask):

    # Example description used to live here, removed 2-20-2017 by MSSG

    ConfigClass = CharacterizeSpotsConfig
    _DefaultName = "characterizeSpots"

    def __init__(self, butler=None, schema=None, **kwargs):
        """!Construct a CharacterizeSpotsTask

        @param[in] butler  A butler object is passed to the refObjLoader constructor in case
            it is needed to load catalogs.  May be None if a catalog-based star selector is
            not used, if the reference object loader constructor does not require a butler,
            or if a reference object loader is passed directly via the refObjLoader argument.
        @param[in,out] schema  initial schema (an lsst.afw.table.SourceTable), or None
        @param[in,out] kwargs  other keyword arguments for lsst.pipe.base.CmdLineTask
        """
        super().__init__(**kwargs)

        if schema is None:
            schema = SourceTable.makeMinimalSchema()
        self.schema = schema
        self.makeSubtask("background")
        self.makeSubtask("installSimplePsf")
        self.makeSubtask("repair")
        self.makeSubtask("measurePsf", schema=self.schema)
        self.algMetadata = dafBase.PropertyList()
        self.makeSubtask('detection', schema=self.schema)
        self.makeSubtask('measurement', schema=self.schema, algMetadata=self.algMetadata)
        self.makeSubtask('catalogCalculation', schema=self.schema)
        self._initialFrame = getDebugFrame(self._display, "frame") or 1
        self._frame = self._initialFrame
        self.schema.checkUnits(parse_strict=self.config.checkUnitsParseStrict)
        self.outputSchema = afwTable.SourceCatalog(self.schema)

    def getInitOutputDatasets(self):
        outputCatSchema = afwTable.SourceCatalog(self.schema)
        outputCatSchema.getTable().setMetadata(self.algMetadata)
        return {'outputSchema': outputCatSchema}

    @pipeBase.timeMethod
    def run(self, exposure, exposureIdInfo=None, background=None):
        """!Characterize a science image

        Peforms the following operations:
        - Iterate the following config.psfIterations times, or once if config.doMeasurePsf false:
            - detect and measure sources and estimate PSF (see detectMeasureAndEstimatePsf for details)
        - interpolate over cosmic rays
        - perform final measurement

        @param[in,out] exposure  exposure to characterize (an lsst.afw.image.ExposureF or similar).
            The following changes are made:
            - update or set psf
            - update detection and cosmic ray mask planes
            - subtract background and interpolate over cosmic rays
        @param[in] exposureIdInfo  ID info for exposure (an lsst.obs.base.ExposureIdInfo).
            If not provided, returned SourceCatalog IDs will not be globally unique.
        @param[in,out] background  initial model of background already subtracted from exposure
            (an lsst.afw.math.BackgroundList). May be None if no background has been subtracted,
            which is typical for image characterization.

        @return pipe_base Struct containing these fields, all from the final iteration
        of detectMeasureAndEstimatePsf:
        - exposure: characterized exposure; image is repaired by interpolating over cosmic rays,
            mask is updated accordingly, and the PSF model is set
        - sourceCat: detected sources (an lsst.afw.table.SourceCatalog)
        - background: model of background subtracted from exposure (an lsst.afw.math.BackgroundList)
        - psfCellSet: spatial cells of PSF candidates (an lsst.afw.math.SpatialCellSet)
        """
        self._frame = self._initialFrame  # reset debug display frame

        if not self.config.doMeasurePsf and not exposure.hasPsf():
            self.log.warn("Source catalog detected and measured with placeholder or default PSF")
            self.installSimplePsf.run(exposure=exposure)

        if exposureIdInfo is None:
            exposureIdInfo = ExposureIdInfo()

        # subtract an initial estimate of background level
        background = self.background.run(exposure).background

        psfIterations = self.config.psfIterations if self.config.doMeasurePsf else 1
        for i in range(psfIterations):
            dmeRes = self.detectMeasureAndEstimatePsf(
                exposure=exposure,
                exposureIdInfo=exposureIdInfo,
                background=background,
            )

            psf = dmeRes.exposure.getPsf()
            psfSigma = psf.computeShape().getDeterminantRadius()
            psfDimensions = psf.computeImage().getDimensions()
            medBackground = np.median(dmeRes.background.getImage().getArray())
            self.log.info("iter %s; PSF sigma=%0.2f, dimensions=%s; median background=%0.2f" %
                          (i + 1, psfSigma, psfDimensions, medBackground))

        self.display("psf", exposure=dmeRes.exposure, sourceCat=dmeRes.sourceCat)

        # perform final repair with final PSF
        self.repair.run(exposure=dmeRes.exposure)
        self.display("repair", exposure=dmeRes.exposure, sourceCat=dmeRes.sourceCat)

        # perform final measurement with final PSF, including measuring and applying aperture correction,
        # if wanted
        self.measurement.run(measCat=dmeRes.sourceCat, exposure=dmeRes.exposure,
                             exposureId=exposureIdInfo.expId)
        self.catalogCalculation.run(dmeRes.sourceCat)

        self.display("measure", exposure=dmeRes.exposure, sourceCat=dmeRes.sourceCat)

        return pipeBase.Struct(
            exposure=dmeRes.exposure,
            sourceCat=dmeRes.sourceCat,
            background=dmeRes.background,
            psfCellSet=dmeRes.psfCellSet,

            characterized=dmeRes.exposure,
            backgroundModel=dmeRes.background
        )

    @pipeBase.timeMethod
    def detectMeasureAndEstimatePsf(self, exposure, exposureIdInfo, background):
        """!Perform one iteration of detect, measure and estimate PSF

        Performs the following operations:
        - if config.doMeasurePsf or not exposure.hasPsf():
            - install a simple PSF model (replacing the existing one, if need be)
        - interpolate over cosmic rays with keepCRs=True
        - estimate background and subtract it from the exposure
        - detect, deblend and measure sources, and subtract a refined background model;
        - if config.doMeasurePsf:
            - measure PSF

        @param[in,out] exposure  exposure to characterize (an lsst.afw.image.ExposureF or similar)
            The following changes are made:
            - update or set psf
            - update detection and cosmic ray mask planes
            - subtract background
        @param[in] exposureIdInfo  ID info for exposure (an lsst.obs_base.ExposureIdInfo)
        @param[in,out] background  initial model of background already subtracted from exposure
            (an lsst.afw.math.BackgroundList).

        @return pipe_base Struct containing these fields, all from the final iteration
        of detect sources, measure sources and estimate PSF:
        - exposure  characterized exposure; image is repaired by interpolating over cosmic rays,
            mask is updated accordingly, and the PSF model is set
        - sourceCat  detected sources (an lsst.afw.table.SourceCatalog)
        - background  model of background subtracted from exposure (an lsst.afw.math.BackgroundList)
        - psfCellSet  spatial cells of PSF candidates (an lsst.afw.math.SpatialCellSet)
        """
        # install a simple PSF model, if needed or wanted
        if not exposure.hasPsf() or (self.config.doMeasurePsf and self.config.useSimplePsf):
            self.log.warn("Source catalog detected and measured with placeholder or default PSF")
            self.installSimplePsf.run(exposure=exposure)

        # run repair, but do not interpolate over cosmic rays (do that elsewhere, with the final PSF model)
        if self.config.requireCrForPsf:
            self.repair.run(exposure=exposure, keepCRs=True)
        else:
            try:
                self.repair.run(exposure=exposure, keepCRs=True)
            except LengthError:
                self.log.warn("Skipping cosmic ray detection: Too many CR pixels (max %0.f)" %
                              self.config.repair.cosmicray.nCrPixelMax)

        self.display("repair_iter", exposure=exposure)

        if background is None:
            background = BackgroundList()

        sourceIdFactory = exposureIdInfo.makeSourceIdFactory()
        table = SourceTable.make(self.schema, sourceIdFactory)
        table.setMetadata(self.algMetadata)

        detRes = self.detection.run(table=table, exposure=exposure, doSmooth=True)
        sourceCat = detRes.sources
        if detRes.fpSets.background:
            for bg in detRes.fpSets.background:
                background.append(bg)

        self.measurement.run(measCat=sourceCat, exposure=exposure, exposureId=exposureIdInfo.expId)

        measPsfRes = pipeBase.Struct(cellSet=None)
        if self.config.doMeasurePsf:
            if self.measurePsf.usesMatches:
                matches = self.ref_match.loadAndMatch(exposure=exposure, sourceCat=sourceCat).matches
            else:
                matches = None
            measPsfRes = self.measurePsf.run(exposure=exposure, sources=sourceCat, matches=matches,
                                             expId=exposureIdInfo.expId)
        self.display("measure_iter", exposure=exposure, sourceCat=sourceCat)

        return pipeBase.Struct(
            exposure=exposure,
            sourceCat=sourceCat,
            background=background,
            psfCellSet=measPsfRes.cellSet,
        )

    def display(self, itemName, exposure, sourceCat=None):
        """Display exposure and sources on next frame, if display of itemName has been requested

        @param[in] itemName  name of item in debugInfo
        @param[in] exposure  exposure to display
        @param[in] sourceCat  source catalog to display
        """
        val = getDebugFrame(self._display, itemName)
        if not val:
            return

        displayAstrometry(exposure=exposure, sourceCat=sourceCat, frame=self._frame, pause=False)
        self._frame += 1
