#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
import math

from lsstDebug import getDebugFrame
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.afw.table as afwTable
from lsst.meas.astrom import AstrometryTask, displayAstrometry, denormalizeMatches
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask, SkyObjectsTask
from lsst.obs.base import ExposureIdInfo
import lsst.daf.base as dafBase
from lsst.afw.math import BackgroundList
from lsst.afw.table import SourceTable
from lsst.meas.algorithms import SourceDetectionTask, ReferenceObjectLoader
from lsst.meas.base import (SingleFrameMeasurementTask,
                            ApplyApCorrTask,
                            CatalogCalculationTask)
from lsst.meas.deblender import SourceDeblendTask
from lsst.pipe.tasks.setPrimaryFlags import SetPrimaryFlagsTask
from .fakes import BaseFakeSourcesTask
from .photoCal import PhotoCalTask

__all__ = ["CalibrateConfig", "CalibrateTask"]


class CalibrateConnections(pipeBase.PipelineTaskConnections, dimensions=("instrument", "exposure", "detector"),
                           defaultTemplates={}):

    icSourceSchema = cT.InitInput(
        doc="Schema produced by characterize image task, used to initialize this task",
        name="icSrc_schema",
        storageClass="SourceCatalog",
    )

    outputSchema = cT.InitOutput(
        doc="Schema after CalibrateTask has been initialized",
        name="src_schema",
        storageClass="SourceCatalog",
    )

    exposure = cT.Input(
        doc="Input image to calibrate",
        name="icExp",
        storageClass="ExposureF",
        dimensions=("instrument", "exposure", "detector"),
    )

    background = cT.Input(
        doc="Backgrounds determined by characterize task",
        name="icExpBackground",
        storageClass="Background",
        dimensions=("instrument", "exposure", "detector"),
    )

    icSourceCat = cT.Input(
        doc="Source catalog created by characterize task",
        name="icSrc",
        storageClass="SourceCatalog",
        dimensions=("instrument", "exposure", "detector"),
    )

    astromRefCat = cT.PrerequisiteInput(
        doc="Reference catalog to use for astrometry",
        name="cal_ref_cat",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )

    photoRefCat = cT.PrerequisiteInput(
        doc="Reference catalog to use for photometric calibration",
        name="cal_ref_cat",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True
    )

    outputExposure = cT.Output(
        doc="Exposure after running calibration task",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "exposure", "detector"),
    )

    outputCat = cT.Output(
        doc="Source catalog produced in calibrate task",
        name="src",
        storageClass="SourceCatalog",
        dimensions=("instrument", "exposure", "detector"),
    )

    outputBackground = cT.Output(
        doc="Background models estimated in calibration task",
        name="calexpBackground",
        storageClass="Background",
        dimensions=("instrument", "exposure", "detector"),
    )

    matches = cT.Output(
        doc="Source/refObj matches from the astrometry solver",
        name="srcMatch",
        storageClass="Catalog",
        dimensions=("instrument", "exposure", "detector"),
    )

    matchesDenormalized = cT.Output(
        doc="Denormalized matches from astrometry solver",
        name="srcMatchFull",
        storageClass="Catalog",
        dimensions=("instrument", "exposure", "detector"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if config.doAstrometry is False:
            self.prerequisiteInputs.remove("astromRefCat")
        if config.doPhotoCal is False:
            self.prerequisiteInputs.remove("photoRefCat")

        if config.doWriteMatches is False or config.doAstrometry is False:
            self.outputs.remove("matches")
        if config.doWriteMatchesDenormalized is False or config.doAstrometry is False:
            self.outputs.remove("matchesDenormalized")


class CalibrateConfig(pipeBase.PipelineTaskConfig, pipelineConnections=CalibrateConnections):
    """Config for CalibrateTask"""
    doWrite = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Save calibration results?",
    )
    doWriteHeavyFootprintsInSources = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Include HeavyFootprint data in source table? If false then heavy "
            "footprints are saved as normal footprints, which saves some space"
    )
    doWriteMatches = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Write reference matches (ignored if doWrite or doAstrometry false)?",
    )
    doWriteMatchesDenormalized = pexConfig.Field(
        dtype=bool,
        default=False,
        doc=("Write reference matches in denormalized format? "
             "This format uses more disk space, but is more convenient to "
             "read. Ignored if doWriteMatches=False or doWrite=False."),
    )
    doAstrometry = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Perform astrometric calibration?",
    )
    astromRefObjLoader = pexConfig.ConfigurableField(
        target=LoadIndexedReferenceObjectsTask,
        doc="reference object loader for astrometric calibration",
    )
    photoRefObjLoader = pexConfig.ConfigurableField(
        target=LoadIndexedReferenceObjectsTask,
        doc="reference object loader for photometric calibration",
    )
    astrometry = pexConfig.ConfigurableField(
        target=AstrometryTask,
        doc="Perform astrometric calibration to refine the WCS",
    )
    requireAstrometry = pexConfig.Field(
        dtype=bool,
        default=True,
        doc=("Raise an exception if astrometry fails? Ignored if doAstrometry "
             "false."),
    )
    doPhotoCal = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Perform phometric calibration?",
    )
    requirePhotoCal = pexConfig.Field(
        dtype=bool,
        default=True,
        doc=("Raise an exception if photoCal fails? Ignored if doPhotoCal "
             "false."),
    )
    photoCal = pexConfig.ConfigurableField(
        target=PhotoCalTask,
        doc="Perform photometric calibration",
    )
    icSourceFieldsToCopy = pexConfig.ListField(
        dtype=str,
        default=("calib_psf_candidate", "calib_psf_used", "calib_psf_reserved"),
        doc=("Fields to copy from the icSource catalog to the output catalog "
             "for matching sources Any missing fields will trigger a "
             "RuntimeError exception. Ignored if icSourceCat is not provided.")
    )
    matchRadiusPix = pexConfig.Field(
        dtype=float,
        default=3,
        doc=("Match radius for matching icSourceCat objects to sourceCat "
             "objects (pixels)"),
    )
    checkUnitsParseStrict = pexConfig.Field(
        doc=("Strictness of Astropy unit compatibility check, can be 'raise', "
             "'warn' or 'silent'"),
        dtype=str,
        default="raise",
    )
    detection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Detect sources"
    )
    doDeblend = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Run deblender input exposure"
    )
    deblend = pexConfig.ConfigurableField(
        target=SourceDeblendTask,
        doc="Split blended sources into their components"
    )
    doSkySources = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Generate sky sources?",
    )
    skySources = pexConfig.ConfigurableField(
        target=SkyObjectsTask,
        doc="Generate sky sources",
    )
    measurement = pexConfig.ConfigurableField(
        target=SingleFrameMeasurementTask,
        doc="Measure sources"
    )
    setPrimaryFlags = pexConfig.ConfigurableField(
        target=SetPrimaryFlagsTask,
        doc=("Set flags for primary source classification in single frame "
             "processing. True if sources are not sky sources and not a parent.")
    )
    doApCorr = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Run subtask to apply aperture correction"
    )
    applyApCorr = pexConfig.ConfigurableField(
        target=ApplyApCorrTask,
        doc="Subtask to apply aperture corrections"
    )
    # If doApCorr is False, and the exposure does not have apcorrections
    # already applied, the active plugins in catalogCalculation almost
    # certainly should not contain the characterization plugin
    catalogCalculation = pexConfig.ConfigurableField(
        target=CatalogCalculationTask,
        doc="Subtask to run catalogCalculation plugins on catalog"
    )
    doInsertFakes = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Run fake sources injection task"
    )
    insertFakes = pexConfig.ConfigurableField(
        target=BaseFakeSourcesTask,
        doc="Injection of fake sources for testing purposes (must be "
            "retargeted)"
    )
    doWriteExposure = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Write the calexp? If fakes have been added then we do not want to write out the calexp as a "
            "normal calexp but as a fakes_calexp."
    )

    def setDefaults(self):
        super().setDefaults()
        self.detection.doTempLocalBackground = False
        self.deblend.maxFootprintSize = 2000
        self.measurement.plugins.names |= ["base_LocalPhotoCalib", "base_LocalWcs"]

    def validate(self):
        super().validate()
        astromRefCatGen2 = getattr(self.astromRefObjLoader, "ref_dataset_name", None)
        if astromRefCatGen2 is not None and astromRefCatGen2 != self.connections.astromRefCat:
            raise ValueError(
                f"Gen2 ({astromRefCatGen2}) and Gen3 ({self.connections.astromRefCat}) astrometry reference "
                f"catalogs are different.  These options must be kept in sync until Gen2 is retired."
            )
        photoRefCatGen2 = getattr(self.photoRefObjLoader, "ref_dataset_name", None)
        if photoRefCatGen2 is not None and photoRefCatGen2 != self.connections.photoRefCat:
            raise ValueError(
                f"Gen2 ({photoRefCatGen2}) and Gen3 ({self.connections.photoRefCat}) photometry reference "
                f"catalogs are different.  These options must be kept in sync until Gen2 is retired."
            )


## \addtogroup LSST_task_documentation
## \{
## \page CalibrateTask
## \ref CalibrateTask_ "CalibrateTask"
## \copybrief CalibrateTask
## \}

class CalibrateTask(pipeBase.PipelineTask, pipeBase.CmdLineTask):
    r"""!Calibrate an exposure: measure sources and perform astrometric and
        photometric calibration

    @anchor CalibrateTask_

    @section pipe_tasks_calibrate_Contents  Contents

     - @ref pipe_tasks_calibrate_Purpose
     - @ref pipe_tasks_calibrate_Initialize
     - @ref pipe_tasks_calibrate_IO
     - @ref pipe_tasks_calibrate_Config
     - @ref pipe_tasks_calibrate_Metadata
     - @ref pipe_tasks_calibrate_Debug


    @section pipe_tasks_calibrate_Purpose  Description

    Given an exposure with a good PSF model and aperture correction map
    (e.g. as provided by @ref CharacterizeImageTask), perform the following
     operations:
    - Run detection and measurement
    - Run astrometry subtask to fit an improved WCS
    - Run photoCal subtask to fit the exposure's photometric zero-point

    @section pipe_tasks_calibrate_Initialize  Task initialisation

    @copydoc \_\_init\_\_

    @section pipe_tasks_calibrate_IO  Invoking the Task

    If you want this task to unpersist inputs or persist outputs, then call
    the `runDataRef` method (a wrapper around the `run` method).

    If you already have the inputs unpersisted and do not want to persist the
    output then it is more direct to call the `run` method:

    @section pipe_tasks_calibrate_Config  Configuration parameters

    See @ref CalibrateConfig

    @section pipe_tasks_calibrate_Metadata  Quantities set in exposure Metadata

    Exposure metadata
    <dl>
        <dt>MAGZERO_RMS  <dd>MAGZERO's RMS == sigma reported by photoCal task
        <dt>MAGZERO_NOBJ <dd>Number of stars used == ngood reported by photoCal
                             task
        <dt>COLORTERM1   <dd>?? (always 0.0)
        <dt>COLORTERM2   <dd>?? (always 0.0)
        <dt>COLORTERM3   <dd>?? (always 0.0)
    </dl>

    @section pipe_tasks_calibrate_Debug  Debug variables

    The @link lsst.pipe.base.cmdLineTask.CmdLineTask command line task@endlink
    interface supports a flag
    `--debug` to import `debug.py` from your `$PYTHONPATH`; see @ref baseDebug
    for more about `debug.py`.

    CalibrateTask has a debug dictionary containing one key:
    <dl>
    <dt>calibrate
    <dd>frame (an int; <= 0 to not display) in which to display the exposure,
        sources and matches. See @ref lsst.meas.astrom.displayAstrometry for
        the meaning of the various symbols.
    </dl>

    For example, put something like:
    @code{.py}
        import lsstDebug
        def DebugInfo(name):
            di = lsstDebug.getInfo(name)  # N.b. lsstDebug.Info(name) would
                                          # call us recursively
            if name == "lsst.pipe.tasks.calibrate":
                di.display = dict(
                    calibrate = 1,
                )

            return di

        lsstDebug.Info = DebugInfo
    @endcode
    into your `debug.py` file and run `calibrateTask.py` with the `--debug`
    flag.

    Some subtasks may have their own debug variables; see individual Task
    documentation.
    """

    # Example description used to live here, removed 2-20-2017 as per
    # https://jira.lsstcorp.org/browse/DM-9520

    ConfigClass = CalibrateConfig
    _DefaultName = "calibrate"
    RunnerClass = pipeBase.ButlerInitializedTaskRunner

    def __init__(self, butler=None, astromRefObjLoader=None,
                 photoRefObjLoader=None, icSourceSchema=None,
                 initInputs=None, **kwargs):
        """!Construct a CalibrateTask

        @param[in] butler  The butler is passed to the refObjLoader constructor
            in case it is needed.  Ignored if the refObjLoader argument
            provides a loader directly.
        @param[in] astromRefObjLoader  An instance of LoadReferenceObjectsTasks
            that supplies an external reference catalog for astrometric
            calibration.  May be None if the desired loader can be constructed
            from the butler argument or all steps requiring a reference catalog
            are disabled.
        @param[in] photoRefObjLoader  An instance of LoadReferenceObjectsTasks
            that supplies an external reference catalog for photometric
            calibration.  May be None if the desired loader can be constructed
            from the butler argument or all steps requiring a reference catalog
            are disabled.
        @param[in] icSourceSchema  schema for icSource catalog, or None.
            Schema values specified in config.icSourceFieldsToCopy will be
            taken from this schema. If set to None, no values will be
            propagated from the icSourceCatalog
        @param[in,out] kwargs  other keyword arguments for
            lsst.pipe.base.CmdLineTask
        """
        super().__init__(**kwargs)

        if icSourceSchema is None and butler is not None:
            # Use butler to read icSourceSchema from disk.
            icSourceSchema = butler.get("icSrc_schema", immediate=True).schema

        if icSourceSchema is None and butler is None and initInputs is not None:
            icSourceSchema = initInputs['icSourceSchema'].schema

        if icSourceSchema is not None:
            # use a schema mapper to avoid copying each field separately
            self.schemaMapper = afwTable.SchemaMapper(icSourceSchema)
            minimumSchema = afwTable.SourceTable.makeMinimalSchema()
            self.schemaMapper.addMinimalSchema(minimumSchema, False)

            # Add fields to copy from an icSource catalog
            # and a field to indicate that the source matched a source in that
            # catalog. If any fields are missing then raise an exception, but
            # first find all missing fields in order to make the error message
            # more useful.
            self.calibSourceKey = self.schemaMapper.addOutputField(
                afwTable.Field["Flag"]("calib_detected",
                                       "Source was detected as an icSource"))
            missingFieldNames = []
            for fieldName in self.config.icSourceFieldsToCopy:
                try:
                    schemaItem = icSourceSchema.find(fieldName)
                except Exception:
                    missingFieldNames.append(fieldName)
                else:
                    # field found; if addMapping fails then raise an exception
                    self.schemaMapper.addMapping(schemaItem.getKey())

            if missingFieldNames:
                raise RuntimeError("isSourceCat is missing fields {} "
                                   "specified in icSourceFieldsToCopy"
                                   .format(missingFieldNames))

            # produce a temporary schema to pass to the subtasks; finalize it
            # later
            self.schema = self.schemaMapper.editOutputSchema()
        else:
            self.schemaMapper = None
            self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.makeSubtask('detection', schema=self.schema)

        self.algMetadata = dafBase.PropertyList()

        # Only create a subtask for fakes if configuration option is set
        # N.B. the config for fake object task must be retargeted to a child
        # of BaseFakeSourcesTask
        if self.config.doInsertFakes:
            self.makeSubtask("insertFakes")

        if self.config.doDeblend:
            self.makeSubtask("deblend", schema=self.schema)
        if self.config.doSkySources:
            self.makeSubtask("skySources")
            self.skySourceKey = self.schema.addField("sky_source", type="Flag", doc="Sky objects.")
        self.makeSubtask('measurement', schema=self.schema,
                         algMetadata=self.algMetadata)
        self.makeSubtask("setPrimaryFlags", schema=self.schema, isSingleFrame=True)
        if self.config.doApCorr:
            self.makeSubtask('applyApCorr', schema=self.schema)
        self.makeSubtask('catalogCalculation', schema=self.schema)

        if self.config.doAstrometry:
            if astromRefObjLoader is None and butler is not None:
                self.makeSubtask('astromRefObjLoader', butler=butler)
                astromRefObjLoader = self.astromRefObjLoader
            self.makeSubtask("astrometry", refObjLoader=astromRefObjLoader,
                             schema=self.schema)
        if self.config.doPhotoCal:
            if photoRefObjLoader is None and butler is not None:
                self.makeSubtask('photoRefObjLoader', butler=butler)
                photoRefObjLoader = self.photoRefObjLoader
            self.makeSubtask("photoCal", refObjLoader=photoRefObjLoader,
                             schema=self.schema)

        if initInputs is not None and (astromRefObjLoader is not None or photoRefObjLoader is not None):
            raise RuntimeError("PipelineTask form of this task should not be initialized with "
                               "reference object loaders.")

        if self.schemaMapper is not None:
            # finalize the schema
            self.schema = self.schemaMapper.getOutputSchema()
        self.schema.checkUnits(parse_strict=self.config.checkUnitsParseStrict)

        sourceCatSchema = afwTable.SourceCatalog(self.schema)
        sourceCatSchema.getTable().setMetadata(self.algMetadata)
        self.outputSchema = sourceCatSchema

    @pipeBase.timeMethod
    def runDataRef(self, dataRef, exposure=None, background=None, icSourceCat=None,
                   doUnpersist=True):
        """!Calibrate an exposure, optionally unpersisting inputs and
            persisting outputs.

        This is a wrapper around the `run` method that unpersists inputs
        (if `doUnpersist` true) and persists outputs (if `config.doWrite` true)

        @param[in] dataRef  butler data reference corresponding to a science
            image
        @param[in,out] exposure  characterized exposure (an
            lsst.afw.image.ExposureF or similar), or None to unpersist existing
            icExp and icBackground. See `run` method for details of what is
            read and written.
        @param[in,out] background  initial model of background already
            subtracted from exposure (an lsst.afw.math.BackgroundList). May be
            None if no background has been subtracted, though that is unusual
            for calibration. A refined background model is output. Ignored if
            exposure is None.
        @param[in] icSourceCat  catalog from which to copy the fields specified
            by icSourceKeys, or None;
        @param[in] doUnpersist  unpersist data:
            - if True, exposure, background and icSourceCat are read from
              dataRef and those three arguments must all be None;
            - if False the exposure must be provided; background and
              icSourceCat are optional. True is intended for running as a
              command-line task, False for running as a subtask
        @return same data as the calibrate method
        """
        self.log.info("Processing %s" % (dataRef.dataId))

        if doUnpersist:
            if any(item is not None for item in (exposure, background,
                                                 icSourceCat)):
                raise RuntimeError("doUnpersist true; exposure, background "
                                   "and icSourceCat must all be None")
            exposure = dataRef.get("icExp", immediate=True)
            background = dataRef.get("icExpBackground", immediate=True)
            icSourceCat = dataRef.get("icSrc", immediate=True)
        elif exposure is None:
            raise RuntimeError("doUnpersist false; exposure must be provided")

        exposureIdInfo = dataRef.get("expIdInfo")

        calRes = self.run(
            exposure=exposure,
            exposureIdInfo=exposureIdInfo,
            background=background,
            icSourceCat=icSourceCat,
        )

        if self.config.doWrite:
            self.writeOutputs(
                dataRef=dataRef,
                exposure=calRes.exposure,
                background=calRes.background,
                sourceCat=calRes.sourceCat,
                astromMatches=calRes.astromMatches,
                matchMeta=calRes.matchMeta,
            )

        return calRes

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs['exposureIdInfo'] = ExposureIdInfo.fromDataId(butlerQC.quantum.dataId, "exposure_detector")

        if self.config.doAstrometry:
            refObjLoader = ReferenceObjectLoader(dataIds=[ref.datasetRef.dataId
                                                          for ref in inputRefs.astromRefCat],
                                                 refCats=inputs.pop('astromRefCat'),
                                                 config=self.config.astromRefObjLoader, log=self.log)
            self.astrometry.setRefObjLoader(refObjLoader)

        if self.config.doPhotoCal:
            photoRefObjLoader = ReferenceObjectLoader(dataIds=[ref.datasetRef.dataId
                                                      for ref in inputRefs.photoRefCat],
                                                      refCats=inputs.pop('photoRefCat'),
                                                      config=self.config.photoRefObjLoader,
                                                      log=self.log)
            self.photoCal.match.setRefObjLoader(photoRefObjLoader)

        outputs = self.run(**inputs)

        if self.config.doWriteMatches and self.config.doAstrometry:
            normalizedMatches = afwTable.packMatches(outputs.astromMatches)
            normalizedMatches.table.setMetadata(outputs.matchMeta)
            if self.config.doWriteMatchesDenormalized:
                denormMatches = denormalizeMatches(outputs.astromMatches, outputs.matchMeta)
                outputs.matchesDenormalized = denormMatches
            outputs.matches = normalizedMatches
        butlerQC.put(outputs, outputRefs)

    @pipeBase.timeMethod
    def run(self, exposure, exposureIdInfo=None, background=None,
            icSourceCat=None):
        """!Calibrate an exposure (science image or coadd)

        @param[in,out] exposure  exposure to calibrate (an
            lsst.afw.image.ExposureF or similar);
            in:
            - MaskedImage
            - Psf
            out:
            - MaskedImage has background subtracted
            - Wcs is replaced
            - PhotoCalib is replaced
        @param[in] exposureIdInfo  ID info for exposure (an
            lsst.obs.base.ExposureIdInfo) If not provided, returned
            SourceCatalog IDs will not be globally unique.
        @param[in,out] background  background model already subtracted from
            exposure (an lsst.afw.math.BackgroundList). May be None if no
            background has been subtracted, though that is unusual for
            calibration. A refined background model is output.
        @param[in] icSourceCat  A SourceCatalog from CharacterizeImageTask
            from which we can copy some fields.

        @return pipe_base Struct containing these fields:
        - exposure  calibrate science exposure with refined WCS and PhotoCalib
        - background  model of background subtracted from exposure (an
          lsst.afw.math.BackgroundList)
        - sourceCat  catalog of measured sources
        - astromMatches  list of source/refObj matches from the astrometry
          solver
        """
        # detect, deblend and measure sources
        if exposureIdInfo is None:
            exposureIdInfo = ExposureIdInfo()

        if background is None:
            background = BackgroundList()
        sourceIdFactory = exposureIdInfo.makeSourceIdFactory()
        table = SourceTable.make(self.schema, sourceIdFactory)
        table.setMetadata(self.algMetadata)

        detRes = self.detection.run(table=table, exposure=exposure,
                                    doSmooth=True)
        sourceCat = detRes.sources
        if detRes.fpSets.background:
            for bg in detRes.fpSets.background:
                background.append(bg)
        if self.config.doSkySources:
            skySourceFootprints = self.skySources.run(mask=exposure.mask, seed=exposureIdInfo.expId)
            if skySourceFootprints:
                for foot in skySourceFootprints:
                    s = sourceCat.addNew()
                    s.setFootprint(foot)
                    s.set(self.skySourceKey, True)
        if self.config.doDeblend:
            self.deblend.run(exposure=exposure, sources=sourceCat)
        self.measurement.run(
            measCat=sourceCat,
            exposure=exposure,
            exposureId=exposureIdInfo.expId
        )
        if self.config.doApCorr:
            self.applyApCorr.run(
                catalog=sourceCat,
                apCorrMap=exposure.getInfo().getApCorrMap()
            )
        self.catalogCalculation.run(sourceCat)

        self.setPrimaryFlags.run(sourceCat)

        if icSourceCat is not None and \
           len(self.config.icSourceFieldsToCopy) > 0:
            self.copyIcSourceFields(icSourceCat=icSourceCat,
                                    sourceCat=sourceCat)

        # TODO DM-11568: this contiguous check-and-copy could go away if we
        # reserve enough space during SourceDetection and/or SourceDeblend.
        # NOTE: sourceSelectors require contiguous catalogs, so ensure
        # contiguity now, so views are preserved from here on.
#        if not sourceCat.isContiguous():
#            sourceCat = sourceCat.copy(deep=True)

        # perform astrometry calibration:
        # fit an improved WCS and update the exposure's WCS in place
        astromMatches = None
        matchMeta = None
        if self.config.doAstrometry:
            try:
                astromRes = self.astrometry.run(
                    exposure=exposure,
                    sourceCat=sourceCat,
                )
                astromMatches = astromRes.matches
                matchMeta = astromRes.matchMeta
            except Exception as e:
                if self.config.requireAstrometry:
                    raise
                self.log.warn("Unable to perform astrometric calibration "
                              "(%s): attempting to proceed" % e)

        # compute photometric calibration
        if self.config.doPhotoCal:
            try:
                photoRes = self.photoCal.run(exposure, sourceCat=sourceCat, expId=exposureIdInfo.expId)
                exposure.setPhotoCalib(photoRes.photoCalib)
                # TODO: reword this to phrase it in terms of the calibration factor?
                self.log.info("Photometric zero-point: %f" %
                              photoRes.photoCalib.instFluxToMagnitude(1.0))
                self.setMetadata(exposure=exposure, photoRes=photoRes)
            except Exception as e:
                if self.config.requirePhotoCal:
                    raise
                self.log.warn("Unable to perform photometric calibration "
                              "(%s): attempting to proceed" % e)
                self.setMetadata(exposure=exposure, photoRes=None)

        if self.config.doInsertFakes:
            self.insertFakes.run(exposure, background=background)

            table = SourceTable.make(self.schema, sourceIdFactory)
            table.setMetadata(self.algMetadata)

            detRes = self.detection.run(table=table, exposure=exposure,
                                        doSmooth=True)
            sourceCat = detRes.sources
            if detRes.fpSets.background:
                for bg in detRes.fpSets.background:
                    background.append(bg)
            if self.config.doDeblend:
                self.deblend.run(exposure=exposure, sources=sourceCat)
            self.measurement.run(
                measCat=sourceCat,
                exposure=exposure,
                exposureId=exposureIdInfo.expId
            )
            if self.config.doApCorr:
                self.applyApCorr.run(
                    catalog=sourceCat,
                    apCorrMap=exposure.getInfo().getApCorrMap()
                )
            self.catalogCalculation.run(sourceCat)

            if icSourceCat is not None and len(self.config.icSourceFieldsToCopy) > 0:
                self.copyIcSourceFields(icSourceCat=icSourceCat,
                                        sourceCat=sourceCat)

        frame = getDebugFrame(self._display, "calibrate")
        if frame:
            displayAstrometry(
                sourceCat=sourceCat,
                exposure=exposure,
                matches=astromMatches,
                frame=frame,
                pause=False,
            )

        return pipeBase.Struct(
            exposure=exposure,
            background=background,
            sourceCat=sourceCat,
            astromMatches=astromMatches,
            matchMeta=matchMeta,
            # These are duplicate entries with different names for use with
            # gen3 middleware
            outputExposure=exposure,
            outputCat=sourceCat,
            outputBackground=background,
        )

    def writeOutputs(self, dataRef, exposure, background, sourceCat,
                     astromMatches, matchMeta):
        """Write output data to the output repository

        @param[in] dataRef  butler data reference corresponding to a science
            image
        @param[in] exposure  exposure to write
        @param[in] background  background model for exposure
        @param[in] sourceCat  catalog of measured sources
        @param[in] astromMatches  list of source/refObj matches from the
            astrometry solver
        """
        dataRef.put(sourceCat, "src")
        if self.config.doWriteMatches and astromMatches is not None:
            normalizedMatches = afwTable.packMatches(astromMatches)
            normalizedMatches.table.setMetadata(matchMeta)
            dataRef.put(normalizedMatches, "srcMatch")
            if self.config.doWriteMatchesDenormalized:
                denormMatches = denormalizeMatches(astromMatches, matchMeta)
                dataRef.put(denormMatches, "srcMatchFull")
        if self.config.doWriteExposure:
            dataRef.put(exposure, "calexp")
            dataRef.put(background, "calexpBackground")

    def getSchemaCatalogs(self):
        """Return a dict of empty catalogs for each catalog dataset produced
        by this task.
        """
        sourceCat = afwTable.SourceCatalog(self.schema)
        sourceCat.getTable().setMetadata(self.algMetadata)
        return {"src": sourceCat}

    def setMetadata(self, exposure, photoRes=None):
        """!Set task and exposure metadata

        Logs a warning and continues if needed data is missing.

        @param[in,out] exposure  exposure whose metadata is to be set
        @param[in]  photoRes  results of running photoCal; if None then it was
            not run
        """
        if photoRes is None:
            return

        metadata = exposure.getMetadata()

        # convert zero-point to (mag/sec/adu) for task MAGZERO metadata
        try:
            exposureTime = exposure.getInfo().getVisitInfo().getExposureTime()
            magZero = photoRes.zp - 2.5*math.log10(exposureTime)
        except Exception:
            self.log.warn("Could not set normalized MAGZERO in header: no "
                          "exposure time")
            magZero = math.nan

        try:
            metadata.set('MAGZERO', magZero)
            metadata.set('MAGZERO_RMS', photoRes.sigma)
            metadata.set('MAGZERO_NOBJ', photoRes.ngood)
            metadata.set('COLORTERM1', 0.0)
            metadata.set('COLORTERM2', 0.0)
            metadata.set('COLORTERM3', 0.0)
        except Exception as e:
            self.log.warn("Could not set exposure metadata: %s" % (e,))

    def copyIcSourceFields(self, icSourceCat, sourceCat):
        """!Match sources in icSourceCat and sourceCat and copy the specified fields

        @param[in] icSourceCat  catalog from which to copy fields
        @param[in,out] sourceCat  catalog to which to copy fields

        The fields copied are those specified by `config.icSourceFieldsToCopy`
        that actually exist in the schema. This was set up by the constructor
        using self.schemaMapper.
        """
        if self.schemaMapper is None:
            raise RuntimeError("To copy icSource fields you must specify "
                               "icSourceSchema nd icSourceKeys when "
                               "constructing this task")
        if icSourceCat is None or sourceCat is None:
            raise RuntimeError("icSourceCat and sourceCat must both be "
                               "specified")
        if len(self.config.icSourceFieldsToCopy) == 0:
            self.log.warn("copyIcSourceFields doing nothing because "
                          "icSourceFieldsToCopy is empty")
            return

        mc = afwTable.MatchControl()
        mc.findOnlyClosest = False  # return all matched objects
        matches = afwTable.matchXy(icSourceCat, sourceCat,
                                   self.config.matchRadiusPix, mc)
        if self.config.doDeblend:
            deblendKey = sourceCat.schema["deblend_nChild"].asKey()
            # if deblended, keep children
            matches = [m for m in matches if m[1].get(deblendKey) == 0]

        # Because we had to allow multiple matches to handle parents, we now
        # need to prune to the best matches
        # closest matches as a dict of icSourceCat source ID:
        # (icSourceCat source, sourceCat source, distance in pixels)
        bestMatches = {}
        for m0, m1, d in matches:
            id0 = m0.getId()
            match = bestMatches.get(id0)
            if match is None or d <= match[2]:
                bestMatches[id0] = (m0, m1, d)
        matches = list(bestMatches.values())

        # Check that no sourceCat sources are listed twice (we already know
        # that each match has a unique icSourceCat source ID, due to using
        # that ID as the key in bestMatches)
        numMatches = len(matches)
        numUniqueSources = len(set(m[1].getId() for m in matches))
        if numUniqueSources != numMatches:
            self.log.warn("{} icSourceCat sources matched only {} sourceCat "
                          "sources".format(numMatches, numUniqueSources))

        self.log.info("Copying flags from icSourceCat to sourceCat for "
                      "%s sources" % (numMatches,))

        # For each match: set the calibSourceKey flag and copy the desired
        # fields
        for icSrc, src, d in matches:
            src.setFlag(self.calibSourceKey, True)
            # src.assign copies the footprint from icSrc, which we don't want
            # (DM-407)
            # so set icSrc's footprint to src's footprint before src.assign,
            # then restore it
            icSrcFootprint = icSrc.getFootprint()
            try:
                icSrc.setFootprint(src.getFootprint())
                src.assign(icSrc, self.schemaMapper)
            finally:
                icSrc.setFootprint(icSrcFootprint)
