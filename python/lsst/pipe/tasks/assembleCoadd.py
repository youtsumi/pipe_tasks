# This file is part of pipe_tasks.
#
# LSST Data Management System
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
# See COPYRIGHT file at the top of the source tree.
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
import os
import copy
import numpy
import warnings
import lsst.pex.config as pexConfig
import lsst.pex.exceptions as pexExceptions
import lsst.geom as geom
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.detection as afwDet
import lsst.coadd.utils as coaddUtils
import lsst.pipe.base as pipeBase
import lsst.meas.algorithms as measAlg
import lsst.log as log
import lsstDebug
import lsst.utils as utils
from lsst.skymap import BaseSkyMap
from .coaddBase import CoaddBaseTask, SelectDataIdContainer, makeSkyInfo, makeCoaddSuffix, reorderAndPadList
from .interpImage import InterpImageTask
from .scaleZeroPoint import ScaleZeroPointTask
from .coaddHelpers import groupPatchExposures, getGroupDataRef
from .scaleVariance import ScaleVarianceTask
from .maskStreaks import MaskStreaksTask
from .healSparseMapping import HealSparseInputMapTask
from lsst.meas.algorithms import SourceDetectionTask
from lsst.daf.butler import DeferredDatasetHandle

__all__ = ["AssembleCoaddTask", "AssembleCoaddConnections", "AssembleCoaddConfig",
           "SafeClipAssembleCoaddTask", "SafeClipAssembleCoaddConfig",
           "CompareWarpAssembleCoaddTask", "CompareWarpAssembleCoaddConfig"]


class AssembleCoaddConnections(pipeBase.PipelineTaskConnections,
                               dimensions=("tract", "patch", "band", "skymap"),
                               defaultTemplates={"inputCoaddName": "deep",
                                                 "outputCoaddName": "deep",
                                                 "warpType": "direct",
                                                 "warpTypeSuffix": "",
                                                 "fakesType": ""}):
    inputWarps = pipeBase.connectionTypes.Input(
        doc=("Input list of warps to be assemebled i.e. stacked."
             "WarpType (e.g. direct, psfMatched) is controlled by the warpType config parameter"),
        name="{inputCoaddName}Coadd_{warpType}Warp",
        storageClass="ExposureF",
        dimensions=("tract", "patch", "skymap", "visit", "instrument"),
        deferLoad=True,
        multiple=True
    )
    skyMap = pipeBase.connectionTypes.Input(
        doc="Input definition of geometry/bbox and projection/wcs for coadded exposures",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap", ),
    )
    selectedVisits = pipeBase.connectionTypes.Input(
        doc="Selected visits to be coadded.",
        name="{outputCoaddName}VisitsDict",
        storageClass="StructuredDataDict",
        dimensions=("instrument", "tract", "patch", "skymap", "band")
    )
    brightObjectMask = pipeBase.connectionTypes.PrerequisiteInput(
        doc=("Input Bright Object Mask mask produced with external catalogs to be applied to the mask plane"
             " BRIGHT_OBJECT."),
        name="brightObjectMask",
        storageClass="ObjectMaskCatalog",
        dimensions=("tract", "patch", "skymap", "band"),
    )
    coaddExposure = pipeBase.connectionTypes.Output(
        doc="Output coadded exposure, produced by stacking input warps",
        name="{fakesType}{outputCoaddName}Coadd{warpTypeSuffix}",
        storageClass="ExposureF",
        dimensions=("tract", "patch", "skymap", "band"),
    )
    nImage = pipeBase.connectionTypes.Output(
        doc="Output image of number of input images per pixel",
        name="{outputCoaddName}Coadd_nImage",
        storageClass="ImageU",
        dimensions=("tract", "patch", "skymap", "band"),
    )
    inputMap = pipeBase.connectionTypes.Output(
        doc="Output healsparse map of input images",
        name="{outputCoaddName}Coadd_inputMap",
        storageClass="HealSparseMap",
        dimensions=("tract", "patch", "skymap", "band"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        # Override the connection's name template with config to replicate Gen2 behavior
        # This duplicates some of the logic in the base class, due to wanting Gen2 and
        # Gen3 configs to stay in sync. This should be removed when gen2 is deprecated
        templateValues = {name: getattr(config.connections, name) for name in self.defaultTemplates}
        templateValues['warpType'] = config.warpType
        templateValues['warpTypeSuffix'] = makeCoaddSuffix(config.warpType)
        if config.hasFakes:
            templateValues['fakesType'] = "_fakes"
        self._nameOverrides = {name: getattr(config.connections, name).format(**templateValues)
                               for name in self.allConnections}
        self._typeNameToVarName = {v: k for k, v in self._nameOverrides.items()}
        # End code to remove after deprecation

        if not config.doMaskBrightObjects:
            self.prerequisiteInputs.remove("brightObjectMask")

        if not config.doSelectVisits:
            self.inputs.remove("selectedVisits")

        if not config.doNImage:
            self.outputs.remove("nImage")

        if not self.config.doInputMap:
            self.outputs.remove("inputMap")


class AssembleCoaddConfig(CoaddBaseTask.ConfigClass, pipeBase.PipelineTaskConfig,
                          pipelineConnections=AssembleCoaddConnections):
    """Configuration parameters for the `AssembleCoaddTask`.

    Notes
    -----
    The `doMaskBrightObjects` and `brightObjectMaskName` configuration options
    only set the bitplane config.brightObjectMaskName. To make this useful you
    *must* also configure the flags.pixel algorithm, for example by adding

    .. code-block:: none

       config.measurement.plugins["base_PixelFlags"].masksFpCenter.append("BRIGHT_OBJECT")
       config.measurement.plugins["base_PixelFlags"].masksFpAnywhere.append("BRIGHT_OBJECT")

    to your measureCoaddSources.py and forcedPhotCoadd.py config overrides.
    """
    warpType = pexConfig.Field(
        doc="Warp name: one of 'direct' or 'psfMatched'",
        dtype=str,
        default="direct",
    )
    subregionSize = pexConfig.ListField(
        dtype=int,
        doc="Width, height of stack subregion size; "
        "make small enough that a full stack of images will fit into memory at once.",
        length=2,
        default=(2000, 2000),
    )
    statistic = pexConfig.Field(
        dtype=str,
        doc="Main stacking statistic for aggregating over the epochs.",
        default="MEANCLIP",
    )
    doSigmaClip = pexConfig.Field(
        dtype=bool,
        doc="Perform sigma clipped outlier rejection with MEANCLIP statistic? (DEPRECATED)",
        default=False,
    )
    sigmaClip = pexConfig.Field(
        dtype=float,
        doc="Sigma for outlier rejection; ignored if non-clipping statistic selected.",
        default=3.0,
    )
    clipIter = pexConfig.Field(
        dtype=int,
        doc="Number of iterations of outlier rejection; ignored if non-clipping statistic selected.",
        default=2,
    )
    calcErrorFromInputVariance = pexConfig.Field(
        dtype=bool,
        doc="Calculate coadd variance from input variance by stacking statistic."
            "Passed to StatisticsControl.setCalcErrorFromInputVariance()",
        default=True,
    )
    scaleZeroPoint = pexConfig.ConfigurableField(
        target=ScaleZeroPointTask,
        doc="Task to adjust the photometric zero point of the coadd temp exposures",
    )
    doInterp = pexConfig.Field(
        doc="Interpolate over NaN pixels? Also extrapolate, if necessary, but the results are ugly.",
        dtype=bool,
        default=True,
    )
    interpImage = pexConfig.ConfigurableField(
        target=InterpImageTask,
        doc="Task to interpolate (and extrapolate) over NaN pixels",
    )
    doWrite = pexConfig.Field(
        doc="Persist coadd?",
        dtype=bool,
        default=True,
    )
    doNImage = pexConfig.Field(
        doc="Create image of number of contributing exposures for each pixel",
        dtype=bool,
        default=False,
    )
    doUsePsfMatchedPolygons = pexConfig.Field(
        doc="Use ValidPolygons from shrunk Psf-Matched Calexps? Should be set to True by CompareWarp only.",
        dtype=bool,
        default=False,
    )
    maskPropagationThresholds = pexConfig.DictField(
        keytype=str,
        itemtype=float,
        doc=("Threshold (in fractional weight) of rejection at which we propagate a mask plane to "
             "the coadd; that is, we set the mask bit on the coadd if the fraction the rejected frames "
             "would have contributed exceeds this value."),
        default={"SAT": 0.1},
    )
    removeMaskPlanes = pexConfig.ListField(dtype=str, default=["NOT_DEBLENDED"],
                                           doc="Mask planes to remove before coadding")
    doMaskBrightObjects = pexConfig.Field(dtype=bool, default=False,
                                          doc="Set mask and flag bits for bright objects?")
    brightObjectMaskName = pexConfig.Field(dtype=str, default="BRIGHT_OBJECT",
                                           doc="Name of mask bit used for bright objects")
    coaddPsf = pexConfig.ConfigField(
        doc="Configuration for CoaddPsf",
        dtype=measAlg.CoaddPsfConfig,
    )
    doAttachTransmissionCurve = pexConfig.Field(
        dtype=bool, default=False, optional=False,
        doc=("Attach a piecewise TransmissionCurve for the coadd? "
             "(requires all input Exposures to have TransmissionCurves).")
    )
    hasFakes = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Should be set to True if fake sources have been inserted into the input data."
    )
    doSelectVisits = pexConfig.Field(
        doc="Coadd only visits selected by a SelectVisitsTask",
        dtype=bool,
        default=False,
    )
    doInputMap = pexConfig.Field(
        doc="Create a bitwise map of coadd inputs",
        dtype=bool,
        default=False,
    )
    inputMapper = pexConfig.ConfigurableField(
        doc="Input map creation subtask.",
        target=HealSparseInputMapTask,
    )

    def setDefaults(self):
        super().setDefaults()
        self.badMaskPlanes = ["NO_DATA", "BAD", "SAT", "EDGE"]

    def validate(self):
        super().validate()
        if self.doPsfMatch:
            # Backwards compatibility.
            # Configs do not have loggers
            log.warn("Config doPsfMatch deprecated. Setting warpType='psfMatched'")
            self.warpType = 'psfMatched'
        if self.doSigmaClip and self.statistic != "MEANCLIP":
            log.warn('doSigmaClip deprecated. To replicate behavior, setting statistic to "MEANCLIP"')
            self.statistic = "MEANCLIP"
        if self.doInterp and self.statistic not in ['MEAN', 'MEDIAN', 'MEANCLIP', 'VARIANCE', 'VARIANCECLIP']:
            raise ValueError("Must set doInterp=False for statistic=%s, which does not "
                             "compute and set a non-zero coadd variance estimate." % (self.statistic))

        unstackableStats = ['NOTHING', 'ERROR', 'ORMASK']
        if not hasattr(afwMath.Property, self.statistic) or self.statistic in unstackableStats:
            stackableStats = [str(k) for k in afwMath.Property.__members__.keys()
                              if str(k) not in unstackableStats]
            raise ValueError("statistic %s is not allowed. Please choose one of %s."
                             % (self.statistic, stackableStats))


class AssembleCoaddTask(CoaddBaseTask, pipeBase.PipelineTask):
    """Assemble a coadded image from a set of warps (coadded temporary exposures).

    We want to assemble a coadded image from a set of Warps (also called
    coadded temporary exposures or ``coaddTempExps``).
    Each input Warp covers a patch on the sky and corresponds to a single
    run/visit/exposure of the covered patch. We provide the task with a list
    of Warps (``selectDataList``) from which it selects Warps that cover the
    specified patch (pointed at by ``dataRef``).
    Each Warp that goes into a coadd will typically have an independent
    photometric zero-point. Therefore, we must scale each Warp to set it to
    a common photometric zeropoint. WarpType may be one of 'direct' or
    'psfMatched', and the boolean configs `config.makeDirect` and
    `config.makePsfMatched` set which of the warp types will be coadded.
    The coadd is computed as a mean with optional outlier rejection.
    Criteria for outlier rejection are set in `AssembleCoaddConfig`.
    Finally, Warps can have bad 'NaN' pixels which received no input from the
    source calExps. We interpolate over these bad (NaN) pixels.

    `AssembleCoaddTask` uses several sub-tasks. These are

    - `ScaleZeroPointTask`
    - create and use an ``imageScaler`` object to scale the photometric zeropoint for each Warp
    - `InterpImageTask`
    - interpolate across bad pixels (NaN) in the final coadd

    You can retarget these subtasks if you wish.

    Notes
    -----
    The `lsst.pipe.base.cmdLineTask.CmdLineTask` interface supports a
    flag ``-d`` to import ``debug.py`` from your ``PYTHONPATH``; see
    `baseDebug` for more about ``debug.py`` files. `AssembleCoaddTask` has
    no debug variables of its own. Some of the subtasks may support debug
    variables. See the documentation for the subtasks for further information.

    Examples
    --------
    `AssembleCoaddTask` assembles a set of warped images into a coadded image.
    The `AssembleCoaddTask` can be invoked by running ``assembleCoadd.py``
    with the flag '--legacyCoadd'. Usage of assembleCoadd.py expects two
    inputs: a data reference to the tract patch and filter to be coadded, and
    a list of Warps to attempt to coadd. These are specified using ``--id`` and
    ``--selectId``, respectively:

    .. code-block:: none

       --id = [KEY=VALUE1[^VALUE2[^VALUE3...] [KEY=VALUE1[^VALUE2[^VALUE3...] ...]]
       --selectId [KEY=VALUE1[^VALUE2[^VALUE3...] [KEY=VALUE1[^VALUE2[^VALUE3...] ...]]

    Only the Warps that cover the specified tract and patch will be coadded.
    A list of the available optional arguments can be obtained by calling
    ``assembleCoadd.py`` with the ``--help`` command line argument:

    .. code-block:: none

       assembleCoadd.py --help

    To demonstrate usage of the `AssembleCoaddTask` in the larger context of
    multi-band processing, we will generate the HSC-I & -R band coadds from
    HSC engineering test data provided in the ``ci_hsc`` package. To begin,
    assuming that the lsst stack has been already set up, we must set up the
    obs_subaru and ``ci_hsc`` packages. This defines the environment variable
    ``$CI_HSC_DIR`` and points at the location of the package. The raw HSC
    data live in the ``$CI_HSC_DIR/raw directory``. To begin assembling the
    coadds, we must first

    - processCcd
    - process the individual ccds in $CI_HSC_RAW to produce calibrated exposures
    - makeSkyMap
    - create a skymap that covers the area of the sky present in the raw exposures
    - makeCoaddTempExp
    - warp the individual calibrated exposures to the tangent plane of the coadd

    We can perform all of these steps by running

    .. code-block:: none

       $CI_HSC_DIR scons warp-903986 warp-904014 warp-903990 warp-904010 warp-903988

    This will produce warped exposures for each visit. To coadd the warped
    data, we call assembleCoadd.py as follows:

    .. code-block:: none

       assembleCoadd.py --legacyCoadd $CI_HSC_DIR/DATA --id patch=5,4 tract=0 filter=HSC-I \
       --selectId visit=903986 ccd=16 --selectId visit=903986 ccd=22 --selectId visit=903986 ccd=23 \
       --selectId visit=903986 ccd=100 --selectId visit=904014 ccd=1 --selectId visit=904014 ccd=6 \
       --selectId visit=904014 ccd=12 --selectId visit=903990 ccd=18 --selectId visit=903990 ccd=25 \
       --selectId visit=904010 ccd=4 --selectId visit=904010 ccd=10 --selectId visit=904010 ccd=100 \
       --selectId visit=903988 ccd=16 --selectId visit=903988 ccd=17 --selectId visit=903988 ccd=23 \
       --selectId visit=903988 ccd=24

    that will process the HSC-I band data. The results are written in
    ``$CI_HSC_DIR/DATA/deepCoadd-results/HSC-I``.

    You may also choose to run:

    .. code-block:: none

       scons warp-903334 warp-903336 warp-903338 warp-903342 warp-903344 warp-903346
       assembleCoadd.py --legacyCoadd $CI_HSC_DIR/DATA --id patch=5,4 tract=0 filter=HSC-R \
       --selectId visit=903334 ccd=16 --selectId visit=903334 ccd=22 --selectId visit=903334 ccd=23 \
       --selectId visit=903334 ccd=100 --selectId visit=903336 ccd=17 --selectId visit=903336 ccd=24 \
       --selectId visit=903338 ccd=18 --selectId visit=903338 ccd=25 --selectId visit=903342 ccd=4 \
       --selectId visit=903342 ccd=10 --selectId visit=903342 ccd=100 --selectId visit=903344 ccd=0 \
       --selectId visit=903344 ccd=5 --selectId visit=903344 ccd=11 --selectId visit=903346 ccd=1 \
       --selectId visit=903346 ccd=6 --selectId visit=903346 ccd=12

    to generate the coadd for the HSC-R band if you are interested in
    following multiBand Coadd processing as discussed in `pipeTasks_multiBand`
    (but note that normally, one would use the `SafeClipAssembleCoaddTask`
    rather than `AssembleCoaddTask` to make the coadd.
    """
    ConfigClass = AssembleCoaddConfig
    _DefaultName = "assembleCoadd"

    def __init__(self, *args, **kwargs):
        # TODO: DM-17415 better way to handle previously allowed passed args e.g.`AssembleCoaddTask(config)`
        if args:
            argNames = ["config", "name", "parentTask", "log"]
            kwargs.update({k: v for k, v in zip(argNames, args)})
            warnings.warn("AssembleCoadd received positional args, and casting them as kwargs: %s. "
                          "PipelineTask will not take positional args" % argNames, FutureWarning)

        super().__init__(**kwargs)
        self.makeSubtask("interpImage")
        self.makeSubtask("scaleZeroPoint")

        if self.config.doMaskBrightObjects:
            mask = afwImage.Mask()
            try:
                self.brightObjectBitmask = 1 << mask.addMaskPlane(self.config.brightObjectMaskName)
            except pexExceptions.LsstCppException:
                raise RuntimeError("Unable to define mask plane for bright objects; planes used are %s" %
                                   mask.getMaskPlaneDict().keys())
            del mask

        if self.config.doInputMap:
            self.makeSubtask("inputMapper")

        self.warpType = self.config.warpType

    @utils.inheritDoc(pipeBase.PipelineTask)
    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Docstring to be formatted with info from PipelineTask.runQuantum
        """
        Notes
        -----
        Assemble a coadd from a set of Warps.

        PipelineTask (Gen3) entry point to Coadd a set of Warps.
        Analogous to `runDataRef`, it prepares all the data products to be
        passed to `run`, and processes the results before returning a struct
        of results to be written out. AssembleCoadd cannot fit all Warps in memory.
        Therefore, its inputs are accessed subregion by subregion
        by the Gen3 `DeferredDatasetHandle` that is analagous to the Gen2
        `lsst.daf.persistence.ButlerDataRef`. Any updates to this method should
        correspond to an update in `runDataRef` while both entry points
        are used.
        """
        inputData = butlerQC.get(inputRefs)

        # Construct skyInfo expected by run
        # Do not remove skyMap from inputData in case makeSupplementaryDataGen3 needs it
        skyMap = inputData["skyMap"]
        outputDataId = butlerQC.quantum.dataId

        inputData['skyInfo'] = makeSkyInfo(skyMap,
                                           tractId=outputDataId['tract'],
                                           patchId=outputDataId['patch'])

        if self.config.doSelectVisits:
            warpRefList = self.filterWarps(inputData['inputWarps'], inputData['selectedVisits'])
        else:
            warpRefList = inputData['inputWarps']

        # Perform same middle steps as `runDataRef` does
        inputs = self.prepareInputs(warpRefList)
        self.log.info("Found %d %s", len(inputs.tempExpRefList),
                      self.getTempExpDatasetName(self.warpType))
        if len(inputs.tempExpRefList) == 0:
            self.log.warn("No coadd temporary exposures found")
            return

        supplementaryData = self.makeSupplementaryDataGen3(butlerQC, inputRefs, outputRefs)
        retStruct = self.run(inputData['skyInfo'], inputs.tempExpRefList, inputs.imageScalerList,
                             inputs.weightList, supplementaryData=supplementaryData)

        inputData.setdefault('brightObjectMask', None)
        self.processResults(retStruct.coaddExposure, inputData['brightObjectMask'], outputDataId)

        if self.config.doWrite:
            butlerQC.put(retStruct, outputRefs)
        return retStruct

    @pipeBase.timeMethod
    def runDataRef(self, dataRef, selectDataList=None, warpRefList=None):
        """Assemble a coadd from a set of Warps.

        Pipebase.CmdlineTask entry point to Coadd a set of Warps.
        Compute weights to be applied to each Warp and
        find scalings to match the photometric zeropoint to a reference Warp.
        Assemble the Warps using `run`. Interpolate over NaNs and
        optionally write the coadd to disk. Return the coadded exposure.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.butlerSubset.ButlerDataRef`
            Data reference defining the patch for coaddition and the
            reference Warp (if ``config.autoReference=False``).
            Used to access the following data products:
            - ``self.config.coaddName + "Coadd_skyMap"``
            - ``self.config.coaddName + "Coadd_ + <warpType> + "Warp"`` (optionally)
            - ``self.config.coaddName + "Coadd"``
        selectDataList : `list`
            List of data references to Calexps. Data to be coadded will be
            selected from this list based on overlap with the patch defined
            by dataRef, grouped by visit, and converted to a list of data
            references to warps.
        warpRefList : `list`
            List of data references to Warps to be coadded.
            Note: `warpRefList` is just the new name for `tempExpRefList`.

        Returns
        -------
        retStruct : `lsst.pipe.base.Struct`
           Result struct with components:

           - ``coaddExposure``: coadded exposure (``Exposure``).
           - ``nImage``: exposure count image (``Image``).
        """
        if selectDataList and warpRefList:
            raise RuntimeError("runDataRef received both a selectDataList and warpRefList, "
                               "and which to use is ambiguous. Please pass only one.")

        skyInfo = self.getSkyInfo(dataRef)
        if warpRefList is None:
            calExpRefList = self.selectExposures(dataRef, skyInfo, selectDataList=selectDataList)
            if len(calExpRefList) == 0:
                self.log.warn("No exposures to coadd")
                return
            self.log.info("Coadding %d exposures", len(calExpRefList))

            warpRefList = self.getTempExpRefList(dataRef, calExpRefList)

        inputData = self.prepareInputs(warpRefList)
        self.log.info("Found %d %s", len(inputData.tempExpRefList),
                      self.getTempExpDatasetName(self.warpType))
        if len(inputData.tempExpRefList) == 0:
            self.log.warn("No coadd temporary exposures found")
            return

        supplementaryData = self.makeSupplementaryData(dataRef, warpRefList=inputData.tempExpRefList)

        retStruct = self.run(skyInfo, inputData.tempExpRefList, inputData.imageScalerList,
                             inputData.weightList, supplementaryData=supplementaryData)

        brightObjects = self.readBrightObjectMasks(dataRef) if self.config.doMaskBrightObjects else None
        self.processResults(retStruct.coaddExposure, brightObjectMasks=brightObjects, dataId=dataRef.dataId)

        if self.config.doWrite:
            if self.getCoaddDatasetName(self.warpType) == "deepCoadd" and self.config.hasFakes:
                coaddDatasetName = "fakes_" + self.getCoaddDatasetName(self.warpType)
            else:
                coaddDatasetName = self.getCoaddDatasetName(self.warpType)
            self.log.info("Persisting %s" % coaddDatasetName)
            dataRef.put(retStruct.coaddExposure, coaddDatasetName)
        if self.config.doNImage and retStruct.nImage is not None:
            dataRef.put(retStruct.nImage, self.getCoaddDatasetName(self.warpType) + '_nImage')

        return retStruct

    def processResults(self, coaddExposure, brightObjectMasks=None, dataId=None):
        """Interpolate over missing data and mask bright stars.

        Parameters
        ----------
        coaddExposure : `lsst.afw.image.Exposure`
            The coadded exposure to process.
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference for supplementary data.
        """
        if self.config.doInterp:
            self.interpImage.run(coaddExposure.getMaskedImage(), planeName="NO_DATA")
            # The variance must be positive; work around for DM-3201.
            varArray = coaddExposure.variance.array
            with numpy.errstate(invalid="ignore"):
                varArray[:] = numpy.where(varArray > 0, varArray, numpy.inf)

        if self.config.doMaskBrightObjects:
            self.setBrightObjectMasks(coaddExposure, brightObjectMasks, dataId)

    def makeSupplementaryData(self, dataRef, selectDataList=None, warpRefList=None):
        """Make additional inputs to run() specific to subclasses (Gen2)

        Duplicates interface of `runDataRef` method
        Available to be implemented by subclasses only if they need the
        coadd dataRef for performing preliminary processing before
        assembling the coadd.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference for supplementary data.
        selectDataList : `list` (optional)
            Optional List of data references to Calexps.
        warpRefList : `list` (optional)
            Optional List of data references to Warps.
        """
        return pipeBase.Struct()

    def makeSupplementaryDataGen3(self, butlerQC, inputRefs, outputRefs):
        """Make additional inputs to run() specific to subclasses (Gen3)

        Duplicates interface of `runQuantum` method.
        Available to be implemented by subclasses only if they need the
        coadd dataRef for performing preliminary processing before
        assembling the coadd.

        Parameters
        ----------
        butlerQC : `lsst.pipe.base.ButlerQuantumContext`
            Gen3 Butler object for fetching additional data products before
            running the Task specialized for quantum being processed
        inputRefs : `lsst.pipe.base.InputQuantizedConnection`
            Attributes are the names of the connections describing input dataset types.
            Values are DatasetRefs that task consumes for corresponding dataset type.
            DataIds are guaranteed to match data objects in ``inputData``.
        outputRefs : `lsst.pipe.base.OutputQuantizedConnection`
            Attributes are the names of the connections describing output dataset types.
            Values are DatasetRefs that task is to produce
            for corresponding dataset type.
        """
        return pipeBase.Struct()

    def getTempExpRefList(self, patchRef, calExpRefList):
        """Generate list data references corresponding to warped exposures
        that lie within the patch to be coadded.

        Parameters
        ----------
        patchRef : `dataRef`
            Data reference for patch.
        calExpRefList : `list`
            List of data references for input calexps.

        Returns
        -------
        tempExpRefList : `list`
            List of Warp/CoaddTempExp data references.
        """
        butler = patchRef.getButler()
        groupData = groupPatchExposures(patchRef, calExpRefList, self.getCoaddDatasetName(self.warpType),
                                        self.getTempExpDatasetName(self.warpType))
        tempExpRefList = [getGroupDataRef(butler, self.getTempExpDatasetName(self.warpType),
                                          g, groupData.keys) for
                          g in groupData.groups.keys()]
        return tempExpRefList

    def prepareInputs(self, refList):
        """Prepare the input warps for coaddition by measuring the weight for
        each warp and the scaling for the photometric zero point.

        Each Warp has its own photometric zeropoint and background variance.
        Before coadding these Warps together, compute a scale factor to
        normalize the photometric zeropoint and compute the weight for each Warp.

        Parameters
        ----------
        refList : `list`
            List of data references to tempExp

        Returns
        -------
        result : `lsst.pipe.base.Struct`
           Result struct with components:

           - ``tempExprefList``: `list` of data references to tempExp.
           - ``weightList``: `list` of weightings.
           - ``imageScalerList``: `list` of image scalers.
        """
        statsCtrl = afwMath.StatisticsControl()
        statsCtrl.setNumSigmaClip(self.config.sigmaClip)
        statsCtrl.setNumIter(self.config.clipIter)
        statsCtrl.setAndMask(self.getBadPixelMask())
        statsCtrl.setNanSafe(True)
        # compute tempExpRefList: a list of tempExpRef that actually exist
        # and weightList: a list of the weight of the associated coadd tempExp
        # and imageScalerList: a list of scale factors for the associated coadd tempExp
        tempExpRefList = []
        weightList = []
        imageScalerList = []
        tempExpName = self.getTempExpDatasetName(self.warpType)
        for tempExpRef in refList:
            # Gen3's DeferredDatasetHandles are guaranteed to exist and
            # therefore have no datasetExists() method
            if not isinstance(tempExpRef, DeferredDatasetHandle):
                if not tempExpRef.datasetExists(tempExpName):
                    self.log.warn("Could not find %s %s; skipping it", tempExpName, tempExpRef.dataId)
                    continue

            tempExp = tempExpRef.get(datasetType=tempExpName, immediate=True)
            # Ignore any input warp that is empty of data
            if numpy.isnan(tempExp.image.array).all():
                continue
            maskedImage = tempExp.getMaskedImage()
            imageScaler = self.scaleZeroPoint.computeImageScaler(
                exposure=tempExp,
                dataRef=tempExpRef,
            )
            try:
                imageScaler.scaleMaskedImage(maskedImage)
            except Exception as e:
                self.log.warn("Scaling failed for %s (skipping it): %s", tempExpRef.dataId, e)
                continue
            statObj = afwMath.makeStatistics(maskedImage.getVariance(), maskedImage.getMask(),
                                             afwMath.MEANCLIP, statsCtrl)
            meanVar, meanVarErr = statObj.getResult(afwMath.MEANCLIP)
            weight = 1.0 / float(meanVar)
            if not numpy.isfinite(weight):
                self.log.warn("Non-finite weight for %s: skipping", tempExpRef.dataId)
                continue
            self.log.info("Weight of %s %s = %0.3f", tempExpName, tempExpRef.dataId, weight)

            del maskedImage
            del tempExp

            tempExpRefList.append(tempExpRef)
            weightList.append(weight)
            imageScalerList.append(imageScaler)

        return pipeBase.Struct(tempExpRefList=tempExpRefList, weightList=weightList,
                               imageScalerList=imageScalerList)

    def prepareStats(self, mask=None):
        """Prepare the statistics for coadding images.

        Parameters
        ----------
        mask : `int`, optional
            Bit mask value to exclude from coaddition.

        Returns
        -------
        stats : `lsst.pipe.base.Struct`
            Statistics structure with the following fields:

            - ``statsCtrl``: Statistics control object for coadd
                (`lsst.afw.math.StatisticsControl`)
            - ``statsFlags``: Statistic for coadd (`lsst.afw.math.Property`)
        """
        if mask is None:
            mask = self.getBadPixelMask()
        statsCtrl = afwMath.StatisticsControl()
        statsCtrl.setNumSigmaClip(self.config.sigmaClip)
        statsCtrl.setNumIter(self.config.clipIter)
        statsCtrl.setAndMask(mask)
        statsCtrl.setNanSafe(True)
        statsCtrl.setWeighted(True)
        statsCtrl.setCalcErrorFromInputVariance(self.config.calcErrorFromInputVariance)
        for plane, threshold in self.config.maskPropagationThresholds.items():
            bit = afwImage.Mask.getMaskPlane(plane)
            statsCtrl.setMaskPropagationThreshold(bit, threshold)
        statsFlags = afwMath.stringToStatisticsProperty(self.config.statistic)
        return pipeBase.Struct(ctrl=statsCtrl, flags=statsFlags)

    @pipeBase.timeMethod
    def run(self, skyInfo, tempExpRefList, imageScalerList, weightList,
            altMaskList=None, mask=None, supplementaryData=None):
        """Assemble a coadd from input warps

        Assemble the coadd using the provided list of coaddTempExps. Since
        the full coadd covers a patch (a large area), the assembly is
        performed over small areas on the image at a time in order to
        conserve memory usage. Iterate over subregions within the outer
        bbox of the patch using `assembleSubregion` to stack the corresponding
        subregions from the coaddTempExps with the statistic specified.
        Set the edge bits the coadd mask based on the weight map.

        Parameters
        ----------
        skyInfo : `lsst.pipe.base.Struct`
            Struct with geometric information about the patch.
        tempExpRefList : `list`
            List of data references to Warps (previously called CoaddTempExps).
        imageScalerList : `list`
            List of image scalers.
        weightList : `list`
            List of weights
        altMaskList : `list`, optional
            List of alternate masks to use rather than those stored with
            tempExp.
        mask : `int`, optional
            Bit mask value to exclude from coaddition.
        supplementaryData : lsst.pipe.base.Struct, optional
            Struct with additional data products needed to assemble coadd.
            Only used by subclasses that implement `makeSupplementaryData`
            and override `run`.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
           Result struct with components:

           - ``coaddExposure``: coadded exposure (``lsst.afw.image.Exposure``).
           - ``nImage``: exposure count image (``lsst.afw.image.Image``), if requested.
           - ``inputMap``: bit-wise map of inputs, if requested.
           - ``warpRefList``: input list of refs to the warps (
                              ``lsst.daf.butler.DeferredDatasetHandle`` or
                              ``lsst.daf.persistence.ButlerDataRef``)
                              (unmodified)
           - ``imageScalerList``: input list of image scalers (unmodified)
           - ``weightList``: input list of weights (unmodified)
        """
        tempExpName = self.getTempExpDatasetName(self.warpType)
        self.log.info("Assembling %s %s", len(tempExpRefList), tempExpName)
        stats = self.prepareStats(mask=mask)

        if altMaskList is None:
            altMaskList = [None]*len(tempExpRefList)

        coaddExposure = afwImage.ExposureF(skyInfo.bbox, skyInfo.wcs)
        coaddExposure.setPhotoCalib(self.scaleZeroPoint.getPhotoCalib())
        coaddExposure.getInfo().setCoaddInputs(self.inputRecorder.makeCoaddInputs())
        self.assembleMetadata(coaddExposure, tempExpRefList, weightList)
        coaddMaskedImage = coaddExposure.getMaskedImage()
        subregionSizeArr = self.config.subregionSize
        subregionSize = geom.Extent2I(subregionSizeArr[0], subregionSizeArr[1])
        # if nImage is requested, create a zero one which can be passed to assembleSubregion
        if self.config.doNImage:
            nImage = afwImage.ImageU(skyInfo.bbox)
        else:
            nImage = None
        # If inputMap is requested, create the initial version that can be masked in
        # assembleSubregion.
        if self.config.doInputMap:
            self.inputMapper.buildCcdInputMap(skyInfo.bbox,
                                              skyInfo.wcs,
                                              coaddExposure.getInfo().getCoaddInputs().ccds)

        for subBBox in self._subBBoxIter(skyInfo.bbox, subregionSize):
            try:
                self.assembleSubregion(coaddExposure, subBBox, tempExpRefList, imageScalerList,
                                       weightList, altMaskList, stats.flags, stats.ctrl,
                                       nImage=nImage)
            except Exception as e:
                self.log.fatal("Cannot compute coadd %s: %s", subBBox, e)

        # If inputMap is requested, we must finalize the map after the accumulation.
        if self.config.doInputMap:
            self.inputMapper.finalizeCcdInputMapMask()
            inputMap = self.inputMapper.ccdInputMap
        else:
            inputMap = None

        self.setInexactPsf(coaddMaskedImage.getMask())
        # Despite the name, the following doesn't really deal with "EDGE" pixels: it identifies
        # pixels that didn't receive any unmasked inputs (as occurs around the edge of the field).
        coaddUtils.setCoaddEdgeBits(coaddMaskedImage.getMask(), coaddMaskedImage.getVariance())
        return pipeBase.Struct(coaddExposure=coaddExposure, nImage=nImage,
                               warpRefList=tempExpRefList, imageScalerList=imageScalerList,
                               weightList=weightList, inputMap=inputMap)

    def assembleMetadata(self, coaddExposure, tempExpRefList, weightList):
        """Set the metadata for the coadd.

        This basic implementation sets the filter from the first input.

        Parameters
        ----------
        coaddExposure : `lsst.afw.image.Exposure`
            The target exposure for the coadd.
        tempExpRefList : `list`
            List of data references to tempExp.
        weightList : `list`
            List of weights.
        """
        assert len(tempExpRefList) == len(weightList), "Length mismatch"
        tempExpName = self.getTempExpDatasetName(self.warpType)
        # We load a single pixel of each coaddTempExp, because we just want to get at the metadata
        # (and we need more than just the PropertySet that contains the header), which is not possible
        # with the current butler (see #2777).
        bbox = geom.Box2I(coaddExposure.getBBox().getMin(), geom.Extent2I(1, 1))

        if isinstance(tempExpRefList[0], DeferredDatasetHandle):
            # Gen 3 API
            tempExpList = [tempExpRef.get(parameters={'bbox': bbox}) for tempExpRef in tempExpRefList]
        else:
            # Gen 2 API. Delete this when Gen 2 retired
            tempExpList = [tempExpRef.get(tempExpName + "_sub", bbox=bbox, immediate=True)
                           for tempExpRef in tempExpRefList]
        numCcds = sum(len(tempExp.getInfo().getCoaddInputs().ccds) for tempExp in tempExpList)

        # Set the coadd FilterLabel to the band of the first input exposure:
        # Coadds are calibrated, so the physical label is now meaningless.
        coaddExposure.setFilterLabel(afwImage.FilterLabel(tempExpList[0].getFilterLabel().bandLabel))
        coaddInputs = coaddExposure.getInfo().getCoaddInputs()
        coaddInputs.ccds.reserve(numCcds)
        coaddInputs.visits.reserve(len(tempExpList))

        for tempExp, weight in zip(tempExpList, weightList):
            self.inputRecorder.addVisitToCoadd(coaddInputs, tempExp, weight)

        if self.config.doUsePsfMatchedPolygons:
            self.shrinkValidPolygons(coaddInputs)

        coaddInputs.visits.sort()
        if self.warpType == "psfMatched":
            # The modelPsf BBox for a psfMatchedWarp/coaddTempExp was dynamically defined by
            # ModelPsfMatchTask as the square box bounding its spatially-variable, pre-matched WarpedPsf.
            # Likewise, set the PSF of a PSF-Matched Coadd to the modelPsf
            # having the maximum width (sufficient because square)
            modelPsfList = [tempExp.getPsf() for tempExp in tempExpList]
            modelPsfWidthList = [modelPsf.computeBBox().getWidth() for modelPsf in modelPsfList]
            psf = modelPsfList[modelPsfWidthList.index(max(modelPsfWidthList))]
        else:
            psf = measAlg.CoaddPsf(coaddInputs.ccds, coaddExposure.getWcs(),
                                   self.config.coaddPsf.makeControl())
        coaddExposure.setPsf(psf)
        apCorrMap = measAlg.makeCoaddApCorrMap(coaddInputs.ccds, coaddExposure.getBBox(afwImage.PARENT),
                                               coaddExposure.getWcs())
        coaddExposure.getInfo().setApCorrMap(apCorrMap)
        if self.config.doAttachTransmissionCurve:
            transmissionCurve = measAlg.makeCoaddTransmissionCurve(coaddExposure.getWcs(), coaddInputs.ccds)
            coaddExposure.getInfo().setTransmissionCurve(transmissionCurve)

    def assembleSubregion(self, coaddExposure, bbox, tempExpRefList, imageScalerList, weightList,
                          altMaskList, statsFlags, statsCtrl, nImage=None):
        """Assemble the coadd for a sub-region.

        For each coaddTempExp, check for (and swap in) an alternative mask
        if one is passed. Remove mask planes listed in
        `config.removeMaskPlanes`. Finally, stack the actual exposures using
        `lsst.afw.math.statisticsStack` with the statistic specified by
        statsFlags. Typically, the statsFlag will be one of lsst.afw.math.MEAN for
        a mean-stack or `lsst.afw.math.MEANCLIP` for outlier rejection using
        an N-sigma clipped mean where N and iterations are specified by
        statsCtrl.  Assign the stacked subregion back to the coadd.

        Parameters
        ----------
        coaddExposure : `lsst.afw.image.Exposure`
            The target exposure for the coadd.
        bbox : `lsst.geom.Box`
            Sub-region to coadd.
        tempExpRefList : `list`
            List of data reference to tempExp.
        imageScalerList : `list`
            List of image scalers.
        weightList : `list`
            List of weights.
        altMaskList : `list`
            List of alternate masks to use rather than those stored with
            tempExp, or None.  Each element is dict with keys = mask plane
            name to which to add the spans.
        statsFlags : `lsst.afw.math.Property`
            Property object for statistic for coadd.
        statsCtrl : `lsst.afw.math.StatisticsControl`
            Statistics control object for coadd.
        nImage : `lsst.afw.image.ImageU`, optional
            Keeps track of exposure count for each pixel.
        """
        self.log.debug("Computing coadd over %s", bbox)
        tempExpName = self.getTempExpDatasetName(self.warpType)
        coaddExposure.mask.addMaskPlane("REJECTED")
        coaddExposure.mask.addMaskPlane("CLIPPED")
        coaddExposure.mask.addMaskPlane("SENSOR_EDGE")
        maskMap = self.setRejectedMaskMapping(statsCtrl)
        clipped = afwImage.Mask.getPlaneBitMask("CLIPPED")
        maskedImageList = []
        if nImage is not None:
            subNImage = afwImage.ImageU(bbox.getWidth(), bbox.getHeight())
        for tempExpRef, imageScaler, altMask in zip(tempExpRefList, imageScalerList, altMaskList):

            if isinstance(tempExpRef, DeferredDatasetHandle):
                # Gen 3 API
                exposure = tempExpRef.get(parameters={'bbox': bbox})
            else:
                # Gen 2 API. Delete this when Gen 2 retired
                exposure = tempExpRef.get(tempExpName + "_sub", bbox=bbox)

            maskedImage = exposure.getMaskedImage()
            mask = maskedImage.getMask()
            if altMask is not None:
                self.applyAltMaskPlanes(mask, altMask)
            imageScaler.scaleMaskedImage(maskedImage)

            # Add 1 for each pixel which is not excluded by the exclude mask.
            # In legacyCoadd, pixels may also be excluded by afwMath.statisticsStack.
            if nImage is not None:
                subNImage.getArray()[maskedImage.getMask().getArray() & statsCtrl.getAndMask() == 0] += 1
            if self.config.removeMaskPlanes:
                self.removeMaskPlanes(maskedImage)
            maskedImageList.append(maskedImage)

            if self.config.doInputMap:
                visit = exposure.getInfo().getCoaddInputs().visits[0].getId()
                self.inputMapper.maskWarpBBox(bbox, visit, mask, statsCtrl.getAndMask())

        with self.timer("stack"):
            coaddSubregion = afwMath.statisticsStack(maskedImageList, statsFlags, statsCtrl, weightList,
                                                     clipped,  # also set output to CLIPPED if sigma-clipped
                                                     maskMap)
        coaddExposure.maskedImage.assign(coaddSubregion, bbox)
        if nImage is not None:
            nImage.assign(subNImage, bbox)

    def removeMaskPlanes(self, maskedImage):
        """Unset the mask of an image for mask planes specified in the config.

        Parameters
        ----------
        maskedImage : `lsst.afw.image.MaskedImage`
            The masked image to be modified.
        """
        mask = maskedImage.getMask()
        for maskPlane in self.config.removeMaskPlanes:
            try:
                mask &= ~mask.getPlaneBitMask(maskPlane)
            except pexExceptions.InvalidParameterError:
                self.log.debug("Unable to remove mask plane %s: no mask plane with that name was found.",
                               maskPlane)

    @staticmethod
    def setRejectedMaskMapping(statsCtrl):
        """Map certain mask planes of the warps to new planes for the coadd.

        If a pixel is rejected due to a mask value other than EDGE, NO_DATA,
        or CLIPPED, set it to REJECTED on the coadd.
        If a pixel is rejected due to EDGE, set the coadd pixel to SENSOR_EDGE.
        If a pixel is rejected due to CLIPPED, set the coadd pixel to CLIPPED.

        Parameters
        ----------
        statsCtrl : `lsst.afw.math.StatisticsControl`
            Statistics control object for coadd

        Returns
        -------
        maskMap : `list` of `tuple` of `int`
            A list of mappings of mask planes of the warped exposures to
            mask planes of the coadd.
        """
        edge = afwImage.Mask.getPlaneBitMask("EDGE")
        noData = afwImage.Mask.getPlaneBitMask("NO_DATA")
        clipped = afwImage.Mask.getPlaneBitMask("CLIPPED")
        toReject = statsCtrl.getAndMask() & (~noData) & (~edge) & (~clipped)
        maskMap = [(toReject, afwImage.Mask.getPlaneBitMask("REJECTED")),
                   (edge, afwImage.Mask.getPlaneBitMask("SENSOR_EDGE")),
                   (clipped, clipped)]
        return maskMap

    def applyAltMaskPlanes(self, mask, altMaskSpans):
        """Apply in place alt mask formatted as SpanSets to a mask.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            Original mask.
        altMaskSpans : `dict`
            SpanSet lists to apply. Each element contains the new mask
            plane name (e.g. "CLIPPED and/or "NO_DATA") as the key,
            and list of SpanSets to apply to the mask.

        Returns
        -------
        mask : `lsst.afw.image.Mask`
            Updated mask.
        """
        if self.config.doUsePsfMatchedPolygons:
            if ("NO_DATA" in altMaskSpans) and ("NO_DATA" in self.config.badMaskPlanes):
                # Clear away any other masks outside the validPolygons. These pixels are no longer
                # contributing to inexact PSFs, and will still be rejected because of NO_DATA
                # self.config.doUsePsfMatchedPolygons should be True only in CompareWarpAssemble
                # This mask-clearing step must only occur *before* applying the new masks below
                for spanSet in altMaskSpans['NO_DATA']:
                    spanSet.clippedTo(mask.getBBox()).clearMask(mask, self.getBadPixelMask())

        for plane, spanSetList in altMaskSpans.items():
            maskClipValue = mask.addMaskPlane(plane)
            for spanSet in spanSetList:
                spanSet.clippedTo(mask.getBBox()).setMask(mask, 2**maskClipValue)
        return mask

    def shrinkValidPolygons(self, coaddInputs):
        """Shrink coaddInputs' ccds' ValidPolygons in place.

        Either modify each ccd's validPolygon in place, or if CoaddInputs
        does not have a validPolygon, create one from its bbox.

        Parameters
        ----------
        coaddInputs : `lsst.afw.image.coaddInputs`
            Original mask.

        """
        for ccd in coaddInputs.ccds:
            polyOrig = ccd.getValidPolygon()
            validPolyBBox = polyOrig.getBBox() if polyOrig else ccd.getBBox()
            validPolyBBox.grow(-self.config.matchingKernelSize//2)
            if polyOrig:
                validPolygon = polyOrig.intersectionSingle(validPolyBBox)
            else:
                validPolygon = afwGeom.polygon.Polygon(geom.Box2D(validPolyBBox))
            ccd.setValidPolygon(validPolygon)

    def readBrightObjectMasks(self, dataRef):
        """Retrieve the bright object masks.

        Returns None on failure.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.butlerSubset.ButlerDataRef`
            A Butler dataRef.

        Returns
        -------
        result : `lsst.daf.persistence.butlerSubset.ButlerDataRef`
            Bright object mask from the Butler object, or None if it cannot
            be retrieved.
        """
        try:
            return dataRef.get(datasetType="brightObjectMask", immediate=True)
        except Exception as e:
            self.log.warn("Unable to read brightObjectMask for %s: %s", dataRef.dataId, e)
            return None

    def setBrightObjectMasks(self, exposure, brightObjectMasks, dataId=None):
        """Set the bright object masks.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure under consideration.
        dataId : `lsst.daf.persistence.dataId`
            Data identifier dict for patch.
        brightObjectMasks : `lsst.afw.table`
            Table of bright objects to mask.
        """

        if brightObjectMasks is None:
            self.log.warn("Unable to apply bright object mask: none supplied")
            return
        self.log.info("Applying %d bright object masks to %s", len(brightObjectMasks), dataId)
        mask = exposure.getMaskedImage().getMask()
        wcs = exposure.getWcs()
        plateScale = wcs.getPixelScale().asArcseconds()

        for rec in brightObjectMasks:
            center = geom.PointI(wcs.skyToPixel(rec.getCoord()))
            if rec["type"] == "box":
                assert rec["angle"] == 0.0, ("Angle != 0 for mask object %s" % rec["id"])
                width = rec["width"].asArcseconds()/plateScale    # convert to pixels
                height = rec["height"].asArcseconds()/plateScale  # convert to pixels

                halfSize = geom.ExtentI(0.5*width, 0.5*height)
                bbox = geom.Box2I(center - halfSize, center + halfSize)

                bbox = geom.BoxI(geom.PointI(int(center[0] - 0.5*width), int(center[1] - 0.5*height)),
                                 geom.PointI(int(center[0] + 0.5*width), int(center[1] + 0.5*height)))
                spans = afwGeom.SpanSet(bbox)
            elif rec["type"] == "circle":
                radius = int(rec["radius"].asArcseconds()/plateScale)   # convert to pixels
                spans = afwGeom.SpanSet.fromShape(radius, offset=center)
            else:
                self.log.warn("Unexpected region type %s at %s" % rec["type"], center)
                continue
            spans.clippedTo(mask.getBBox()).setMask(mask, self.brightObjectBitmask)

    def setInexactPsf(self, mask):
        """Set INEXACT_PSF mask plane.

        If any of the input images isn't represented in the coadd (due to
        clipped pixels or chip gaps), the `CoaddPsf` will be inexact. Flag
        these pixels.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            Coadded exposure's mask, modified in-place.
        """
        mask.addMaskPlane("INEXACT_PSF")
        inexactPsf = mask.getPlaneBitMask("INEXACT_PSF")
        sensorEdge = mask.getPlaneBitMask("SENSOR_EDGE")  # chip edges (so PSF is discontinuous)
        clipped = mask.getPlaneBitMask("CLIPPED")  # pixels clipped from coadd
        rejected = mask.getPlaneBitMask("REJECTED")  # pixels rejected from coadd due to masks
        array = mask.getArray()
        selected = array & (sensorEdge | clipped | rejected) > 0
        array[selected] |= inexactPsf

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser.
        """
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", cls.ConfigClass().coaddName + "Coadd_"
                               + cls.ConfigClass().warpType + "Warp",
                               help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=AssembleCoaddDataIdContainer)
        parser.add_id_argument("--selectId", "calexp", help="data ID, e.g. --selectId visit=6789 ccd=0..9",
                               ContainerClass=SelectDataIdContainer)
        return parser

    @staticmethod
    def _subBBoxIter(bbox, subregionSize):
        """Iterate over subregions of a bbox.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I`
            Bounding box over which to iterate.
        subregionSize: `lsst.geom.Extent2I`
            Size of sub-bboxes.

        Yields
        ------
        subBBox : `lsst.geom.Box2I`
            Next sub-bounding box of size ``subregionSize`` or smaller; each ``subBBox``
            is contained within ``bbox``, so it may be smaller than ``subregionSize`` at
            the edges of ``bbox``, but it will never be empty.
        """
        if bbox.isEmpty():
            raise RuntimeError("bbox %s is empty" % (bbox,))
        if subregionSize[0] < 1 or subregionSize[1] < 1:
            raise RuntimeError("subregionSize %s must be nonzero" % (subregionSize,))

        for rowShift in range(0, bbox.getHeight(), subregionSize[1]):
            for colShift in range(0, bbox.getWidth(), subregionSize[0]):
                subBBox = geom.Box2I(bbox.getMin() + geom.Extent2I(colShift, rowShift), subregionSize)
                subBBox.clip(bbox)
                if subBBox.isEmpty():
                    raise RuntimeError("Bug: empty bbox! bbox=%s, subregionSize=%s, "
                                       "colShift=%s, rowShift=%s" %
                                       (bbox, subregionSize, colShift, rowShift))
                yield subBBox

    def filterWarps(self, inputs, goodVisits):
        """Return list of only inputRefs with visitId in goodVisits ordered by goodVisit

        Parameters
        ----------
        inputs : list
            List of `lsst.pipe.base.connections.DeferredDatasetRef` with dataId containing visit
        goodVisit : `dict`
            Dictionary with good visitIds as the keys. Value ignored.

        Returns:
        --------
        filteredInputs : `list`
            Filtered and sorted list of `lsst.pipe.base.connections.DeferredDatasetRef`
        """
        inputWarpDict = {inputRef.ref.dataId['visit']: inputRef for inputRef in inputs}
        filteredInputs = []
        for visit in goodVisits.keys():
            filteredInputs.append(inputWarpDict[visit])
        return filteredInputs


class AssembleCoaddDataIdContainer(pipeBase.DataIdContainer):
    """A version of `lsst.pipe.base.DataIdContainer` specialized for assembleCoadd.
    """

    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList.

        Parameters
        ----------
        namespace
            Results of parsing command-line (with ``butler`` and ``log`` elements).
        """
        datasetType = namespace.config.coaddName + "Coadd"
        keysCoadd = namespace.butler.getKeys(datasetType=datasetType, level=self.level)

        for dataId in self.idList:
            # tract and patch are required
            for key in keysCoadd:
                if key not in dataId:
                    raise RuntimeError("--id must include " + key)

            dataRef = namespace.butler.dataRef(
                datasetType=datasetType,
                dataId=dataId,
            )
            self.refList.append(dataRef)


def countMaskFromFootprint(mask, footprint, bitmask, ignoreMask):
    """Function to count the number of pixels with a specific mask in a
    footprint.

    Find the intersection of mask & footprint. Count all pixels in the mask
    that are in the intersection that have bitmask set but do not have
    ignoreMask set. Return the count.

    Parameters
    ----------
    mask : `lsst.afw.image.Mask`
        Mask to define intersection region by.
    footprint : `lsst.afw.detection.Footprint`
        Footprint to define the intersection region by.
    bitmask
        Specific mask that we wish to count the number of occurances of.
    ignoreMask
        Pixels to not consider.

    Returns
    -------
    result : `int`
        Count of number of pixels in footprint with specified mask.
    """
    bbox = footprint.getBBox()
    bbox.clip(mask.getBBox(afwImage.PARENT))
    fp = afwImage.Mask(bbox)
    subMask = mask.Factory(mask, bbox, afwImage.PARENT)
    footprint.spans.setMask(fp, bitmask)
    return numpy.logical_and((subMask.getArray() & fp.getArray()) > 0,
                             (subMask.getArray() & ignoreMask) == 0).sum()


class SafeClipAssembleCoaddConfig(AssembleCoaddConfig, pipelineConnections=AssembleCoaddConnections):
    """Configuration parameters for the SafeClipAssembleCoaddTask.
    """
    clipDetection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Detect sources on difference between unclipped and clipped coadd")
    minClipFootOverlap = pexConfig.Field(
        doc="Minimum fractional overlap of clipped footprint with visit DETECTED to be clipped",
        dtype=float,
        default=0.6
    )
    minClipFootOverlapSingle = pexConfig.Field(
        doc="Minimum fractional overlap of clipped footprint with visit DETECTED to be "
        "clipped when only one visit overlaps",
        dtype=float,
        default=0.5
    )
    minClipFootOverlapDouble = pexConfig.Field(
        doc="Minimum fractional overlap of clipped footprints with visit DETECTED to be "
        "clipped when two visits overlap",
        dtype=float,
        default=0.45
    )
    maxClipFootOverlapDouble = pexConfig.Field(
        doc="Maximum fractional overlap of clipped footprints with visit DETECTED when "
        "considering two visits",
        dtype=float,
        default=0.15
    )
    minBigOverlap = pexConfig.Field(
        doc="Minimum number of pixels in footprint to use DETECTED mask from the single visits "
        "when labeling clipped footprints",
        dtype=int,
        default=100
    )

    def setDefaults(self):
        """Set default values for clipDetection.

        Notes
        -----
        The numeric values for these configuration parameters were
        empirically determined, future work may further refine them.
        """
        AssembleCoaddConfig.setDefaults(self)
        self.clipDetection.doTempLocalBackground = False
        self.clipDetection.reEstimateBackground = False
        self.clipDetection.returnOriginalFootprints = False
        self.clipDetection.thresholdPolarity = "both"
        self.clipDetection.thresholdValue = 2
        self.clipDetection.nSigmaToGrow = 2
        self.clipDetection.minPixels = 4
        self.clipDetection.isotropicGrow = True
        self.clipDetection.thresholdType = "pixel_stdev"
        self.sigmaClip = 1.5
        self.clipIter = 3
        self.statistic = "MEAN"

    def validate(self):
        if self.doSigmaClip:
            log.warn("Additional Sigma-clipping not allowed in Safe-clipped Coadds. "
                     "Ignoring doSigmaClip.")
            self.doSigmaClip = False
        if self.statistic != "MEAN":
            raise ValueError("Only MEAN statistic allowed for final stacking in SafeClipAssembleCoadd "
                             "(%s chosen). Please set statistic to MEAN."
                             % (self.statistic))
        AssembleCoaddTask.ConfigClass.validate(self)


class SafeClipAssembleCoaddTask(AssembleCoaddTask):
    """Assemble a coadded image from a set of coadded temporary exposures,
    being careful to clip & flag areas with potential artifacts.

    In ``AssembleCoaddTask``, we compute the coadd as an clipped mean (i.e.,
    we clip outliers). The problem with doing this is that when computing the
    coadd PSF at a given location, individual visit PSFs from visits with
    outlier pixels contribute to the coadd PSF and cannot be treated correctly.
    In this task, we correct for this behavior by creating a new
    ``badMaskPlane`` 'CLIPPED'. We populate this plane on the input
    coaddTempExps and the final coadd where

        i. difference imaging suggests that there is an outlier and
        ii. this outlier appears on only one or two images.

    Such regions will not contribute to the final coadd. Furthermore, any
    routine to determine the coadd PSF can now be cognizant of clipped regions.
    Note that the algorithm implemented by this task is preliminary and works
    correctly for HSC data. Parameter modifications and or considerable
    redesigning of the algorithm is likley required for other surveys.

    ``SafeClipAssembleCoaddTask`` uses a ``SourceDetectionTask``
    "clipDetection" subtask and also sub-classes ``AssembleCoaddTask``.
    You can retarget the ``SourceDetectionTask`` "clipDetection" subtask
    if you wish.

    Notes
    -----
    The `lsst.pipe.base.cmdLineTask.CmdLineTask` interface supports a
    flag ``-d`` to import ``debug.py`` from your ``PYTHONPATH``;
    see `baseDebug` for more about ``debug.py`` files.
    `SafeClipAssembleCoaddTask` has no debug variables of its own.
    The ``SourceDetectionTask`` "clipDetection" subtasks may support debug
    variables. See the documetation for `SourceDetectionTask` "clipDetection"
    for further information.

    Examples
    --------
    `SafeClipAssembleCoaddTask` assembles a set of warped ``coaddTempExp``
    images into a coadded image. The `SafeClipAssembleCoaddTask` is invoked by
    running assembleCoadd.py *without* the flag '--legacyCoadd'.

    Usage of ``assembleCoadd.py`` expects a data reference to the tract patch
    and filter to be coadded (specified using
    '--id = [KEY=VALUE1[^VALUE2[^VALUE3...] [KEY=VALUE1[^VALUE2[^VALUE3...] ...]]')
    along with a list of coaddTempExps to attempt to coadd (specified using
    '--selectId [KEY=VALUE1[^VALUE2[^VALUE3...] [KEY=VALUE1[^VALUE2[^VALUE3...] ...]]').
    Only the coaddTempExps that cover the specified tract and patch will be
    coadded. A list of the available optional arguments can be obtained by
    calling assembleCoadd.py with the --help command line argument:

    .. code-block:: none

       assembleCoadd.py --help

    To demonstrate usage of the `SafeClipAssembleCoaddTask` in the larger
    context of multi-band processing, we will generate the HSC-I & -R band
    coadds from HSC engineering test data provided in the ci_hsc package.
    To begin, assuming that the lsst stack has been already set up, we must
    set up the obs_subaru and ci_hsc packages. This defines the environment
    variable $CI_HSC_DIR and points at the location of the package. The raw
    HSC data live in the ``$CI_HSC_DIR/raw`` directory. To begin assembling
    the coadds, we must first

    - ``processCcd``
        process the individual ccds in $CI_HSC_RAW to produce calibrated exposures
    - ``makeSkyMap``
        create a skymap that covers the area of the sky present in the raw exposures
    - ``makeCoaddTempExp``
        warp the individual calibrated exposures to the tangent plane of the coadd</DD>

    We can perform all of these steps by running

    .. code-block:: none

       $CI_HSC_DIR scons warp-903986 warp-904014 warp-903990 warp-904010 warp-903988

    This will produce warped coaddTempExps for each visit. To coadd the
    warped data, we call ``assembleCoadd.py`` as follows:

    .. code-block:: none

       assembleCoadd.py $CI_HSC_DIR/DATA --id patch=5,4 tract=0 filter=HSC-I \
       --selectId visit=903986 ccd=16 --selectId visit=903986 ccd=22 --selectId visit=903986 ccd=23 \
       --selectId visit=903986 ccd=100--selectId visit=904014 ccd=1 --selectId visit=904014 ccd=6 \
       --selectId visit=904014 ccd=12 --selectId visit=903990 ccd=18 --selectId visit=903990 ccd=25 \
       --selectId visit=904010 ccd=4 --selectId visit=904010 ccd=10 --selectId visit=904010 ccd=100 \
       --selectId visit=903988 ccd=16 --selectId visit=903988 ccd=17 --selectId visit=903988 ccd=23 \
       --selectId visit=903988 ccd=24

    This will process the HSC-I band data. The results are written in
    ``$CI_HSC_DIR/DATA/deepCoadd-results/HSC-I``.

    You may also choose to run:

    .. code-block:: none

       scons warp-903334 warp-903336 warp-903338 warp-903342 warp-903344 warp-903346 nnn
       assembleCoadd.py $CI_HSC_DIR/DATA --id patch=5,4 tract=0 filter=HSC-R --selectId visit=903334 ccd=16 \
       --selectId visit=903334 ccd=22 --selectId visit=903334 ccd=23 --selectId visit=903334 ccd=100 \
       --selectId visit=903336 ccd=17 --selectId visit=903336 ccd=24 --selectId visit=903338 ccd=18 \
       --selectId visit=903338 ccd=25 --selectId visit=903342 ccd=4 --selectId visit=903342 ccd=10 \
       --selectId visit=903342 ccd=100 --selectId visit=903344 ccd=0 --selectId visit=903344 ccd=5 \
       --selectId visit=903344 ccd=11 --selectId visit=903346 ccd=1 --selectId visit=903346 ccd=6 \
       --selectId visit=903346 ccd=12

    to generate the coadd for the HSC-R band if you are interested in following
    multiBand Coadd processing as discussed in ``pipeTasks_multiBand``.
    """
    ConfigClass = SafeClipAssembleCoaddConfig
    _DefaultName = "safeClipAssembleCoadd"

    def __init__(self, *args, **kwargs):
        AssembleCoaddTask.__init__(self, *args, **kwargs)
        schema = afwTable.SourceTable.makeMinimalSchema()
        self.makeSubtask("clipDetection", schema=schema)

    @utils.inheritDoc(AssembleCoaddTask)
    @pipeBase.timeMethod
    def run(self, skyInfo, tempExpRefList, imageScalerList, weightList, *args, **kwargs):
        """Assemble the coadd for a region.

        Compute the difference of coadds created with and without outlier
        rejection to identify coadd pixels that have outlier values in some
        individual visits.
        Detect clipped regions on the difference image and mark these regions
        on the one or two individual coaddTempExps where they occur if there
        is significant overlap between the clipped region and a source. This
        leaves us with a set of footprints from the difference image that have
        been identified as having occured on just one or two individual visits.
        However, these footprints were generated from a difference image. It
        is conceivable for a large diffuse source to have become broken up
        into multiple footprints acrosss the coadd difference in this process.
        Determine the clipped region from all overlapping footprints from the
        detected sources in each visit - these are big footprints.
        Combine the small and big clipped footprints and mark them on a new
        bad mask plane.
        Generate the coadd using `AssembleCoaddTask.run` without outlier
        removal. Clipped footprints will no longer make it into the coadd
        because they are marked in the new bad mask plane.

        Notes
        -----
        args and kwargs are passed but ignored in order to match the call
        signature expected by the parent task.
        """
        exp = self.buildDifferenceImage(skyInfo, tempExpRefList, imageScalerList, weightList)
        mask = exp.getMaskedImage().getMask()
        mask.addMaskPlane("CLIPPED")

        result = self.detectClip(exp, tempExpRefList)

        self.log.info('Found %d clipped objects', len(result.clipFootprints))

        maskClipValue = mask.getPlaneBitMask("CLIPPED")
        maskDetValue = mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE")
        # Append big footprints from individual Warps to result.clipSpans
        bigFootprints = self.detectClipBig(result.clipSpans, result.clipFootprints, result.clipIndices,
                                           result.detectionFootprints, maskClipValue, maskDetValue,
                                           exp.getBBox())
        # Create mask of the current clipped footprints
        maskClip = mask.Factory(mask.getBBox(afwImage.PARENT))
        afwDet.setMaskFromFootprintList(maskClip, result.clipFootprints, maskClipValue)

        maskClipBig = maskClip.Factory(mask.getBBox(afwImage.PARENT))
        afwDet.setMaskFromFootprintList(maskClipBig, bigFootprints, maskClipValue)
        maskClip |= maskClipBig

        # Assemble coadd from base class, but ignoring CLIPPED pixels
        badMaskPlanes = self.config.badMaskPlanes[:]
        badMaskPlanes.append("CLIPPED")
        badPixelMask = afwImage.Mask.getPlaneBitMask(badMaskPlanes)
        return AssembleCoaddTask.run(self, skyInfo, tempExpRefList, imageScalerList, weightList,
                                     result.clipSpans, mask=badPixelMask)

    def buildDifferenceImage(self, skyInfo, tempExpRefList, imageScalerList, weightList):
        """Return an exposure that contains the difference between unclipped
        and clipped coadds.

        Generate a difference image between clipped and unclipped coadds.
        Compute the difference image by subtracting an outlier-clipped coadd
        from an outlier-unclipped coadd. Return the difference image.

        Parameters
        ----------
        skyInfo : `lsst.pipe.base.Struct`
            Patch geometry information, from getSkyInfo
        tempExpRefList : `list`
            List of data reference to tempExp
        imageScalerList : `list`
            List of image scalers
        weightList : `list`
            List of weights

        Returns
        -------
        exp : `lsst.afw.image.Exposure`
            Difference image of unclipped and clipped coadd wrapped in an Exposure
        """
        config = AssembleCoaddConfig()
        # getattr necessary because subtasks do not survive Config.toDict()
        # exclude connections because the class of self.config.connections is not
        # the same as AssembleCoaddConfig.connections, and the connections are not
        # needed to run this task anyway.
        configIntersection = {k: getattr(self.config, k)
                              for k, v in self.config.toDict().items()
                              if (k in config.keys() and k != "connections")}
        config.update(**configIntersection)

        # statistic MEAN copied from self.config.statistic, but for clarity explicitly assign
        config.statistic = 'MEAN'
        task = AssembleCoaddTask(config=config)
        coaddMean = task.run(skyInfo, tempExpRefList, imageScalerList, weightList).coaddExposure

        config.statistic = 'MEANCLIP'
        task = AssembleCoaddTask(config=config)
        coaddClip = task.run(skyInfo, tempExpRefList, imageScalerList, weightList).coaddExposure

        coaddDiff = coaddMean.getMaskedImage().Factory(coaddMean.getMaskedImage())
        coaddDiff -= coaddClip.getMaskedImage()
        exp = afwImage.ExposureF(coaddDiff)
        exp.setPsf(coaddMean.getPsf())
        return exp

    def detectClip(self, exp, tempExpRefList):
        """Detect clipped regions on an exposure and set the mask on the
        individual tempExp masks.

        Detect footprints in the difference image after smoothing the
        difference image with a Gaussian kernal. Identify footprints that
        overlap with one or two input ``coaddTempExps`` by comparing the
        computed overlap fraction to thresholds set in the config. A different
        threshold is applied depending on the number of overlapping visits
        (restricted to one or two). If the overlap exceeds the thresholds,
        the footprint is considered "CLIPPED" and is marked as such on the
        coaddTempExp. Return a struct with the clipped footprints, the indices
        of the ``coaddTempExps`` that end up overlapping with the clipped
        footprints, and a list of new masks for the ``coaddTempExps``.

        Parameters
        ----------
        exp : `lsst.afw.image.Exposure`
            Exposure to run detection on.
        tempExpRefList : `list`
            List of data reference to tempExp.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
           Result struct with components:

           - ``clipFootprints``: list of clipped footprints.
           - ``clipIndices``: indices for each ``clippedFootprint`` in
                ``tempExpRefList``.
           - ``clipSpans``: List of dictionaries containing spanSet lists
                to clip. Each element contains the new maskplane name
                ("CLIPPED") as the key and list of ``SpanSets`` as the value.
           - ``detectionFootprints``: List of DETECTED/DETECTED_NEGATIVE plane
                compressed into footprints.
        """
        mask = exp.getMaskedImage().getMask()
        maskDetValue = mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE")
        fpSet = self.clipDetection.detectFootprints(exp, doSmooth=True, clearMask=True)
        # Merge positive and negative together footprints together
        fpSet.positive.merge(fpSet.negative)
        footprints = fpSet.positive
        self.log.info('Found %d potential clipped objects', len(footprints.getFootprints()))
        ignoreMask = self.getBadPixelMask()

        clipFootprints = []
        clipIndices = []
        artifactSpanSets = [{'CLIPPED': list()} for _ in tempExpRefList]

        # for use by detectClipBig
        visitDetectionFootprints = []

        dims = [len(tempExpRefList), len(footprints.getFootprints())]
        overlapDetArr = numpy.zeros(dims, dtype=numpy.uint16)
        ignoreArr = numpy.zeros(dims, dtype=numpy.uint16)

        # Loop over masks once and extract/store only relevant overlap metrics and detection footprints
        for i, warpRef in enumerate(tempExpRefList):
            tmpExpMask = warpRef.get(datasetType=self.getTempExpDatasetName(self.warpType),
                                     immediate=True).getMaskedImage().getMask()
            maskVisitDet = tmpExpMask.Factory(tmpExpMask, tmpExpMask.getBBox(afwImage.PARENT),
                                              afwImage.PARENT, True)
            maskVisitDet &= maskDetValue
            visitFootprints = afwDet.FootprintSet(maskVisitDet, afwDet.Threshold(1))
            visitDetectionFootprints.append(visitFootprints)

            for j, footprint in enumerate(footprints.getFootprints()):
                ignoreArr[i, j] = countMaskFromFootprint(tmpExpMask, footprint, ignoreMask, 0x0)
                overlapDetArr[i, j] = countMaskFromFootprint(tmpExpMask, footprint, maskDetValue, ignoreMask)

        # build a list of clipped spans for each visit
        for j, footprint in enumerate(footprints.getFootprints()):
            nPixel = footprint.getArea()
            overlap = []  # hold the overlap with each visit
            indexList = []  # index of visit in global list
            for i in range(len(tempExpRefList)):
                ignore = ignoreArr[i, j]
                overlapDet = overlapDetArr[i, j]
                totPixel = nPixel - ignore

                # If we have more bad pixels than detection skip
                if ignore > overlapDet or totPixel <= 0.5*nPixel or overlapDet == 0:
                    continue
                overlap.append(overlapDet/float(totPixel))
                indexList.append(i)

            overlap = numpy.array(overlap)
            if not len(overlap):
                continue

            keep = False   # Should this footprint be marked as clipped?
            keepIndex = []  # Which tempExps does the clipped footprint belong to

            # If footprint only has one overlap use a lower threshold
            if len(overlap) == 1:
                if overlap[0] > self.config.minClipFootOverlapSingle:
                    keep = True
                    keepIndex = [0]
            else:
                # This is the general case where only visit should be clipped
                clipIndex = numpy.where(overlap > self.config.minClipFootOverlap)[0]
                if len(clipIndex) == 1:
                    keep = True
                    keepIndex = [clipIndex[0]]

                # Test if there are clipped objects that overlap two different visits
                clipIndex = numpy.where(overlap > self.config.minClipFootOverlapDouble)[0]
                if len(clipIndex) == 2 and len(overlap) > 3:
                    clipIndexComp = numpy.where(overlap <= self.config.minClipFootOverlapDouble)[0]
                    if numpy.max(overlap[clipIndexComp]) <= self.config.maxClipFootOverlapDouble:
                        keep = True
                        keepIndex = clipIndex

            if not keep:
                continue

            for index in keepIndex:
                globalIndex = indexList[index]
                artifactSpanSets[globalIndex]['CLIPPED'].append(footprint.spans)

            clipIndices.append(numpy.array(indexList)[keepIndex])
            clipFootprints.append(footprint)

        return pipeBase.Struct(clipFootprints=clipFootprints, clipIndices=clipIndices,
                               clipSpans=artifactSpanSets, detectionFootprints=visitDetectionFootprints)

    def detectClipBig(self, clipList, clipFootprints, clipIndices, detectionFootprints,
                      maskClipValue, maskDetValue, coaddBBox):
        """Return individual warp footprints for large artifacts and append
        them to ``clipList`` in place.

        Identify big footprints composed of many sources in the coadd
        difference that may have originated in a large diffuse source in the
        coadd. We do this by indentifying all clipped footprints that overlap
        significantly with each source in all the coaddTempExps.

        Parameters
        ----------
        clipList : `list`
            List of alt mask SpanSets with clipping information. Modified.
        clipFootprints : `list`
            List of clipped footprints.
        clipIndices : `list`
            List of which entries in tempExpClipList each footprint belongs to.
        maskClipValue
            Mask value of clipped pixels.
        maskDetValue
            Mask value of detected pixels.
        coaddBBox : `lsst.geom.Box`
            BBox of the coadd and warps.

        Returns
        -------
        bigFootprintsCoadd : `list`
            List of big footprints
        """
        bigFootprintsCoadd = []
        ignoreMask = self.getBadPixelMask()
        for index, (clippedSpans, visitFootprints) in enumerate(zip(clipList, detectionFootprints)):
            maskVisitDet = afwImage.MaskX(coaddBBox, 0x0)
            for footprint in visitFootprints.getFootprints():
                footprint.spans.setMask(maskVisitDet, maskDetValue)

            # build a mask of clipped footprints that are in this visit
            clippedFootprintsVisit = []
            for foot, clipIndex in zip(clipFootprints, clipIndices):
                if index not in clipIndex:
                    continue
                clippedFootprintsVisit.append(foot)
            maskVisitClip = maskVisitDet.Factory(maskVisitDet.getBBox(afwImage.PARENT))
            afwDet.setMaskFromFootprintList(maskVisitClip, clippedFootprintsVisit, maskClipValue)

            bigFootprintsVisit = []
            for foot in visitFootprints.getFootprints():
                if foot.getArea() < self.config.minBigOverlap:
                    continue
                nCount = countMaskFromFootprint(maskVisitClip, foot, maskClipValue, ignoreMask)
                if nCount > self.config.minBigOverlap:
                    bigFootprintsVisit.append(foot)
                    bigFootprintsCoadd.append(foot)

            for footprint in bigFootprintsVisit:
                clippedSpans["CLIPPED"].append(footprint.spans)

        return bigFootprintsCoadd


class CompareWarpAssembleCoaddConnections(AssembleCoaddConnections):
    psfMatchedWarps = pipeBase.connectionTypes.Input(
        doc=("PSF-Matched Warps are required by CompareWarp regardless of the coadd type requested. "
             "Only PSF-Matched Warps make sense for image subtraction. "
             "Therefore, they must be an additional declared input."),
        name="{inputCoaddName}Coadd_psfMatchedWarp",
        storageClass="ExposureF",
        dimensions=("tract", "patch", "skymap", "visit"),
        deferLoad=True,
        multiple=True
    )
    templateCoadd = pipeBase.connectionTypes.Output(
        doc=("Model of the static sky, used to find temporal artifacts. Typically a PSF-Matched, "
             "sigma-clipped coadd. Written if and only if assembleStaticSkyModel.doWrite=True"),
        name="{fakesType}{outputCoaddName}CoaddPsfMatched",
        storageClass="ExposureF",
        dimensions=("tract", "patch", "skymap", "band"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        if not config.assembleStaticSkyModel.doWrite:
            self.outputs.remove("templateCoadd")
        config.validate()


class CompareWarpAssembleCoaddConfig(AssembleCoaddConfig,
                                     pipelineConnections=CompareWarpAssembleCoaddConnections):
    assembleStaticSkyModel = pexConfig.ConfigurableField(
        target=AssembleCoaddTask,
        doc="Task to assemble an artifact-free, PSF-matched Coadd to serve as a"
            " naive/first-iteration model of the static sky.",
    )
    detect = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Detect outlier sources on difference between each psfMatched warp and static sky model"
    )
    detectTemplate = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Detect sources on static sky model. Only used if doPreserveContainedBySource is True"
    )
    maskStreaks = pexConfig.ConfigurableField(
        target=MaskStreaksTask,
        doc="Detect streaks on difference between each psfMatched warp and static sky model. Only used if "
            "doFilterMorphological is True. Adds a mask plane to an exposure, with the mask plane name set by"
            "streakMaskName"
    )
    streakMaskName = pexConfig.Field(
        dtype=str,
        default="STREAK",
        doc="Name of mask bit used for streaks"
    )
    maxNumEpochs = pexConfig.Field(
        doc="Charactistic maximum local number of epochs/visits in which an artifact candidate can appear  "
            "and still be masked.  The effective maxNumEpochs is a broken linear function of local "
            "number of epochs (N): min(maxFractionEpochsLow*N, maxNumEpochs + maxFractionEpochsHigh*N). "
            "For each footprint detected on the image difference between the psfMatched warp and static sky "
            "model, if a significant fraction of pixels (defined by spatialThreshold) are residuals in more "
            "than the computed effective maxNumEpochs, the artifact candidate is deemed persistant rather "
            "than transient and not masked.",
        dtype=int,
        default=2
    )
    maxFractionEpochsLow = pexConfig.RangeField(
        doc="Fraction of local number of epochs (N) to use as effective maxNumEpochs for low N. "
            "Effective maxNumEpochs = "
            "min(maxFractionEpochsLow * N, maxNumEpochs + maxFractionEpochsHigh * N)",
        dtype=float,
        default=0.4,
        min=0., max=1.,
    )
    maxFractionEpochsHigh = pexConfig.RangeField(
        doc="Fraction of local number of epochs (N) to use as effective maxNumEpochs for high N. "
            "Effective maxNumEpochs = "
            "min(maxFractionEpochsLow * N, maxNumEpochs + maxFractionEpochsHigh * N)",
        dtype=float,
        default=0.03,
        min=0., max=1.,
    )
    spatialThreshold = pexConfig.RangeField(
        doc="Unitless fraction of pixels defining how much of the outlier region has to meet the "
            "temporal criteria. If 0, clip all. If 1, clip none.",
        dtype=float,
        default=0.5,
        min=0., max=1.,
        inclusiveMin=True, inclusiveMax=True
    )
    doScaleWarpVariance = pexConfig.Field(
        doc="Rescale Warp variance plane using empirical noise?",
        dtype=bool,
        default=True,
    )
    scaleWarpVariance = pexConfig.ConfigurableField(
        target=ScaleVarianceTask,
        doc="Rescale variance on warps",
    )
    doPreserveContainedBySource = pexConfig.Field(
        doc="Rescue artifacts from clipping that completely lie within a footprint detected"
            "on the PsfMatched Template Coadd. Replicates a behavior of SafeClip.",
        dtype=bool,
        default=True,
    )
    doPrefilterArtifacts = pexConfig.Field(
        doc="Ignore artifact candidates that are mostly covered by the bad pixel mask, "
            "because they will be excluded anyway. This prevents them from contributing "
            "to the outlier epoch count image and potentially being labeled as persistant."
            "'Mostly' is defined by the config 'prefilterArtifactsRatio'.",
        dtype=bool,
        default=True
    )
    prefilterArtifactsMaskPlanes = pexConfig.ListField(
        doc="Prefilter artifact candidates that are mostly covered by these bad mask planes.",
        dtype=str,
        default=('NO_DATA', 'BAD', 'SAT', 'SUSPECT'),
    )
    prefilterArtifactsRatio = pexConfig.Field(
        doc="Prefilter artifact candidates with less than this fraction overlapping good pixels",
        dtype=float,
        default=0.05
    )
    doFilterMorphological = pexConfig.Field(
        doc="Filter artifact candidates based on morphological criteria, i.g. those that appear to "
            "be streaks.",
        dtype=bool,
        default=False
    )

    def setDefaults(self):
        AssembleCoaddConfig.setDefaults(self)
        self.statistic = 'MEAN'
        self.doUsePsfMatchedPolygons = True

        # Real EDGE removed by psfMatched NO_DATA border half the width of the matching kernel
        # CompareWarp applies psfMatched EDGE pixels to directWarps before assembling
        if "EDGE" in self.badMaskPlanes:
            self.badMaskPlanes.remove('EDGE')
        self.removeMaskPlanes.append('EDGE')
        self.assembleStaticSkyModel.badMaskPlanes = ["NO_DATA", ]
        self.assembleStaticSkyModel.warpType = 'psfMatched'
        self.assembleStaticSkyModel.connections.warpType = 'psfMatched'
        self.assembleStaticSkyModel.statistic = 'MEANCLIP'
        self.assembleStaticSkyModel.sigmaClip = 2.5
        self.assembleStaticSkyModel.clipIter = 3
        self.assembleStaticSkyModel.calcErrorFromInputVariance = False
        self.assembleStaticSkyModel.doWrite = False
        self.detect.doTempLocalBackground = False
        self.detect.reEstimateBackground = False
        self.detect.returnOriginalFootprints = False
        self.detect.thresholdPolarity = "both"
        self.detect.thresholdValue = 5
        self.detect.minPixels = 4
        self.detect.isotropicGrow = True
        self.detect.thresholdType = "pixel_stdev"
        self.detect.nSigmaToGrow = 0.4
        # The default nSigmaToGrow for SourceDetectionTask is already 2.4,
        # Explicitly restating because ratio with detect.nSigmaToGrow matters
        self.detectTemplate.nSigmaToGrow = 2.4
        self.detectTemplate.doTempLocalBackground = False
        self.detectTemplate.reEstimateBackground = False
        self.detectTemplate.returnOriginalFootprints = False

    def validate(self):
        super().validate()
        if self.assembleStaticSkyModel.doNImage:
            raise ValueError("No dataset type exists for a PSF-Matched Template N Image."
                             "Please set assembleStaticSkyModel.doNImage=False")

        if self.assembleStaticSkyModel.doWrite and (self.warpType == self.assembleStaticSkyModel.warpType):
            raise ValueError("warpType (%s) == assembleStaticSkyModel.warpType (%s) and will compete for "
                             "the same dataset name. Please set assembleStaticSkyModel.doWrite to False "
                             "or warpType to 'direct'. assembleStaticSkyModel.warpType should ways be "
                             "'PsfMatched'" % (self.warpType, self.assembleStaticSkyModel.warpType))


class CompareWarpAssembleCoaddTask(AssembleCoaddTask):
    """Assemble a compareWarp coadded image from a set of warps
    by masking artifacts detected by comparing PSF-matched warps.

    In ``AssembleCoaddTask``, we compute the coadd as an clipped mean (i.e.,
    we clip outliers). The problem with doing this is that when computing the
    coadd PSF at a given location, individual visit PSFs from visits with
    outlier pixels contribute to the coadd PSF and cannot be treated correctly.
    In this task, we correct for this behavior by creating a new badMaskPlane
    'CLIPPED' which marks pixels in the individual warps suspected to contain
    an artifact. We populate this plane on the input warps by comparing
    PSF-matched warps with a PSF-matched median coadd which serves as a
    model of the static sky. Any group of pixels that deviates from the
    PSF-matched template coadd by more than config.detect.threshold sigma,
    is an artifact candidate. The candidates are then filtered to remove
    variable sources and sources that are difficult to subtract such as
    bright stars. This filter is configured using the config parameters
    ``temporalThreshold`` and ``spatialThreshold``. The temporalThreshold is
    the maximum fraction of epochs that the deviation can appear in and still
    be considered an artifact. The spatialThreshold is the maximum fraction of
    pixels in the footprint of the deviation that appear in other epochs
    (where other epochs is defined by the temporalThreshold). If the deviant
    region meets this criteria of having a significant percentage of pixels
    that deviate in only a few epochs, these pixels have the 'CLIPPED' bit
    set in the mask. These regions will not contribute to the final coadd.
    Furthermore, any routine to determine the coadd PSF can now be cognizant
    of clipped regions. Note that the algorithm implemented by this task is
    preliminary and works correctly for HSC data. Parameter modifications and
    or considerable redesigning of the algorithm is likley required for other
    surveys.

    ``CompareWarpAssembleCoaddTask`` sub-classes
    ``AssembleCoaddTask`` and instantiates ``AssembleCoaddTask``
    as a subtask to generate the TemplateCoadd (the model of the static sky).

    Notes
    -----
    The `lsst.pipe.base.cmdLineTask.CmdLineTask` interface supports a
    flag ``-d`` to import ``debug.py`` from your ``PYTHONPATH``; see
    ``baseDebug`` for more about ``debug.py`` files.

    This task supports the following debug variables:

    - ``saveCountIm``
        If True then save the Epoch Count Image as a fits file in the `figPath`
    - ``figPath``
        Path to save the debug fits images and figures

    For example, put something like:

    .. code-block:: python

       import lsstDebug
       def DebugInfo(name):
           di = lsstDebug.getInfo(name)
           if name == "lsst.pipe.tasks.assembleCoadd":
               di.saveCountIm = True
               di.figPath = "/desired/path/to/debugging/output/images"
           return di
       lsstDebug.Info = DebugInfo

    into your ``debug.py`` file and run ``assemebleCoadd.py`` with the
    ``--debug`` flag. Some subtasks may have their own debug variables;
    see individual Task documentation.

    Examples
    --------
    ``CompareWarpAssembleCoaddTask`` assembles a set of warped images into a
    coadded image. The ``CompareWarpAssembleCoaddTask`` is invoked by running
    ``assembleCoadd.py`` with the flag ``--compareWarpCoadd``.
    Usage of ``assembleCoadd.py`` expects a data reference to the tract patch
    and filter to be coadded (specified using
    '--id = [KEY=VALUE1[^VALUE2[^VALUE3...] [KEY=VALUE1[^VALUE2[^VALUE3...] ...]]')
    along with a list of coaddTempExps to attempt to coadd (specified using
    '--selectId [KEY=VALUE1[^VALUE2[^VALUE3...] [KEY=VALUE1[^VALUE2[^VALUE3...] ...]]').
    Only the warps that cover the specified tract and patch will be coadded.
    A list of the available optional arguments can be obtained by calling
    ``assembleCoadd.py`` with the ``--help`` command line argument:

    .. code-block:: none

       assembleCoadd.py --help

    To demonstrate usage of the ``CompareWarpAssembleCoaddTask`` in the larger
    context of multi-band processing, we will generate the HSC-I & -R band
    oadds from HSC engineering test data provided in the ``ci_hsc`` package.
    To begin, assuming that the lsst stack has been already set up, we must
    set up the ``obs_subaru`` and ``ci_hsc`` packages.
    This defines the environment variable ``$CI_HSC_DIR`` and points at the
    location of the package. The raw HSC data live in the ``$CI_HSC_DIR/raw``
    directory. To begin assembling the coadds, we must first

      - processCcd
        process the individual ccds in $CI_HSC_RAW to produce calibrated exposures
      - makeSkyMap
        create a skymap that covers the area of the sky present in the raw exposures
      - makeCoaddTempExp
        warp the individual calibrated exposures to the tangent plane of the coadd

    We can perform all of these steps by running

    .. code-block:: none

       $CI_HSC_DIR scons warp-903986 warp-904014 warp-903990 warp-904010 warp-903988

    This will produce warped ``coaddTempExps`` for each visit. To coadd the
    warped data, we call ``assembleCoadd.py`` as follows:

    .. code-block:: none

       assembleCoadd.py --compareWarpCoadd $CI_HSC_DIR/DATA --id patch=5,4 tract=0 filter=HSC-I \
       --selectId visit=903986 ccd=16 --selectId visit=903986 ccd=22 --selectId visit=903986 ccd=23 \
       --selectId visit=903986 ccd=100 --selectId visit=904014 ccd=1 --selectId visit=904014 ccd=6 \
       --selectId visit=904014 ccd=12 --selectId visit=903990 ccd=18 --selectId visit=903990 ccd=25 \
       --selectId visit=904010 ccd=4 --selectId visit=904010 ccd=10 --selectId visit=904010 ccd=100 \
       --selectId visit=903988 ccd=16 --selectId visit=903988 ccd=17 --selectId visit=903988 ccd=23 \
       --selectId visit=903988 ccd=24

    This will process the HSC-I band data. The results are written in
    ``$CI_HSC_DIR/DATA/deepCoadd-results/HSC-I``.
    """
    ConfigClass = CompareWarpAssembleCoaddConfig
    _DefaultName = "compareWarpAssembleCoadd"

    def __init__(self, *args, **kwargs):
        AssembleCoaddTask.__init__(self, *args, **kwargs)
        self.makeSubtask("assembleStaticSkyModel")
        detectionSchema = afwTable.SourceTable.makeMinimalSchema()
        self.makeSubtask("detect", schema=detectionSchema)
        if self.config.doPreserveContainedBySource:
            self.makeSubtask("detectTemplate", schema=afwTable.SourceTable.makeMinimalSchema())
        if self.config.doScaleWarpVariance:
            self.makeSubtask("scaleWarpVariance")
        if self.config.doFilterMorphological:
            self.makeSubtask("maskStreaks")

    @utils.inheritDoc(AssembleCoaddTask)
    def makeSupplementaryDataGen3(self, butlerQC, inputRefs, outputRefs):
        """
        Generate a templateCoadd to use as a naive model of static sky to
        subtract from PSF-Matched warps.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
           Result struct with components:

           - ``templateCoadd`` : coadded exposure (``lsst.afw.image.Exposure``)
           - ``nImage`` : N Image (``lsst.afw.image.Image``)
        """
        # Ensure that psfMatchedWarps are used as input warps for template generation
        staticSkyModelInputRefs = copy.deepcopy(inputRefs)
        staticSkyModelInputRefs.inputWarps = inputRefs.psfMatchedWarps

        # Because subtasks don't have connections we have to make one.
        # The main task's `templateCoadd` is the subtask's `coaddExposure`
        staticSkyModelOutputRefs = copy.deepcopy(outputRefs)
        if self.config.assembleStaticSkyModel.doWrite:
            staticSkyModelOutputRefs.coaddExposure = staticSkyModelOutputRefs.templateCoadd
            # Remove template coadd from both subtask's and main tasks outputs,
            # because it is handled by the subtask as `coaddExposure`
            del outputRefs.templateCoadd
            del staticSkyModelOutputRefs.templateCoadd

        # A PSF-Matched nImage does not exist as a dataset type
        if 'nImage' in staticSkyModelOutputRefs.keys():
            del staticSkyModelOutputRefs.nImage

        templateCoadd = self.assembleStaticSkyModel.runQuantum(butlerQC, staticSkyModelInputRefs,
                                                               staticSkyModelOutputRefs)
        if templateCoadd is None:
            raise RuntimeError(self._noTemplateMessage(self.assembleStaticSkyModel.warpType))

        return pipeBase.Struct(templateCoadd=templateCoadd.coaddExposure,
                               nImage=templateCoadd.nImage,
                               warpRefList=templateCoadd.warpRefList,
                               imageScalerList=templateCoadd.imageScalerList,
                               weightList=templateCoadd.weightList)

    @utils.inheritDoc(AssembleCoaddTask)
    def makeSupplementaryData(self, dataRef, selectDataList=None, warpRefList=None):
        """
        Generate a templateCoadd to use as a naive model of static sky to
        subtract from PSF-Matched warps.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
           Result struct with components:

           - ``templateCoadd``: coadded exposure (``lsst.afw.image.Exposure``)
           - ``nImage``: N Image (``lsst.afw.image.Image``)
        """
        templateCoadd = self.assembleStaticSkyModel.runDataRef(dataRef, selectDataList, warpRefList)
        if templateCoadd is None:
            raise RuntimeError(self._noTemplateMessage(self.assembleStaticSkyModel.warpType))

        return pipeBase.Struct(templateCoadd=templateCoadd.coaddExposure,
                               nImage=templateCoadd.nImage,
                               warpRefList=templateCoadd.warpRefList,
                               imageScalerList=templateCoadd.imageScalerList,
                               weightList=templateCoadd.weightList)

    def _noTemplateMessage(self, warpType):
        warpName = (warpType[0].upper() + warpType[1:])
        message = """No %(warpName)s warps were found to build the template coadd which is
            required to run CompareWarpAssembleCoaddTask. To continue assembling this type of coadd,
            first either rerun makeCoaddTempExp with config.make%(warpName)s=True or
            coaddDriver with config.makeCoadTempExp.make%(warpName)s=True, before assembleCoadd.

            Alternatively, to use another algorithm with existing warps, retarget the CoaddDriverConfig to
            another algorithm like:

                from lsst.pipe.tasks.assembleCoadd import SafeClipAssembleCoaddTask
                config.assemble.retarget(SafeClipAssembleCoaddTask)
        """ % {"warpName": warpName}
        return message

    @utils.inheritDoc(AssembleCoaddTask)
    @pipeBase.timeMethod
    def run(self, skyInfo, tempExpRefList, imageScalerList, weightList,
            supplementaryData, *args, **kwargs):
        """Assemble the coadd.

        Find artifacts and apply them to the warps' masks creating a list of
        alternative masks with a new "CLIPPED" plane and updated "NO_DATA"
        plane. Then pass these alternative masks to the base class's `run`
        method.

        The input parameters ``supplementaryData`` is a `lsst.pipe.base.Struct`
        that must contain a ``templateCoadd`` that serves as the
        model of the static sky.
        """

        # Check and match the order of the supplementaryData
        # (PSF-matched) inputs to the order of the direct inputs,
        # so that the artifact mask is applied to the right warp
        dataIds = [ref.dataId for ref in tempExpRefList]
        psfMatchedDataIds = [ref.dataId for ref in supplementaryData.warpRefList]

        if dataIds != psfMatchedDataIds:
            self.log.info("Reordering and or/padding PSF-matched visit input list")
            supplementaryData.warpRefList = reorderAndPadList(supplementaryData.warpRefList,
                                                              psfMatchedDataIds, dataIds)
            supplementaryData.imageScalerList = reorderAndPadList(supplementaryData.imageScalerList,
                                                                  psfMatchedDataIds, dataIds)

        # Use PSF-Matched Warps (and corresponding scalers) and coadd to find artifacts
        spanSetMaskList = self.findArtifacts(supplementaryData.templateCoadd,
                                             supplementaryData.warpRefList,
                                             supplementaryData.imageScalerList)

        badMaskPlanes = self.config.badMaskPlanes[:]
        badMaskPlanes.append("CLIPPED")
        badPixelMask = afwImage.Mask.getPlaneBitMask(badMaskPlanes)

        result = AssembleCoaddTask.run(self, skyInfo, tempExpRefList, imageScalerList, weightList,
                                       spanSetMaskList, mask=badPixelMask)

        # Propagate PSF-matched EDGE pixels to coadd SENSOR_EDGE and INEXACT_PSF
        # Psf-Matching moves the real edge inwards
        self.applyAltEdgeMask(result.coaddExposure.maskedImage.mask, spanSetMaskList)
        return result

    def applyAltEdgeMask(self, mask, altMaskList):
        """Propagate alt EDGE mask to SENSOR_EDGE AND INEXACT_PSF planes.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            Original mask.
        altMaskList : `list`
            List of Dicts containing ``spanSet`` lists.
            Each element contains the new mask plane name (e.g. "CLIPPED
            and/or "NO_DATA") as the key, and list of ``SpanSets`` to apply to
            the mask.
        """
        maskValue = mask.getPlaneBitMask(["SENSOR_EDGE", "INEXACT_PSF"])
        for visitMask in altMaskList:
            if "EDGE" in visitMask:
                for spanSet in visitMask['EDGE']:
                    spanSet.clippedTo(mask.getBBox()).setMask(mask, maskValue)

    def findArtifacts(self, templateCoadd, tempExpRefList, imageScalerList):
        """Find artifacts.

        Loop through warps twice. The first loop builds a map with the count
        of how many epochs each pixel deviates from the templateCoadd by more
        than ``config.chiThreshold`` sigma. The second loop takes each
        difference image and filters the artifacts detected in each using
        count map to filter out variable sources and sources that are
        difficult to subtract cleanly.

        Parameters
        ----------
        templateCoadd : `lsst.afw.image.Exposure`
            Exposure to serve as model of static sky.
        tempExpRefList : `list`
            List of data references to warps.
        imageScalerList : `list`
            List of image scalers.

        Returns
        -------
        altMasks : `list`
            List of dicts containing information about CLIPPED
            (i.e., artifacts), NO_DATA, and EDGE pixels.
        """

        self.log.debug("Generating Count Image, and mask lists.")
        coaddBBox = templateCoadd.getBBox()
        slateIm = afwImage.ImageU(coaddBBox)
        epochCountImage = afwImage.ImageU(coaddBBox)
        nImage = afwImage.ImageU(coaddBBox)
        spanSetArtifactList = []
        spanSetNoDataMaskList = []
        spanSetEdgeList = []
        spanSetBadMorphoList = []
        badPixelMask = self.getBadPixelMask()

        # mask of the warp diffs should = that of only the warp
        templateCoadd.mask.clearAllMaskPlanes()

        if self.config.doPreserveContainedBySource:
            templateFootprints = self.detectTemplate.detectFootprints(templateCoadd)
        else:
            templateFootprints = None

        for warpRef, imageScaler in zip(tempExpRefList, imageScalerList):
            warpDiffExp = self._readAndComputeWarpDiff(warpRef, imageScaler, templateCoadd)
            if warpDiffExp is not None:
                # This nImage only approximates the final nImage because it uses the PSF-matched mask
                nImage.array += (numpy.isfinite(warpDiffExp.image.array)
                                 * ((warpDiffExp.mask.array & badPixelMask) == 0)).astype(numpy.uint16)
                fpSet = self.detect.detectFootprints(warpDiffExp, doSmooth=False, clearMask=True)
                fpSet.positive.merge(fpSet.negative)
                footprints = fpSet.positive
                slateIm.set(0)
                spanSetList = [footprint.spans for footprint in footprints.getFootprints()]

                # Remove artifacts due to defects before they contribute to the epochCountImage
                if self.config.doPrefilterArtifacts:
                    spanSetList = self.prefilterArtifacts(spanSetList, warpDiffExp)

                # Clear mask before adding prefiltered spanSets
                self.detect.clearMask(warpDiffExp.mask)
                for spans in spanSetList:
                    spans.setImage(slateIm, 1, doClip=True)
                    spans.setMask(warpDiffExp.mask, warpDiffExp.mask.getPlaneBitMask("DETECTED"))
                epochCountImage += slateIm

                if self.config.doFilterMorphological:
                    maskName = self.config.streakMaskName
                    _ = self.maskStreaks.run(warpDiffExp)
                    streakMask = warpDiffExp.mask
                    spanSetStreak = afwGeom.SpanSet.fromMask(streakMask,
                                                             streakMask.getPlaneBitMask(maskName)).split()

                # PSF-Matched warps have less available area (~the matching kernel) because the calexps
                # undergo a second convolution. Pixels with data in the direct warp
                # but not in the PSF-matched warp will not have their artifacts detected.
                # NaNs from the PSF-matched warp therefore must be masked in the direct warp
                nans = numpy.where(numpy.isnan(warpDiffExp.maskedImage.image.array), 1, 0)
                nansMask = afwImage.makeMaskFromArray(nans.astype(afwImage.MaskPixel))
                nansMask.setXY0(warpDiffExp.getXY0())
                edgeMask = warpDiffExp.mask
                spanSetEdgeMask = afwGeom.SpanSet.fromMask(edgeMask,
                                                           edgeMask.getPlaneBitMask("EDGE")).split()
            else:
                # If the directWarp has <1% coverage, the psfMatchedWarp can have 0% and not exist
                # In this case, mask the whole epoch
                nansMask = afwImage.MaskX(coaddBBox, 1)
                spanSetList = []
                spanSetEdgeMask = []
                spanSetStreak = []

            spanSetNoDataMask = afwGeom.SpanSet.fromMask(nansMask).split()

            spanSetNoDataMaskList.append(spanSetNoDataMask)
            spanSetArtifactList.append(spanSetList)
            spanSetEdgeList.append(spanSetEdgeMask)
            if self.config.doFilterMorphological:
                spanSetBadMorphoList.append(spanSetStreak)

        if lsstDebug.Info(__name__).saveCountIm:
            path = self._dataRef2DebugPath("epochCountIm", tempExpRefList[0], coaddLevel=True)
            epochCountImage.writeFits(path)

        for i, spanSetList in enumerate(spanSetArtifactList):
            if spanSetList:
                filteredSpanSetList = self.filterArtifacts(spanSetList, epochCountImage, nImage,
                                                           templateFootprints)
                spanSetArtifactList[i] = filteredSpanSetList
            if self.config.doFilterMorphological:
                spanSetArtifactList[i] += spanSetBadMorphoList[i]

        altMasks = []
        for artifacts, noData, edge in zip(spanSetArtifactList, spanSetNoDataMaskList, spanSetEdgeList):
            altMasks.append({'CLIPPED': artifacts,
                             'NO_DATA': noData,
                             'EDGE': edge})
        return altMasks

    def prefilterArtifacts(self, spanSetList, exp):
        """Remove artifact candidates covered by bad mask plane.

        Any future editing of the candidate list that does not depend on
        temporal information should go in this method.

        Parameters
        ----------
        spanSetList : `list`
            List of SpanSets representing artifact candidates.
        exp : `lsst.afw.image.Exposure`
            Exposure containing mask planes used to prefilter.

        Returns
        -------
        returnSpanSetList : `list`
            List of SpanSets with artifacts.
        """
        badPixelMask = exp.mask.getPlaneBitMask(self.config.prefilterArtifactsMaskPlanes)
        goodArr = (exp.mask.array & badPixelMask) == 0
        returnSpanSetList = []
        bbox = exp.getBBox()
        x0, y0 = exp.getXY0()
        for i, span in enumerate(spanSetList):
            y, x = span.clippedTo(bbox).indices()
            yIndexLocal = numpy.array(y) - y0
            xIndexLocal = numpy.array(x) - x0
            goodRatio = numpy.count_nonzero(goodArr[yIndexLocal, xIndexLocal])/span.getArea()
            if goodRatio > self.config.prefilterArtifactsRatio:
                returnSpanSetList.append(span)
        return returnSpanSetList

    def filterArtifacts(self, spanSetList, epochCountImage, nImage, footprintsToExclude=None):
        """Filter artifact candidates.

        Parameters
        ----------
        spanSetList : `list`
            List of SpanSets representing artifact candidates.
        epochCountImage : `lsst.afw.image.Image`
            Image of accumulated number of warpDiff detections.
        nImage : `lsst.afw.image.Image`
            Image of the accumulated number of total epochs contributing.

        Returns
        -------
        maskSpanSetList : `list`
            List of SpanSets with artifacts.
        """

        maskSpanSetList = []
        x0, y0 = epochCountImage.getXY0()
        for i, span in enumerate(spanSetList):
            y, x = span.indices()
            yIdxLocal = [y1 - y0 for y1 in y]
            xIdxLocal = [x1 - x0 for x1 in x]
            outlierN = epochCountImage.array[yIdxLocal, xIdxLocal]
            totalN = nImage.array[yIdxLocal, xIdxLocal]

            # effectiveMaxNumEpochs is broken line (fraction of N) with characteristic config.maxNumEpochs
            effMaxNumEpochsHighN = (self.config.maxNumEpochs
                                    + self.config.maxFractionEpochsHigh*numpy.mean(totalN))
            effMaxNumEpochsLowN = self.config.maxFractionEpochsLow * numpy.mean(totalN)
            effectiveMaxNumEpochs = int(min(effMaxNumEpochsLowN, effMaxNumEpochsHighN))
            nPixelsBelowThreshold = numpy.count_nonzero((outlierN > 0)
                                                        & (outlierN <= effectiveMaxNumEpochs))
            percentBelowThreshold = nPixelsBelowThreshold / len(outlierN)
            if percentBelowThreshold > self.config.spatialThreshold:
                maskSpanSetList.append(span)

        if self.config.doPreserveContainedBySource and footprintsToExclude is not None:
            # If a candidate is contained by a footprint on the template coadd, do not clip
            filteredMaskSpanSetList = []
            for span in maskSpanSetList:
                doKeep = True
                for footprint in footprintsToExclude.positive.getFootprints():
                    if footprint.spans.contains(span):
                        doKeep = False
                        break
                if doKeep:
                    filteredMaskSpanSetList.append(span)
            maskSpanSetList = filteredMaskSpanSetList

        return maskSpanSetList

    def _readAndComputeWarpDiff(self, warpRef, imageScaler, templateCoadd):
        """Fetch a warp from the butler and return a warpDiff.

        Parameters
        ----------
        warpRef : `lsst.daf.persistence.butlerSubset.ButlerDataRef`
            Butler dataRef for the warp.
        imageScaler : `lsst.pipe.tasks.scaleZeroPoint.ImageScaler`
            An image scaler object.
        templateCoadd : `lsst.afw.image.Exposure`
            Exposure to be substracted from the scaled warp.

        Returns
        -------
        warp : `lsst.afw.image.Exposure`
            Exposure of the image difference between the warp and template.
        """

        # If the PSF-Matched warp did not exist for this direct warp
        # None is holding its place to maintain order in Gen 3
        if warpRef is None:
            return None
        # Warp comparison must use PSF-Matched Warps regardless of requested coadd warp type
        warpName = self.getTempExpDatasetName('psfMatched')
        if not isinstance(warpRef, DeferredDatasetHandle):
            if not warpRef.datasetExists(warpName):
                self.log.warn("Could not find %s %s; skipping it", warpName, warpRef.dataId)
                return None
        warp = warpRef.get(datasetType=warpName, immediate=True)
        # direct image scaler OK for PSF-matched Warp
        imageScaler.scaleMaskedImage(warp.getMaskedImage())
        mi = warp.getMaskedImage()
        if self.config.doScaleWarpVariance:
            try:
                self.scaleWarpVariance.run(mi)
            except Exception as exc:
                self.log.warn("Unable to rescale variance of warp (%s); leaving it as-is" % (exc,))
        mi -= templateCoadd.getMaskedImage()
        return warp

    def _dataRef2DebugPath(self, prefix, warpRef, coaddLevel=False):
        """Return a path to which to write debugging output.

        Creates a hyphen-delimited string of dataId values for simple filenames.

        Parameters
        ----------
        prefix : `str`
            Prefix for filename.
        warpRef : `lsst.daf.persistence.butlerSubset.ButlerDataRef`
            Butler dataRef to make the path from.
        coaddLevel : `bool`, optional.
            If True, include only coadd-level keys (e.g., 'tract', 'patch',
            'filter', but no 'visit').

        Returns
        -------
        result : `str`
            Path for debugging output.
        """
        if coaddLevel:
            keys = warpRef.getButler().getKeys(self.getCoaddDatasetName(self.warpType))
        else:
            keys = warpRef.dataId.keys()
        keyList = sorted(keys, reverse=True)
        directory = lsstDebug.Info(__name__).figPath if lsstDebug.Info(__name__).figPath else "."
        filename = "%s-%s.fits" % (prefix, '-'.join([str(warpRef.dataId[k]) for k in keyList]))
        return os.path.join(directory, filename)
