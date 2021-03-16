#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import numpy as np
import lsst.sphgeom
import lsst.pex.config as pexConfig
import lsst.pex.exceptions as pexExceptions
import lsst.geom as geom
import lsst.pipe.base as pipeBase

__all__ = ["BaseSelectImagesTask", "BaseExposureInfo", "WcsSelectImagesTask", "PsfWcsSelectImagesTask",
           "DatabaseSelectImagesConfig", "BestSeeingWcsSelectImagesTask", "BestSeeingSelectVisitsTask",
           "BestSeeingPercentileSelectVisitsTask"]


class DatabaseSelectImagesConfig(pexConfig.Config):
    """Base configuration for subclasses of BaseSelectImagesTask that use a database"""
    host = pexConfig.Field(
        doc="Database server host name",
        dtype=str,
    )
    port = pexConfig.Field(
        doc="Database server port",
        dtype=int,
    )
    database = pexConfig.Field(
        doc="Name of database",
        dtype=str,
    )
    maxExposures = pexConfig.Field(
        doc="maximum exposures to select; intended for debugging; ignored if None",
        dtype=int,
        optional=True,
    )


class BaseExposureInfo(pipeBase.Struct):
    """Data about a selected exposure
    """

    def __init__(self, dataId, coordList):
        """Create exposure information that can be used to generate data references

        The object has the following fields:
        - dataId: data ID of exposure (a dict)
        - coordList: ICRS coordinates of the corners of the exposure (list of lsst.geom.SpherePoint)
        plus any others items that are desired
        """
        super(BaseExposureInfo, self).__init__(dataId=dataId, coordList=coordList)


class BaseSelectImagesTask(pipeBase.Task):
    """Base task for selecting images suitable for coaddition
    """
    ConfigClass = pexConfig.Config
    _DefaultName = "selectImages"

    @pipeBase.timeMethod
    def run(self, coordList):
        """Select images suitable for coaddition in a particular region

        @param[in] coordList: list of coordinates defining region of interest; if None then select all images
        subclasses may add additional keyword arguments, as required

        @return a pipeBase Struct containing:
        - exposureInfoList: a list of exposure information objects (subclasses of BaseExposureInfo),
            which have at least the following fields:
            - dataId: data ID dictionary
            - coordList: ICRS coordinates of the corners of the exposure (list of lsst.geom.SpherePoint)
        """
        raise NotImplementedError()

    def _runArgDictFromDataId(self, dataId):
        """Extract keyword arguments for run (other than coordList) from a data ID

        @return keyword arguments for run (other than coordList), as a dict
        """
        raise NotImplementedError()

    def runDataRef(self, dataRef, coordList, makeDataRefList=True, selectDataList=[]):
        """Run based on a data reference

        This delegates to run() and _runArgDictFromDataId() to do the actual
        selection. In the event that the selectDataList is non-empty, this will
        be used to further restrict the selection, providing the user with
        additional control over the selection.

        @param[in] dataRef: data reference; must contain any extra keys needed by the subclass
        @param[in] coordList: list of coordinates defining region of interest; if None, search the whole sky
        @param[in] makeDataRefList: if True, return dataRefList
        @param[in] selectDataList: List of SelectStruct with dataRefs to consider for selection
        @return a pipeBase Struct containing:
        - exposureInfoList: a list of objects derived from ExposureInfo
        - dataRefList: a list of data references (None if makeDataRefList False)
        """
        runArgDict = self._runArgDictFromDataId(dataRef.dataId)
        exposureInfoList = self.run(coordList, **runArgDict).exposureInfoList

        if len(selectDataList) > 0 and len(exposureInfoList) > 0:
            # Restrict the exposure selection further
            ccdKeys, ccdValues = _extractKeyValue(exposureInfoList)
            inKeys, inValues = _extractKeyValue([s.dataRef for s in selectDataList], keys=ccdKeys)
            inValues = set(inValues)
            newExposureInfoList = []
            for info, ccdVal in zip(exposureInfoList, ccdValues):
                if ccdVal in inValues:
                    newExposureInfoList.append(info)
                else:
                    self.log.info("De-selecting exposure %s: not in selectDataList" % info.dataId)
            exposureInfoList = newExposureInfoList

        if makeDataRefList:
            butler = dataRef.butlerSubset.butler
            dataRefList = [butler.dataRef(datasetType="calexp",
                                          dataId=expInfo.dataId,
                                          ) for expInfo in exposureInfoList]
        else:
            dataRefList = None

        return pipeBase.Struct(
            dataRefList=dataRefList,
            exposureInfoList=exposureInfoList,
        )


def _extractKeyValue(dataList, keys=None):
    """Extract the keys and values from a list of dataIds

    The input dataList is a list of objects that have 'dataId' members.
    This allows it to be used for both a list of data references and a
    list of ExposureInfo
    """
    assert len(dataList) > 0
    if keys is None:
        keys = sorted(dataList[0].dataId.keys())
    keySet = set(keys)
    values = list()
    for data in dataList:
        thisKeys = set(data.dataId.keys())
        if thisKeys != keySet:
            raise RuntimeError("DataId keys inconsistent: %s vs %s" % (keySet, thisKeys))
        values.append(tuple(data.dataId[k] for k in keys))
    return keys, values


class SelectStruct(pipeBase.Struct):
    """A container for data to be passed to the WcsSelectImagesTask"""

    def __init__(self, dataRef, wcs, bbox):
        super(SelectStruct, self).__init__(dataRef=dataRef, wcs=wcs, bbox=bbox)


class WcsSelectImagesTask(BaseSelectImagesTask):
    """Select images using their Wcs

        We use the "convexHull" method of lsst.sphgeom.ConvexPolygon to define
        polygons on the celestial sphere, and test the polygon of the
        patch for overlap with the polygon of the image.

        We use "convexHull" instead of generating a ConvexPolygon
        directly because the standard for the inputs to ConvexPolygon
        are pretty high and we don't want to be responsible for reaching them.
        """

    def runDataRef(self, dataRef, coordList, makeDataRefList=True, selectDataList=[]):
        """Select images in the selectDataList that overlap the patch

        This method is the old entry point for the Gen2 commandline tasks and drivers
        Will be deprecated in v22.

        @param dataRef: Data reference for coadd/tempExp (with tract, patch)
        @param coordList: List of ICRS coordinates (lsst.geom.SpherePoint) specifying boundary of patch
        @param makeDataRefList: Construct a list of data references?
        @param selectDataList: List of SelectStruct, to consider for selection
        """
        dataRefList = []
        exposureInfoList = []

        patchVertices = [coord.getVector() for coord in coordList]
        patchPoly = lsst.sphgeom.ConvexPolygon.convexHull(patchVertices)

        for data in selectDataList:
            dataRef = data.dataRef
            imageWcs = data.wcs
            imageBox = data.bbox

            imageCorners = self.getValidImageCorners(imageWcs, imageBox, patchPoly, dataId=None)
            if imageCorners:
                dataRefList.append(dataRef)
                exposureInfoList.append(BaseExposureInfo(dataRef.dataId, imageCorners))

        return pipeBase.Struct(
            dataRefList=dataRefList if makeDataRefList else None,
            exposureInfoList=exposureInfoList,
        )

    def run(self, wcsList, bboxList, coordList, dataIds=None, **kwargs):
        """Return indices of provided lists that meet the selection criteria

        Parameters:
        -----------
        wcsList : `list` of `lsst.afw.geom.SkyWcs`
            specifying the WCS's of the input ccds to be selected
        bboxList : `list` of `lsst.geom.Box2I`
            specifying the bounding boxes of the input ccds to be selected
        coordList : `list` of `lsst.geom.SpherePoint`
            ICRS coordinates specifying boundary of the patch.

        Returns:
        --------
        result: `list` of `int`
            of indices of selected ccds
        """
        if dataIds is None:
            dataIds = [None] * len(wcsList)
        patchVertices = [coord.getVector() for coord in coordList]
        patchPoly = lsst.sphgeom.ConvexPolygon.convexHull(patchVertices)
        result = []
        for i, (imageWcs, imageBox, dataId) in enumerate(zip(wcsList, bboxList, dataIds)):
            imageCorners = self.getValidImageCorners(imageWcs, imageBox, patchPoly, dataId)
            if imageCorners:
                result.append(i)
        return result

    def getValidImageCorners(self, imageWcs, imageBox, patchPoly, dataId=None):
        "Return corners or None if bad"
        try:
            imageCorners = [imageWcs.pixelToSky(pix) for pix in geom.Box2D(imageBox).getCorners()]
        except (pexExceptions.DomainError, pexExceptions.RuntimeError) as e:
            # Protecting ourselves from awful Wcs solutions in input images
            self.log.debug("WCS error in testing calexp %s (%s): deselecting", dataId, e)
            return

        imagePoly = lsst.sphgeom.ConvexPolygon.convexHull([coord.getVector() for coord in imageCorners])
        if imagePoly is None:
            self.log.debug("Unable to create polygon from image %s: deselecting", dataId)
            return

        if patchPoly.intersects(imagePoly):
            # "intersects" also covers "contains" or "is contained by"
            self.log.info("Selecting calexp %s" % dataId)
            return imageCorners


def sigmaMad(array):
    "Return median absolute deviation scaled to normally distributed data"
    return 1.4826*np.median(np.abs(array - np.median(array)))


class PsfWcsSelectImagesConnections(pipeBase.PipelineTaskConnections,
                                    dimensions=("tract", "patch", "skymap", "instrument", "visit"),
                                    defaultTemplates={"coaddName": "deep"}):
    pass


class PsfWcsSelectImagesConfig(pipeBase.PipelineTaskConfig,
                               pipelineConnections=PsfWcsSelectImagesConnections):
    maxEllipResidual = pexConfig.Field(
        doc="Maximum median ellipticity residual",
        dtype=float,
        default=0.007,
        optional=True,
    )
    maxSizeScatter = pexConfig.Field(
        doc="Maximum scatter in the size residuals",
        dtype=float,
        optional=True,
    )
    maxScaledSizeScatter = pexConfig.Field(
        doc="Maximum scatter in the size residuals, scaled by the median size",
        dtype=float,
        default=0.009,
        optional=True,
    )
    starSelection = pexConfig.Field(
        doc="select star with this field",
        dtype=str,
        default='calib_psf_used'
    )
    starShape = pexConfig.Field(
        doc="name of star shape",
        dtype=str,
        default='base_SdssShape'
    )
    psfShape = pexConfig.Field(
        doc="name of psf shape",
        dtype=str,
        default='base_SdssShape_psf'
    )


class PsfWcsSelectImagesTask(WcsSelectImagesTask):
    """Select images using their Wcs and cuts on the PSF properties

        The PSF quality criteria are based on the size and ellipticity residuals from the
        adaptive second moments of the star and the PSF.

        The criteria are:
          - the median of the ellipticty residuals
          - the robust scatter of the size residuals (using the median absolute deviation)
          - the robust scatter of the size residuals scaled by the square of
            the median size
    """

    ConfigClass = PsfWcsSelectImagesConfig
    _DefaultName = "PsfWcsSelectImages"

    def runDataRef(self, dataRef, coordList, makeDataRefList=True, selectDataList=[]):
        """Select images in the selectDataList that overlap the patch and satisfy PSF quality critera.

        This method is the old entry point for the Gen2 commandline tasks and drivers
        Will be deprecated in v22.

        @param dataRef: Data reference for coadd/tempExp (with tract, patch)
        @param coordList: List of ICRS coordinates (lsst.geom.SpherePoint) specifying boundary of patch
        @param makeDataRefList: Construct a list of data references?
        @param selectDataList: List of SelectStruct, to consider for selection
        """
        result = super(PsfWcsSelectImagesTask, self).runDataRef(dataRef, coordList, makeDataRefList,
                                                                selectDataList)

        dataRefList = []
        exposureInfoList = []
        for dataRef, exposureInfo in zip(result.dataRefList, result.exposureInfoList):
            butler = dataRef.butlerSubset.butler
            srcCatalog = butler.get('src', dataRef.dataId)
            valid = self.isValid(srcCatalog, dataRef.dataId)
            if valid is False:
                continue

            dataRefList.append(dataRef)
            exposureInfoList.append(exposureInfo)

        return pipeBase.Struct(
            dataRefList=dataRefList,
            exposureInfoList=exposureInfoList,
        )

    def run(self, wcsList, bboxList, coordList, srcList, dataIds=None, **kwargs):
        """Return indices of provided lists that meet the selection criteria

        Parameters:
        -----------
        wcsList : `list` of `lsst.afw.geom.SkyWcs`
            specifying the WCS's of the input ccds to be selected
        bboxList : `list` of `lsst.geom.Box2I`
            specifying the bounding boxes of the input ccds to be selected
        coordList : `list` of `lsst.geom.SpherePoint`
            ICRS coordinates specifying boundary of the patch.
        srcList : `list` of `lsst.afw.table.SourceCatalog`
            containing the PSF shape information for the input ccds to be selected

        Returns:
        --------
        goodPsf: `list` of `int`
            of indices of selected ccds
        """
        goodWcs = super(PsfWcsSelectImagesTask, self).run(wcsList=wcsList, bboxList=bboxList,
                                                          coordList=coordList, dataIds=dataIds)

        goodPsf = []
        if dataIds is None:
            dataIds = [None] * len(srcList)
        for i, (srcCatalog, dataId) in enumerate(zip(srcList, dataIds)):
            if i not in goodWcs:
                continue
            if self.isValid(srcCatalog, dataId):
                goodPsf.append(i)

        return goodPsf

    def isValid(self, srcCatalog, dataId=None):
        """Should this ccd be selected based on its PSF shape information

        Parameters
        ----------
        srcCatalog : `lsst.afw.table.SourceCatalog`
        dataId : `dict` of dataId keys, optional.
            Used only for logging. Defaults to None.

        Returns
        -------
        valid : `bool`
            True if selected.
        """
        mask = srcCatalog[self.config.starSelection]

        starXX = srcCatalog[self.config.starShape+'_xx'][mask]
        starYY = srcCatalog[self.config.starShape+'_yy'][mask]
        starXY = srcCatalog[self.config.starShape+'_xy'][mask]
        psfXX = srcCatalog[self.config.psfShape+'_xx'][mask]
        psfYY = srcCatalog[self.config.psfShape+'_yy'][mask]
        psfXY = srcCatalog[self.config.psfShape+'_xy'][mask]

        starSize = np.power(starXX*starYY - starXY**2, 0.25)
        starE1 = (starXX - starYY)/(starXX + starYY)
        starE2 = 2*starXY/(starXX + starYY)
        medianSize = np.median(starSize)

        psfSize = np.power(psfXX*psfYY - psfXY**2, 0.25)
        psfE1 = (psfXX - psfYY)/(psfXX + psfYY)
        psfE2 = 2*psfXY/(psfXX + psfYY)

        medianE1 = np.abs(np.median(starE1 - psfE1))
        medianE2 = np.abs(np.median(starE2 - psfE2))
        medianE = np.sqrt(medianE1**2 + medianE2**2)

        scatterSize = sigmaMad(starSize - psfSize)
        scaledScatterSize = scatterSize/medianSize**2

        valid = True
        if self.config.maxEllipResidual and medianE > self.config.maxEllipResidual:
            self.log.info("Removing visit %s because median e residual too large: %f vs %f" %
                          (dataId, medianE, self.config.maxEllipResidual))
            valid = False
        elif self.config.maxSizeScatter and scatterSize > self.config.maxSizeScatter:
            self.log.info("Removing visit %s because size scatter is too large: %f vs %f" %
                          (dataId, scatterSize, self.config.maxSizeScatter))
            valid = False
        elif self.config.maxScaledSizeScatter and scaledScatterSize > self.config.maxScaledSizeScatter:
            self.log.info("Removing visit %s because scaled size scatter is too large: %f vs %f" %
                          (dataId, scaledScatterSize, self.config.maxScaledSizeScatter))
            valid = False

        return valid


class BestSeeingWcsSelectImageConfig(WcsSelectImagesTask.ConfigClass):
    """Base configuration for BestSeeingSelectImagesTask.
    """
    nImagesMax = pexConfig.RangeField(
        dtype=int,
        doc="Maximum number of images to select",
        default=5,
        min=0)
    maxPsfFwhm = pexConfig.Field(
        dtype=float,
        doc="Maximum PSF FWHM (in arcseconds) to select",
        default=1.5,
        optional=True)
    minPsfFwhm = pexConfig.Field(
        dtype=float,
        doc="Minimum PSF FWHM (in arcseconds) to select",
        default=0.,
        optional=True)


class BestSeeingWcsSelectImagesTask(WcsSelectImagesTask):
    """Select up to a maximum number of the best-seeing images using their Wcs.
    """
    ConfigClass = BestSeeingWcsSelectImageConfig

    def runDataRef(self, dataRef, coordList, makeDataRefList=True,
                   selectDataList=None):
        """Select the best-seeing images in the selectDataList that overlap the patch.

        This method is the old entry point for the Gen2 commandline tasks and drivers
        Will be deprecated in v22.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Data reference for coadd/tempExp (with tract, patch)
        coordList : `list` of `lsst.geom.SpherePoint`
            List of ICRS sky coordinates specifying boundary of patch
        makeDataRefList : `boolean`, optional
            Construct a list of data references?
        selectDataList : `list` of `SelectStruct`
            List of SelectStruct, to consider for selection

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:
            - ``exposureList``: the selected exposures
                (`list` of `lsst.pipe.tasks.selectImages.BaseExposureInfo`).
            - ``dataRefList``: the optional data references corresponding to
                each element of ``exposureList``
                (`list` of `lsst.daf.persistence.ButlerDataRef`, or `None`).
        """
        psfSizes = []
        dataRefList = []
        exposureInfoList = []

        if selectDataList is None:
            selectDataList = []

        result = super().runDataRef(dataRef, coordList, makeDataRefList=True, selectDataList=selectDataList)

        for dataRef, exposureInfo in zip(result.dataRefList, result.exposureInfoList):
            cal = dataRef.get("calexp", immediate=True)

            # if min/max PSF values are defined, remove images out of bounds
            pixToArcseconds = cal.getWcs().getPixelScale().asArcseconds()
            psfSize = cal.getPsf().computeShape().getDeterminantRadius()*pixToArcseconds
            sizeFwhm = psfSize * np.sqrt(8.*np.log(2.))
            if self.config.maxPsfFwhm and sizeFwhm > self.config.maxPsfFwhm:
                continue
            if self.config.minPsfFwhm and sizeFwhm < self.config.minPsfFwhm:
                continue
            psfSizes.append(sizeFwhm)
            dataRefList.append(dataRef)
            exposureInfoList.append(exposureInfo)

        if len(psfSizes) > self.config.nImagesMax:
            sortedIndices = np.argsort(psfSizes)[:self.config.nImagesMax]
            filteredDataRefList = [dataRefList[i] for i in sortedIndices]
            filteredExposureInfoList = [exposureInfoList[i] for i in sortedIndices]
            self.log.info(f"{len(sortedIndices)} images selected with FWHM "
                          f"range of {psfSizes[sortedIndices[0]]}--{psfSizes[sortedIndices[-1]]} arcseconds")

        else:
            if len(psfSizes) == 0:
                self.log.warn("0 images selected.")
            else:
                self.log.debug(f"{len(psfSizes)} images selected with FWHM range "
                               f"of {psfSizes[0]}--{psfSizes[-1]} arcseconds")
            filteredDataRefList = dataRefList
            filteredExposureInfoList = exposureInfoList

        return pipeBase.Struct(
            dataRefList=filteredDataRefList if makeDataRefList else None,
            exposureInfoList=filteredExposureInfoList,
        )

    def run(self, wcsList, bboxList, coordList, psfList, dataIds, **kwargs):
        """Return indices of good calexps from a list.

        This task does not make sense for use with makeWarp where there quanta
        are per-visit rather than per-patch. This task selectes the best ccds
        of ONE VISIT that overlap the patch.

        This includes some code duplication with runDataRef,
        but runDataRef will be deprecated as of v22.

        Parameters:
        -----------
        wcsList : `list` of `lsst.afw.geom.SkyWcs`
            specifying the WCS's of the input ccds to be selected
        bboxList : `list` of `lsst.geom.Box2I`
            specifying the bounding boxes of the input ccds to be selected
        coordList : `list` of `lsst.geom.SpherePoint`
            ICRS coordinates specifying boundary of the patch.
        psfList : `list` of `lsst.afw.detection.Psf`
            specifying the PSF model of the input ccds to be selected

        Returns:
        --------
        output: `list` of `int`
            of indices of selected ccds sorted by seeing
        """
        goodWcs = super().run(wcsList=wcsList, bboxList=bboxList, coordList=coordList, dataIds=dataIds)

        psfSizes = []
        indices = []
        for i, (wcs, psf) in enumerate(wcsList, psfList):
            if i not in goodWcs:
                continue
            # if min/max PSF values are defined, remove images out of bounds
            pixToArcseconds = wcs.getPixelScale().asArcseconds()
            psfSize = psf.computeShape().getDeterminantRadius()*pixToArcseconds
            sizeFwhm = psfSize * np.sqrt(8.*np.log(2.))
            if self.config.maxPsfFwhm and sizeFwhm > self.config.maxPsfFwhm:
                continue
            if self.config.minPsfFwhm and sizeFwhm < self.config.minPsfFwhm:
                continue
            psfSizes.append(sizeFwhm)
            indices.append(i)

        sortedIndices = [ind for (_, ind) in sorted(zip(psfSizes, indices))]
        output = sortedIndices[:self.config.nImagesMax]
        self.log.info(f"{len(output)} images selected with FWHM "
                      f"range of {psfSizes[indices.index(output[0])]}"
                      f"--{psfSizes[indices.index(output[-1])]} arcseconds")
        return output


class BestSeeingSelectVisitsConnections(pipeBase.PipelineTaskConnections,
                                        dimensions=("tract", "patch", "skymap", "band", "instrument"),
                                        defaultTemplates={"coaddName": "bestSeeing"}):
    sourceTables = pipeBase.connectionTypes.Input(
        doc="Source table with PSF size information",
        name="sourceTable_visit",
        storageClass="DataFrame",
        dimensions=("instrument", "visit"),
        multiple=True,
        deferLoad=True
    )
    goodVisits = pipeBase.connectionTypes.Output(
        doc="Selected visits to be coadded.",
        name="{coaddName}VisitsDict",
        storageClass="StructuredDataDict",
        dimensions=("instrument", "tract", "patch", "skymap", "band"),
    )


class BestSeeingSelectVisitsConfig(pipeBase.PipelineTaskConfig,
                                   pipelineConnections=BestSeeingSelectVisitsConnections):
    nImagesMax = pexConfig.RangeField(
        dtype=int,
        doc="Maximum number of images to select",
        default=5,
        min=0)
    maxPsfFwhm = pexConfig.Field(
        dtype=float,
        doc="Maximum PSF FWHM (in arcseconds) to select",
        default=1.5,
        optional=True)
    minPsfFwhm = pexConfig.Field(
        dtype=float,
        doc="Minimum PSF FWHM (in arcseconds) to select",
        default=0.,
        optional=True)


class BestSeeingSelectVisitsTask(pipeBase.PipelineTask):
    ConfigClass = BestSeeingSelectVisitsConfig
    _DefaultName = 'bestSeeingSelectVisits'

    def run(self, sourceTables):
        inputVisits = [sourceTable.ref.dataId['visit'] for sourceTable in sourceTables]
        psfSizes = []
        visits = []
        for visit, sourceTable in zip(inputVisits, sourceTables):
            df = sourceTable.get(parameters={"columns": ['IxxPsf', 'IyyPsf', 'IxyPsf',
                                                         'LocalWcs_CDMatrix_2_1', 'LocalWcs_CDMatrix_1_1',
                                                         'LocalWcs_CDMatrix_1_2', 'LocalWcs_CDMatrix_2_2']})
            wcsDet = (df.LocalWcs_CDMatrix_1_1 * df.LocalWcs_CDMatrix_2_2
                      - df.LocalWcs_CDMatrix_1_2 * df.LocalWcs_CDMatrix_1_2)
            # if min/max PSF values are defined, remove images out of bounds
            pixToArcseconds = np.nanmedian(3600*np.degrees(np.sqrt(np.fabs(wcsDet))))
            determinantRadius = np.power(np.nanmedian(df.IxxPsf * df.IxxPsf - df.IxyPsf**2), 0.25)
            psfSize = determinantRadius*pixToArcseconds
            sizeFwhm = psfSize * np.sqrt(8.*np.log(2.))
            if self.config.maxPsfFwhm and sizeFwhm > self.config.maxPsfFwhm:
                continue
            if self.config.minPsfFwhm and sizeFwhm < self.config.minPsfFwhm:
                continue
            psfSizes.append(sizeFwhm)
            visits.append(visit)

        sortedVisits = [ind for (_, ind) in sorted(zip(psfSizes, visits))]
        output = sortedVisits[:self.config.nImagesMax]
        self.log.info(f"{len(output)} images selected with FWHM "
                      f"range of {psfSizes[visits.index(output[0])]}"
                      f"--{psfSizes[visits.index(output[-1])]} arcseconds")
        # In order to store as a StructuredDataDict, convert list to dict
        goodVisits = {key: True for key in output}
        print(goodVisits)
        return pipeBase.Struct(goodVisits=goodVisits)


class BestSeeingPercentileSelectVisitsConfig(pipeBase.PipelineTaskConfig,
                                             pipelineConnections=BestSeeingSelectVisitsConnections):
    percentile = pexConfig.RangeField(
        doc="Select top percentile of BestPercentile seeing visits. 1 selects all visits. 0 selects None",
        dtype=float,
        default=0.33,
        min=0,
        max=1,
    )
    nMin = pexConfig.Field(
        doc="At least this number of visits selected and supercedes percentile. For example, if 10 visits "
            "cover this patch, percentile=0.33, and nMin=5, the best 5 visits will be selected.",
        dtype=int,
        default=3,
    )


class BestSeeingPercentileSelectVisitsTask(pipeBase.PipelineTask):
    ConfigClass = BestSeeingPercentileSelectVisitsConfig
    _DefaultName = 'bestSeeingPercentileSelectVisits'

    def run(self, sourceTables):

        visits = [sourceTable.ref.dataId['visit'] for sourceTable in sourceTables]

        radius = []
        for sourceTable in sourceTables:
            df = sourceTable.get(parameters={"columns": ['IxxPsf', 'IyyPsf']})
            traceRadius = np.sqrt(0.5) * np.sqrt(df.IxxPsf + df.IyyPsf)
            radius.append(traceRadius.median())

        sortedVisits = [v for rad, v in sorted(zip(radius, visits))]
        cutoff = max(int(np.round(self.config.percentile*len(visits))), self.config.nMin)
        goodVisits = sortedVisits[:cutoff]

        # In order to store as a StructuredDataDict, convert list to dict
        goodVisits = {visit: True for visit in goodVisits}
        print(goodVisits)
        return pipeBase.Struct(goodVisits=goodVisits)
