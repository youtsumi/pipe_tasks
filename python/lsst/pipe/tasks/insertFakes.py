# This file is part of pipe tasks
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Insert fakes into deepCoadds
"""
import galsim
from astropy.table import Table
import numpy as np

import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.pipe.base import CmdLineTask, PipelineTask, PipelineTaskConfig, PipelineTaskConnections
import lsst.pipe.base.connectionTypes as cT
from lsst.pex.exceptions import LogicError
from lsst.coadd.utils.coaddDataIdContainer import ExistingCoaddDataIdContainer
from lsst.geom import SpherePoint, radians, Box2D

__all__ = ["InsertFakesConfig", "InsertFakesTask"]


def addFakeSources(image, objects, calibFluxRadius=12.0, logger=None):
    """Add fake sources to the given image

    Parameters
    ----------
    image : `lsst.afw.image.exposure.exposure.ExposureF`
        The image into which the fake sources should be added
    objects : `typing.Iterator` [`tuple` ['lsst.geom.SpherePoint`, `galsim.GSObject`]]
        An iterator of tuples that contains (or generates) locations and object
        surface brightness profiles to inject.
    calibFluxRadius : `float`, optional
        Aperture radius (in pixels) used to define the calibration for this
        image+catalog.  This is used to produce the correct instrumental fluxes
        within the radius.  The value should match that of the field defined in
        slot_CalibFlux_instFlux.
    logger : `lsst.log.log.log.Log` or `logging.Logger`, optional
        Logger.
    """
    image.mask.addMaskPlane("FAKE")
    bitmask = image.mask.getPlaneBitMask("FAKE")
    if logger:
        logger.info(f"Adding mask plane with bitmask {bitmask}")

    wcs = image.getWcs()
    psf = image.getPsf()

    bbox = image.getBBox()
    fullBounds = galsim.BoundsI(bbox.minX, bbox.maxX, bbox.minY, bbox.maxY)
    gsImg = galsim.Image(image.image.array, bounds=fullBounds)

    for spt, gsObj in objects:
        pt = wcs.skyToPixel(spt)
        posd = galsim.PositionD(pt.x, pt.y)
        posi = galsim.PositionI(pt.x//1, pt.y//1)
        if logger:
            logger.debug(f"Adding fake source at {pt}")

        mat = wcs.linearizePixelToSky(spt, geom.arcseconds).getMatrix()
        gsWCS = galsim.JacobianWCS(mat[0, 0], mat[0, 1], mat[1, 0], mat[1, 1])

        psfArr = psf.computeKernelImage(pt).array
        apCorr = psf.computeApertureFlux(calibFluxRadius)
        psfArr /= apCorr
        gsPSF = galsim.InterpolatedImage(galsim.Image(psfArr), wcs=gsWCS)

        conv = galsim.Convolve(gsObj, gsPSF)
        stampSize = conv.getGoodImageSize(gsWCS.minLinearScale())
        subBounds = galsim.BoundsI(posi).withBorder(stampSize//2)
        subBounds &= fullBounds

        if subBounds.area() > 0:
            subImg = gsImg[subBounds]
            offset = posd - subBounds.true_center
            # Note, for calexp injection, pixel is already part of the PSF and
            # for coadd injection, it's incorrect to include the output pixel.
            # So for both cases, we draw using method='no_pixel'.
            conv.drawImage(
                subImg,
                add_to_image=True,
                offset=offset,
                wcs=gsWCS,
                method='no_pixel'
            )

            subBox = geom.Box2I(
                geom.Point2I(subBounds.xmin, subBounds.ymin),
                geom.Point2I(subBounds.xmax, subBounds.ymax)
            )
            image[subBox].mask.array |= bitmask


class InsertFakesConnections(PipelineTaskConnections,
                             defaultTemplates={"coaddName": "deep",
                                               "fakesType": "fakes_"},
                             dimensions=("tract", "patch", "band", "skymap")):

    image = cT.Input(
        doc="Image into which fakes are to be added.",
        name="{coaddName}Coadd",
        storageClass="ExposureF",
        dimensions=("tract", "patch", "band", "skymap")
    )

    fakeCat = cT.Input(
        doc="Catalog of fake sources to draw inputs from.",
        name="{fakesType}fakeSourceCat",
        storageClass="DataFrame",
        dimensions=("tract", "skymap")
    )

    imageWithFakes = cT.Output(
        doc="Image with fake sources added.",
        name="{fakesType}{coaddName}Coadd",
        storageClass="ExposureF",
        dimensions=("tract", "patch", "band", "skymap")
    )


class InsertFakesConfig(PipelineTaskConfig,
                        pipelineConnections=InsertFakesConnections):
    """Config for inserting fake sources

    Notes
    -----
    The default column names are those from the University of Washington sims database.
    """

    raColName = pexConfig.Field(
        doc="RA column name used in the fake source catalog.",
        dtype=str,
        default="raJ2000",
    )

    decColName = pexConfig.Field(
        doc="Dec. column name used in the fake source catalog.",
        dtype=str,
        default="decJ2000",
    )

    doCleanCat = pexConfig.Field(
        doc="If true removes bad sources from the catalog.",
        dtype=bool,
        default=True,
    )

    diskHLR = pexConfig.Field(
        doc="Column name for the disk half light radius used in the fake source catalog.",
        dtype=str,
        default="DiskHalfLightRadius",
    )

    bulgeHLR = pexConfig.Field(
        doc="Column name for the bulge half light radius used in the fake source catalog.",
        dtype=str,
        default="BulgeHalfLightRadius",
    )

    magVar = pexConfig.Field(
        doc="The column name for the magnitude calculated taking variability into account. In the format "
            "``filter name``magVar, e.g. imagVar for the magnitude in the i band.",
        dtype=str,
        default="%smagVar",
    )

    nDisk = pexConfig.Field(
        doc="The column name for the sersic index of the disk component used in the fake source catalog.",
        dtype=str,
        default="disk_n",
    )

    nBulge = pexConfig.Field(
        doc="The column name for the sersic index of the bulge component used in the fake source catalog.",
        dtype=str,
        default="bulge_n",
    )

    aDisk = pexConfig.Field(
        doc="The column name for the semi major axis length of the disk component used in the fake source"
            "catalog.",
        dtype=str,
        default="a_d",
    )

    aBulge = pexConfig.Field(
        doc="The column name for the semi major axis length of the bulge component.",
        dtype=str,
        default="a_b",
    )

    bDisk = pexConfig.Field(
        doc="The column name for the semi minor axis length of the disk component.",
        dtype=str,
        default="b_d",
    )

    bBulge = pexConfig.Field(
        doc="The column name for the semi minor axis length of the bulge component used in the fake source "
            "catalog.",
        dtype=str,
        default="b_b",
    )

    paDisk = pexConfig.Field(
        doc="The column name for the PA of the disk component used in the fake source catalog.",
        dtype=str,
        default="pa_disk",
    )

    paBulge = pexConfig.Field(
        doc="The column name for the PA of the bulge component used in the fake source catalog.",
        dtype=str,
        default="pa_bulge",
    )

    sourceType = pexConfig.Field(
        doc="The column name for the source type used in the fake source catalog.",
        dtype=str,
        default="sourceType",
    )

    fakeType = pexConfig.Field(
        doc="What type of fake catalog to use, snapshot (includes variability in the magnitudes calculated "
            "from the MJD of the image), static (no variability) or filename for a user defined fits"
            "catalog.",
        dtype=str,
        default="static",
    )

    calibFluxRadius = pexConfig.Field(
        doc="Aperture radius (in pixels) that was used to define the calibration for this image+catalog. "
        "This will be used to produce the correct instrumental fluxes within the radius. "
        "This value should match that of the field defined in slot_CalibFlux_instFlux.",
        dtype=float,
        default=12.0,
    )

    coaddName = pexConfig.Field(
        doc="The name of the type of coadd used",
        dtype=str,
        default="deep",
    )

    doSubSelectSources = pexConfig.Field(
        doc="Set to True if you wish to sub select sources to be input based on the value in the column"
            "set in the sourceSelectionColName config option.",
        dtype=bool,
        default=False
    )

    sourceSelectionColName = pexConfig.Field(
        doc="The name of the column in the input fakes catalogue to be used to determine which sources to"
            "add, default is none and when this is used all sources are added.",
        dtype=str,
        default="templateSource"
    )

    insertImages = pexConfig.Field(
        doc="Insert images directly? True or False.",
        dtype=bool,
        default=False,
    )

    doProcessAllDataIds = pexConfig.Field(
        doc="If True, all input data IDs will be processed, even those containing no fake sources.",
        dtype=bool,
        default=False,
    )

    trimBuffer = pexConfig.Field(
        doc="Size of the pixel buffer surrounding the image. Only those fake sources with a centroid"
        "falling within the image+buffer region will be considered for fake source injection.",
        dtype=int,
        default=100,
    )


class InsertFakesTask(PipelineTask, CmdLineTask):
    """Insert fake objects into images.

    Add fake stars and galaxies to the given image, read in through the dataRef. Galaxy parameters are read in
    from the specified file and then modelled using galsim.
    """

    _DefaultName = "insertFakes"
    ConfigClass = InsertFakesConfig

    def runDataRef(self, dataRef):
        """Read in/write out the required data products and add fake sources to the deepCoadd.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.butlerSubset.ButlerDataRef`
            Data reference defining the image to have fakes added to it
            Used to access the following data products:
                deepCoadd
                deepCoadd_fakeSourceCat
        """

        infoStr = "Adding fakes to: tract: %d, patch: %s, filter: %s" % (dataRef.dataId["tract"],
                                                                         dataRef.dataId["patch"],
                                                                         dataRef.dataId["filter"])
        self.log.info(infoStr)

        # To do: should it warn when asked to insert variable sources into the coadd

        if self.config.fakeType == "static":
            fakeCat = dataRef.get("deepCoadd_fakeSourceCat").toDataFrame()
            # To do: DM-16254, the read and write of the fake catalogs will be changed once the new pipeline
            # task structure for ref cats is in place.
            self.fakeSourceCatType = "deepCoadd_fakeSourceCat"
        else:
            fakeCat = Table.read(self.config.fakeType).to_pandas()

        coadd = dataRef.get("deepCoadd")
        wcs = coadd.getWcs()
        photoCalib = coadd.getPhotoCalib()

        imageWithFakes = self.run(fakeCat, coadd, wcs, photoCalib)

        dataRef.put(imageWithFakes.imageWithFakes, "fakes_deepCoadd")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs["wcs"] = inputs["image"].getWcs()
        inputs["photoCalib"] = inputs["image"].getPhotoCalib()

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    @classmethod
    def _makeArgumentParser(cls):
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="deepCoadd",
                               help="data IDs for the deepCoadd, e.g. --id tract=12345 patch=1,2 filter=r",
                               ContainerClass=ExistingCoaddDataIdContainer)
        return parser

    def run(self, fakeCat, image, wcs, photoCalib):
        """Add fake sources to an image.

        Parameters
        ----------
        fakeCat : `pandas.core.frame.DataFrame`
                    The catalog of fake sources to be input
        image : `lsst.afw.image.exposure.exposure.ExposureF`
                    The image into which the fake sources should be added
        wcs : `lsst.afw.geom.SkyWcs`
                    WCS to use to add fake sources
        photoCalib : `lsst.afw.image.photoCalib.PhotoCalib`
                    Photometric calibration to be used to calibrate the fake sources

        Returns
        -------
        resultStruct : `lsst.pipe.base.struct.Struct`
            contains : image : `lsst.afw.image.exposure.exposure.ExposureF`

        Notes
        -----
        Adds pixel coordinates for each source to the fakeCat and removes objects with bulge or disk half
        light radius = 0 (if ``config.doCleanCat = True``).

        Adds the ``Fake`` mask plane to the image which is then set by `addFakeSources` to mark where fake
        sources have been added. Uses the information in the ``fakeCat`` to make fake galaxies (using galsim)
        and fake stars, using the PSF models from the PSF information for the image. These are then added to
        the image and the image with fakes included returned.

        The galsim galaxies are made using a double sersic profile, one for the bulge and one for the disk,
        this is then convolved with the PSF at that point.
        """
        # Attach overriding wcs and photoCalib to image, but retain originals
        # so we can reset at the end.
        origWcs = image.getWcs()
        origPhotoCalib = image.getPhotoCalib()
        image.setWcs(wcs)
        image.setPhotoCalib(photoCalib)

        fakeCat = self.addPixCoords(fakeCat, image)
        fakeCat = self.trimFakeCat(fakeCat, image)

        if len(fakeCat) > 0:
            if isinstance(fakeCat[self.config.sourceType].iloc[0], str):
                galCheckVal = "galaxy"
                starCheckVal = "star"
            elif isinstance(fakeCat[self.config.sourceType].iloc[0], bytes):
                galCheckVal = b"galaxy"
                starCheckVal = b"star"
            elif isinstance(fakeCat[self.config.sourceType].iloc[0], (int, float)):
                galCheckVal = 1
                starCheckVal = 0
            else:
                raise TypeError("sourceType column does not have required type, should be str, bytes or int")

            if not self.config.insertImages:
                if self.config.doCleanCat:
                    fakeCat = self.cleanCat(fakeCat, starCheckVal)

                generator = self.generateGSObjectsFromCatalog(image, fakeCat, galCheckVal, starCheckVal)
            else:
                generator = self.generateGSObjectsFromImages(image, fakeCat)
            addFakeSources(image, generator, calibFluxRadius=self.config.calibFluxRadius, logger=self.log)
        elif len(fakeCat) == 0 and self.config.doProcessAllDataIds:
            self.log.warn("No fakes found for this dataRef; processing anyway.")
        else:
            raise RuntimeError("No fakes found for this dataRef.")

        # restore originals
        image.setWcs(origWcs)
        image.setPhotoCalib(origPhotoCalib)

        resultStruct = pipeBase.Struct(imageWithFakes=image)

        return resultStruct

    def generateGSObjectsFromCatalog(self, image, fakeCat, galCheckVal, starCheckVal):
        """Process catalog to generate `galsim.GSObject` s.

        Parameters
        ----------
        image : `lsst.afw.image.exposure.exposure.ExposureF`
            The image into which the fake sources should be added
        fakeCat : `pandas.core.frame.DataFrame`
            The catalog of fake sources to be input
        galCheckVal : `str`, `bytes` or `int`
            The value that is set in the sourceType column to specifiy an object is a galaxy.
        starCheckVal : `str`, `bytes` or `int`
            The value that is set in the sourceType column to specifiy an object is a star.

        Yields
        ------
        gsObjects : `generator`
            A generator of tuples of `lsst.geom.SpherePoint` and `galsim.GSObject`.
        """
        band = image.getFilterLabel().bandLabel
        wcs = image.getWcs()
        photoCalib = image.getPhotoCalib()

        self.log.info(f"Making {len(fakeCat)} objects for insertion")

        for (index, row) in fakeCat.iterrows():
            ra = row[self.config.raColName]
            dec = row[self.config.decColName]
            skyCoord = SpherePoint(ra, dec, radians)
            xy = wcs.skyToPixel(skyCoord)

            try:
                flux = photoCalib.magnitudeToInstFlux(row[self.config.magVar % band], xy)
            except LogicError:
                continue

            sourceType = row[self.config.sourceType]
            if sourceType == galCheckVal:
                bulge = galsim.Sersic(n=row[self.config.nBulge], half_light_radius=row[self.config.bulgeHLR])
                axisRatioBulge = row[self.config.bBulge]/row[self.config.aBulge]
                bulge = bulge.shear(q=axisRatioBulge, beta=((90 - row[self.config.paBulge])*galsim.degrees))

                disk = galsim.Sersic(n=row[self.config.nDisk], half_light_radius=row[self.config.diskHLR])
                axisRatioDisk = row[self.config.bDisk]/row[self.config.aDisk]
                disk = disk.shear(q=axisRatioDisk, beta=((90 - row[self.config.paDisk])*galsim.degrees))

                gal = bulge + disk
                gal = gal.withFlux(flux)

                yield skyCoord, gal
            elif sourceType == starCheckVal:
                star = galsim.DeltaFunction()
                star = star.withFlux(flux)
                yield skyCoord, star
            else:
                raise TypeError(f"Unknown sourceType {sourceType}")

    def generateGSObjectsFromImages(self, image, fakeCat):
        """Process catalog to generate `galsim.GSObject` s.

        Parameters
        ----------
        image : `lsst.afw.image.exposure.exposure.ExposureF`
            The image into which the fake sources should be added
        fakeCat : `pandas.core.frame.DataFrame`
            The catalog of fake sources to be input

        Yields
        ------
        gsObjects : `generator`
            A generator of tuples of `lsst.geom.SpherePoint` and `galsim.GSObject`.
        """
        band = image.getFilterLabel().bandLabel
        wcs = image.getWcs()
        photoCalib = image.getPhotoCalib()

        self.log.info(f"Processing {len(fakeCat)} fake images")

        for (index, row) in fakeCat.iterrows():
            ra = row[self.config.raColName]
            dec = row[self.config.decColName]
            skyCoord = SpherePoint(ra, dec, radians)
            xy = wcs.skyToPixel(skyCoord)

            try:
                flux = photoCalib.magnitudeToInstFlux(row[self.config.magVar % band], xy)
            except LogicError:
                continue

            imFile = row[band+"imFilename"]
            # Following implicitly sets image WCS from header
            im = galsim.fits.read(imFile)
            obj = galsim.InterpolatedImage(im)
            obj = obj.withFlux(flux)
            yield skyCoord, obj

    def addPixCoords(self, fakeCat, image):

        """Add pixel coordinates to the catalog of fakes.

        Parameters
        ----------
        fakeCat : `pandas.core.frame.DataFrame`
                    The catalog of fake sources to be input
        image : `lsst.afw.image.exposure.exposure.ExposureF`
                    The image into which the fake sources should be added

        Returns
        -------
        fakeCat : `pandas.core.frame.DataFrame`
        """
        wcs = image.getWcs()
        ras = fakeCat[self.config.raColName].values
        decs = fakeCat[self.config.decColName].values
        xs, ys = wcs.skyToPixelArray(ras, decs)
        fakeCat["x"] = xs
        fakeCat["y"] = ys

        return fakeCat

    def trimFakeCat(self, fakeCat, image):
        """Trim the fake cat to about the size of the input image.

        `fakeCat` must be processed with addPixCoords before using this method.

        Parameters
        ----------
        fakeCat : `pandas.core.frame.DataFrame`
                    The catalog of fake sources to be input
        image : `lsst.afw.image.exposure.exposure.ExposureF`
                    The image into which the fake sources should be added

        Returns
        -------
        fakeCat : `pandas.core.frame.DataFrame`
                    The original fakeCat trimmed to the area of the image
        """

        bbox = Box2D(image.getBBox()).dilatedBy(self.config.trimBuffer)
        xs = fakeCat["x"].values
        ys = fakeCat["y"].values

        isContained = xs >= bbox.minX
        isContained &= xs <= bbox.maxX
        isContained &= ys >= bbox.minY
        isContained &= ys <= bbox.maxY

        return fakeCat[isContained]

    def cleanCat(self, fakeCat, starCheckVal):
        """Remove rows from the fakes catalog which have HLR = 0 for either the buldge or disk component,
           also remove galaxies that have Sersic index outside the galsim min and max
           allowed (0.3 <= n <= 6.2).

        Parameters
        ----------
        fakeCat : `pandas.core.frame.DataFrame`
                    The catalog of fake sources to be input
        starCheckVal : `str`, `bytes` or `int`
                    The value that is set in the sourceType column to specifiy an object is a star.

        Returns
        -------
        fakeCat : `pandas.core.frame.DataFrame`
                    The input catalog of fake sources but with the bad objects removed

        Notes
        -----
        If the config option sourceSelectionColName is set then only objects with this column set to True
        will be added.
        """

        rowsToKeep = (((fakeCat[self.config.bulgeHLR] != 0.0) & (fakeCat[self.config.diskHLR] != 0.0))
                      | (fakeCat[self.config.sourceType] == starCheckVal))
        numRowsNotUsed = len(fakeCat) - len(np.where(rowsToKeep)[0])
        self.log.info("Removing %d rows with HLR = 0 for either the bulge or disk" % numRowsNotUsed)
        fakeCat = fakeCat[rowsToKeep]

        minN = galsim.Sersic._minimum_n
        maxN = galsim.Sersic._maximum_n
        rowsWithGoodSersic = (((fakeCat[self.config.nBulge] >= minN) & (fakeCat[self.config.nBulge] <= maxN)
                              & (fakeCat[self.config.nDisk] >= minN) & (fakeCat[self.config.nDisk] <= maxN))
                              | (fakeCat[self.config.sourceType] == starCheckVal))
        numRowsNotUsed = len(fakeCat) - len(np.where(rowsWithGoodSersic)[0])
        self.log.info("Removing %d rows of galaxies with nBulge or nDisk outside of %0.2f <= n <= %0.2f" %
                      (numRowsNotUsed, minN, maxN))
        fakeCat = fakeCat[rowsWithGoodSersic]

        if self.config.doSubSelectSources:
            try:
                rowsSelected = (fakeCat[self.config.sourceSelectionColName])
            except KeyError:
                raise KeyError("Given column, %s, for source selection not found." %
                               self.config.sourceSelectionColName)
            numRowsNotUsed = len(fakeCat) - len(rowsSelected)
            self.log.info("Removing %d rows which were not designated as template sources" % numRowsNotUsed)
            fakeCat = fakeCat[rowsSelected]

        return fakeCat

    def _getMetadataName(self):
        """Disable metadata writing"""
        return None
