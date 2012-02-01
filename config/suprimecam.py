root.doWriteIsr = False
root.isr.methodList=["doConversionForIsr", "doSaturationDetection",
                     "doOverscanCorrection", "doVariance", "doFlatCorrection"]
root.isr.doWrite = False

root.calibrate.repair.doCosmicRay = True
root.calibrate.repair.cosmicray.nCrPixelMax = 100000

root.calibrate.measurePsf.starSelector.name = "secondMoment"
root.calibrate.measurePsf.psfDeterminer.name = "pca"

root.calibrate.photometry.detect.thresholdValue = 50.0
root.calibrate.photometry.measure.source.astrom = "NAIVE"
root.calibrate.photometry.measure.source.apFlux = "NAIVE"
root.calibrate.photometry.measure.source.modelFlux = "GAUSSIAN"
root.calibrate.photometry.measure.source.psfFlux = "PSF"
root.calibrate.photometry.measure.source.shape = "SDSS"
root.calibrate.photometry.measure.astrometry.names = ["GAUSSIAN", "NAIVE", "SDSS"]
root.calibrate.photometry.measure.shape.names = ["SDSS"]
root.calibrate.photometry.measure.photometry.names = ["NAIVE", "GAUSSIAN", "PSF", "SINC"]
root.calibrate.photometry.measure.photometry["NAIVE"].radius = 7.0
root.calibrate.photometry.measure.photometry["GAUSSIAN"].shiftmax = 10
root.calibrate.photometry.measure.photometry["SINC"].radius = 7.0

root.calibrate.apCorr.alg1.name = "PSF"
root.calibrate.apCorr.alg2.name = "SINC"
root.calibrate.apCorr.alg1[root.calibrate.apCorr.alg1.name] = root.calibrate.photometry.measure.photometry[root.calibrate.apCorr.alg1.name]
root.calibrate.apCorr.alg2[root.calibrate.apCorr.alg2.name] = root.calibrate.photometry.measure.photometry[root.calibrate.apCorr.alg2.name]

root.photometry.measure.source.astrom = "NAIVE"
root.photometry.measure.source.apFlux = "NAIVE"
root.photometry.measure.source.modelFlux = "GAUSSIAN"
root.photometry.measure.source.psfFlux = "PSF"
root.photometry.measure.source.shape = "SDSS"
root.photometry.measure.astrometry.names = ["GAUSSIAN", "NAIVE", "SDSS"]
root.photometry.measure.shape.names = ["SDSS"]
root.photometry.measure.photometry.names = ["NAIVE", "GAUSSIAN", "PSF", "SINC"]
root.photometry.measure.photometry["NAIVE"].radius = 3.0
root.photometry.measure.photometry["GAUSSIAN"].shiftmax = 10
root.photometry.measure.photometry["SINC"].radius = 3.0

