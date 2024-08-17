# Description of current tool will be updated #

--------------------------------------------------------------------------------------

## Background

NASAs James Webb Space Telescope (JWST) Parked at the L2 Legrange point, is a state of the art Infrared space observatory,
the most advanced ever built and launched into space (as of August 2024). the JWST will use its 4 instruments to collect
5 years of science observations (based on mission requirements). The most relevant instrument for our research is the
Mid Infrared Instrument (MIRI) which has three detectors, two used as part of the Medium Resolution Spectrometer
(MRS) and one used for imaging (MIRIM: MIRI Imaging Module). While science is being conducted using the MRS, MIRIM has the
capability to simultaneously capture the field of view (by default and as the recommended mode) and could inadvertently
detect a considerable number of asteroids. This data could be accessible without burdening the instrument time and perhaps
without any extra costs. Husárová’s 2018 MSc thesis [23] estimates a 4.1% chance of interception (chance of asteroids
intercepting the field of view) and a nearly guaranteed (96%) detection in case of an intercept, translating to 916,570
detections (733,256 unknown) over the lifetime of the JWST. This would amount to a considerable addition to any asteroid
catalogue.

One of the advantages that the JWST IR detectors offer, compared to a CCD camera is that these instruments record the
build-up of charge in each pixel with time. These non-destructive readouts are eventually translated to an image which
depicts the average charge build-up rate (slope) for each each pixel, by excluding frames with outliers to produce the calibrated
outputs (extreme simplification of the calibration process). This image is a 2D representation of the higher
dimensional data that we collect in each exposure and is produced by the JWST calibration pipeline in 3 stages. The tool 
demonstrated here interacts with intermediate optional outputs in the the stage 1 pipeline. The Stage 1 calibrations are
detector-level corrections, they are applied in nearly all cases.

The "jump" step flags outliers in the linear charge build-up in a pixel, which are higher than a specified sigma threshold
value (default is 4.0). This check is conducted for each integration within an exposure and the identified jumps are indicated by
setting the pixel’s "JUMP_DET" flag in the corresponding GROUPDQ array. The detection of jumps uses the two-point
difference method described by Anderson and Gordon 2011. The "ramp_fitting" step uses flags from the "jump"
and "saturation" steps to calculate a cosmic ray and saturation free count rate (average counts per second) for each
pixel. This is determined by applying a linear ordinary least squares fit to the remaining data after removing the flagged
segments. The fit uses the ’optimal’ weighing scheme described by Fixsen et al [34]. The module produces two outputs,
a "rate" file with the slope at each pixel averaged over all integrations, and a "rateints" file with slope images for each
integration. A third optional output is available by setting the "save_opt" keyword of the "ramp_fitting" step to "True"
(default is False). These optional outputs and the flexibility afforded to us by the pipeline are necessary tools for any
bespoke analysis of JWST data.

We initially tested our automated asteroid detection product using simulations of asteroids using a combination of the 
MIRISim package and translation and interpolation tool developed by Tsilia in her 2020 MSc thesis. With JWST science beginning
to be released to public now, it is time to test the tool with real JWST MIRI exposures and comment on it's utility / efficacy.

-------------------------------------------------------------

## Useful Links
RUNNING THE PIPELINE: https://jwst-pipeline.readthedocs.io/en/latest/jwst/user_documentation/running_pipeline_python.html
PIPELINE PRODUCTS: https://jwst-pipeline.readthedocs.io/en/stable/jwst/data_products/science_products.html#ramp-data-ramp
DQ FLAG BITS: https://jwst-pipeline.readthedocs.io/en/latest/jwst/references_general/references_general.html#data-quality-flags
Others:
https://jwst-pipeline.readthedocs.io/en/stable/jwst/pipeline/main.html#pipelines
https://jwst-pipeline.readthedocs.io/en/latest/jwst/user_documentation/available_pipelines.html
https://jwst-pipeline.readthedocs.io/en/latest/jwst/dq_init/description.html
