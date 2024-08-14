# IMPORTS #################################################################
import os
from astropy.io import fits
import numpy as np
##########################################################################

# setup o/p #############################################################
verbosity = 2
save = 0
plot = 1
##########################################################################

# HANDLE MIRISIM ENV VARIABLES ############################################
os.environ["CRDS_PATH"] = "/home/aqwork/miriroot/crds"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
os.environ["CRDS_CONTEXT"] = "jwst_1256.pmap"

# printing environment variables
if verbosity:
    print("MIRISIM ROOT : ",os.environ['MIRISIM_ROOT'])
    print("PYSYN CDBS : ",os.environ['PYSYN_CDBS'])
    print("CDP DIR : ",os.environ['CDP_DIR'])
    print("CRDS PATH : ",os.environ['CRDS_PATH'])
    print("SRDS SERVER URL : ",os.environ['CRDS_SERVER_URL'])
    print("SRDS CONTEXT : ",os.environ['CRDS_CONTEXT'])
##########################################################################

# setup directories and paths ############################################
mydir = os.getcwd()
mirisimopsource = '/20210830_144816_mirisim/det_images/'       # using output in /20210830_144816_mirisim , change mirisimopsource variable to switch folder
ipname = "det_image_seq1_MIRIMAGE_F1280Wexp1.fits"             # using "det_image_seq1_MIRIMAGE_F1280exp1.fits" input file for msc, change ipname if other
mirisimopdir = mydir + mirisimopsource
if verbosity:
    print("using simulated op stored in : ",mirisimopdir)

# input file path
ippath = mirisimopdir + ipname
if verbosity:
    print("using simulated file named : ",ippath)

# # path creation moved to runpip
# # output directory setup
# my_output_dir = "pipeline_outputs/"
# if not os.path.exists(my_output_dir): 
#     os.mkdir(my_output_dir)
# print("moved op stored in : ",mirisimopdir)
##########################################################################

# mask reference file exploration ########################################

# get pdq data
maskrefpath = '/home/aqwork/miriroot/crds/references/jwst/miri/jwst_miri_mask_0023.fits'
if verbosity:
    print("using mask named : ", maskrefpath)
hdumask = fits.open(maskrefpath)
maskpdq = hdumask[1].data
if verbosity:
    print("The mask from ref files has size : ",maskpdq.shape)
hdumask.close()

# extract dnu pixels and plot
dnupix = np.where(np.bitwise_and(maskpdq,int(2**0)),1,0)
# fig,axs = plt.subplots(figsize = (10,10))
# axs.imshow(dnupix,origin='lower',cmap='Greys',interpolation='nearest')

# apply refmask and plot first frame
refmask = np.tile(dnupix,(1,40,1,1))
hdut = fits.open(ippath)
data = hdut[1].data
ndata = np.where(refmask == 1, np.nan, data)
##########################################################################

# exec setup #############################################################
moverate = 10
wasteroid = 0.010996
t_1 = 2.775
hdu = fits.open(ippath)
frames = hdu[1].shape[1]
hdu.close()
pixelscale = 0.11
base = 416

# runpip setup ##########################################################

# reuse = 1
# ep = 5
# core = 50
# span = 1

# reuse = 1
# ep = 17
# core = 150
# span = 1

reuse = 1
ep = 5
core = np.around(2.5 * np.pi * (ep**2) * 0.7).astype(int)
span = 1