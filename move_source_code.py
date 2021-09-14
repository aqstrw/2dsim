# tsilia@mail.strw.leidenuniv.nl  
# 2020/07/27

import astropy
from astropy.io import fits

import numpy as np
import math

import function_move_source as fl

def move(inputdir, outputdir, velin, refmask, stattomoving = True):
    pixel_scale=0.11 # pixel scale of the instrument in arcsec/pixel
    t_1=2.775 # time between individual frames in sec

    flag_move_source=stattomoving # True if script is used to generate a moving source out of a stationary one
    input1=inputdir # path of the input file
    output1=outputdir # path of the output file
    vel_arc=0.010996 * velin # velocity in arcsec/sec 
    theta=math.pi * 1.25 # angle (rads) between source's direction and positive x-axis  

    vel=vel_arc*(t_1/pixel_scale) # convert velocity from arcsec/sec to pixels/frame

    # open FITS file 
    hdul=fits.open(input1)
    data=hdul[1].data

    if flag_move_source: 
        theta=theta+math.pi

    # run function
    mv_source_final=fl.move_source(data,vel,theta,refmask)
    hdul[1].data=mv_source_final

    # save results
    hdul.writeto(output1, overwrite = True)
    
    hdul.close()

