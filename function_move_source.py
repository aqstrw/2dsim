# tsilia@mail.strw.leidenuniv.nl  
# 2020/07/27

import astropy
from astropy.io import fits

import numpy as np
import math

from scipy.ndimage.interpolation import shift

def move_source(data,vel,theta,refmask):
    
	"""
	Extract a fixed source of light from a moving one (or the other way around). Calculate the differences between consequent frames, shift them, 		and add them up to effectively collapse a source back on itself (or smear a source across the field of view).
	"""
	data = np.where(refmask == 1, np.nan, data)
	nz=data.shape[1]-1 # number of frames within the ramp - 1
	ny=data.shape[2] # number of pixels along y axis (1024 by default)
	nx=data.shape[3] # number of pixels along x axis (1032 by default)
	exp=data.shape[0] # number of exposures 
	
	velx=vel*math.cos(theta) # velocity along x-axis
	vely=vel*math.sin(theta) # velocity along y-axis
	
	diff=np.empty([exp,nz,ny,nx]) # create empty array of differences
	
	# calculate differences between frames
	for w in range(exp):
		for k in range(nz): 
			diff[w,k,:,:]=data[w,k+1,:,:]-data[w,k,:,:] 

	# shift array of differences - use interpolation
	for w in range(exp):
		for k in range(nz):
			diff[w,k,:,:]=shift(diff[w,k,:,:], (vely*(k+1),velx*(k+1)), order=3, mode='constant', cval=math.nan, prefilter=False)
	
	toy=np.empty([exp,nz+1,ny,nx]) # create empty array of model 
	
	# create first frame 
	for w in range(exp):
		toy[w,0,:,:]=data[w,0,:,:] 

	# create subsequent frames 
	for w in range(exp):    
		for k in range(1,nz+1): 
			toy[w,k,:,:]=toy[w,k-1,:,:]+diff[w,k-1,:,:]
	
	mv_source_final=np.array(toy[:,:,:])
	
	return mv_source_final # return results
	

