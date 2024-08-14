# imports ###############################################################

from myglobals import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from sklearn.cluster import DBSCAN
from scipy.optimize import curve_fit
from scipy.stats import mode
import plotly.express as px

# move source imports ###########################################################

import move_source_code as msc
import change_pixeldq as pixel
from jwst.pipeline import Detector1Pipeline
from jwst import datamodels

# import the dqinit step
from jwst.dq_init import dq_init_step
# print the description and options
print(dq_init_step.DQInitStep.__doc__)
print(dq_init_step.DQInitStep.spec)

# import the jump step
from jwst.jump import jump_step

# print the description and options
print(jump_step.JumpStep.__doc__)
print(jump_step.JumpStep.spec)

# import the rampfitting step
from jwst.ramp_fitting import ramp_fit_step
# print the description and options
print(ramp_fit_step.RampFitStep.__doc__)
print(ramp_fit_step.RampFitStep.spec)

########################################################################

def sparseGDQ_to_3Dplot(data,plottitle,zscale = 10, fig = None):
    '''
    Takes sparse gdq data and plots all clusters
    '''
    # get axes ranges
    _,zr,yr,xr = data.shape
    # get coords from sparse GDQ data
    coords = np.nonzero(data[0])
    
    if fig is None:
        # plot
        fig = px.scatter_3d(
            x = coords[2],
            y=coords[1],
            z=coords[0]*zscale,
            range_x = [0,xr],
            range_y = [0,yr],
            # range_z = [0,zr]*zscale,
            title = plottitle,
            size_max=2,
            opacity=0.005,
            height = 800,
            width = 1200,
        )
        fig.update_scenes(
            zaxis = dict(
                tickmode = 'array',
                tickvals = np.arange(0,(zr*zscale)+1,5*zscale),
                ticktext = np.arange(0,zr+1,5).astype(str),
                ),
            aspectmode='data',
        )
    else:
        # add trace to existing figs
        print("not implemented yet")
    fig.show()


def runpip(use_existing_file,ip, msc_output_prefix = "pipeline_outputs"):
    '''
        Takes in locations to the mirisim simulated scene and produces a moving target based on moverate and saves based in the provided_directory+xmoverate
        ip - the input for msc
    '''

    # moverate = mvx # sets the rate of motion
    
    msc_op_dir = msc_output_prefix+'x'+str(moverate)+'/'        # op directory for msc
    msc_op_file = msc_op_dir+ipname[0:-5]+'_moved.fits'        # op file for msc
    
    # if op directory does not exist, create it
    if not os.path.exists(msc_op_dir):
        if verbosity:
            print('msc_op_dir does not exist, creating')
        os.mkdir(msc_op_dir)

    # if path exists and using old file, then pull data
    if use_existing_file:
        try:
            if verbosity:
                print ('Trying to use existing file')
            dm_slopes =  datamodels.open(msc_op_dir+ipname[0:-5]+'_rate.fits')
            if verbosity:
                print('first output contains datamodel')

            if verbosity:
                print('second output contains path to\n')
            hdumoved = fits.open(msc_op_file)
            if verbosity == 2:
                hdumoved.info()
            hdumoved.close()

            if verbosity:
                print('third output contains path to\n')
            opramp_path = msc_op_dir+ipname[0:-5]+'_ramp.fits'
            hduramp = fits.open(opramp_path)
            if verbosity == 2:
                hduramp.info()
            hduramp.close()

            if verbosity:
                print('fourth output contains moverate: ',moverate)
            return(dm_slopes,msc_op_file,opramp_path)
        except OSError as e:
            print("caught exception : ",e)
            print("missing file(s), choose among existing simulation outputs (moverate = 0.05, 5, 10, 20) or fix mirisim and move_source_code workflows")
            print("exiting prematurely")
            # print('missing file(s), running custom pipeline') # old workflow
            return(0,0,0)

    print('############################## running msc ##################################\n\n\n')

    # get refmask 
    refmask = np.tile(dnupix,(1,40,1,1))
    
    # run the msc code
    msc.move(ip, msc_op_file, refmask = refmask, velin = moverate)

    # use same input filename as above op file name
    '''If individual steps are executed without an output file name specified via the output_file parameter, 
    the stpipe infrastructure automatically uses the input file name as the root of the output file name and 
    appends the name of the step as an additional suffix to the input file name. If the input file name already 
    has a known suffix, that suffix will be replaced'''
    dqi_op_file = msc_op_dir+ipname[0:-5]+'_moved_dqinitstep.fits'        # op file for msc

    # os.environ["STPIPE_DISABLE_CRDS_STEPPARS"] = 'True'
    # run the dqinit step
    dm_dqinit = dq_init_step.DQInitStep.call(msc_op_file, output_use_model=True, output_dir = msc_op_dir, save_results=True,)

    cpdq_op_file = msc_op_dir+ipname        # op file for cpdq
    pixel.cpdq(dqi_op_file, cpdq_op_file)

    # remaining pipeline
    dm_slopes = Detector1Pipeline.call(cpdq_op_file, output_use_model=True, save_results=True, 
                                output_dir=msc_op_dir, save_calibrated_ramp = True, steps={'dq_init': {'skip':True},'ipc': {'skip': True},
                                                                'refpix': {'skip': True}, 'ramp_fit': {'save_opt': True}})

    opramp_path = msc_op_dir+ipname[0:-5]+'_ramp.fits'
    # os.environ["STPIPE_DISABLE_CRDS_STEPPARS"] = 'False'
    return dm_slopes, msc_op_file, opramp_path

def excol(base = 416):
    '''The asteroids movement simulation process causes the size of the frame to be extended. The op nancols is the new size. '''
    extracols = np.ceil(wasteroid * moverate * np.cos((np.pi)/4) * t_1 * frames / pixelscale)
    nancols = int(base + extracols)
    print('extracols added and nancols :',extracols, nancols)
    return(nancols)

# apply masks
def maskandextract(ramps, gdq, bit = 2):
    '''
        Takes in ramps and gdq, masks them and returns crhits (val = 2) flattened and 3d
        ip: ramps
            gdp
            moverate
            val - the bit to extract (2 for jump), more avalable at https://jwst-pipeline.readthedocs.io/en/latest/jwst/references_general/references_general.html#data-quality-flags
        op: mcrramps - masked cosmic ray hits
            mramps - masked ramps
            mgdq - masked gdq
            mcrhitscum - flattened
            xdata - flattened
            ydata - flattened
            '''
    cols = excol()
    print('nancols and base :',cols,base)
    rows = cols - base
    mramps = ramps[:,:,5:-5,cols:-10]
    mgdq = gdq[:,:,5:-5,cols:-10]
    mcrramps = np.bitwise_and(mgdq,int(2**bit))
    mcrramps = np.bitwise_and(mgdq,int(2**bit))
    if save == 2:
        hduds9 = fits.open(ippath)
        hduds9[1].data = mcrramps
        hduds9.writeto("fords9x"+str(moverate)+".fits", overwrite=True)
        hduds9.close()
    mcrhitscum = np.sum(mcrramps,axis = 1)[0]
    xdata = np.sum(mcrhitscum,axis = 0)
    ydata = np.sum(mcrhitscum,axis = 1)
    return(mcrramps,mramps,mgdq,mcrhitscum,xdata,ydata)

def build_arrays_for_DBSCAN(crramp,eps = 5):
    '''
    Gets non-zero coordinates and scales using eps 
    '''
    if verbosity:
        print('The cosmic ray hits contain the following unique values', np.unique(crramp))

    mcrhitramps_nonzero_coords = np.nonzero(crramp)

    if plot == 2:
        sparseGDQ_to_3Dplot(crramp, "All detected clusters")

    sampst = np.vstack((mcrhitramps_nonzero_coords[0],mcrhitramps_nonzero_coords[1],mcrhitramps_nonzero_coords[2],mcrhitramps_nonzero_coords[3]))
    sampsst = np.vstack((mcrhitramps_nonzero_coords[0],mcrhitramps_nonzero_coords[1]*(eps/2),mcrhitramps_nonzero_coords[2],mcrhitramps_nonzero_coords[3]))
    samps_scaled = sampsst.T
    samps = sampst.T
    if verbosity == 2:
        print('shape of samps array: ',samps.shape)
        print('shape of samps_scaled array: ',samps_scaled.shape)
    return(samps,samps_scaled)

def apply_DBSCAN(crramps, scaled_samples, eps = 5, core = 180):

    # cluster the crhits
    if verbosity:
        print("running DBSCAN")
    db = DBSCAN(eps=eps, min_samples = core).fit(scaled_samples)
    _,zr,yr,xr = crramps.shape
    clusters = np.unique(db.labels_)
    clusters = clusters[clusters != -1]
    if verbosity ==2 :
        print("array of cluster labels : ",clusters)
    no_clusters = len(np.unique(db.labels_))
    no_noise = np.sum(np.array(db.labels_) == -1, axis=0)
    if verbosity:
        print('Estimated no. of clusters: %d' % no_clusters,' with labels :', np.unique(db.labels_))
        print('Estimated no. of noise points: %d' % no_noise)
        print('number of samples categorised: ',db.labels_.shape)
    return db, clusters, no_clusters, no_noise

def identify_asteroid_clusters(db,clusters,unscaled_samples):
    # if span == 1: # flag that prompts checks for continuity throughout all frames

    # clusters_per_frame = [for frame in unscaled_samples]

    # the 'll' variable ideally should be automatically determined based on simulation speed
    if moverate == 30:
        ll = 9
    else:
        ll = 4
    labels_filtered_for_asteroids = db.labels_
    for cname in clusters:
        locs = unscaled_samples[np.where(db.labels_ == cname)]
        frame_range_min, frame_range_max = min(locs[:,1]), max(locs[:,1])
        print('range of frames for cluster',cname,':',frame_range_min, frame_range_max )
        if frame_range_min >= ll or frame_range_max <= 34:
            if verbosity:
                print('cluster being dropped due to min/max frame criteria:', cname)
            labels_filtered_for_asteroids = np.where(labels_filtered_for_asteroids == cname, -1, labels_filtered_for_asteroids)
        else:
            if verbosity:
                print('cluster passed asteroid filter', cname)
    astclusters = np.unique(labels_filtered_for_asteroids)
    astclusters = np.delete(astclusters,np.where(astclusters==-1))
    return (astclusters, labels_filtered_for_asteroids)

def plot_3dclusters(coords, labels, box_shape, isolate= 0, **plot_kwargs):
    '''
    coords: unscaled coords
    labels: same length as coords and 4D
    '''
    _,zr,yr,xr = box_shape
    dfcoords = pd.DataFrame(np.hstack([coords,labels[:,np.newaxis]]),columns = ["integration","z","y","x","clusters"])

    dfcoords["clusters"] =  dfcoords["clusters"].astype(str)
    if isolate:
        fig = px.scatter_3d(dfcoords, x = "x", y = "y", z = "z", color = 'clusters',**plot_kwargs)
        fig.update_scenes(
            zaxis = dict(
                tickmode = 'array',
                tickvals = np.arange(0,zr+1,5),
                ticktext = np.arange(0,zr+1,5).astype(str),
                ),
            aspectmode='cube',
        )
    else:
        fig = px.scatter_3d(dfcoords, x = "x", y = "y", z = "z", color = 'clusters', range_x = [0,xr], range_y = [0,yr],**plot_kwargs)
        fig.update_scenes(
            zaxis = dict(
                tickmode = 'array',
                tickvals = np.arange(0,zr+1,5),
                ticktext = np.arange(0,zr+1,5).astype(str),
                ),
            aspectmode='manual',
            aspectratio = {'x':1,'y':2,'z':0.7}
        )
    return fig

def highlight_core(coords,db, asteroid_clusters, ready_fig):
    '''
    coords: unscaled, will be converted to core coords using model
    asteroid_clusters: clusters identified as asteroids
    ready_fig: fig where plot needs to be added
    '''
    core_labels = db.labels_[db.core_sample_indices_]
    asteroidcluster_core_labels = [label for label in core_labels if label in asteroid_clusters]
    asteroids_core_indices = db.core_sample_indices_[[label in asteroid_clusters for label in core_labels]]
    core_coords = np.around(coords[asteroids_core_indices]).astype(int)
    dfcoords = pd.DataFrame(np.hstack([core_coords,np.array(asteroidcluster_core_labels)[:,np.newaxis]]), columns = ["integration","z", "y", "x","clusters"])
    dfcoords["clusters"] =  dfcoords["clusters"].astype(str)

    if verbosity:
        print("core index array has length : ",core_coords.shape)
    ready_fig.add_traces(list(px.scatter_3d(dfcoords, x = "x", y = "y", z = "z", color = "clusters",color_discrete_sequence= ['black'], opacity = 1).select_traces()))
    return ready_fig,core_coords

def fig_update_traceprops(fig_for_update, marker_size = 3):
    return fig_for_update.update_traces(
        marker=dict(
            size=marker_size,
        )
    )

def prep_isolated_asteroid_fig_params(sampss, ast_labels):
    isolated_asteroid_labels = ast_labels[np.invert(ast_labels == -1)]
    isolated_asteroid_sampss = sampss[np.invert(ast_labels == -1)]
    return isolated_asteroid_labels, isolated_asteroid_sampss

def get_leastsquarefit_from_svd(core_coords):
    # Calculate the mean of the points, i.e. the 'center' of the cloud
    meancoords = core_coords[:,1:].mean(axis=0)
    # Do an SVD on the mean-centered data.
    uu, dd, vv = np.linalg.svd(core_coords[:,1:] - meancoords)
    return vv, meancoords

def line_3d_coords(slopes,means,distance_from_center = 100):
    linepts = slopes[0] * np.mgrid[-distance_from_center:distance_from_center:2j][:, np.newaxis]
    linepts += means
    return linepts

def add_line_to_plot(fig_with_asteroid, df_with_coords_and_color):
    fig_with_asteroid.add_traces(list(px.line_3d(df_with_coords_and_color, x = 'x',y = 'y', z='z',color = 'bestfit', height = 500, width = 500).select_traces()))
    return fig_with_asteroid

def calculate_angular_velocity(lsf_params, pixelscale, frametime, round_to = 10):
    # optimal slope multiplied by pixel scale and divided by frame time
    calculated_v = round((np.linalg.norm(lsf_params[0][1:])/lsf_params[0][0])*(pixelscale/frametime),round_to)
    input_v = wasteroid*moverate
    return calculated_v,input_v,calculated_v-input_v,(calculated_v-input_v)/input_v

########################## from old file ########################################
# def curve(x,m,c):
#     return((m*x)+c)

# def fit2dline(func,f,x, iclust, loop = 'na:'):

#     popt, pcov = curve_fit(curve, f,x)

#     print(loop+' for cluster'+str(iclust)+' optimal slope multiplied by pixel scale and divided by frame time :', round(popt[0]*0.11/2.775,10))
    
#     return (popt,pcov)

# def fit2dlinecent(func,f,x, iclust, loop = 'na:', mean = 1):
#     lx = []
#     lf = []
#     for fs in np.unique(f):
#         xs = x[f==fs]
#         if mean:
#             xcentarr=np.mean(xs)
#             lx.append(xcentarr)
#         else:
#             xcentarr=mode(xs)[0]
#             lx.append(xcentarr[0])
#         lf.append(fs)
        
# #     for i,x in enumerate(lx):
# #         if x == 272:
# #             lx[i] = 268
            
# #     print(lx,lf)
    
#     popt, pcov = curve_fit(curve, lf,lx)

#     print(loop+' for cluster'+str(iclust)+' optimal slope multiplied by pixel scale and divided by frame time :', round(popt[0]*0.11/2.775,10))
    
#     return (popt,pcov,lx,lf)

# def plot3dgdq(data,cluster = 99,useaxes = (3,2,1), figdims =(10,10), fname = 'demo', az = -60, el = 30, at =  (0,0,0), xd = 0, setat = 0):
#     # pass flag 99 to print multiple groups
#     fig = plt.figure(figsize=figdims)
#     ax = fig.add_subplot(1, 1, 1, projection='3d', azim = az, elev = el)
#     if at == (0,0,0):
#         at = (1,1.5,0.5)
#     if cluster == 99:
#         for cluster in np.unique(data):
#             # if cluster is not nan, set
#             if np.isnan(cluster) == 0:
#                 locs = np.where(data  == cluster)
#                 if cluster == -1:
#                     plotlab = 'Noise'
#                 else:
#                     plotlab = 'cluster '+str(int(cluster))
#                 ax.scatter(locs[useaxes[0]],locs[useaxes[1]],locs[useaxes[2]], label = plotlab,marker = 'o',alpha = 0.1)
#         ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
#         ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
#         ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
#         ax.set_box_aspect(at)
#         print ('99 aspect used is :', at)
#         ax.zaxis.set_rotate_label(False)
#         if az == -90:
            
#             ax.set_xlabel('\n \nColumns\n \n', fontsize=14)
# #             ax.set_ylabel('Rows', fontsize=14, rotation = 45)
#             ax.set_zlabel('\n \nFrames', fontsize=14, rotation = 90)
#             for line in ax.yaxis.get_ticklabels():
#                 line.set_visible(False)
#             for line in ax.xaxis.get_ticklabels():
#                 line.set_visible(False)
# #             fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0)
#             fig.tight_layout()
#             fig.savefig(fname+".png",format='png')
#             fig.savefig(fname+".pdf",format='pdf')
#             fig.show()
#         elif az == 0:
            
# #             ax.set_xlabel('Columns', fontsize=14)
#             ax.set_ylabel('\n \nRows', fontsize=14, rotation = 45)
#             ax.set_zlabel('Frames\n \n', fontsize=14, rotation = 90)
#             for line in ax.xaxis.get_ticklabels():
#                 line.set_visible(False)
# #             fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0)
#             fig.tight_layout()
#             fig.savefig(fname+".png",format='png')
#             fig.savefig(fname+".pdf",format='pdf')
#             fig.show()
#         else:
#             ax.set_xlabel('\n \nColumns\n \n', fontsize=14)
#             ax.set_ylabel('\n \nRows\n \n', fontsize=14, rotation = 45)
#             ax.set_zlabel('\n \nFrames\n \n', fontsize=14, rotation = 90)
# #             for line in ax.yaxis.get_ticklines():
# #                 line.set_visible(False)
# #             fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0)
#             fig.tight_layout()
#             fig.savefig(fname+".png",format='png')
#             fig.savefig(fname+".pdf",format='pdf')
#             fig.show()
#     else:
#         locs = np.where(data == cluster)
# #         print(locs)
#         # if setat:
#         #     at = (max(locs[useaxes[0]]), max(locs[useaxes[1]]), 0.5*(max(max(locs[useaxes[0]]), max(locs[useaxes[1]]))))
#         # print ('aspect used is :', at)
#         # ax.set_box_aspect(at)
#         # ax.scatter(locs[useaxes[0]],locs[useaxes[1]],locs[useaxes[2]], label = 'cluster_'+str(cluster), c = 'red', marker = 'o',alpha = 0.1)
#         _, xr,yr,zr = data.shape
#         fig_pl = px.scatter_3d(x = locs[useaxes[0]], y = locs[useaxes[1]], z = locs[useaxes[2]],
#                                range_x = [0,xr], range_y = [0,yr], range_z = [0,zr],
#                                size_max=18, opacity=0.2, height = 800, width = 800)
#         fig_pl.show()
#         # ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
#         # ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
#         # ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
#         # ax.zaxis.set_rotate_label(False)
# #         if az == -90:
# #             ax.set_xlabel('\n \nColumns', fontsize=14)
# # #             ax.set_ylabel('\n \nRows', fontsize=14, rotation = 45)
# #             ax.set_zlabel('\n \nFrames', fontsize=14, rotation = 90)
# # #             fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0)
# #             for line in ax.yaxis.get_ticklabels():
# #                         line.set_visible(False)
# #             fig.tight_layout()
# #             fig.savefig(fname+".png",format='png')
# #             fig.savefig(fname+".pdf",format='pdf')
# #             fig.show()
# #         elif az == 0:
# # #             ax.set_xlabel('\n \nColumns', fontsize=14)
# #             ax.set_ylabel('\n \nRows', fontsize=14, rotation = 45)
# #             ax.set_zlabel('Frames\n \n', fontsize=14, rotation = 90)
# # #             fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0)
# #             fig.tight_layout()
# #             for line in ax.xaxis.get_ticklabels():
# #                 line.set_visible(False)
# #             fig.savefig(fname+".png",format='png')
# #             fig.savefig(fname+".pdf",format='pdf')
# #             fig.show()
# #         else:
# #             ax.set_xlabel('\n \nColumns\n \n', fontsize=14)
# #             ax.set_ylabel('\n \nRows\n \n', fontsize=14, rotation = 45)
# #             ax.set_zlabel('\n \nFrames\n \n', fontsize=14, rotation = 90)
# # #             fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0)
# #             fig.tight_layout()
# # #             for line in ax.xaxis.get_ticklines():
# # #                 line.set_visible(False)
# #             fig.savefig(fname+".png",format='png')
# #             fig.savefig(fname+".pdf",format='pdf')
# #             fig.show()
            
#     plt.legend()
#     # return(fig,ax)


# def DBSCAN_and_checks(crramps, unscaled_samples, scaled_samples, eps = 5, core = 180, grpno = 0,span = 0, vel = 0):

#     # cluster the crhits
#     db = DBSCAN(eps=eps, min_samples = core).fit(scaled_samples)
#     _,zr,yr,xr = crramps.shape
#     labels = db.labels_
#     clusters = np.unique(labels)
#     clusters = clusters[clusters != -1]
#     if verbosity == 2 :
#         print("array of cluster labels : ",clusters)
#     no_clusters = len(np.unique(labels))
#     no_noise = np.sum(np.array(labels) == -1, axis=0)
#     if verbosity:
#         print('Estimated no. of clusters: %d' % no_clusters,' with labels :', np.unique(labels))
#         print('Estimated no. of noise points: %d' % no_noise)
#         print('number of samples categorised: ',labels.shape)
#     if span == 1: # flag that prompts checks for continuity throughout all frames, most likely not required for real data
#         # this variable ideally should be automatically determined based on simulation speed
#         if vel == 30:
#             ll = 9
#         else:
#             ll = 4
#         for cname in clusters:
#             locs = unscaled_samples[np.where(labels == cname)]
#             print('range of frames for cluster',cname,':',min(locs[:,1]), max(locs[:,1]))
#             if min(locs[:,1]) >= ll:
#                 if verbosity:
#                     print('cluster being dropped due to min frame criteria:', cname)
#                 labels = np.where(labels == cname, -1, labels)
#             elif max(locs[:,1]) <= 34: # probably needs to be automatically determined based on speed
#                 if verbosity:
#                     print('cluster being dropped due to max frame criteria:', cname)
#                 labels = np.where(labels == cname, -1, labels)
#             else:
#                 if verbosity:
#                     print('asteroid detected in cluster', cname)
#     astclusters = np.unique(labels)

#     # make sure that all clusters weren't dropped before moving forward with plots etc.
#     if len(astclusters) == 1:
#         if verbosity:
#                     print('no asteroid found')
#         return(0,0,0)
#     else:
#         print('clusters representing asteroids :', astclusters[1:])
#         lcls = []
#         lsls = []
#         for grp in astclusters:
#             complabelledset,speclabelledset = lsets(crramps, unscaled_samples, scaled_samples, labels, grpno = grp)
#             lcls.append(complabelledset)
#             lsls.append(speclabelledset)
#     return(lcls,lsls,astclusters)

# def lsets(crramps, unscaled_samples, scaled_samples, labels, grpno = 0):
#     complabelledset = np.copy(crramps).astype('float64')
#     complabelledset = np.where(complabelledset == 0, np.nan, complabelledset)
#     speclabelledset = np.copy(crramps).astype('float64')
#     for num,val in enumerate(unscaled_samples):
# #         print(val)
#         complabelledset[val[0],val[1],val[2],val[3]] = labels[num]
#     speclabels = np.where(labels == grpno,1,np.nan)
#     for num,val in enumerate(unscaled_samples):
# #         print(val)
#         speclabelledset[val[0],val[1],val[2],val[3]] = speclabels[num]
#     return(complabelledset,speclabelledset)










