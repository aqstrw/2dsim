# Ambar Qadeer
# ambar.qadeer@gmail.com, qadeer@strw.leidenuniv.nl
# Michael Migo Mueller
# migo@strw.leidenuniv.nl, m.mueller@astro.rug.nl


import os
from utilfuncs import *
# from combinepipe import *
from myglobals import *

def workflow(reuse = 1):

    # Takes in locations to the mirisim simulated scene and produces a moving target based on mvx and saves based in the provided_directory+xmoverate
    dm,moved,jumpdet = runpip(reuse,ippath) # return(dm_slopes,msc_op_file,opramp_path,moverate)
    # at this point the asteroid has been simulated and the outputs are available in the opdirectory : msc_op_dir = msc_output_prefix+'x'+str(moverate)+'/'
    #####################################################################

    hduramp = fits.open(jumpdet)
    ramp = hduramp[1].data # ramps
    gdq = hduramp[3].data # dat quality flags for each group
    hduramp.close()

    # Takes in ramps and gdq, masks them and returns crhits (val = 2) flattened and 3d
    mcrhitramps,mramp,mgdq,mcrhitscum,xdata,ydata  = maskandextract(ramp, gdq)

    # build samples for DBSCAN
    samps,sampss = build_arrays_for_DBSCAN(mcrhitramps, eps = ep)

    model, clusters, no_clusters, no_noise = apply_DBSCAN(mcrhitramps, sampss, eps = ep , core = core)

    asteroid_clusters, asteroid_labels = identify_asteroid_clusters(model,clusters,samps)

    # plotting ##########################################################

    f = fig_update_traceprops(
        highlight_core(
            sampss,
            model,
            asteroid_clusters,
            plot_3dclusters(
                sampss,
                model,
                asteroid_labels, # controls what is displayed
                box_shape = mgdq.shape, # used to get plot ranges and descaling zaxis labels
                zscale = 3,
                opacity = 0.05,
                height = 800,
                width = 800,
                title = 'Asteroid vs. Noise',
            ),
            zscale = 3,
        )[0]
    )
    f.update_layout(showlegend=True)
    f.show()

    iso_labs, iso_samps = prep_isolated_asteroid_fig_params(sampss, asteroid_labels)
    i = fig_update_traceprops(
        highlight_core(
            sampss,
            model,
            plot_3dclusters(
                iso_samps,
                model,
                iso_labs, # controls what is displayed
                isolate = 1,
                box_shape = mgdq.shape,
                zscale = 3,
                opacity = 0.05,
                height = 800,
                width = 800,
                title = 'Asteroid',
            ),
            zscale = 3,
        )[0]
    )
    i.show()

    # cls, sls, clusters = DBSCAN_and_checks(mcrhitramps, samps, sampss, eps = ep , core = core,
    #                                  span = span, vel = vel)

    # if (cls, sls, clusters) == (0,0,0):
    #     print('premature exit')

if __name__ == "__main__":
    print("running workflow with :\n\
        wasteroid = {}\n\
        t_1 = {}\n\
        frames = {}\n\
        pixelscale = {}\n\
        base = {}\n\
        verbosity = {}\n\
        save = {}\n\
        plot = {}\n\
        moverate = {}\n\
        reuse = {}\n\
        ep = {}\n\
        core = {}\n\
        span = {}\n".format(wasteroid,t_1,frames,pixelscale,base,verbosity, save, plot, moverate,reuse, ep,core,span))
    workflow()