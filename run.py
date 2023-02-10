#!/usr/bin/env python

import yaml,sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.special import erf

def desi_lbg_redshift_success(efftime_hours,rmag,template) :
    # 50% efficiency in 1.3 hours at rmag = 24

    if template == 0 :
        nsig=2. # controls the derivative
        efftime_50 = 1.25
        rmag_50 = 24.
    elif template == 1 :
        nsig=1.6 # controls the derivative
        efftime_50 = 1.1
        rmag_50 = 24.
    elif template == 2 :
        nsig=1.4 # controls the derivative
        efftime_50 = 0.8
        rmag_50 = 24.
    elif template == 3 :
        nsig=1.15 # controls the derivative
        efftime_50 = 0.3
        rmag_50 = 24.
    else :
        print("Error unknown template",template)
        sys.exit(12)

    snr = np.sqrt(efftime_hours/efftime_50 )*10**(0.4*(rmag_50-rmag))
    return 0.5*(1 + erf(nsig*(snr-1)/np.sqrt(2)))


def get_target_density(rmag_min=0,rmag_max=24.) :
    # for u_ext selection (fit of distribution from C. Yeche)
    return 600.*(np.exp(2.5*(rmag_max-24))-np.exp(2.5*(rmag_min-24)))

def read_config_params(param_filename) :
    with open(param_filename, 'r') as file :
        params=yaml.safe_load(file)
    return params

def fiber_assignement_prob(target_density_deg2,params) :
    ntargets              = int(params["FOCAL_PLANE_AREA_DEG2"]*target_density_deg2)
    nfibers               = params["NUMBER_OF_FIBERS"]-params["NUMBER_OF_SKY_FIBERS"]-params["NUMBER_OF_STAR_FIBERS"]
    fiber_to_focal_plane  = params["FIBER_PATROL_AREA_DEG2"]/params["FOCAL_PLANE_AREA_DEG2"]
    nassigned=0
    for target in range(ntargets) :
        prob = max(0,min(1,(nfibers-nassigned)*fiber_to_focal_plane))
        nassigned += prob
        if nassigned >= nfibers : break
    nfree = max(0,nfibers-nassigned)
    return nassigned/ntargets , nfree

def fiber_assignement_prob_rnd(target_density_deg2,params) :
    ntargets              = int(params["FOCAL_PLANE_AREA_DEG2"]*target_density_deg2)
    nfibers               = params["NUMBER_OF_FIBERS"]-params["NUMBER_OF_SKY_FIBERS"]-params["NUMBER_OF_STAR_FIBERS"]
    fiber_to_focal_plane  = params["FIBER_PATROL_AREA_DEG2"]/params["FOCAL_PLANE_AREA_DEG2"]
    nassigned=0
    probs=[]
    for target in range(ntargets) :
        prob = max(0,min(1,(nfibers-nassigned)*fiber_to_focal_plane))
        #print(prob)
        probs.append(prob)
        nassigned += int(prob>np.random.uniform())
        if nassigned >= nfibers : break
    nfree = max(0,nfibers-nassigned)
    return probs , nfree


def main() :

    parser = argparse.ArgumentParser(
        description="Survey emulator",
        epilog="""Simple emulator to estimate a survey redshift distribution.
                  example: ./run.py -c desi-refurbished-config.yaml""")
    parser.add_argument('-c','--config', type=str, required=True,
                        help = 'path to survey config yaml file')
    parser.add_argument('-d','--debug', action = 'store_true',
                        help = 'verbose and debugging tests')

    args = parser.parse_args()

    params=read_config_params(args.config)

    print("Inputs:")
    for p in params :
        print(f"{p} = {params[p]}")

    rmag_max=params["RMAG_MAX"]

    computed_total_target_density = get_target_density(rmag_min=0,rmag_max=rmag_max)
    print("With rmag<{} the computed target density is {:.1f} per deg2".format(rmag_max,computed_total_target_density))
    if "TOTAL_TARGET_DENSITY" in params.keys() :
        actual_total_target_density = params["TOTAL_TARGET_DENSITY"]
        print("Will scale total target density from {:.1f} to {:.1f} per deg2".format(computed_total_target_density,actual_total_target_density))
    else :
        actual_total_target_density = computed_total_target_density

    print("")
    print("Outputs:")
    fprob , nfree = fiber_assignement_prob(actual_total_target_density,params)
    print("fiber assignment probability = {:.3f}".format(fprob))

    if args.debug :
        x=[]
        for r in range(40) :
            x.append(fiber_assignement_prob_rnd(args.target_density,params))
        plt.figure("fprob")
        plt.hist(np.hstack(x),bins=50)
        plt.axvline(fprob,color="k")
        plt.xlabel("Fiber assignement probability")
        print("fiber assignment probability (from rdm) = {:.3f}".format(np.mean(x)))

    ntiles=params["SURVEY_DURATION_YEARS"]*params["EFFECTIVE_TIME_PER_YEAR_HOURS"]*3600./params["EFFECTIVE_TIME_PER_TILE_SECONDS"]
    print("number of tiles = {:.1f}".format(ntiles))
    nobs=(ntiles*params["FOCAL_PLANE_AREA_DEG2"])/params["SURVEY_AREA_DEG2"]
    print("number of observations of each pos in the survey = {:.1f}".format(nobs))

    print("number of free fibers per exposure  = {:.1f}".format(nfree))
    print("density of free fibers per exposure = {:.1f} per deg2".format(nfree/params["FOCAL_PLANE_AREA_DEG2"]))
    print("density of free fibers on sky       = {:.1f} per deg2".format(nfree/params["FOCAL_PLANE_AREA_DEG2"]*nobs))



    sum_frac = 0
    lbg_template_frac={}
    for t in range(4) :
        k=f"TARGET_FRAC_LBG_TEMPLATE_{t}"
        if k not in params.keys() :
            print(f"warning, missing keyword {k} in config, assuming 0 of those")
            continue
        lbg_template_frac[t]=float(params[k])
        #print(k,lbg_template_frac[t])
        sum_frac += lbg_template_frac[t]
    k="TARGET_FRAC_NON_LBG"
    if k not in params.keys() :
        print(f"error, missing keyword {k} in config (need it for cross-check)")
        sys.exit(12)
    sum_frac += float(params[k])
    if np.abs(sum_frac-1)>0.01 :
        print("Error, not all targets are described. The sum of TARGET_FRAC_LBG_TEMPLATE_X + TARGET_FRAC_NON_LBG should be 1, and I have {:.3f}".format(sum_frac))
        sys.exit(12)

    efftime_per_target = nobs*fprob*params["EFFECTIVE_TIME_PER_TILE_SECONDS"]
    print("mean effective time per target = {:.2f} hours".format(efftime_per_target/3600.))

    # random distribution of efftime per target
    nt=100000
    nobs_per_target=np.zeros(nt)
    for o in range(int(nobs)) :
        nobs_per_target += (np.random.uniform(size=nt)<fprob)
    efftime_per_target = nobs_per_target*params["EFFECTIVE_TIME_PER_TILE_SECONDS"]
    #print("mean effective time per target = {:.2f} hours".format(np.mean(efftime_per_target)/3600.))

    #print("Compute redshift efficiency for known magnitude distribution")

    mag_bin_edges   = np.linspace(22,rmag_max,int((rmag_max-22+1e-12)/0.1)+1)
    mag_bin_centers = (mag_bin_edges[1:]+mag_bin_edges[:-1])/2.
    mag_bin_size    = (mag_bin_edges[1:]-mag_bin_edges[:-1])

    target_density=[]
    redshift_density=[]


    target_density_per_template={}
    redshift_density_per_template={}
    for template in lbg_template_frac :
        target_density_per_template[template]=[]
        redshift_density_per_template[template]=[]


    for rmag,dmag in zip(mag_bin_centers,mag_bin_size) :

        density = get_target_density(rmag-dmag/2,rmag+dmag/2) * actual_total_target_density / computed_total_target_density
        density_with_redshifts = 0.

        for template in lbg_template_frac :
            prob = np.mean(desi_lbg_redshift_success(efftime_per_target/3600.,rmag,template=template))
            density_with_redshifts += lbg_template_frac[template] * density * prob

            target_density_per_template[template].append(lbg_template_frac[template] * density)
            redshift_density_per_template[template].append(lbg_template_frac[template] * density * prob)

        target_density.append(density)
        redshift_density.append(density_with_redshifts)

    target_density=np.array(target_density)
    redshift_density=np.array(redshift_density)


    plt.figure("redshift success")
    ax=plt.subplot(121)
    plt.plot(mag_bin_centers,target_density,"--",c="k")
    plt.plot(mag_bin_centers,redshift_density,"-",c="k")
    for template in lbg_template_frac :
        target_density_per_template[template]=np.array(target_density_per_template[template])
        redshift_density_per_template[template]=np.array(redshift_density_per_template[template])

        plt.plot(mag_bin_centers,target_density_per_template[template],"--",color=f"C{template}")
        plt.plot(mag_bin_centers,redshift_density_per_template[template],"-",color=f"C{template}")

    total_target_density = np.sum(target_density)
    total_redshift_density = np.sum(redshift_density)
    plt.text(0.05,0.95,'Target density = {:.1f}\nRedshift density = {:.1f} per deg2'.format(total_target_density,total_redshift_density),transform=ax.transAxes,verticalalignment="top")

    print("Target density   = {:.1f} per deg2".format(total_target_density))
    print("Redshift density = {:.1f} per deg2".format(total_redshift_density))
    print("Redshift success rate for all targets      = {:.1f}%".format(100*total_redshift_density/total_target_density))
    print("Redshift success rate for true LBG targets = {:.1f}%".format(100*total_redshift_density/(total_target_density*(1-params["TARGET_FRAC_NON_LBG"]))))


    plt.xlabel("rmag")
    plt.ylabel("Number density per 0.1 mag per deg2")
    #plt.legend(loc="lower left")
    plt.grid()

    plt.subplot(122)
    plt.plot(mag_bin_centers,redshift_density/target_density,"-",c="k",label="All LBGs")
    for template in lbg_template_frac :
        plt.plot(mag_bin_centers,redshift_density_per_template[template]/target_density_per_template[template],"-",color=f"C{template}",label=f"LBG template {template}")

    plt.xlabel("rmag")
    plt.ylabel("Redshift success")
    plt.legend(loc="lower left")
    plt.grid()
    plt.show()






    plt.figure("efftime")
    plt.hist(efftime_per_target/3600.,bins=50)
    plt.xlabel("Total effective time (hours)")
    plt.grid()


    if args.debug:
        plt.show()

main()
