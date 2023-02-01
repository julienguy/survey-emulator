#!/usr/bin/env python

import yaml
import numpy as np
import matplotlib.pyplot as plt
import argparse

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
    return nassigned/ntargets

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
    return probs


def main() :

    parser = argparse.ArgumentParser(
        description="Survey emulator",
        epilog="""Simple emulator to estimate a survey redshift distribution.
                  example: ./run.py -c desi-refurbished-config.yaml""")
    parser.add_argument('-c','--config', type=str, required=True,
                        help = 'path to survey config yaml file')
    parser.add_argument('-t','--target-density', type=float, required=False,
                        default=1400,
                        help = 'target density per deg2')
    parser.add_argument('-d','--debug', action = 'store_true',
                        help = 'verbose and debugging tests')

    args = parser.parse_args()

    params=read_config_params(args.config)
    print("Inputs:")
    for p in params :
        print(f"{p} = {params[p]}")
    print("TARGET DENSITY =",args.target_density)
    print("")
    print("Outputs:")
    fprob=fiber_assignement_prob(args.target_density,params)
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
    efftime_per_target = nobs*fprob*params["EFFECTIVE_TIME_PER_TILE_SECONDS"]
    print("mean effective time per target = {:.2f} hours".format(efftime_per_target/3600.))

    # distribution
    nt=10000
    nobs_per_target=np.zeros(nt)
    for o in range(int(nobs)) :
        nobs_per_target += (np.random.uniform(size=nt)<fprob)
    efftime_per_target = nobs_per_target*params["EFFECTIVE_TIME_PER_TILE_SECONDS"]

    plt.figure("efftime")
    plt.hist(efftime_per_target/3600.)
    plt.xlabel("Total effective time (hours)")
    plt.grid()


    if args.debug:
        plt.show()

main()
