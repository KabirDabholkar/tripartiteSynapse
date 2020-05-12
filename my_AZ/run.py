import os, sys
import numpy as np
import pickle as pk
from random import randint
from itertools import product, chain
import scipy.interpolate as itp
from multiprocessing import Pool, Process

from analysis import *
from peaks import *


"""
resultPath = "./"
dir = 'nVDCC_4_dVDCC_60_nAZ_10_RRP_5_ISI_20_nAP_10
seeds = range(1, 201) # Make sure seeds are not repeated
resample = 1000
nAP = 2
isi = 20  # ms
tc = 0.02 # s 
ts = 1.5  # ms
#"""

########################################
######## Running MFB Simulations #######

from hits_MFB_model import *
mdl, sim, r, binding_reacs = get_MFB_model()


def runAZTrials(time,hit_times, seed=False):
    vesData = []
    if seed:
        resAZ, vesRel, bReac = simAZ(time, hit_times, sim,  r, binding_reacs, seed=seed)
    else:
        resAZ, vesRel, bReac = simAZ(time, hit_times, sim, r, binding_reacs)

    vesRelTot = np.sum(vesRel, axis=1)
    pks = detect_peaks(vesRelTot, edge='rising', show=False)
    vesData = time[pks]
    
    return vesData, resAZ, bReac

def run_and_average(time,hit_times,resultPath="./", dir='nVDCC_4_dVDCC_60_nAZ_10_RRP_5_ISI_20_nAP_10', isi=20, seeds=range(1, 201)):
    
    
    resample = 1000
    nAP = 2
    #isi = 20  # ms
    tc = 0.02 # s 
    ts = 1.5  # ms

    
    #fname = 'CaConc.dat' # the name of the calcium avg data file
    #CaFile = os.path.join(resultPath, dir, fname)
    #time = np.genfromtxt(CaFile, usecols=0)
    
    #fname = 'hit_times.dat' # the name of the calcium avg data file
    #hit_times_file = os.path.join(resultPath, dir, fname)
    #hit_times = np.genfromtxt(hit_times_file, usecols=0)
    
    print(hit_times.shape)
    
    ### Run Simulations
    p = Pool()
    info = [(time,hit_times,s) for s in seeds]
    vesRelTime, AZstates, bReac = np.array(p.starmap(runAZTrials, info)).T
    vesRelTime = [a.tolist() for a in vesRelTime]
    p.close()


    AZstates = np.mean(AZstates, axis=0)
    AZdata = np.hstack((np.array([time]).T, AZstates))
    bReac = np.mean(bReac, axis=0)
    bReacdata = np.hstack((np.array([time]).T, bReac))
    
    fAZ = os.path.join(resultPath, dir, f'az.dat')
    bReac_loc = os.path.join(resultPath, dir, f'bReac.dat')
    
    prevSeeds = 0

    if os.path.exists(os.path.join(resultPath, dir, 'vesData.dat')):
        #print('vesData.dat exists. Appending!')

        with open(os.path.join(resultPath, dir, 'vesData.dat'), "rb") as infile:
            vesData = pk.load(infile)

        #vesRelTime = np.concatenate((vesData, vesRelTime)).tolist()
        vesRelTime = vesData+vesRelTime
        
        with open(fAZ, 'r') as f:
            prevSeeds = int(f.readline().split('=')[1])

        prevAZdata = np.genfromtxt(fAZ)
        prevbReacdata = np.genfromtxt(bReac_loc)

        seedFrac = len(seeds)/(len(seeds) + prevSeeds)
        AZdata = seedFrac*AZdata + (1-seedFrac)*prevAZdata
        bReacdata = seedFrac*bReacdata + (1-seedFrac)*prevbReacdata
        
    with open(os.path.join(resultPath, dir, 'vesData.dat'),"wb") as outfile:
        pk.dump(vesRelTime, outfile)

    header = f'seeds={len(seeds)+prevSeeds}\nt (s)\tAZ00\tAZ10\tAZ20\tAZ30\tAZ40\tAZ50\tAZ01\tAZ11\tAZ21\tAZ31\tAZ41\tAZ51\tAZ02\tAZ12\tAZ22\tAZ32\tAZ42\tAZ52\tdAZ00\tdAZ10\tdAZ20\tdAZ30\tdAZ40\tdAZ50\tdAZ01\tdAZ11\tdAZ21\tdAZ31\tdAZ41\tdAZ51\tdAZ02\tdAZ12\tdAZ22\tdAZ32\tdAZ42\tdAZ52'
    np.savetxt(fAZ, AZdata, fmt=['%0.5f']+['%0.4f']*36, header=header, delimiter="\t")

    header=f'seeds={len(seeds)+prevSeeds}\nt (s)\tAZ0010\tAZ1020\tAZ2030\tAZ3040\tAZ4050\tAZ0111\tAZ1121\tAZ2131\tAZ3141\tAZ4151\tAZ0212\tAZ1222\tAZ2232\tAZ3242\tAZ4252\tAZ0001\tAZ0102\tAZ1011\tAZ1112\tAZ2021\tAZ2122\tAZ3031\tAZ3132\tAZ4041\tAZ4142\tAZ5051\tAZ5152\tdAZ0010\tdAZ1020\tdAZ2030\tdAZ3040\tdAZ4050\tdAZ0111\tdAZ1121\tdAZ2131\tdAZ3141\tdAZ4151\tdAZ0212\tdAZ1222\tdAZ2232\tdAZ3242\tdAZ4252\tdAZ0001\tdAZ0102\tdAZ1011\tdAZ1112\tdAZ2021\tdAZ2122\tdAZ3031\tdAZ3132\tdAZ4041\tdAZ4142\tdAZ5051\tdAZ5152'
    np.savetxt(bReac_loc, bReacdata, fmt=['%0.5f']+['%0.4f']*54, header=header, delimiter="\t")

    ########################################
    ############# get Pr Stat ##############

    #with open(os.path.join(resultPath, dir, 'vesData.dat'), "rb") as infile:
    #    vesRelTime = pk.load(infile)

    ts = [(i*isi+ts)/1000.0 for i in range(nAP)]

    p = Pool()
    seeds = len(vesRelTime)

    prs, rels, prratios = [], [], []
    #for _ in tqdm(range(resample), desc='Resampling'):
    for _ in range(resample):
        trials = [randint(0,seeds-1) for p in range(seeds)]

        timesInfo = list(product([vesRelTime[i] for i in trials[:]], [ts], [tc]))
        relData = np.array(p.starmap(getTimesData, timesInfo))
        relData = np.sum(relData, axis=0)

        nRel = relData[0]
        allRel = relData[1]
        prratio = nRel/nRel[0]

        prs.append(np.array(nRel))
        rels.append(np.array(allRel))
        prratios.append(np.array(prratio))

    prs = np.array(prs)/seeds
    rels = np.array(rels)/seeds

    mpr = np.mean(prs, axis=0)
    spr = np.std(prs, axis=0)

    mrel = np.mean(rels, axis=0)
    srel = np.std(rels, axis=0)

    mprratios = np.mean(prratios, axis=0)
    sprratios = np.std(prratios, axis=0)

    result = np.array([mpr, spr, mrel, srel, mprratios, sprratios]).T
    header = 'Pr\tePr\t\tRel\t\teRel\t\tFac\t\teFac'
    print(result)

    np.savetxt(os.path.join(resultPath, dir, 'result.dat'), 
               result, header=header, fmt=['%0.4f']*6, delimiter='\t')

    p.close()
    #"""