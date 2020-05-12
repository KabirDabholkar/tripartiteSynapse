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

from MFB_model import *
mdl, sim, r = get_MFB_model()


def runAZTrials(CaData, seed=False):
    vesData = []
    if seed:
        resAZ, vesRel = simAZ(CaData[0], CaData[1], sim, r, seed=seed)
    else:
        resAZ, vesRel = simAZ(CaData[0], CaData[1], sim, r)

    vesRelTot = np.sum(vesRel, axis=1)
    pks = detect_peaks(vesRelTot, edge='rising', show=False)
    vesData = CaData[0,pks]
    
    return vesData, resAZ

def run_and_average(resultPath="./", dir='nVDCC_4_dVDCC_60_nAZ_10_RRP_5_ISI_20_nAP_10', seeds=range(1, 201)):
    
    
    resample = 1000
    nAP = 2
    isi = 20  # ms
    tc = 0.02 # s 
    ts = 1.5  # ms

    
    fname = 'CaConc.dat' # the name of the calcium avg data file
    CaFile = os.path.join(resultPath, dir, fname)
    CaData = np.genfromtxt(CaFile, unpack=True, usecols=(0,1)) # in uM

    ### Run Simulations
    p = Pool()
    info = [(CaData, s) for s in seeds]
    vesRelTime, AZstates = np.array(p.starmap(runAZTrials, info)).T
    vesRelTime = [a.tolist() for a in vesRelTime]
    p.close()


    AZstates = np.mean(AZstates, axis=0)
    AZdata = np.hstack((np.array([CaData[0]]).T, AZstates))

    fAZ = os.path.join(resultPath, dir, f'az.dat')
    prevSeeds = 0

    if os.path.exists(os.path.join(resultPath, dir, 'vesData.dat')):
        #print('vesData.dat exists. Appending!')

        with open(os.path.join(resultPath, dir, 'vesData.dat'), "rb") as infile:
            vesData = pk.load(infile)

        vesRelTime = np.concatenate((vesData, vesRelTime)).tolist()

        with open(fAZ, 'r') as f:
            prevSeeds = int(f.readline().split('=')[1])

        prevAZdata = np.genfromtxt(fAZ)

        seedFrac = len(seeds)/(len(seeds) + prevSeeds)
        AZdata = seedFrac*AZdata + (1-seedFrac)*prevAZdata

    with open(os.path.join(resultPath, dir, 'vesData.dat'),"wb") as outfile:
        pk.dump(vesRelTime, outfile)

    header = f'seeds={len(seeds)+prevSeeds}\nt (s)\tAZ00\tAZ10\tAZ20\tAZ30\tAZ40\tAZ50\tAZ01\tAZ11\tAZ21\tAZ31\tAZ41\tAZ51\tAZ02\tAZ12\tAZ22\tAZ32\tAZ42\tAZ52'
    np.savetxt(fAZ, AZdata, fmt=['%0.5f']+['%0.4f']*18, header=header, delimiter="\t")

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