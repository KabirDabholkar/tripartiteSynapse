import os, sys
import numpy as np
import pickle as pk
from random import randint
from itertools import product, chain
import scipy.interpolate as itp
from multiprocessing import Pool, Process

from analysis import *
from peaks import *


resultPath = "/home/kabir/Project/tripartiteSynapse/results/STEPS_AZ"


########################################
######## Running MFB Simulations #######

from MFB_model import *
mdl, sim, r = get_MFB_model()


def runAZTrials(CaData):
    vesData = []

    resAZ, vesRel = simAZ(CaData[0], CaData[1], sim, r)
    vesRelTot = np.sum(vesRel, axis=1)
    pks = detect_peaks(vesRelTot, edge='rising', show=False)
    vesData = CaData[0,pks]
    
    #return np.array([vesData, resAZ])
    return [vesData, resAZ]


#dir = 'nVDCC_9_dVDCC_100_nAZ_9_RRP_10_ISI_20_nAP_10'
sims=["1apr2020/R150control","1apr2020/R150ER2x","1apr2020/R300ER2x","1apr2020/R150ER3x","1apr2020/R300ER3x"]
out_vdcc_range=np.arange(60,201,20)
for v,sim_type in product(out_vdcc_range,sims):

    fname = 'ca_trace.dat' # the name of the calcium avg data file
    CaFile = os.path.join(resultPath, sim_type, "vdcc%d" %v, fname)
    CaData = np.genfromtxt(CaFile, unpack=True, usecols=(0,1)) # in uM

    seeds = 1000
    p = Pool()
    vesRelTime = np.array(p.map(runAZTrials, [CaData for _ in range(seeds)]))
    p.close()

    AZstates = vesRelTime[:,1]
    vesRelTime = vesRelTime[:,0]
    """
    vesRelTime is a list of lists
    each list is for a single trial/seed 
    and contains the times at which a vesicle was released
    """

    AZstates = np.mean(AZstates, axis=0)


    fAZ = os.path.join(resultPath, sim_type, "vdcc%d" %v, f'az.dat')
    AZdata = np.hstack((np.array([CaData[0]]).T, AZstates))

    #np.savetxt(fAZ, AZdata, fmt=['%0.5f']+['%0.4f']*18, delimiter="\t")

    #with open(os.path.join(resultPath, dir, 'vesData.dat'),"wb") as outfile:
    #    pk.dump(vesRelTime, outfile)


    ########################################
    ############# get Pr Stat ##############

    resample = 1000
    nAP = 2
    isi = 20  # ms
    tc = 0.02 # s 
    ts = 1.5  # ms

    ts = [(i*isi+ts)/1000.0 for i in range(nAP)]

    p = Pool()

    seeds = len(vesRelTime)

    prs, rels = [], []
    for _ in tqdm(range(resample), desc='Resampling'):
        trials = [randint(0,seeds-1) for p in range(seeds)]

        timesInfo = list(product([vesRelTime[i] for i in trials[:]], [ts], [tc]))
        relData = np.array(p.starmap(getTimesData, timesInfo))
        relData = np.sum(relData, axis=0)

        nRel = relData[0]
        allRel = relData[1]

        prs.append(np.array(nRel))
        rels.append(np.array(allRel))

    prs = np.array(prs)/seeds
    rels = np.array(rels)/seeds

    mpr = np.mean(prs, axis=0)
    spr = np.std(prs, axis=0)

    mrel = np.mean(rels, axis=0)
    srel = np.std(rels, axis=0)

    result = np.array([mpr, spr, mrel, srel]).T
    header = 'Pr\tePr\t\tRel\t\teRel'
    print(result)

    np.savetxt(os.path.join(resultPath, sim_type, "vdcc%d" %v, 'result.dat'), \
               result, header=header, fmt=['%0.4f']*4, delimiter='\t')

    p.close()