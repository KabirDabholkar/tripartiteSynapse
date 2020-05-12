import os, sys
import numpy as np
import pickle as pk
from random import randint
from itertools import product, chain
import scipy.interpolate as itp
from multiprocessing import Pool, Process

from analysis import *
from peaks import *


resultPath = "/home/kabir/Project/tripartiteSynapse/results/ca_traces"


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


#dir = 'nVDCC_9_dVDCC_100_nAZ_9_RRP_10_ISI_100_nAP_10'
#save_loc="../results/ca_traces"
out_vdcc_range=np.arange(60,201,10)
for v,sim_type in product(out_vdcc_range,sims):

    fname = 'ca_trace.dat' # the name of the calcium avg data file
    CaFile = os.path.join(resultPath, sim_type, "vdcc"+v, fname)
    CaData = np.genfromtxt(CaFile, unpack=True, usecols=(0,1)) # in uM


    seeds = 10
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
    #print(AZstates)

    fAZ = os.path.join(resultPath, sim_type, "vdcc"+v, f'az.dat')
    AZdata = np.hstack((np.array([CaData[0]]).T, AZstates))
    #print(AZdata, AZdata.shape)
    np.savetxt(fAZ, AZdata, fmt=['%0.5f']+['%0.4f']*18, delimiter="\t")

