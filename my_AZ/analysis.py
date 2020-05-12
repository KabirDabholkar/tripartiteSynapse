import os, sys, re
import numpy as np
from peaks import *
import pickle as pk
from time import time, sleep
from random import randint
from itertools import product
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from multiprocessing import Pool, Process
from tqdm import tqdm, trange,  tnrange
from tqdm import tqdm_notebook as tqdmn

sys.path.append('/home/kabir/tripartiteSynapse/nishant_STEPS_AZ')


##################################################
### Get list of directories in the path starting with 's_'
def getDirs(path, sstr='s_'):
    dirs = [d for d in os.listdir(path) if os.path.isdir(path + '/' + d) and sstr in d]
    dirs.sort()
    return(dirs)

##################################################
### Get data from dataFile in array
def _getData(file):
    return np.genfromtxt(file, unpack=True)

def getData(dataPath, dir, fname, cores=None):
    path = os.path.join(dataPath,dir)
    allDirs = getDirs(path)
    files = [os.path.join(path,d,'dat',fname) for d in allDirs]
    
    if cores:
        p = Pool(cores)
        #data = p.map(np.genfromtxt, files, )
        data = p.map(_getData, files)
        p.close()
        #p.join()

    else:
        data = []
        for d in tqdmn(allDirs, desc=dir):
            file = os.path.join(path,d,'dat',fname)
            temp = np.genfromtxt(file, unpack=True)
            data.append(temp)
    
    return np.array(data)

##################################################
### Average Over all Seeds
def avg_dat(dataPath, resultPath, dir, fname, header='', std=False, 
			fmt=['%.7f'], ret=False, write=True, cores=None):
    data = getData(dataPath, dir, fname, cores=cores)

    seeds = len(data)
    avg = np.mean(data, axis=0)
    ncols = len(avg)

    if not os.path.exists(resultPath):
        #print('made directory!')
        os.makedirs(resultPath)

    if not os.path.exists(os.path.join(resultPath,dir)):
        os.makedirs(os.path.join(resultPath,dir))

    if std:
        std = np.std(data, axis=0)

        temp = np.zeros((avg.shape[0]*2-1, avg.shape[1]))
        temp[0,:] = avg[0]
        temp[1::2,:] = avg[1:]
        temp[2::2,:] = std[1:]

        if write:
            np.savetxt(os.path.join(resultPath,dir,fname), temp.T,
                       fmt=['%.5f']+fmt*(ncols-1)*2, header=header, delimiter='\t')

    else:
        temp = avg
        if write:
            np.savetxt(os.path.join(resultPath,dir,fname), temp.T,
                       fmt=['%.5f']+fmt*(ncols-1), header=header, delimiter='\t')
    #print("Writing averaged data to:\t" + fname)

    if ret:
        return temp
    
##################################################
### Ca Concentration Calculation
def _getConc(data, step):
    #sdata = savgol_filter(data[1], 51, 3)
    sdata = data[1]
    
    #plt.plot(data[0],sdata)
    #plt.ylim(0,0.2)
    
    c_tc = np.multiply(data[0],sdata)
    dt = step*(data[0][1]-data[0][0])

    c_out = []
    for i in range(0,len(data[0])-step-1,step):
        c_out.append((c_tc[i+step]-c_tc[i])/dt)

    return np.array(c_out)

def CaConc(resultPath, dir, fname='ca.dat', step=5, skip=[1]):
    try:
        data = np.genfromtxt(os.path.join(resultPath,dir,fname), unpack=True)

        CaData = [data[0,range(0,len(data[0])-step-1,step)]]
        #print(CaData)
        for i in range(2,data.shape[0]):
            if i in skip:
                continue
            CaData.append(_getConc([data[0],data[i]], step=step))

        CaData = np.array(CaData)

        np.savetxt(os.path.join(resultPath,dir,'CaConc.dat'), CaData.T, delimiter='\t',
                   fmt=['%.5f']+['%.5f']*(CaData.shape[0]-1))
        print(f'{dir} : Done!')
    except:
        print(f"Error in {dir}")
        
        
##################################################
### Get Pr 
## isi in ms
def PrStat(dir, resultPath, resample=1000, tsi=1.5e-3, tc=0.02):
    nVDCC, dVDCC, nAZ = [int(getSimInfo(dir, key)) for key in ['nVDCC', 'dVDCC', 'nAZ']]
    #print(nVDCC, dVDCC, nAZ)
    
    with open(os.path.join(resultPath, dir, 'vesData.dat'), "rb") as infile:
        vesData = pk.load(infile)
    
    print(f'Vesicle release stats:')
    result = []
    for rrp, timedata in tqdmn(vesData.items(), desc=dir):
        seeds = len(timedata)
        #print(seeds)

        prs = []
        for _ in range(resample):
            x = [randint(0,seeds-1) for p in range(seeds)]
            #print(x)

            nRel, nAllRel = 0, 0
            for times in [timedata[i] for i in x]:
                #   Get no. of vesicles released after AP specified by ts
                '''
                ifRel, allRel = 0,0
                for time in times:
                    if (time>ts and time<ts+tc):
                        ifRel = 1
                        allRel += 1
                '''        
                
                ifRel = int(len(list(filter(lambda a: a>tsi and a<tsi+tc, times)))>0)
                allRel = sum(map(lambda a: a>tsi and a<tsi+tc, times))
                
                if ifRel == 1: nRel += 1
                nAllRel += allRel

            #print(nRel, nAllRel)
            prs.append(np.array([nRel, nAllRel])/seeds)

        m = np.mean(prs, axis=0)
        s = np.std(prs, axis=0)

        result.append(np.concatenate((np.array(list(zip(m,s))).flatten(),[nAZ, nVDCC, dVDCC, int(rrp)]), axis=0).tolist())
        
    for res in result:
        print([f'{i:.4f}' for i in res[:4]] + [f'{i:.2f}' for i in res[4:]])
        
    header = 'p1\tep1\t\ttRel\tetRel\tnAZ\tnV\tdV\tRRP'
    np.savetxt(os.path.join(resultPath, dir, 'result.dat'), result, header=header, 
              fmt=['%0.4f']*4+['%d']*4, delimiter="\t")

    return(result)

### Get release stats
def getTimesData(times, ts, tc):
    nRelTrial = [int(len(list(filter(lambda a: a>tsi and a<tsi+tc, times)))>0) for tsi in ts]
    allRelTrial = [sum(map(lambda a: a>tsi and a<tsi+tc, times)) for tsi in ts]

    return np.array([nRelTrial, allRelTrial])

def relStat(dir, resultPath, resample=1000, ts=1.5e-3, tc=0.02):
    with open(os.path.join(resultPath, dir, 'vesData.dat'), "rb") as infile:
        vesData = pk.load(infile)
    
    try:
        n = int(getSimInfo(dir, 'nAP')) # number of AP
    except:
        n = 1

    try:
        isi = int(getSimInfo(dir, 'ISI')) # isi in ms
    except:
        isi = 0
    
    ts = [(i*isi+1.5)/1000.0 for i in range(n)]
    
    result = []
    p = Pool()
    desc = str([int(a) for a in re.findall(r'\d+', dir)[:-1]])
    for rrp, timedata in vesData.items():
        seeds = len(timedata)

        prs, rels = [], []
        for _ in tqdmn(range(resample), desc=f'Resampling {desc}'):
            trials = [randint(0,seeds-1) for p in range(seeds)]
            '''
            nRel = np.zeros(n)
            allRel = np.zeros(n)
            for times in [timedata[i] for i in trials[:]]:

                # Get no. of vesicles released after APs specified by ts
                nRelTrial = [int(len(list(filter(lambda a: a>tsi and a<tsi+tc, times)))>0) for tsi in ts]
                allRelTrial = [sum(map(lambda a: a>tsi and a<tsi+tc, times)) for tsi in ts]
                
                #print(nRelTrial, allRelTrial)
                nRel = [sum(x) for x in zip(nRel, nRelTrial)]
                allRel = [sum(x) for x in zip(allRel, allRelTrial)]
            '''
            timesInfo = list(product([timedata[i] for i in trials[:]], [ts], [tc]))
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
        print(f'RRP: {rrp}\n', result)
        header = 'Pr\tePr\t\tRel\t\teRel'
        np.savetxt(os.path.join(resultPath, dir, f'result_RRP_{rrp}.dat'), 
                   result, header=header, fmt=['%0.4f']*4, delimiter='\t')
    
    p.close()
    return(result)