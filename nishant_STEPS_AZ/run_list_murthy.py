from run import *
import os
from shutil import copyfile
from analysis import *
import pickle as pkl


def seed_num(vdcc_num):
    #vdcc_num=int(fname.split('V')[1].replace('.mdl',''))
    if vdcc_num>=130 and vdcc_num<=200:
        return 500
    elif vdcc_num>=90 and vdcc_num<=120:
        return 1000
    elif vdcc_num>=70 and vdcc_num<=80:
        return 3000
    elif vdcc_num>=40 and vdcc_num<=60:
        return 5000
    else:
        return 10


def old_seed_num(vdcc_num):
    #vdcc_num=int(fname.split('V')[1].replace('.mdl',''))
    if vdcc_num>=130 and vdcc_num<=200:
        return 100
    elif vdcc_num>=90 and vdcc_num<=100:
        return 1
    elif vdcc_num>=70 and vdcc_num<=80:
        return 1
    elif vdcc_num>=40 and vdcc_num<=60:
        return 1
    else:
        return 10
    

sims=["stores_blocked","R150control","R150ER2x","R300ER2x","R150ER3x","R300ER3x"]#, "ryr_old",
oldresultPath="/home/kabir/Project/tripartiteSynapse/results/STEPS_AZ/16apr2020"
resultPath=   "/home/kabir/Project/tripartiteSynapse/results/STEPS_AZ/16apr2020"
print("resultPath: ",resultPath)
isis=list(range(60,181,40))
VDCC=[66, 76, 82, 87, 92, 96, 100, 104, 108, 112, 116, 120, 124, 128, 132, 137, 143, 149, 158, 198]

vdcc_dist={66: 0.053833605220228384, 76: 0.16476345840130505, 82: 0.19575856443719414, 87: 0.1631321370309951,\
           92:0.1370309951060359, 96: 0.08319738988580751, 100: 0.06035889070146819, 104: 0.03915171288743882, 108: \
           0.02936378466557912, 112:0.024469820554649267, 116: 0.011419249592169658, 120: 0.01468189233278956, 124: \
           0.0032626427406199023, 128:0.0032626427406199023, 132: 0.0065252854812398045, 137: 0.0016313213703099511, 143: \
           0.004893964110929853, 149: 0.0, 158:0.0032626427406199023, 198: 0.0}

#print('isi = '+str(isi),'\nvdcc = '+str(VDCC))
"""
folders=[]

for (i,sim) in product(isi,sims):
    folders.append(sim+"/isi%d" %i)
    #seeds.append((1,int(5000*vdcc_dist(v)) #seed_num(v)
print(folders)
#"""

for isi,sim in product(isis,sims):
    folder=os.path.join(sim,"isi%d" %isi)
    
    resample = 1000
    nAP = 2
    #isi = isi #20  # ms
    tc = 0.02 # s 
    ts = 1.5  # ms

    ts = [(i*isi+ts)/1000.0 for i in range(nAP)]
    
    #PrStat(folder, resultPath, resample=1000, tsi=ts, tc=tc)
    print(os.path.join(resultPath,folder,"vesData.dat"))
    with open(os.path.join(resultPath,folder,"vesData.dat"),'rb') as f:
        vesRelTime = pkl.load(f)
    
    
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

    np.savetxt(os.path.join(resultPath, folder, 'result.dat'), 
               result, header=header, fmt=['%0.4f']*6, delimiter='\t')

    p.close()
    
    #"""
    #if not os.path.exists(os.path.join(resultPath,folder)):
    #    os.makedirs(os.path.join(resultPath,folder))
    """
    for v in VDCC:
        copyfile(os.path.join(oldresultPath,folder,"vdcc%d" %v,"CaConc.dat"),os.path.join(resultPath,folder,"CaConc.dat"))
        seed_range_upper=int(10000*vdcc_dist[v])+1
        print(v,seed_range_upper)
        if seed_range_upper>1:
            run_and_average(resultPath=resultPath,dir=folder,seeds=range(1,seed_range_upper))
    #"""