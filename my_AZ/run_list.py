from run import *
import data_collector as dc
from itertools import product
from analysis import *
import pickle as pkl
import os

def seed_num(vdcc_num):
    #vdcc_num=int(fname.split('V')[1].replace('.mdl',''))
    if vdcc_num>=130 and vdcc_num<=200:
        return 500
    elif vdcc_num>=90 and vdcc_num<=120:
        return 500
    elif vdcc_num>=70 and vdcc_num<=80:
        return 1000
    elif vdcc_num>=40 and vdcc_num<=60:
        return 3000
    else:
        return 10


def old_seed_num(vdcc_num):
    #vdcc_num=int(fname.split('V')[1].replace('.mdl',''))
    if vdcc_num>=130 and vdcc_num<=200:
        return 100
    elif vdcc_num>=90 and vdcc_num<=100:
        return 100
    elif vdcc_num>=70 and vdcc_num<=80:
        return 100
    elif vdcc_num>=40 and vdcc_num<=60:
        return 100
    else:
        return 10
    


    
sims=["R150control","R150ER2x","R300ER2x","R150ER3x","R300ER3x","stores_blocked"] #"ryr_old",
resultPath="/home/kabir/Project/tripartiteSynapse/results/STEPS_AZ/hits_6May2020"#+'/'+sim_type
load_cal_path="/media/kabir/ca_binding/"
isis=range(20,181,40)
VDCCs=range(60,201,20)
fnames=["RSI%d"%i+"V%d"%v for i,v in product(isis,VDCCs)]

D=dc.data_collector(load_cal_path,sims,fnames)

print("resultPath: ",resultPath,"Calcium traces path: ")


print('isi = '+str(isis),'\nvdcc = '+str(VDCCs))

"""
folders=[]
seeds=[]
for (i,v,sim) in product(isi,VDCC,sims):
    folders.append(sim+"/vdcc"+str(v))
    seeds.append((1,seed_num(v)+1)) #seed_num(v)
print(folders,seeds)
#"""

for v,isi,sim in product(VDCCs[:1],isis[:1],sims[:1]):
    
    folder=os.path.join(sim,"isi%d" %isi,"vdcc%d"%v)
    
    #loading calcium traces
    fname="RSI%d"%isi+"V%d"%v
    
    #copyfile(os.path.join(oldresultPath,folder,"vdcc%d" %v,"CaConc.dat"),os.path.join(resultPath,folder,"CaConc.dat"))
    loc=os.path.join(load_cal_path,sim,fname)#,,"dat","ca_hits.dat")
    hit_times=[]
    for seed in os.listdir(loc):
        file_loc=os.path.join(loc,seed,'dat','ca_hits.dat')
        hit_times.append(np.sort(np.loadtxt(file_loc,delimiter=' ',usecols=0)))
    ca_loc=file_loc=os.path.join(loc,seed,'dat','ca.dat')
    time=np.loadtxt(ca_loc,delimiter=' ',usecols=0)
    #print(time.shape,conc.shape,int(time.shape[0]/conc.shape[0]))
    #conc=np.repeat(conc,int(time.shape[0]/conc.shape[0]),axis=0)
    
    
    #resampling from conc
    r_sample=np.random.randint(0,len(hit_times),seed_num(v))
    count=np.bincount(r_sample)
    
    path=os.path.join(resultPath,folder)
    if not os.path.exists(path):
            os.makedirs(path)
    for val,freq in zip(r_sample,count):
        if freq>0:
            #print(time[:conc.shape[0]].shape,conc[:,val].shape)
            #out=np.stack((time[:conc.shape[0]],conc[:,val])).T    
            #print(folder,out.shape)
            #np.savetxt(os.path.join(path,"CaConc.dat"),out,fmt=['%0.5f']*2, delimiter='\t')
                run_and_average(time,hit_times[val],resultPath=resultPath,dir=folder,seeds=range(1,freq+1))
    #"""
    """
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
    """
    for v in VDCC:
        copyfile(os.path.join(oldresultPath,folder,"vdcc%d" %v,"CaConc.dat"),os.path.join(resultPath,folder,"CaConc.dat"))
        seed_range_upper=int(10000*vdcc_dist[v])+1
        print(v,seed_range_upper)
        if seed_range_upper>1:
            run_and_average(resultPath=resultPath,dir=folder,seeds=range(1,seed_range_upper))
    #"""
