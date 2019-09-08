
# coding: utf-8

# In[2]:


import numpy as np
from scipy.integrate import *
from random import randint
import matplotlib.pyplot as plt
import os, sys
from time import time
from peaks import *


# In[15]:


class analysis:
    #   Get list of directories in the path starting with 's_'
    def __init__(self, dataDirName = "RSI40V90", dataType = "ppf/AD/"):
        self.dataDirName = dataDirName
        self.dataType = dataType

        self.dataPath = "/data/kabir/output/" + self.dataType + self.dataDirName
        self.resultPath = "/home/kabir/Project/tripartiteSynapse/results/" + self.dataType + self.dataDirName

    def getDirs(self, path, sstr='s_'):
        dirs = [d for d in os.listdir(path) if os.path.isdir(path + '/' + d) and sstr in d]
        dirs.sort()
        return(dirs)
    
    def makeDirs(self):
        if not os.path.exists(self.resultPath):
            os.makedirs(self.resultPath)
    
    #   Average Over all Seeds
    def avg_dat(self, inFile="/dat/ca.dat", outFile="/ca.dat"):
        print("\nCalculating Average of", inFile)
        print(self.dataDirName)

        # Get list of directories in data path
        dirs = self.getDirs(self.dataPath)
        seeds = len(dirs)
        print('seeds: ', seeds)

        j=1
        avg = np.genfromtxt(self.dataPath+'/'+dirs[0]+inFile)
        l = len(avg)
        for i in range(1, seeds):
            temp = np.genfromtxt(self.dataPath+'/'+dirs[i]+inFile, invalid_raise=False)
            if len(temp) != l:
                l = len(temp)
                print(i, dirs[i], l, temp[-1])
            else:
                avg += temp #genfromtxt(dataPath+'/'+dirs[i]+inFile)
                j += 1
        avg = avg/j #seeds
        #print('j = ', j)
        
        self.makeDirs()
        
        of = self.resultPath + outFile
        print("Writing average to: " + of)
        #print(avg)
        np.savetxt(of, avg, fmt='%.6f')

    #   Ca Concentration Calculation
    def conc_calc(self, step=5, inFile="/ca.dat", outFile="/CaConc"):
        data = np.genfromtxt(self.resultPath + inFile, usecols=(0,1), unpack=True)
        print("Calculating Calcium Concentration...")

        c_tc = np.multiply(data[0],data[1])
        dt=step*(data[0][1]-data[0][0])
        c_out = []
        for i in range(0,len(data[0])-step-1,step):
            c_out.append([data[0][i], (c_tc[i+step]-c_tc[i])/dt])
        self.makeDirs()
        print("Writing Ca Conc. to file:" + outFile)
        np.savetxt(self.resultPath + outFile, c_out, fmt='%.6f')



    #avg_dat(inFile="/dat/ca.dat", outFile="/ca.dat")
    #conc_calc(inFile="/ca.dat", outFile="/CaConc")

    #   Get Vesicle Release Statistics for PPF
    def relppf(self, isi, vdcc, resample=1000, tc=0.02): # isi in ms
        n = 2 # number of AP
        ts = [(i*isi+2.0)/1000.0 for i in range(n)]
        alldirs = self.getDirs(self.dataPath)
        ndirs = len(alldirs)
        print('seeds: ', ndirs)

        for d in [(self.dataPath + '/' + dir + '/dat/') for dir in alldirs]:
            os.system("cd " + d + "; cat vdcc.* > rel.dat")

        prs = []
        for r in range(resample):
            if (r+1)%100==0: print('resampling:', r+1)
            x=[randint(0,ndirs-1) for p in range(0,ndirs)]
            dirs = [alldirs[i] for i in x]

            nRel = np.zeros(n) # [Rel1, Rel2,... Reln]
            pr = []
            cp = np.zeros(4) # [P00, P01, P10, P11]
            for d in [(self.dataPath + '/' + dir + '/dat/') for dir in dirs]:
                fpath = (d + "/rel.dat")
                f = open(fpath, 'r')

                #   Get no. of vesicles released after AP specified by ts
                temp = np.zeros(n)
                p = [0, 0]
                time = 0
                for line in f:
                    time = float(line.strip("\n").split(" ")[0])
                    for i in range(len(ts)):
                        if (time>ts[i] and time<ts[i]+tc):
                            temp[i] = 1

                    for i in range(len(ts)):
                        if (time>ts[i] and time<ts[i]+tc):
                            p[i] = 1

                for i in range(n):
                    if temp[i] == 1: nRel[i] += 1


                # Calculate Conditional Release Probabilities
                if(p[0]==0 and p[1]==0): cp[0] += 1
                if(p[0]==0 and p[1]==1): cp[1] += 1
                if(p[0]==1 and p[1]==0): cp[2] += 1
                if(p[0]==1 and p[1]==1): cp[3] += 1


            #print('num of rel', cp)
            cp = [float(i)/ndirs for i in cp]
            #print('fraction of rel', cp)
            print('cp:',cp)
            pp = [cp[0]/(cp[0]+cp[1]), cp[1]/(cp[0]+cp[1]), cp[2]/(cp[2]+cp[3]), cp[3]/(cp[2]+cp[3])]
            #print('pp: ', pp)

            for i in range(n):
                pr.append(nRel[i]/float(len(dirs)))
            for i in range(1,n):
                pr.append(pr[i]/pr[0])


            for i in range(4):
                pr.append(pp[i])

            #print('pr: ', pr)
            prs.append(pr)

        m = np.mean(prs, axis=0)
        s = np.std(prs, axis=0)

        self.makeDirs()

        result = np.concatenate((np.array(list(zip(m,s))).flatten(),[isi, vdcc]), axis=0)
        result = list(result)
        print('Vesicle release stats:\n', result)
        header = 'p1\tep1\t\tp2\t\tep2\t\tppr\t\teppr\tP00\t\teP00\tP01\t\teP01\tP10\t\teP10\tP11\t\teP11\tISI\tVDCC\n'
        np.savetxt(self.resultPath + '/result', [result], header=header, fmt=['%0.4f']*14+['%d']*2, delimiter="\t")

        os.system("cat " + self.dataPath + "/*/dat/rel.dat > " + self.resultPath + "/vesRel")
        os.system("cat " + self.dataPath + "/*/dat/vdcc.async_*.dat > " + self.resultPath + "/asyncRel")
        os.system("cat " + self.dataPath + "/*/dat/vdcc.sync_*.dat > " + self.resultPath + "/syncRel")


    #isi = int(dataDirName.split("I")[1].split("V")[0])
    #vdcc = int(dataDirName.split("V")[1])
    #print('isi: ', isi, '\nvdcc: ', vdcc)

    #relppf(isi, vdcc, resample=1000)

    def caStat(self,outFile='/caStat.dat', showFig=False):
        data = np.genfromtxt(self.resultPath + '/CaConc', unpack=True)
        #print data

        pk = detect_peaks(data[1], mph=2, mpd=360, threshold=0, show=showFig)
        pkValue = [data[1][i] for i in pk]
        pkTime = [data[0][i] for i in pk]

        cumVal = []
        dt = data[0][1]-data[0][0]
        for p in [p-30 for p in pk]:
            dataChunk = data[1][p:p+400]
            cumVal.append(simps(dataChunk, dx=dt))
            #plt.plot(data[0][p:p+400], dataChunk)
        #plt.show()

        isi = int(self.dataDirName.split("I")[1].split("V")[0])
        vdcc = int(self.dataDirName.split("V")[1])
        print('isi: ', isi, '\nvdcc: ', vdcc)
        #print(range(1,len(pk)+1), pkTime, pkValue)
        data = np.array(list(zip(pkTime, pkValue, cumVal))).flatten()
        print('Ca stats:\n', data)
        data = np.concatenate([data, [isi, vdcc]])
        print('Ca stats:\n', data)

        header = 't_pk1\tpk1\t\tc_pk1\tt_pk2\tpk2\t\tc_pk2\tISI\tVDCC\n'
        self.makeDirs()
        np.savetxt(self.resultPath + outFile, [data], fmt=['%0.4f','%0.2f','%0.4f']*2+['%d']*2, header=header, delimiter='\t')


    #caStat(showFig=True)

