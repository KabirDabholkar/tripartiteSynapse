#!/usr/bin/python

from numpy import *
from scipy.integrate import *
from random import randint
#from pylab import *
import os, sys
from time import time
from peaks import *

class analysis:

    #   Set dataPath and resultPath
    def __init__(self, dp, rp):
        self.dataPath = dp
        self.resultPath = rp

    #   Get list of directories in the path starting with 's_'
    def getDirs(self, path, sstr='s_'):
        dirs = [d for d in os.listdir(path) if os.path.isdir(path + '/' + d) and sstr in d]
        dirs.sort()
        return(dirs)

    #   Get data from dataFile in array
    def getData(self, dataFile):
        data = []
        f = open(dataFile, "r")
        for l in f:
            if not l.startswith('#'): data.append([float(x) for x in l.strip('\n').strip(" ").split(" ")])
        return (data)

    #   Average Over all Seeds
    def avg_dat(self, inFile="/dat/ca.dat", outFile="/ca.dat"):
        print("\nCalculating Average of", inFile)
        # Get list of directories in data path
        dirs = self.getDirs(self.dataPath)
        self.seeds = len(dirs)
        print(self.seeds)

        j=1
        avg = genfromtxt(self.dataPath+'/'+dirs[0]+inFile)
        l = len(avg)
        for i in range(1,self.seeds):
            temp = genfromtxt(self.dataPath+'/'+dirs[i]+inFile, invalid_raise=False)
            if len(temp) != l:
                l = len(temp)
                print(i, dirs[i], l, temp[-1])    
            else:
                avg += temp #genfromtxt(self.dataPath+'/'+dirs[i]+inFile)
                j += 1
        avg = avg/j #self.seeds
        print('j = ',j)

        of = self.resultPath + outFile
        print("Writing average to: " + of)
        savetxt(of, avg, fmt='%.6f')

    #   Ca Concentration Calculation
    def conc_calc(self, step=5, inFile="/ca.dat", outFile="/CaConc"):
        data = genfromtxt(self.resultPath + inFile, usecols=(0,1), unpack=True)
        print("Calculating Calcium Concentration...")

        c_tc = multiply(data[0],data[1])
        dt=step*(data[0][1]-data[0][0])
        c_out = []
        for i in range(0,len(data[0])-step-1,step):
            c_out.append([data[0][i], (c_tc[i+step]-c_tc[i])/dt])

        print("Writing Ca Conc. to file:" + outFile)
        savetxt(self.resultPath + outFile, c_out, fmt='%.6f')

    #   Get Vesicle Release Statistics for PPF
    def relppf(self, isi, vdcc, resample=1000, tc=0.02): # isi in ms
        n = 2
        self.ts = [(i*isi+2.0)/1000.0 for i in range(n)]
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

            nRel = [0]*n # [Rel1, Rel2,... Reln]
            pr = []
            cp = [0]*4 # [P00, P01, P10, P11]
            for d in [(self.dataPath + '/' + dir + '/dat/') for dir in dirs]:
                fpath = (d + "/rel.dat")
                f = open(fpath, 'r')

                #   Get no. of vesicles released after AP specified by self.ts
                temp = [0]*n
                p = [0, 0]
                time = 0
                for line in f:
                    time = float(line.strip("\n").split(" ")[0])
                    for i in range(len(self.ts)):
                        if (time>self.ts[i] and time<self.ts[i]+tc):
                            temp[i] = 1

                    '''
                    for i in range(len(self.ts)):
                        if (time>self.ts[i] and time<self.ts[i]+tc):
                            p[i] = 1
                    '''

                for i in range(n):
                    if temp[i] == 1: nRel[i] += 1

                '''
                # Calculate Conditional Release Probabilities
                if(p[0]==0 and p[1]==0): cp[0] += 1
                if(p[0]==0 and p[1]==1): cp[1] += 1
                if(p[0]==1 and p[1]==0): cp[2] += 1
                if(p[0]==1 and p[1]==1): cp[3] += 1
                '''
            '''
            print cp
            cp = [float(i)/ndirs for i in cp]
            print cp
            pp = [cp[0]/(cp[0]+cp[1]), cp[1]/(cp[0]+cp[1]), cp[2]/(cp[2]+cp[3]), cp[3]/(cp[2]+cp[3])]
            pp1 = concatenate((cp, pp), axis=0)
            '''
            for i in range(n):
                pr.append(nRel[i]/float(len(dirs)))
            for i in range(1,n):
                pr.append(pr[i]/pr[0])

            '''
            for i in range(8):
                pr.append(pp1[i])
            '''

            prs.append(pr)

        m = mean(prs, axis=0)
        s = std(prs, axis=0)

        #r = open(self.resultPath + '/result', 'w')
        result = list(concatenate((array(list(zip(m,s))).flatten(),[isi, vdcc]), axis=0))
        print(result)
        savetxt(self.resultPath + '/result', [result], fmt=['%0.4f']*6+['%d']*2, delimiter="\t")

        os.system("cat " + self.dataPath + "/*/dat/rel.dat > " + self.resultPath + "/vesRel")
        os.system("cat " + self.dataPath + "/*/dat/vdcc.async_*.dat > " + self.resultPath + "/asyncRel")
        os.system("cat " + self.dataPath + "/*/dat/vdcc.sync_*.dat > " + self.resultPath + "/syncRel")

    #   Get Vesicle Release Statistics for PTP
    def relptp(self, n, isi, resample=1000, tc=0.02): # isi in ms

        self.ts = [(i*isi+2.0)/1000.0 for i in range(n)]
        alldirs = self.getDirs(self.dataPath)
        ndirs = len(alldirs)
        print('seeds: ', ndirs)

        for d in [(self.dataPath + '/' + dir + '/dat/') for dir in alldirs]:
            if not os.path.exists(d+'/rel.dat'):
                os.system("cd " + d + "; cat vdcc.* > rel.dat")

        prs = []
        for r in range(resample):
            if (r+1)%100==0: print('resampling:', r+1)
            x=[randint(0,ndirs-1) for p in range(0,ndirs)]
            dirs = [alldirs[i] for i in x]

            nRel = [0]*len(self.ts) # [Rel1, Rel2,... Reln]
            pr = []
            for d in [(self.dataPath + '/' + dir + '/dat/') for dir in dirs]:
                fpath = (d + "/rel.dat")
                f = open(fpath,'r')

                #   Get no. of vesicles released after AP specified by self.ts
                temp = [0]*n
                for line in f:
                    time = float(line.strip("\n").split(" ")[0])
                    for i in range(len(self.ts)):
                        if (time>self.ts[i] and time<self.ts[i]+tc):
                            #nRel[i] += 1
                            temp[i] = 1
                for i in range(n):
                    if temp[i] == 1: nRel[i] += 1

            for i in range(n):
                pr.append(nRel[i]/float(len(dirs)))
            for i in range(n):
                pr.append(pr[i]/pr[0])

            prs.append(pr)

        m = mean(prs, axis=0)
        s = std(prs, axis=0)

        os.system("cat " + self.dataPath + "/*/dat/rel.dat > " + self.resultPath + "/vesRel")
        os.system("cat " + self.dataPath + "/*/dat/vdcc.async_*.dat > " + self.resultPath + "/asyncRel")
        os.system("cat " + self.dataPath + "/*/dat/vdcc.sync_*.dat > " + self.resultPath + "/syncRel")

        #   Print Results
        r = open(self.resultPath + '/result', 'w')
        savetxt(r,array(zip(range(1,n+1),m[:n],s[:n],m[n:],s[n:])),fmt=['%d','%0.4f','%0.4f','%0.4f','%0.4f'], delimiter="\t")

     #   Calculate Current (pA)
    def fluxCurrent(self, inFile, outFile, step=10, ncharge=1, line=1):
        dataFile = self.resultPath + inFile
        data = self.getData(dataFile)

        charge = 1#1.602e-7 #pico Coulomb
        c_out = []
        dt=step*(data[1][0]-data[0][0])
        for i in range(0,len(data)-step-1,step):
            temp = [data[i+step][0]-dt/2]
            for l in range(1,line+1):
                temp.append(ncharge*charge*(data[i+step][l]-data[i][l])/dt)
            #c_out.append([data[i+step][0]-dt/2, ncharge*charge*(data[i+step][line]-data[i][line])/dt])
            c_out.append(temp)

        of = self.resultPath + outFile
        print("\nWriting average to file: " + of)
        outfile = open(of,'w')
        for l in c_out:
            s = ""
            for d in l: s += str(d) + " "
            outfile.write(s + '\n')
        outfile.close()

    def caStat(self, outFile='/caStat.dat', showFig=False):
        data = genfromtxt(self.resultPath + '/CaConc', unpack=True)
        #print data

        pk = detect_peaks(data[1], mph=2, mpd=360, threshold=0, show=showFig)
        pkValue = [data[1][i] for i in pk]
        pkTime = [data[0][i] for i in pk]

        cumVal = []
        dt = data[0][1]-data[0][0]
        for p in [p-30 for p in pk]:
            dataChunk = data[1][p:p+400]
            cumVal.append(simps(dataChunk, dx=dt))
            #plot(data[0][p:p+400], dataChunk)

        #show()
        print(range(1,len(pk)+1), pkTime, pkValue)
        savetxt(self.resultPath + outFile,zip(range(1,len(pk)+1), pkTime, pkValue, cumVal), \
        fmt=('%d','%0.6f','%0.2f','%0.4f'), delimiter='\t')
