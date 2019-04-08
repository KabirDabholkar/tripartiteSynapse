from analysis import *
import sys

name=sys.argv[1]
sim=sys.argv[2]

def all_analysis(name,sim):
    print(name)
    dataDirName = name.replace(".mdl","")
    dataType = "ppf/"+sim+"/"
    
    M=analysis(dataDirName,dataType)

    M.avg_dat(inFile="/dat/ca.dat", outFile="/ca.dat")
    M.conc_calc(inFile="/ca.dat", outFile="/CaConc")

    isi = int(M.dataDirName.split("I")[1].split("V")[0])
    vdcc = int(M.dataDirName.split("V")[1])
    print('isi: ', isi, '\nvdcc: ', vdcc)

    M.relppf(isi, vdcc, resample=1000)

    M.caStat(showFig=False)

all_analysis(name,sim)
