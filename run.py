import subprocess
import os

initial_conc = range(350,750,50)
serca_f = [2,3,4,5,6]
fnames=[]

for ic in initial_conc:
    for sf in serca_f:
        fname="/home/kabir/Project/tripartiteSynapse/brute750/RSnostim_"+str(sf)+"_"+str(ic)+".mdl"
        fnames.append(fname)
        #with open(fname,'w') as wfile:
        #    mytext=ftext[0]+str(ic)+ftext[1]+str(sf)+ftext[2]
        #    wfile.write(mytext)

p=subprocess.call(["parallel","mcell",":::","/home/kabir/Project/tripartiteSynapse/ryr_750/RSnostim.mdl"])
