from pylab import *

# Note: stimuli are 5 ms wide.

n = 20  # n = number of pules in tetanic pulse
nBurst = 2  # number of bursts
isi = 10e-3 # isi = Inter-Spike Interval in burst
ibi = 5  # inter-burst internal (sec)


infile=[0]*9
outfile=[0]*9

infile[0]="VDCC_PQ_C01.dat"
infile[1]="VDCC_PQ_C10.dat" 
infile[2]="VDCC_PQ_C12.dat"
infile[3]="VDCC_PQ_C21.dat"
infile[4]="VDCC_PQ_C23.dat"
infile[5]="VDCC_PQ_C32.dat"
infile[6]="VDCC_PQ_C34.dat"
infile[7]="VDCC_PQ_C43.dat"
infile[8]="VDCC_PQ_Ca.dat"

# Generate outfile
for (i,f) in enumerate(infile):
    f = f.split(".")[0]+"_"
    #outfile[i] = "ptp/"+f+str(n)+"_"+str(int(1/isi))+"hz_train.dat"
    outfile[i] = f+"n_"+str(n)+"_nB_"+str(nBurst)+"_isi_"+str(int(1000*isi))+"ms_ibi_"+str(ibi)+"s_burstTrain.dat"
    print outfile[i]


# Generate Stim for PTP
for i in range(len(infile)):
    tdelay = 0
    ofile = open(outfile[i], 'w')
    for nb in range(nBurst): #4, 6, 10, 15, 20, 30]:
        tdelay += nb*ibi
        
        # Read The Stim
        idata = genfromtxt(infile[i])
        
        # Post-Tetanic Pulse with number of spikes n and interspike interval isi
        for j in range(n):
            
            #print tdelay
            for k in range(len(idata)):
                ofile.write(str(idata[k][0]+tdelay)+"\t"+str(idata[k][1])+"\n")
            tdelay += isi if j<n-1 else 0

