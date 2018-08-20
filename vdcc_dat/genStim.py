from pylab import *

# Note: stimuli are 5 ms wide.

#n = 30         # n = number of pules in tetanic pulse
isi = 400e-3 # isi = Inter-Spike Interval

for n in [2]: #4, 6, 10, 15, 20, 30]:
    d = 0e-3    # d = delay in first pulse

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
        outfile[i] = f+str(int(isi*1000))+"ms_ppf.dat"
        print outfile[i]

    # Generate Stim for PTP
    for i in range(len(infile)):

        # Read The Stim
        idata = genfromtxt(infile[i])

                # Make The Stim
        ofile = open(outfile[i], 'w')
        duration = idata[-1][0]

        # Post-Tetanic Pulse with number of spikes n and interspike interval isi
        for i in range(n):
            tdelay = d + i*isi
            for j in range(len(idata)):
                ofile.write(str(idata[j][0]+tdelay)+"\t"+str(idata[j][1])+"\n")
