from pylab import *
from numpy import *

# Note: stimuli are 5 ms wide.

n = 50		# n = number of pules
f = 10		# Hz
d = 0		# d = delay in first pulse (ms)

p = random.poisson(1000/f, n)
#print 1.0*sum(p)/len(p)	#mean

P = [p[0]]
for i in range(1,len(p)):
	P.append(p[i]+P[i-1]+d)
	
print P

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
for (i,fn) in enumerate(infile):
	fn = fn.split(".")[0]+"_"
	
	outfile[i] = "poisson/"+fn+str(n)+"_"+str(f)+"hz_poisson.dat"
	print outfile[i]

# Generate Stim for PTP
for i in range(len(infile)):

	# Read The Stim
	idata = genfromtxt(infile[i])

	# Make The Stim
	ofile = open(outfile[i], 'w')
	duration = idata[-1][0]

	# Post-Tetanic Pulse with number of spikes n and interspike interval isi
	for t in P:
		tdelay = d+t
		for j in range(len(idata)):
			ofile.write(str(idata[j][0]+tdelay)+"\t"+str(idata[j][1])+"\n")

