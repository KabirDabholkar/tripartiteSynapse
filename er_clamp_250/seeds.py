import multiprocessing as mp, os, sys

# MDL File and seeds
if len(sys.argv) == 2: # If no seeds mentioned => seeds = 1
	mdlFile = sys.argv[1]
	seeds = 1
	print "seeds = ", seeds

else:
	if len(sys.argv) == 3:
		mdlFile = sys.argv[1]
		seeds = int(sys.argv[2])
		if seeds > 0:
			print "seeds = ", seeds
		else:
			print "Please enter valid number of seeds"
			sys.exit()
	else:
		print "Please enter valid arguments \n example: \n python seeds.py mdlFile.mdl 5"
		sys.exit()

# Define an output queue
output = mp.Queue()

# define function to be parallelised
def runmdl(i, mdl):
	os.system("mcell -seed %i -logfreq 10000 %s" % (i, mdl))

# Setup a list of processes that we want to run
processes = [mp.Process(target=runmdl, args=(i, mdlFile)) for i in range(1,1+seeds)]

# Run processes
for p in processes:
	p.start()

# Exit the completed processes
for p in processes:
	p.join()

# Get process results from the output queue
results = [output.get() for p in processes]

print(results)
sys.exit()
