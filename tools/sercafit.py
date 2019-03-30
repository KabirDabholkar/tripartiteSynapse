import numpy as np
from scipy.optimize import curve_fit
sm=[1,2,3,5,10]
data={}

for m in sm:
	data[str(m)]=np.loadtxt('/data/kabir/output/ppf/RSnostim_750_emptyER_sm/sm'+str(m)+'/s_00001/dat/ca.dat')
	
