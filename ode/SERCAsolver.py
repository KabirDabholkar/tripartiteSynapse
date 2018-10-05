from scipy.integrate import *
from scipy.optimize import *
from matplotlib.pyplot import *
from math import *
import numpy as np

cai = 100.0e-9
cae = 1.0e-6

# Reaction Rates
k_orig={
	'kx1_x1a' : 2*1.0e8,
	'kx1a_x2' : 1.0e8,
	'kx1a_x1' : 83.666,
	'kx2_x1a' : 2*83.666,
	'kx2_y2' : 0.6,
	'ky2_x2' : 4.118,
	'ky2_y1a' : 2*30.015,
	'ky1a_y1' : 30.015,
	'ky1a_y2' : 1.0e5,
	'ky1_y1a' : 2*1.0e5,
	'ky1_x1' : 0.4,
	'kx1_y1' : 1.20e-3,
	'kca' : 18.0 #update
}

#########for 750
for par in ['kx1_x1a','kx1a_x2','kx2_y2']:
	k_orig[par]=k_orig[par]*2.0**(2.0/3.0) #multiplier

k=k_orig.copy()

#modification in reaction rates
#kx2_y2=0.6
#ky2_x2=4.118

# Volume Normalisations
vol_er = (3.9*0.1*0.1)
vol_cyt = (4.0*0.5*0.5-vol_er)
vol_tot = vol_cyt + vol_er
ntot = 1
ctot = ntot/vol_tot
fvol_cyt = vol_cyt/vol_tot
fvol_er = vol_er/vol_tot

# SERCA ODE model
def serca_ode(v,t):
	fx1, fx1a, fx2,	fy1, fy1a, fy2, cae = v 
	for key in k.keys():
		exec(key + " = "+str(k[key]))
	'''
	dx1 = fx1*(-kx1_x1a*cai-kx1_y1)+fx1a*kx1a_x1+fy1*ky1_x1
	dx1a = fx1a*(-kx1a_x2*cai-kx1a_x1)+fx1*cai*kx1_x1a+fx2*kx2_x1a
	dx2 = fx2*(-kx2_y2-kx2_x1a)+fx1a*cai*kx1a_x2+fy2*ky2_x2

	dy1 = fy1*(-ky1_y1a*cae-ky1_x1)+fy1a*ky1a_y1+fx1*kx1_y1
	dy1a = fy1a*(-ky1a_y2*cae-ky1a_y1)+fy1*cae*ky1_y1a+fy2*ky2_y1a
	dy2 = fy2*(-ky2_x2-ky2_y1a)+fy1a*cae*ky1a_y2+fx2*kx2_y2
	'''
	dx1 = fx1*(-kx1_x1a*cai-kx1_y1)+fx1a*kx1a_x1+fy1*ky1_x1
	dx1a = fx1a*(-kx1a_x2*cai-kx1a_x1)+fx1*cai*kx1_x1a+fx2*kx2_x1a
	dx2 = fx2*(-kx2_y2-kx2_x1a)+fx1a*cai*kx1a_x2+fy2*ky2_x2

	dy1 = fy1*(-ky1_y1a*cae-ky1_x1)+fy1a*ky1a_y1+fx1*kx1_y1
	dy1a = fy1a*(-ky1a_y2*cae-ky1a_y1)+fy1*cae*ky1_y1a+fy2*ky2_y1a
	dy2 = fy2*(-ky2_x2-ky2_y1a)+fy1a*cae*ky1a_y2+fx2*kx2_y2
	
	dcae = -cae*(fy1a*ky1a_y2 + fy1*ky1_y1a) + (fy1a*ky1a_y1 + fy2*ky2_y1a) + kca*(cai - cae)
	#print dcae

	return [dx1, dx1a, dx2, dy1, dy1a, dy2, dcae]
	
# Initial Conditions
#v0 = [0.75554529, 0.18054556, 0.010757505, 0.015718635, 0.026363904, 0.011069105, cae]
#v0 = [0.79150804, 0.18920662, 0.011307258, 0.0023747239, 0.0039558979, 0.0016474684, cae]
v0 = [0.4, 0.1, 0.0, 0.4, 0.1, 0.0, cae]

tstep = 1e-2
tf = 100
t = np.arange(0, tf, tstep)

# Solve ODE
'''
f_range = np.arange(0,50,5)
eq_cae=[]
for f in f_range:
	k['kca']=f#*k_orig['kca']
	sol = odeint(serca_ode, v0, t)
	print f,sol[-1,-1]
	print k
	eq_cae.append(sol[-1,-1])

grid(True)
plot(f_range,eq_cae)
show()	
'''	

#main stuff 

def func_tbs(f,eq_value):
	print f[0]
	k['kca']=f[0]#*k_orig['kca']
	sol = odeint(serca_ode, v0, t)
	return sol[-1,-1]-eq_value

SOL=fsolve(func_tbs,1,args=(250.0e-6,))
#'''

'''	
k['kx1a_x2']=3.0*k_orig['kx1a_x2']
sol = odeint(serca_ode, v0, t)
print 3.0,sol[-1,-1]
print k
plot(sol[:,-1])
show()	
	
# calculate cae_ss, steady-state ER calcium at zero flux through pump
# from Higgins et al., 2006 p.155
K1 = sqrt((kx2_x1a*kx1a_x1)/(kx1_x1a*kx1a_x2))
K2 = ky2_x2/kx2_y2
K3 = sqrt((ky1_y1a*ky1a_y2)/(ky2_y1a*ky1a_y1))
K4 = kx1_y1/ky1_x1

#print(K1**2*K2*K3**2*K4)
cae_ss = cai/(K1*K3*sqrt(K2*K4))
print("cae_ss = %g" %(cae_ss))
print 'aaa'

'''
'''
for i in range(6):
	plot(t,sol[:,i])
show()
close()
'''
#sol = odeint(serca_ode, v0, t)
#grid(True)
#plot(t,sol[:,6])
#show()


