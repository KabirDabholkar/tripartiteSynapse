from scipy.integrate import *
from matplotlib.pyplot import *
from scipy.optimize import curve_fit
from math import *
import numpy as np

cai = 100.0e-9
cae = 0.0

f=3.0**(2.0/3.0) #multiplier 

# Reaction Rates
k_orig={
	'kx1_x1a' : 2*1.0e8,
	'kx1a_x2' : 1.0e8,
	'kx1a_x1' : 83.666,
	'kx2_x1a' : 2*83.666,
	'kx2_y2' : 0.6,
	'ky2_x2' : 4.118, #corrected
	'ky2_y1a' : 2*30.015,
	'ky1a_y1' : 30.015,
	'ky1a_y2' : 1.0e5,
	'ky1_y1a' : 2*1.0e5,
	'ky1_x1' : 0.4,
	'kx1_y1' : 1.20e-3,
	'kca' : 0#17.13 #40.45 #update 
}

k_over2=k_orig.copy() 
k_over3=k_orig.copy()
k_over2_leak=k_orig.copy()
k_over3_leak=k_orig.copy()
for par in ['kx1_x1a','kx1a_x2','kx2_y2']:
	k_over2_leak[par]=k_over2[par]=k_orig[par]*2.0**(2.0/3.0) #multiplier
	
for par in ['kx1_x1a','kx1a_x2','kx2_y2']:
	k_over3_leak[par]=k_over3[par]=k_orig[par]*3.0**(2.0/3.0) #multiplier
	
k_over2_leak['kca']=17.13
k_over3_leak['kca']=40.45
k=k_orig.copy()


# Volume Normalisations
vol_er = (3.9*0.1*0.1)
vol_cyt = (4.0*0.5*0.5-vol_er)
vol_tot = vol_cyt + vol_er
ntot = 1
ctot = ntot/vol_tot
fvol_cyt = vol_cyt/vol_tot
fvol_er = vol_er/vol_tot

# SERCA ODE model
def serca_ode(v ,t):
	fx1, fx1a, fx2,	fy1, fy1a, fy2, cae = v 
	for key in k.keys():
		exec(key + " = "+str(k[key]))

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
#v0=[0.4, 0.1, 0.0, 0.4, 0.1, 0.0,cae]
v0 = [0, 0.0, 0.0, 1, 0.0, 0.0, cae]



tstep = 1e-3
tf = 100
t = np.linspace(0, tf, tf/tstep+1)

# Solve ODE
#sol = odeint(serca_ode, v0, t)

def func(x,tau):
	A=250e-6
	return A*(1-np.exp(-x/tau))

#Fit (1-exp(-t/T))
#'''

def ref_time(k_option):
	for par in k.keys():
		k[par]=k_option[par]
	sol = odeint(serca_ode, v0, t)	
	popt, pcov = curve_fit(func,t,sol[:,6])
	print "refilling time: "+str(popt[0])+"s"
	plot(t,sol[:,6],label='Cae')
	plot(t,func(t,popt[0]),label='A*(1-np.exp(-x/tau))')
	legend()
	#show()
	return popt[0]

#ref_time(k_over3_leak)

with open('Refilling_time.txt','w') as outfile:
	outfile.write("Normal refilling time: "+str(ref_time(k_orig))+"\n")
	outfile.write("SERCA Overload(2x) + leak, refilling time: "+str(ref_time(k_over2_leak))+"\n")
	outfile.write("SERCA Overload(3x) + leak, refilling time: "+str(ref_time(k_over3_leak))+"\n")
'''
tau=[]
time_factor_range=np.arange(0.5,3.0,0.1)
for tf in time_factor_range:
	for key in k.keys():
		k[key]=k_orig[key]/tf
	sol = odeint(serca_ode, v0, t)
	popt, pcov = curve_fit(func,t,sol[:,6])
	tau.append(popt[0])	
print tau
print time_factor_range

plot(time_factor_range,tau)
show()
'''
'''

state_label=['x0','x1','x2','y0','y1','y2']
for i in range(6):
	plot(t,sol[:,i],label=state_label[i])
legend()
show()

#grid(True)
plot(t,sol[:,6])
plot(t,func(t,popt[0]))
show()
'''
