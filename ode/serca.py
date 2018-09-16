from scipy.integrate import *
from matplotlib.pyplot import *
from math import *
import numpy as np

cai = 100.0e-9
cae = 1.0e-6

# Reaction Rates
kx1_x1a=2*1.0e8
kx1a_x2=1.0e8
kx1a_x1=83.666
kx2_x1a=2*83.666
kx2_y2=0.6
ky2_x2=0.097
ky2_y1a=2*30.015
ky1a_y1=30.015
ky1a_y2=1.0e5
ky1_y1a=2*1.0e5
ky1_x1=0.4
kx1_y1=1.20e-3

f = 1
#modification in reaction rates
kx2_y2=0.6
ky2_x2=4.118


ky1_x1=f*0.4
kx1_y1=f*1.20e-3

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
	
	dx1 = fx1*(-kx1_x1a*cai-kx1_y1)+fx1a*kx1a_x1+fy1*ky1_x1
	dx1a = fx1a*(-kx1a_x2*cai-kx1a_x1)+fx1*cai*kx1_x1a+fx2*kx2_x1a
	dx2 = fx2*(-kx2_y2-kx2_x1a)+fx1a*cai*kx1a_x2+fy2*ky2_x2

	dy1 = fy1*(-ky1_y1a*cae-ky1_x1)+fy1a*ky1a_y1+fx1*kx1_y1
	dy1a = fy1a*(-ky1a_y2*cae-ky1a_y1)+fy1*cae*ky1_y1a+fy2*ky2_y1a
	dy2 = fy2*(-ky2_x2-ky2_y1a)+fy1a*cae*ky1a_y2+fx2*kx2_y2
	
	dcae = -cae*(fy1a*ky1a_y2 + fy1*ky1_y1a) + (fy1a*ky1a_y1 + fy2*ky2_y1a)
	#print dcae

	return [dx1, dx1a, dx2, dy1, dy1a, dy2, dcae]
	
# Initial Conditions
#v0 = [0.75554529, 0.18054556, 0.010757505, 0.015718635, 0.026363904, 0.011069105, cae]
#v0 = [0.79150804, 0.18920662, 0.011307258, 0.0023747239, 0.0039558979, 0.0016474684, cae]
v0 = [1, 0, 0.0, 0, 0, 0.0, cae]

tstep = 1e-4
tf = 50
t = np.linspace(0, tf, tf/tstep+1)

# Solve ODE
sol = odeint(serca_ode, v0, t)

# calculate cae_ss, steady-state ER calcium at zero flux through pump
# from Higgins et al., 2006 p.155
K1 = sqrt((kx2_x1a*kx1a_x1)/(kx1_x1a*kx1a_x2))
K2 = ky2_x2/kx2_y2
K3 = sqrt((ky1_y1a*ky1a_y2)/(ky2_y1a*ky1a_y1))
K4 = kx1_y1/ky1_x1

#print(K1**2*K2*K3**2*K4)
cae_ss = cai/(K1*K3*sqrt(K2*K4))
print("cae_ss = %g" %(cae_ss))

'''
for i in range(6):
	plot(t,sol[:,i])
show()
#close()
'''
grid(True)
plot(t,sol[:,6])
show()
