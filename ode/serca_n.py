from scipy.integrate import *
from matplotlib.pyplot import *
import numpy as np

cai = 100.0e-9
cae = 0
nserca = 8678.0
N_avo=6.022e23

# Reaction Rates
kx0_x1=2*1.0e8
kx1_x2=1.0e8
kx1_x0=83.666
kx2_x1=2*83.666
kx2_y2=0.6
ky2_x2=0.097
ky2_y1=2*30.015
ky1_y0=30.015
ky1_y2=1.0e5
ky0_y1=2*1.0e5
ky0_x0=0.4
kx0_y0=1.20e-3

'''
f=1

#modification in reaction rates
kx2_y2=0.6
ky2_x2=4.118

ky0_x0=f*0.4
kx0_y0=f*1.20e-3
'''

# Volume Normalisations
vol_er = (3.9*0.1*0.1) # all units in um
vol_cyt = (4.0*0.5*0.5-vol_er)
vol_tot = vol_cyt + vol_er
ntot = 1
ctot = ntot/vol_tot
fvol_cyt = vol_cyt/vol_tot
fvol_er = vol_er/vol_tot

ncai = cai*vol_cyt*6.022*1e8
ncae = cae*vol_er*6.022*1e8

# SERCA ODE model
def serca_ode(v ,t):
	x0, x1, x2, y0, y1, y2, ca_e = v 
	
	dx0 = x0*(-kx0_x1*cai-kx0_y0)+x1*kx1_x0+y0*ky0_x0
	dx1 = x1*(-kx1_x2*cai-kx1_x0)+x0*cai*kx0_x1+x2*kx2_x1
	dx2 = x2*(-kx2_y2-kx2_x1)+x1*ncai*kx1_x2+y2*ky2_x2

	dy0 = y0*(-ky0_y1*ca_e-ky0_x0)+y1*ky1_y0+x0*kx0_y0
	dy1 = y1*(-ky1_y2*ca_e-ky1_y0)+y0*ca_e*ky0_y1+y2*ky2_y1
	dy2 = y2*(-ky2_x2-ky2_y1)+y1*ca_e*ky1_y2+x2*kx2_y2
	
	dcae = -ca_e*(y1*ky1_y2 + y0*ky0_y1) + (y1*ky1_y0 + y2*ky2_y1)

	return [dx0, dx1, dx2, dy0, dy1, dy2, dcae]

#constraints


#just comparing if the serca_ode() function in serca_n1.py serca_n.py are the same. I changed only variable names.
#serca_ode([4,5,3,6,4,3,2],5)
# Initial Conditions
v0 = np.array([0.75554529, 0.18054556, 0.010757505, 0.015718635, 0.026363904, 0.011069105,1])
v0 = v0*nserca/(vol_er/1e18)/N_avo
v0[6]=cae

tstep = 10e-9
tf = 10-6
'''
v0 = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
v0 = [i*nserca for i in v0]
v0.append(ncae)

#tstep = 1e-6
#tf = 1e-3
'''
t = np.linspace(0, tf, tf/tstep+1)

# Solve ODE
sol = odeint(serca_ode, v0, t)
print sol[-1,6]
label = ['x0', 'x1', 'x2', 'y0', 'y1', 'y2', 'cae']
a_ca=subplot(3,1,3)
a_ca.plot(t,sol[:,6], label=label[6])
a_ca.set_title(label[6])
a_ca.legend()

ax=[]
for i in range(6):
	ax.append(subplot(3,3,i+1))
	ax[-1].plot(t,sol[:,i], label=label[i])
	ax[-1].set_title(label[i])
	ax[-1].legend()
show()

#close()
#plot(t,sol[:,6]/23.45, label='Cae')
