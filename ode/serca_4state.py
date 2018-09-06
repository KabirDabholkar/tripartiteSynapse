#modelling the 4 state serca Higgins model
from scipy import optimize
from scipy.integrate import *
from matplotlib.pyplot import *
import numpy as np

# Volume Normalisations
vol_er = (3.9*0.1*0.1) # all units in um
vol_cyt = (4.0*0.5*0.5-vol_er)
vol_tot = vol_cyt + vol_er
Pt=15e-6
gamma=vol_er/vol_cyt
c=100e-9
t0_cond=[Pt,0,0,0,0]
initials=[1,1,1,1,1,1,1,1]
eq_value=[Pt/4,Pt/4,gamma*Pt/4,gamma*Pt/4,100e-9,250e-6]

def ode(x,t,k1,k_1,k2,k_2,k3,k_3,k4,k_4):
    x1,x2,y1,y2,ce=x
    dx1 = k4*y1/gamma-k_4*x1-k1*c**2*x1+k_1*x2
    dx2 = k1*c**2*x1-k_1*x2+k_2*y2/gamma-k2*x2
    dy1 = k2*y2-k_3*ce**2*y1-k4
    dy2 = k2*x2/gamma-k_2*y2-k3*y2+k_3*ce**2*y1
    dce = 2*k3*y2-2*k_3*ce**2*y1
    return [dx1,dx2,dy1,dy2,dce]


t=np.linspace(0,0.1,100)

def func_tbs(params, value):
	X=odeint(ode,t0_cond,t,args=tuple(params))
	return X[-1]-value

#print func_tbs([5,10],value)

sol=optimize.root(func_tbs,initials,args=(eq_value,),method='hybr',jac=False,tol=1e-3)
#plot(odeint(ode,t0_cond,t,args=tuple(sol['x'])))
#plot(odeint(ode,t0_cond,t,args=tuple([1,1,1,1,1,1,1,1])))
#show()
#plot(func_tbs())