#!/usr/bin/env python

from math import *
from scipy import linalg, matrix

#From Higgins et al., 2006:
#K1^2 = 0.7e-12 (M-1s-1)^2
#K1^2 = Keff^2
#Keff^2 = (k-1/k+1)*(k-2/k+2) but k-1 = k-2 and k+1 = k+2
#Keff^2 = (k-1/k+1)^2
# so Keff = sqrt(K1^2)
# Keff = sqrt(0.7e-12) = 8.3666002653407553e-07
# k-1 = Keff*k+1
# k-1 = 8.3666002653407553e-07*1.0e8 = 83.666

cai = 100.0e-9
cae = 250e-6

f = 1
p = 1


kx1_x1a=2*1.0e8
kx1a_x2=1.0e8
kx1a_x1=f*(83.666)
kx2_x1a=f*(2*83.666)
kx2_y2=p*0.6
ky2_x2=p*0.097
ky2_y1a=(1/f)*(2*30.015)
ky1a_y1=(1/f)*(30.015)
ky1a_y2=1.0e5
ky1_y1a=2*1.0e5
ky1_x1=0.4
kx1_y1=1.20e-3

ky2_x2=4.118

# Calculate K1^2*K2*K^3*K4 = exp(dG_ATP/RT)
# where dG_ATP = -50 KJ/mol, R = 8.314e-3 KJ/(mol*K) and T = 310 K
# forward loop rates:
rx1_x1a = kx1_x1a
rx1a_x2 = kx1a_x2
rx2_y2 = kx2_y2
ry2_y1a = ky2_y1a
ry1a_y1 = ky1a_y1
ry1_x1 = ky1_x1


prod_fwd = rx1_x1a*rx1a_x2*rx2_y2*ry2_y1a*ry1a_y1*ry1_x1

# backward loop rate constants:
ry1_y1a = ky1_y1a
ry1a_y2 = ky1a_y2
ry2_x2 = ky2_x2
rx2_x1a = kx2_x1a
rx1a_x1 = kx1a_x1
rx1_y1 = kx1_y1

prod_bwrd = ry1_y1a*ry1a_y2*ry2_x2*rx2_x1a*rx1a_x1*rx1_y1

e_constraint = prod_bwrd/prod_fwd

# calculate cae_ss, steady-state ER calcium at zero flux through pump
# from Higgins et al., 2006 p.155
K1 = sqrt((kx2_x1a*kx1a_x1)/(kx1_x1a*kx1a_x2))
K2 = ky2_x2/kx2_y2
K3 = sqrt((ky1_y1a*ky1a_y2)/(ky2_y1a*ky1a_y1))
K4 = kx1_y1/ky1_x1

print(K1**2*K2*K3**2*K4)

cae_ss = cai/(K1*K3*sqrt(K2*K4))
print("energy constraint = %g" %(e_constraint))
print("cae_ss = %g" %(cae_ss))


# differential equations:
# N.B.! use equation in "number of molecules" form:
#dx1 = fx1*(-kx1_x1a*cai-kx1_y1)+fx1a*kx1a_x1+fy1*ky1_x1
#dx1a = fx1a*(-kx1a_x2*cai-kx1a_x1)+fx1*cai*kx1_x1a+fx2*kx2_x1a
#dx2 = fx2*(-kx2_y2-kx2_x1a)+fx1a*cai*kx1a_x2+fy2*ky2_x2

#dy1 = fy1*(-ky1_y1a*cae-ky1_x1)+fy1a*ky1a_y1+fx1*kx1_y1
#dy1a = fy1a*(-ky1a_y2*cae-ky1a_y1)+fy1*cae*ky1_y1a+fy2*ky2_y1a
#dy2 = fy2*(-ky2_x2-ky2_y1a)+fy1a*cae*ky1a_y2+fx2*kx2_y2

# concentration to total number relations: 
# cx1 = nx1/vol_cyt     nx1 = cx1*vol_cyt     fx1 = nx1/ntot
# cx1a = nx1a/vol_cyt   nx1a = cx1a*vol_cyt   fx1a = nx1a/ntot
# cx2 = nx2/vol_cyt     nx2 = cx2*vol_cyt     fx2 = nx2/ntot

# cy1 = ny1/vol_er      ny1 = cy1*vol_er      fy1 = ny1/ntot
# cy1a = ny1a/vol_er    ny1a = cy1a*vol_er    fy1a = ny1a/ntot
# cy2 = ny2/vol_er      ny2 = cy2*vol_er      fy2 = ny2/ntot

# constraint:
# ntot = nx1 + nx1a + nx2 + ny1 + ny1a + ny2
# ntot = vol_cyt*(cx1 + cx1a + cx2) + vol_er*(cy1 + cy1a + cy2)
# vol_tot = vol_cyt+vol_er
# fvol_cyt = vol_cyt/vol_tot
# fvol_er = vol_er/vol_tot
# ctot = ntot/vol_tot
# ctot = fvol_cyt*(cx1 + cx1a + cx2) + fvol_er*(cy1 + cy1a + cy2)

vol_er = (3.9*0.1*0.1)
vol_cyt = (4.0*0.5*0.5-vol_er)
vol_tot = vol_cyt + vol_er
ntot = 1
ctot = ntot/vol_tot
fvol_cyt = vol_cyt/vol_tot
fvol_er = vol_er/vol_tot


#dx1 = fx1*(-kx1_x1a*cai-kx1_y1)+fx1a*kx1a_x1+fy1*ky1_x1
#dx1a = fx1a*(-kx1a_x2*cai-kx1a_x1)+fx1*cai*kx1_x1a+fx2*kx2_x1a
#dx2 = fx2*(-kx2_y2-kx2_x1a)+fx1a*cai*kx1a_x2+fy2*ky2_x2

#dy1 = fy1*(-ky1_y1a*cae-ky1_x1)+fy1a*ky1a_y1+fx1*kx1_y1
#dy1a = fy1a*(-ky1a_y2*cae-ky1a_y1)+fy1*cae*ky1_y1a+fy2*ky2_y1a
#dy2 = fy2*(-ky2_x2-ky2_y1a)+fy1a*cae*ky1a_y2+fx2*kx2_y2

# constraints:
# fx1+fx1a+fx2+fy1+fy1a+fy2 = 1
# fx2 = 1 - (fx1+fx1a+fy1+fy1a+fy2)
# but: dx2 = fx2*(-kx2_y2-kx2_x1a)+fx1a*kx1a_x2+fy2*ky2_x2
# so: dx2 = (1-(fx1+fx1a+fy1+fy1a+fy2))*(-kx2_y2-kx2_x1a)+fx1a*kx1a_x2+fy2*ky2_x2
#     dx2 = ((fx1+fx1a+fy1+fy1a+fy2)-1)*(kx2_y2+kx2_x1a)+fx1a*kx1a_x2+fy2*ky2_x2
#     dx2 = (fx1*(kx2_y2+kx2_x1a) +fx1a*(kx2_y2+kx2_x1a) +fy1*(kx2_y2+kx2_x1a) +fy1a*(kx2_y2+kx2_x1a) +fy2*(kx2_y2+kx2_x1a) -(kx2_y2+kx2_x1a)+fx1a*kx1a_x2+fy2*ky2_x2
#     dx2+(kx2_y2+kx2_x1a) = fx1*(kx2_y2+kx2_x1a) +fx1a*((kx2_y2+kx2_x1a)+kx1a_x2)+fy1*(kx2_y2+kx2_x1a) +fy1a*(kx2_y2+kx2_x1a) +fy2*((kx2_y2+kx2_x1a)+ky2_x2)

vx1_fx1 = (-kx1_x1a*cai-kx1_y1)
vx1_fx1a = kx1a_x1
vx1_fy1 = ky1_x1

vx1a_fx1a = (-kx1a_x2*cai-kx1a_x1)
vx1a_fx1 = cai*kx1_x1a
vx1a_fx2 = kx2_x1a

vx2n_fx2 = (-kx2_y2-kx2_x1a)
vx2n_fx1a = cai*kx1a_x2
vx2n_fy2 = ky2_x2

vx2c_fx1 = (kx2_y2+kx2_x1a)
vx2c_fx1a = ((kx2_y2+kx2_x1a)+kx1a_x2)
vx2c_fy1 = (kx2_y2+kx2_x1a)
vx2c_fy1a = (kx2_y2+kx2_x1a)
vx2c_fy2 = ((kx2_y2+kx2_x1a)+ky2_x2)

vy1_fy1 = (-ky1_y1a*cae-ky1_x1)
vy1_fy1a = ky1a_y1
vy1_fx1 = kx1_y1

vy1a_fy1a = (-ky1a_y2*cae-ky1a_y1)
vy1a_fy1 = cae*ky1_y1a
vy1a_fy2 = ky2_y1a

vy2_fy2 = (-ky2_x2-ky2_y1a)
vy2_fy1a = cae*ky1a_y2
vy2_fx2 = kx2_y2

# matrix A without contraints:
# cx1        cx1a       cx2        cy1        cy1a       cy2

# vx1_fx1   vx1_fx1a  0         vx1_fy1   0         0
# vx1a_fx1  vx1a_fx1a vx1a_fx2  0         0         0
# 0         vx2n_fx1a vx2n_fx2  0         0         vx2n_fy2
# vy1_fx1   0         0         vy1_fy1   vy1_fy1a  0
# 0         0         0         vy1a_fy1  vy1a_fy1a vy1a_fy2
# 0         0         vy2_fx2   0         vy2_fy1a  vy2_fy2
# molecule number constraint:
# 1         1         1         1         1         1
# concentration constraint:
#            [fvol_cyt,  fvol_cyt,  fvol_cyt,  fvol_er,   fvol_er,   fvol_er]])


A = matrix([[vx1_fx1,   vx1_fx1a,  0.0,       vx1_fy1,   0.0,       0.0],
            [vx1a_fx1,  vx1a_fx1a, vx1a_fx2,  0.0,       0.0,       0.0],
            [0.0,       vx2n_fx1a, vx2n_fx2,  0.0,       0.0,       vx2n_fy2],
            [vy1_fx1,   0.0,       0.0,       vy1_fy1,   vy1_fy1a,  0.0],
            [0.0,       0.0,       0.0,       vy1a_fy1,  vy1a_fy1a, vy1a_fy2],
            [1,         1,         1,         1,         1,        1]])

# molecule number constraint:
#            [1,         1,         1,         1,         1,        1]])
# concentration constraint:
#            [fvol_cyt,  fvol_cyt,  fvol_cyt,  fvol_er,   fvol_er,   fvol_er]])

# matrix B:
B = matrix([[0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [1.0]])

det = linalg.det(A)
X = linalg.solve(A,B)


#### modified #####
fx1, fx1a, fx2, fy1, fy1a, fy2 = X
print "d[Cae]/dt = ",-cae_ss*(fy1a*ky1a_y2 + fy1*ky1_y1a) + (fy1a*ky1a_y1 + fy2*ky2_y1a),'\n'
#### end modify ########


print("Det(A):")
print(det)
print("X:")
print(X)
print("A*X:")
print(A*X)

fx1=X[0,0]
fx1a=X[1,0]
fx2=X[2,0]
fy1=X[3,0]
fy1a=X[4,0]
fy2=X[5,0]

print("serca_x1_feq = %.8g" % (fx1))
print("serca_x1a_feq = %.8g" % (fx1a))
print("serca_x2_feq = %.8g" % (fx2))
print("serca_y1_feq = %.8g" % (fy1))
print("serca_y1a_feq = %.8g" % (fy1a))
print("serca_y2_feq = %.8g" % (fy2))

# forward loop rates:
rx1_x1a = fx1*kx1_x1a*cai
rx1a_x2 = fx1a*kx1a_x2*cai
rx2_y2 = fx2*kx2_y2
ry2_y1a = fy2*ky2_y1a
ry1a_y1 = fy1a*ky1a_y1
ry1_x1 = fy1*ky1_x1

prod_fwd = rx1_x1a*rx1a_x2*rx2_y2*ry2_y1a*ry1a_y1*ry1_x1

# backward loop rates:
ry1_y1a = fy1*ky1_y1a*cae
ry1a_y2 = fy1a*ky1a_y2*cae
ry2_x2 = fy2*ky2_x2
rx2_x1a = fx2*kx2_x1a
rx1a_x1 = fx1a*kx1a_x1
rx1_y1 = fx1*kx1_y1

prod_bwrd = ry1_y1a*ry1a_y2*ry2_x2*rx2_x1a*rx1a_x1*rx1_y1

print("product of forward rates = %g" %(prod_fwd))
print("product of backward rates = %g" %(prod_bwrd))


# forward loop rate constants:
rx1_x1a = kx1_x1a*cai
rx1a_x2 = kx1a_x2*cai
rx2_y2 = kx2_y2
ry2_y1a = ky2_y1a
ry1a_y1 = ky1a_y1
ry1_x1 = ky1_x1

prod_fwd = rx1_x1a*rx1a_x2*rx2_y2*ry2_y1a*ry1a_y1*ry1_x1

# backward loop rate constants:
ry1_y1a = ky1_y1a*cae
ry1a_y2 = ky1a_y2*cae
ry2_x2 = ky2_x2
rx2_x1a = kx2_x1a
rx1a_x1 = kx1a_x1
rx1_y1 = kx1_y1

prod_bwrd = ry1_y1a*ry1a_y2*ry2_x2*rx2_x1a*rx1a_x1*rx1_y1

print("product of forward rate constants = %g" %(prod_fwd))
print("product of backward rate constants = %g" %(prod_bwrd))

# Calculate K1^2*K2*K^3*K4 = exp(dG_ATP/RT)
# where dG_ATP = -50 KJ/mol, R = 8.314e-3 KJ/(mol*K) and T = 310 K
# forward loop rates:
rx1_x1a = kx1_x1a
rx1a_x2 = kx1a_x2
rx2_y2 = kx2_y2
ry2_y1a = ky2_y1a
ry1a_y1 = ky1a_y1
ry1_x1 = ky1_x1


prod_fwd = rx1_x1a*rx1a_x2*rx2_y2*ry2_y1a*ry1a_y1*ry1_x1

# backward loop rate constants:
ry1_y1a = ky1_y1a
ry1a_y2 = ky1a_y2
ry2_x2 = ky2_x2
rx2_x1a = kx2_x1a
rx1a_x1 = kx1a_x1
rx1_y1 = kx1_y1

prod_bwrd = ry1_y1a*ry1a_y2*ry2_x2*rx2_x1a*rx1a_x1*rx1_y1

e_constraint = prod_bwrd/prod_fwd

# calculate cae_ss, steady-state ER calcium at zero flux through pump
# from Higgins et al., 2006 p.155
K1 = sqrt((kx2_x1a*kx1a_x1)/(kx1_x1a*kx1a_x2))
K2 = ky2_x2/kx2_y2
K3 = sqrt((ky1_y1a*ky1a_y2)/(ky2_y1a*ky1a_y1))
K4 = kx1_y1/ky1_x1

cae_ss = cai/(K1*K3*sqrt(K2*K4))

print("energy constraint = %g" %(e_constraint))
print("cae_ss = %g" %(cae_ss))

cx1 =  fx1/vol_cyt
cx1a = fx1a/vol_cyt
cx2 =  fx2/vol_cyt
cy1 =  fy1/vol_er
cy1a = fy1a/vol_er
cy2 =  fy2/vol_er

Xp = matrix([[fx1],
             [fx1a],
             [fx2],
             [fy1],
             [fy1a],
             [fy2]])

print("MCell's A*X:")
print(A*Xp)



# matrix A with contraints:
# x1        x1a       x2        y1        y1a       y2
# vx1_fx1   vx1_fx1a  0         vx1_fy1   0         0
# vx1a_fx1  vx1a_fx1a vx1a_fx2  0         0         0
# vx2c_fx1  vx2c_fx1a  0        vx2c_fy1  vx2c_fy1a vx2c_fy2
# vy1_fx1   0         0         vy1_fy1   vy1_fy1a  0
# 0         0         0         vy1a_fy1  vy1a_fy1a vy1a_fy2
# 0         0         vy2_fx2   0         vy2_fy1a  vy2_fy2

A = matrix([[vx1_fx1,   vx1_fx1a,  0,         vx1_fy1,   0,         0],
            [vx1a_fx1,  vx1a_fx1a, vx1a_fx2,  0,         0,         0],
            [vx2c_fx1,  vx2c_fx1a, 0,         vx2c_fy1,  vx2c_fy1a, vx2c_fy2],
            [vy1_fx1,   0,         0,         vy1_fy1,   vy1_fy1a,  0],
            [0,         0,         0,         vy1a_fy1,  vy1a_fy1a, vy1a_fy2],
            [0,         0,         vy2_fx2,   0,         vy2_fy1a,  vy2_fy2]])

# matrix B with contraints:
B = matrix([[0],
            [0],
            [(kx2_y2+kx2_x1a)],
            [0],
            [0],
            [0]])

#det = linalg.det(A)
#X = linalg.solve(A,B)

#print(det)
#print(X)
