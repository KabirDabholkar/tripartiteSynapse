import subprocess
from multiprocessing import Pool,cpu_count
import os,sys
from itertools import product,repeat
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline,interp2d,RectBivariateSpline,SmoothBivariateSpline
from scipy.ndimage.filters import gaussian_filter1d
from scipy import stats
import fileinput as fi
from peaks import *
import seaborn as sns
import pickle as pkl
#probably dont't need most of these but im leaving them for now


class data_collector:
    def __init__(self,res_loc,sims,fnames):
        self.res_loc=res_loc
        self.sims=sims
        self.fnames=fnames
    
    def load(self,sim,fname,file_type='ca.dat',cols=1,step=2,contains=None,return_original_time=False):
        loc=os.path.join(self.res_loc,sim,fname)
        seeds=os.listdir(loc)
        temp_list=[]
        if not contains==None:
            seeds=[s for s in seeds if contains in s]
        for seed in seeds:
            file_loc=os.path.join(loc,seed,'dat',file_type)
            temp_list.append(np.genfromtxt(file_loc, unpack=True)[cols,:])
        
        time = np.genfromtxt(file_loc, unpack=True)[0,:]
        data = np.stack(temp_list,axis=1)
        print(time.shape,data.shape)
        conc=(np.tile(time[step::step],(data.shape[1],1)).T*data[step::step]-np.tile(time[0:-step:step],(data.shape[1],1)).T*data[0:-step:step])/(time[step]-time[0])
        conc=np.logical_not(np.sign(conc)==-1)*conc
        
        if not return_original_time:
            time = time[:-step:step]    
        
        return time,conc

    def load_all(self,file_type='ca.dat',cols=1):
        self.data=dict(zip(self.sims,[{}]*len(self.sims)))
        self.time=dict(zip(self.sims,[{}]*len(self.sims)))
        
        for sim,fname in product(sims,fnames):
            self.time[sim][fname],self.data[sim][fname]=self.load(sim,fname,file_type=file_type,cols=cols)
        