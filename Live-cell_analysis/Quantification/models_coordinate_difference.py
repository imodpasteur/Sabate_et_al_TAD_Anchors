import numpy as np
from os import walk
import sys
import random
from scipy.optimize import curve_fit
from numpy import sqrt, pi, exp, linspace
from time import process_time
from scipy import integrate
from scipy import ndimage
import pandas as pd
import scipy.stats as stats


class Models(object):

    def __init__(self,varloc_x,varloc_y,varloc_z):
        self.varloc_x=varloc_x
        self.varloc_y=varloc_y
        self.varloc_z=varloc_z


    def gaussian(self,x,amplitude,var):
        lenvect=int(len(x)//3)
        g=np.zeros(len(x))
        f=1/((2*np.pi)**.5)
        A=(f/((var+self.varloc_x)**0.5))
        B=(np.divide(np.power(x[:lenvect],2) ,var+self.varloc_x))
        g[:lenvect]=(f/((var+self.varloc_x)**0.5))*amplitude*exp(-0.5* (np.divide(np.power(x[:lenvect],2), var+self.varloc_x)))
        g[lenvect:2*lenvect]=(f/((var+self.varloc_y)**0.5))*amplitude*exp(-0.5* (np.divide(np.power(x[lenvect:2*lenvect],2) ,var+self.varloc_y)))
        g[2*lenvect:]=(f/((var+self.varloc_z)**0.5))*amplitude*exp(-0.5* (np.divide(np.power(x[2*lenvect:],2) ,var+self.varloc_z)))
        return(g.ravel())


    def gaussian_integrated_tmp0(self,x,amplitude,var_loc,var1,var2):
        f=1/((2*np.pi)**.5)
        gauss = lambda var : (f/((var+var_loc)**0.5))*amplitude*exp(-(0.5*(np.divide(x**2 ,var+var_loc))))
        gaussInteg= integrate.quad(gauss, var1, var2)[0]
        return(gaussInteg)


    def gaussian_integrated_tmp0_weightedRandomLoading(self,x,amplitude,var_loc,var1,var2):
        f=1/((2*np.pi)**.5)
        var_l=min(var1,var2)
        var_h=max(var1,var2)
        gauss = lambda var : ((var/(var_l-var_h))+(var_l/(var_h-var_l))+2)*(f/((var+var_loc)**0.5))*amplitude*exp(-(0.5 * (np.divide(x**2 ,var+var_loc))))
        gaussInteg= integrate.quad(gauss, var1, var2)[0]
        return(gaussInteg)


    def gaussian_integrated_fast(self, x, amplitude, var1, var2, weightRandomLoading=False):
        # weightRandomLoading = False assumes a single extrusion speed. If True, assumes that cohesin can switch from V0 to V0/2 upon blocking at one CTCF anchor.
        lenvect=int(len(x)//3)
        g=np.zeros(len(x))
        if weightRandomLoading:
            gaussian_integrated_tmp1_weightedRandomLoading=np.vectorize(self.gaussian_integrated_tmp0_weightedRandomLoading, otypes=[float])
            g[:lenvect]=gaussian_integrated_tmp1_weightedRandomLoading(x[:lenvect],amplitude,self.varloc_x,var1,var2)
            g[lenvect:2*lenvect]=gaussian_integrated_tmp1_weightedRandomLoading(x[lenvect:2*lenvect], amplitude,self.varloc_y, var1, var2)
            g[2*lenvect:]=gaussian_integrated_tmp1_weightedRandomLoading(x[2*lenvect:], amplitude, self.varloc_z, var1, var2)
        else:
            gaussian_integrated_tmp1=np.vectorize(self.gaussian_integrated_tmp0,otypes=[float])
            g[:lenvect]=gaussian_integrated_tmp1(x[:lenvect],amplitude,self.varloc_x,var1,var2)
            g[lenvect:2*lenvect]=gaussian_integrated_tmp1(x[lenvect:2*lenvect],amplitude,self.varloc_y,var1,var2)
            g[2*lenvect:]=gaussian_integrated_tmp1(x[2*lenvect:],amplitude,self.varloc_z,var1,var2)
        return(g)


    def gaussian_integrated_full(self, x, ampl_loop,var_loop, ampl_free,var_free, ampl_extr, weightRandomLoading=False):
        return(self.gaussian(x,ampl_loop,var_loop)   + self.gaussian(x,ampl_free,var_free)   + self.gaussian_integrated_fast(x,ampl_extr,var_loop,var_free,weightRandomLoading) )
