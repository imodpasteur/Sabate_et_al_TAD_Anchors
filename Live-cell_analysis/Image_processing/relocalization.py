

import skimage as sk
from skimage.util import random_noise
import numpy as np
import pandas as pd
import scipy.ndimage as scnd
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit,minimize
#import gc
import csv
import sys
import os
from pathlib import Path
#import matplotlib.pyplot as plt

from statsmodels.nonparametric.kernel_regression import KernelReg

from lib.ChromaticAberrations import *
from lib.Utils import *
from lib.DriftCorrection import *




######################################### PATH:



#original image with DC (no CHAB)
imagepath=sys.argv[1]


chromagnonpath=sys.argv[2]

locpath=sys.argv[3]

print('loc path: ',locpath)

driftpath=sys.argv[4]

print('drift path: ',driftpath)










######################################### PARAMETERS :
#size of pixel [X,Y,Z] (in µm)
voxelsize=[0.1205,0.1205,0.2897]#um

delimiter="," # even ',' or ';'


#image_photons_count=(image+intercept)*rescale
photon_rescale=0.00604
photon_intercept=-507


######################################### HYPER PARAMETERS :

shiftFit=3 #radius around spot during fit (twice standard deviation is largely enough)

#standard deviation of red and green spots (in µm):
fixedStdRed=[.15,.15,.3] #µm
fixedStdGreen=[.15,.15,.3] #µm

#use the following to fit spot size instead of fixing it
#fixedStdRed=None
#fixedStdGreen=None


refit_with_averaged_intensity=False










######################################### METHODS :





#here, we create mesh grids
def make_XYZ_axis(patchsizeZ=2,patchsizeY=None,patchsizeX=None):
    if patchsizeY is None:
        patchsizeY=patchsizeZ
    if patchsizeX is None:
        patchsizeX=patchsizeZ




    ''' #slow computation
    binsZ=np.arange(0,patchsizeZ)
    binsY=np.arange(0,patchsizeY)
    binsX=np.arange(0,patchsizeX)

    totalsize=len(binsZ)*len(binsY)*len(binsX)

    oo=np.ones((len(binsZ),len(binsY),len(binsX)))
    xyz=np.zeros((3,totalsize))
    index=0
    for z in binsZ:
        for y in binsY:
            for x in binsX:
                xyz[0,index]=z
                xyz[1,index]=y
                xyz[2,index]=x
                index+=1

    print(np.shape(xyz))'''

    #fast computation
    binsZ=np.arange(0,patchsizeZ).reshape(patchsizeZ, 1, 1)
    binsY=np.arange(0,patchsizeY).reshape(1, patchsizeY, 1)
    binsX=np.arange(0,patchsizeX).reshape(1, 1, patchsizeX)

    oo=np.ones((len(binsZ.ravel()),len(binsY.ravel()),len(binsX.ravel())))


    bins_x=oo*binsX
    bins_z=(oo*binsZ)#.transpose()
    bins_y=oo*binsY#np.rot90(oo*binsY, axes=(0,1))

    bins_x_flat=bins_x.ravel()
    bins_y_flat=bins_y.ravel()
    bins_z_flat=bins_z.ravel()

    xyz=np.array([bins_z_flat,bins_y_flat,bins_x_flat])



    return xyz



######################################### METHODS :





def readLoc(locpath,voxelsize,imscale=1,offsetX=0,offsetY=0,delimiter=','):

    loc=pd.read_csv(locpath,delimiter=delimiter)
    #print(loc.keys())
    x1=np.array(loc['Spot_1_X'].tolist())
    y1=np.array(loc['Spot_1_Y'].tolist())
    z1=np.array(loc['Spot_1_Z'].tolist())
    x2=np.array(loc['Spot_2_X'].tolist())
    y2=np.array(loc['Spot_2_Y'].tolist())
    z2=np.array(loc['Spot_2_Z'].tolist())
    frame=np.array(loc['Frame'].tolist())
    track1=np.array(loc['Track_pair'].tolist())
    #print("WARNING: id track2 same as track1")
    #track2=np.array(loc['Track_1_id.1'].tolist())

    x1/=voxelsize[0]
    y1/=voxelsize[1]
    z1/=voxelsize[2]

    x2/=voxelsize[0]
    y2/=voxelsize[1]
    z2/=voxelsize[2]

    x1*=imscale
    y1*=imscale
    x2*=imscale
    y2*=imscale

    x1-=offsetX
    x2-=offsetX
    y1-=offsetY
    y2-=offsetY

    X1=[]
    Y1=[]
    Z1=[]
    X2=[]
    Y2=[]
    Z2=[]
    Frame=[]
    Track1=[]

    for t in np.unique(track1):
        X1.append(x1[track1==t].tolist())
        Y1.append(y1[track1==t].tolist())
        Z1.append(z1[track1==t].tolist())
        X2.append(x2[track1==t].tolist())
        Y2.append(y2[track1==t].tolist())
        Z2.append(z2[track1==t].tolist())
        Frame.append(frame[track1==t].tolist())
        Track1.append(track1[track1==t].tolist())
    return Frame,X1,Y1,Z1,X2,Y2,Z2,Track1



def saveLoc(locpathold,locpathnew,voxelsize,X1,Y1,Z1,precX1,precY1,precZ1,stdX1,stdY1,stdZ1,X2,Y2,Z2,precX2,precY2,precZ2,stdX2,stdY2,stdZ2,photons1,photons2,bckg1,bckg2,distancePixelToInit1,distancePixelToInit2,distancePixelToEdge1,distancePixelToEdge2,goodnessOfFit1,goodnessOfFit2,imscale=1,offsetX=0,offsetY=0,delimiter=','):
    loc=pd.read_csv(locpathold,delimiter=delimiter)



    #print(loc.keys())

    #x1=np.array(loc['Spot_1_X'].tolist())
    #y1=np.array(loc['Spot_1_Y'].tolist())
    #z1=np.array(loc['Spot_1_Z'].tolist())
    #x2=np.array(loc['Spot_2_X'].tolist())
    #y2=np.array(loc['Spot_2_Y'].tolist())
    #z2=np.array(loc['Spot_2_Z'].tolist())
    #frame=np.array(loc['Frame'].tolist())
    track1=np.array(loc['Track_pair'].tolist())
    ind=np.arange(len(track1))
    #print("WARNING: id track2 same as track1")
    #track2=np.array(loc['Track_1_id.1'].tolist())
    precX1list=np.zeros(np.shape(track1))
    precY1list=np.zeros(np.shape(track1))
    precZ1list=np.zeros(np.shape(track1))

    precX2list=np.zeros(np.shape(track1))
    precY2list=np.zeros(np.shape(track1))
    precZ2list=np.zeros(np.shape(track1))


    stdX1list=np.zeros(np.shape(track1))
    stdY1list=np.zeros(np.shape(track1))
    stdZ1list=np.zeros(np.shape(track1))

    stdX2list=np.zeros(np.shape(track1))
    stdY2list=np.zeros(np.shape(track1))
    stdZ2list=np.zeros(np.shape(track1))


    X1list=np.zeros(np.shape(track1))
    Y1list=np.zeros(np.shape(track1))
    Z1list=np.zeros(np.shape(track1))

    X2list=np.zeros(np.shape(track1))
    Y2list=np.zeros(np.shape(track1))
    Z2list=np.zeros(np.shape(track1))

    photons1list=np.zeros(np.shape(track1))
    photons2list=np.zeros(np.shape(track1))
    bckg1list=np.zeros(np.shape(track1))
    bckg2list=np.zeros(np.shape(track1))
    precisionDistlist=np.zeros(np.shape(track1))

    distancelist=np.zeros(np.shape(track1))


    distancePixelToInit1list=np.zeros(np.shape(track1))
    distancePixelToInit2list=np.zeros(np.shape(track1))
    distancePixelToEdge1list=np.zeros(np.shape(track1))
    distancePixelToEdge2list=np.zeros(np.shape(track1))

    goodnessOfFit1list=np.zeros(np.shape(track1))
    goodnessOfFit2list=np.zeros(np.shape(track1))

    for ti,t in enumerate(np.unique(track1)):
        indsub=ind[track1==t].tolist()

        for i,li in enumerate(indsub):
            X1list[li]=((X1[ti][i]+offsetX)/imscale)*voxelsize[0]
            Y1list[li]=((Y1[ti][i]+offsetY)/imscale)*voxelsize[1]
            Z1list[li]=(Z1[ti][i])*voxelsize[2]

            X2list[li]=((X2[ti][i]+offsetX)/imscale)*voxelsize[0]
            Y2list[li]=((Y2[ti][i]+offsetY)/imscale)*voxelsize[1]
            Z2list[li]=(Z2[ti][i])*voxelsize[2]

            precX1list[li]=(precX1[ti][i]/imscale)*voxelsize[0]
            precY1list[li]=(precY1[ti][i]/imscale)*voxelsize[1]
            precZ1list[li]=precZ1[ti][i]*voxelsize[2]

            precX2list[li]=(precX2[ti][i]/imscale)*voxelsize[0]
            precY2list[li]=(precY2[ti][i]/imscale)*voxelsize[1]
            precZ2list[li]=precZ2[ti][i]*voxelsize[2]

            stdX1list[li]=(stdX1[ti][i]/imscale)*voxelsize[0]
            stdY1list[li]=(stdY1[ti][i]/imscale)*voxelsize[1]
            stdZ1list[li]=stdZ1[ti][i]*voxelsize[2]

            stdX2list[li]=(stdX2[ti][i]/imscale)*voxelsize[0]
            stdY2list[li]=(stdY2[ti][i]/imscale)*voxelsize[1]
            stdZ2list[li]=stdZ2[ti][i]*voxelsize[2]

            photons1list[li]=photons1[ti][i]
            photons2list[li]=photons2[ti][i]
            bckg1list[li]=bckg1[ti][i]
            bckg2list[li]=bckg2[ti][i]

            distancePixelToInit1list[li]=distancePixelToInit1[ti][i]
            distancePixelToInit2list[li]=distancePixelToInit2[ti][i]
            distancePixelToEdge1list[li]=distancePixelToEdge1[ti][i]
            distancePixelToEdge2list[li]=distancePixelToEdge2[ti][i]


            goodnessOfFit1list[li]=goodnessOfFit1[ti][i]
            goodnessOfFit2list[li]=goodnessOfFit2[ti][i]

            precisionDistlist[li]=np.sqrt(precX1list[li]**2+precY1list[li]**2+precZ1list[li]**2+precX2list[li]**2+precY2list[li]**2+precZ2list[li]**2)



            distancelist[li]=np.sqrt((X1list[li]-X2list[li])**2+(Y1list[li]-Y2list[li])**2+(Z1list[li]-Z2list[li])**2)




    loc['Spot_1_X']=X1list
    loc['Spot_1_Y']=Y1list
    loc['Spot_1_Z']=Z1list

    loc['Spot_2_X']=X2list
    loc['Spot_2_Y']=Y2list
    loc['Spot_2_Z']=Z2list

    loc['Precision_1_X']=precX1list
    loc['Precision_1_Y']=precY1list
    loc['Precision_1_Z']=precZ1list

    loc['Precision_2_X']=precX2list
    loc['Precision_2_Y']=precY2list
    loc['Precision_2_Z']=precZ2list

    loc['SigmaFit_1_X']=stdX1list
    loc['SigmaFit_1_Y']=stdY1list
    loc['SigmaFit_1_Z']=stdZ1list

    loc['SigmaFit_2_X']=stdX2list
    loc['SigmaFit_2_Y']=stdY2list
    loc['SigmaFit_2_Z']=stdZ2list


    loc['Photons_1']=photons1list
    loc['Photons_2']=photons2list
    loc['Background_1']=bckg1list
    loc['Background_2']=bckg2list

    loc['Distance']=distancelist



    loc['precision_Distance']=precisionDistlist

    loc['distanceToEdge_1 (pixels)']=distancePixelToEdge1list
    loc['distanceToEdge_2 (pixels)']=distancePixelToEdge2list

    loc['distanceToInit_1 (pixels)']=distancePixelToInit1list
    loc['distanceToInit_2 (pixels)']=distancePixelToInit2list

    loc['goodness_of_fit_1']=goodnessOfFit1list
    loc['goodness_of_fit_2']=goodnessOfFit2list

    loc.astype(str)
    loc.to_csv(locpathnew,sep=delimiter,index=False,quoting=csv.QUOTE_ALL,doublequote=True)
    print('file saved ',locpathnew)

def gaussian(x,amplitude,muX,muY,muZ,varX,varY,varZ,offset):#gaussian model 3d WITHOUT COVARIANCE !!!
    if amplitude<0:
        amplitude=1
    if offset<0:
        offset=.000001
    amp=(1/(((2*3.141592)**1.5)*((varX*varY*varZ)**(0.5))))*amplitude
    xx=np.divide(np.power(x[2,:]-muX,2) ,varX)
    yy=np.divide(np.power(x[1,:]-muY,2) ,varY)
    zz=np.divide(np.power(x[0,:]-muZ,2) ,varZ)
    g=(amp*np.exp(-0.5*(xx+yy+zz)))+offset

    return(g.ravel())

def gaussianLinBckg(x,amplitude,muX,muY,muZ,varX,varY,varZ,offset,interceptX,interceptY,interceptZ):#gaussian model 3d WITHOUT COVARIANCE !!!
    amp=(1/(((2*3.141592)**1.5)*((varX*varY*varZ)**(0.5))))*amplitude
    xx=np.divide(np.power(x[2,:]-muX,2) ,varX)
    yy=np.divide(np.power(x[1,:]-muY,2) ,varY)
    zz=np.divide(np.power(x[0,:]-muZ,2) ,varZ)
    g=(amp*np.exp(-0.5*(xx+yy+zz)))+offset+x[2,:]*interceptX+x[1,:]*interceptY+x[0,:]*interceptZ

    return(g.ravel())

def fitLocMLEGaussian(x0,y0,z0,xyz,data3D,shiftFit=2,sigma=None,photonoffset=None,withUncertainty=True,eps=.001,show=False):
    #if sigma is none, then it is fitted, otherwise, it is fixed






    #3D gaussian kernet without covariance


    #corner position of the crop
    shiftZYX=[np.max((int(z0+.5)-shiftFit,0)),np.max((int(y0+.5)-shiftFit,0)),np.max((int(x0+.5)-shiftFit,0))]
    crop=data3D[np.max((int(z0+.5)-shiftFit,0)):int(z0+.5)+shiftFit+1,np.max((int(y0+.5)-shiftFit,0)):int(y0+.5)+shiftFit+1,np.max((int(x0+.5)-shiftFit,0)):int(x0+.5)+shiftFit+1]



    data=crop.ravel()

    if len(data)<=0:
        if withUncertainty:
            return None,None,None,None,None,None,None,None,None
        else:
            return None,None,None,None,None,None,None,None


    muz=shiftFit+z0+.5-int(z0+.5)-.5
    muy=shiftFit+y0+.5-int(y0+.5)-.5
    mux=shiftFit+x0+.5-int(x0+.5)-.5

    sh=np.shape(crop)

    patchsizeZ=sh[0]
    patchsizeY=sh[1]
    patchsizeX=sh[2]

    locIsOnImageEdge=False


    #print('shift',sh[0],sh[1],sh[2],(shiftFit*2+1))
    if sh[0]!=shiftFit*2+1 or sh[1]!=shiftFit*2+1 or sh[2]!=shiftFit*2+1 :
        locIsOnImageEdge=True
        xyz=make_XYZ_axis(sh[0],sh[1],sh[2])

        #if crop is on negative edge of image --> shift center
        if int(z0+.5)-shiftFit-.5<0:
            muz+=int(z0+.5)-shiftFit-.5
            #print('muz=',muz)
        if int(y0)-shiftFit-.5<0:
            muy+=int(y0+.5)-shiftFit-.5
            #print('muy=',muy)
        if int(x0+.5)-shiftFit-.5<0:
            mux+=int(x0+.5)-shiftFit-.5
            #print('mux=',mux)





    if sigma is None:
        varz=1
        vary=1
        varx=1
    else:
        varz=sigma[2]**2
        vary=sigma[1]**2
        varx=sigma[0]**2

    if photonoffset is None:
        amplitude=None
        offset=None

    else:
        amplitude=photonoffset[0]
        offset=photonoffset[1]

    offsetinit=np.mean(crop)
    amplitudeinit=np.sum(crop[int(muz),:,:])-offsetinit*sh[2]*sh[1]


    gaussInit = gaussian(xyz, amplitudeinit,mux,muy,muz,varx,vary,varz,offsetinit)

    #print("XXX ",len(data),np.shape(data),np.shape(gaussInit),sh)
    serrInit=np.power(np.subtract(data,gaussInit),2)
    mseInit=np.mean(serrInit)
    #fitting of free polymers only:


    def loss_func_amplitude(x, xobs, yobs):

        # Evaluate the fit function with the current parameter estimates
        ynew = my_gaussian_amplitude(xobs, *x)
        yerr= np.sum(ynew-yobs*np.log(ynew))#MLE (POISSON)
        #yerr = np.sum((ynew - yobs) ** 2) #MSE

        return yerr

    # Define function
    def my_gaussian_amplitude(x,amplitude,offset):#gaussian model 3d WITHOUT COVARIANCE !!!
        if amplitude<0:
            amplitude=1
        if offset<0:
            offset=.000001
        amp=(1/(((2*3.141592)**1.5)*((varx*vary*varz)**(0.5))))*amplitude
        xx=np.divide(np.power(x[2,:]-mux,2) ,varx)
        yy=np.divide(np.power(x[1,:]-muy,2) ,vary)
        zz=np.divide(np.power(x[0,:]-muz,2) ,varz)
        g=(amp*np.exp(-0.5*(xx+yy+zz)))+offset

        return(g.ravel())







    def loss_func_amplitude_position(x, xobs, yobs):

        # Evaluate the fit function with the current parameter estimates
        ynew = my_gaussian_amplitude_position(xobs, *x)
        yerr= np.sum(ynew-yobs*np.log(ynew))#MLE (POISSON)
        #yerr = np.sum((ynew - yobs) ** 2) #MSE

        return yerr

    # Define function
    def my_gaussian_amplitude_position(x,amplitude,muX,muY,muZ,offset):#gaussian model 3d WITHOUT COVARIANCE !!!
        if amplitude<0:
            amplitude=1
        if offset<0:
            offset=.000001
        amp=(1/(((2*3.141592)**1.5)*((varx*vary*varz)**(0.5))))*amplitude
        xx=np.divide(np.power(x[2,:]-muX,2) ,varx)
        yy=np.divide(np.power(x[1,:]-muY,2) ,vary)
        zz=np.divide(np.power(x[0,:]-muZ,2) ,varz)
        g=(amp*np.exp(-0.5*(xx+yy+zz)))+offset

        return(g.ravel())








    def loss_func_position(x, xobs, yobs):

        # Evaluate the fit function with the current parameter estimates
        ynew = my_gaussian_position(xobs, *x)
        yerr= np.sum(ynew-yobs*np.log(ynew))#MLE (POISSON)
        #yerr = np.sum((ynew - yobs) ** 2) #MSE

        return yerr

    # Define function
    def my_gaussian_position(x,muX,muY,muZ):#gaussian model 3d WITHOUT COVARIANCE !!!

        amp=(1/(((2*3.141592)**1.5)*((varx*vary*varz)**(0.5))))*amplitude
        xx=np.divide(np.power(x[2,:]-muX,2) ,varx)
        yy=np.divide(np.power(x[1,:]-muY,2) ,vary)
        zz=np.divide(np.power(x[0,:]-muZ,2) ,varz)
        g=(amp*np.exp(-0.5*(xx+yy+zz)))+offset

        return(g.ravel())










    def loss_func_position_var(x, xobs, yobs):

        # Evaluate the fit function with the current parameter estimates
        ynew = my_gaussian_position_var(xobs, *x)
        yerr= np.sum(ynew-yobs*np.log(ynew))#MLE (POISSON)
        #yerr = np.sum((ynew - yobs) ** 2) #MSE

        return yerr

    # Define function
    def my_gaussian_position_var(x,muX,muY,muZ,varx,vary,varz):#gaussian model 3d WITHOUT COVARIANCE !!!
        if varx<=0:
            varx=.000001
        if vary<=0:
            vary=.000001
        if varz<=0:
            varz=.000001
        amp=(1/(((2*3.141592)**1.5)*((varx*vary*varz)**(0.5))))*amplitude
        xx=np.divide(np.power(x[2,:]-muX,2) ,varx)
        yy=np.divide(np.power(x[1,:]-muY,2) ,vary)
        zz=np.divide(np.power(x[0,:]-muZ,2) ,varz)
        g=(amp*np.exp(-0.5*(xx+yy+zz)))+offset

        return(g.ravel())









    def loss_func_amplitude_position_var(x, xobs, yobs):

        # Evaluate the fit function with the current parameter estimates
        ynew = my_gaussian_amplitude_position_var(xobs, *x)
        yerr= np.sum(ynew-yobs*np.log(ynew))#MLE (POISSON)
        #yerr = np.sum((ynew - yobs) ** 2) #MSE

        return yerr

    # Define function
    def my_gaussian_amplitude_position_var(x,amplitude,muX,muY,muZ,varx,vary,varz,offset):#gaussian model 3d WITHOUT COVARIANCE !!!
        if amplitude<0:
            amplitude=1
        if offset<0:
            offset=.000001
        if varx<=0:
            varx=.000001
        if vary<=0:
            vary=.000001
        if varz<=0:
            varz=.000001
        amp=(1/(((2*3.141592)**1.5)*((varx*vary*varz)**(0.5))))*amplitude
        xx=np.divide(np.power(x[2,:]-muX,2) ,varx)
        yy=np.divide(np.power(x[1,:]-muY,2) ,vary)
        zz=np.divide(np.power(x[0,:]-muZ,2) ,varz)
        g=(amp*np.exp(-0.5*(xx+yy+zz)))+offset

        return(g.ravel())



    #print('init free',init_vals)
    try:
        if photonoffset is None:

            #first: estimate photons in spot and bckg only:
            init_vals = [amplitudeinit,offsetinit]
            #best_vals, covar = curve_fit(lambda xyz,aa,oo: gaussian(xyz,aa,mux,muy,muz,varx,vary,varz,oo), xyz,data,p0=init_vals,maxfev = 500)
            bnds = ((.000001, None), (0, None))
            res = minimize(loss_func_amplitude, x0=init_vals, args=(xyz, data),bounds=bnds)
            best_vals=res.x



            [amplitude_init,offset_init]=best_vals[:]



            #print('amplitude',res.success , amplitude_init,mux,muy,muz,varx,vary,varz,offset_init)

            if sigma is None:
                init_vals = [amplitude_init,mux,muy,muz,varx,vary,varz,offset_init]     # for [amp, cen, wid]
                #best_vals, covar = curve_fit(gaussian, xyz, data, p0=init_vals)
                bnds = ((.000001, None), (-.5, shiftFit*2+.5), (-.5, shiftFit*2+.5), (-.5, shiftFit*2+.5), (.000001, None), (.000001, None), (.000001, None), (0, None))
                res = minimize(loss_func_amplitude_position_var, x0=init_vals, args=(xyz, data),bounds=bnds,method='Nelder-Mead')
                #print('to readd')
                if res.success == False:
                    if withUncertainty:
                        return None,None,None,None,None,None,None,None,None
                    else:
                        return None,None,None,None,None,None,None,None
                best_vals=res.x

                [amplitude_fitted,mux_fitted,muy_fitted,muz_fitted,varx_fitted,vary_fitted,varz_fitted,offset_fitted]=best_vals[:]

            else:
                init_vals = [amplitude_init,mux,muy,muz,offset_init]     # for [amp, cen, wid]
                #best_vals, covar = curve_fit(lambda xyz,aa,xx,yy,zz,oo: gaussian(xyz,aa,xx,yy,zz,varx,vary,varz,oo), xyz,data,p0=init_vals,maxfev = 500)
                bnds = ((.000001, None), (-.5, shiftFit*2+.5), (-.5, shiftFit*2+.5), (-.5, shiftFit*2+.5),  (0, None))
                res = minimize(loss_func_amplitude_position, x0=init_vals, args=(xyz, data),bounds=bnds,method='Nelder-Mead')

                if res.success == False:
                    if withUncertainty:
                        return None,None,None,None,None,None,None,None,None
                    else:
                        return None,None,None,None,None,None,None,None
                best_vals=res.x


                [amplitude_fitted,mux_fitted,muy_fitted,muz_fitted,offset_fitted]=best_vals[:]


                varx_fitted=varx
                vary_fitted=vary
                varz_fitted=varz

        else:



            #print('amplitude',res.success , amplitude_init,mux,muy,muz,varx,vary,varz,offset_init)

            if sigma is None:
                init_vals = [mux,muy,muz,varx,vary,varz]     # for [amp, cen, wid]
                #best_vals, covar = curve_fit(gaussian, xyz, data, p0=init_vals)
                bnds = ((-.5, shiftFit*2+.5), (-.5, shiftFit*2+.5), (-.5, shiftFit*2+.5), (.000001, None), (.000001, None), (.000001, None))
                res = minimize(loss_func_position_var, x0=init_vals, args=(xyz, data),bounds=bnds,method='Nelder-Mead')
                #print('to readd')
                if res.success == False:
                    if withUncertainty:
                        return None,None,None,None,None,None,None,None,None
                    else:
                        return None,None,None,None,None,None,None,None
                best_vals=res.x

                [mux_fitted,muy_fitted,muz_fitted,varx_fitted,vary_fitted,varz_fitted]=best_vals[:]
                offset_fitted=offset
                amplitude_fitted=amplitude
            else:
                init_vals = [mux,muy,muz]     # for [amp, cen, wid]
                #best_vals, covar = curve_fit(lambda xyz,aa,xx,yy,zz,oo: gaussian(xyz,aa,xx,yy,zz,varx,vary,varz,oo), xyz,data,p0=init_vals,maxfev = 500)
                bnds = ((-.5, shiftFit*2+.5), (-.5, shiftFit*2+.5), (-.5, shiftFit*2+.5))
                res = minimize(loss_func_position, x0=init_vals, args=(xyz, data),bounds=bnds,method='Nelder-Mead')

                if res.success == False:
                    if withUncertainty:
                        return None,None,None,None,None,None,None,None,None
                    else:
                        return None,None,None,None,None,None,None,None
                best_vals=res.x


                [mux_fitted,muy_fitted,muz_fitted]=best_vals[:]


                varx_fitted=varx
                vary_fitted=vary
                varz_fitted=varz
                offset_fitted=offset
                amplitude_fitted=amplitude

    except RuntimeError:
        print('RuntimeError is raised')
        if withUncertainty:
            return None,None,None,None,None,None,None,None,None
        else:
            return None,None,None,None,None,None,None,None




    gaussFitted = gaussian(xyz, amplitude_fitted,mux_fitted,muy_fitted,muz_fitted,varx_fitted,vary_fitted,varz_fitted,offset_fitted)



    serr=np.power(np.subtract(data,gaussFitted),2)
    mse=np.mean(serr)

    #print('init_vals',init_vals,mseInit)
    #print('best_vals',best_vals,mse)

    if show:
        gaussFittedForPlot=gaussFitted.reshape((patchsizeZ,patchsizeY,patchsizeX))

        gaussInit = gaussian(xyz, amplitude_fitted,mux,muy,muz,varx,vary,varz,offset_fitted)
        gaussInitForPlot=gaussInit.reshape((patchsizeZ,patchsizeY,patchsizeX))

        rangeplot=[np.min(crop),np.max(crop)]

        fig=plt.figure()
        for k,im in enumerate(crop):
            plt.subplot(1,len(crop),k+1)
            plt.tick_params(left = False, right = False , labelleft = False , labelbottom = False, bottom = False)
            plt.imshow(crop[k], vmin=rangeplot[0], vmax=rangeplot[1])



        fig=plt.figure()
        for k,im in enumerate(gaussInitForPlot):
            plt.subplot(1,len(gaussInitForPlot),k+1)
            plt.tick_params(left = False, right = False , labelleft = False , labelbottom = False, bottom = False)
            plt.imshow(gaussInitForPlot[k], vmin=rangeplot[0], vmax=rangeplot[1])


        fig=plt.figure()
        for k,im in enumerate(gaussFittedForPlot):
            plt.subplot(1,len(gaussFittedForPlot),k+1)
            plt.tick_params(left = False, right = False , labelleft = False , labelbottom = False, bottom = False)
            plt.imshow(gaussFittedForPlot[k], vmin=rangeplot[0], vmax=rangeplot[1])


    res=[mux_fitted+shiftZYX[2],muy_fitted+shiftZYX[1],muz_fitted+shiftZYX[0]]

    distancePixelToEdge=shiftFit*2#init to high value
    distancePixelToEdge=np.min((distancePixelToEdge,np.abs(mux_fitted+.5-0)))
    distancePixelToEdge=np.min((distancePixelToEdge,np.abs(mux_fitted+.5-sh[2])))
    distancePixelToEdge=np.min((distancePixelToEdge,np.abs(muy_fitted+.5-0)))
    distancePixelToEdge=np.min((distancePixelToEdge,np.abs(muy_fitted+.5-sh[1])))
    distancePixelToEdge=np.min((distancePixelToEdge,np.abs(muz_fitted+.5-0)))
    distancePixelToEdge=np.min((distancePixelToEdge,np.abs(muz_fitted+.5-sh[0])))



    distancePixelToInit=np.sqrt((z0-res[2])**2+(y0-res[1])**2+(x0-res[0])**2)

    #goodnessOfFit=np.nanmean(np.divide(np.power(np.subtract(data,gaussFitted),2),gaussFitted))##############################################################""
    gaussFitted3D=gaussFitted.reshape((patchsizeZ,patchsizeY,patchsizeX))
    pzint=int(muz_fitted+.5)
    pyint=int(muy_fitted+.5)
    pxint=int(mux_fitted+.5)
    goodnessOfFit=np.nanmean(np.divide(np.power(np.subtract(crop[np.max((pzint-1,0)):pzint+2,np.max((pyint-1,0)):pyint+2,np.max((pxint-1,0)):pxint+2],gaussFitted3D[np.max((pzint-1,0)):pzint+2,np.max((pyint-1,0)):pyint+2,np.max((pxint-1,0)):pxint+2]),2),gaussFitted3D[np.max((pzint-1,0)):pzint+2,np.max((pyint-1,0)):pyint+2,np.max((pxint-1,0)):pxint+2]))##############################################################""


    if withUncertainty:
        #print('c',coef)
        #compute parabola coefs for slightly shifted parabola:


        modelXshift=gaussian(xyz, amplitude_fitted,mux_fitted+eps,muy_fitted,muz_fitted,varx_fitted,vary_fitted,varz_fitted,offset_fitted)
        modelYshift=gaussian(xyz, amplitude_fitted,mux_fitted,muy_fitted+eps,muz_fitted,varx_fitted,vary_fitted,varz_fitted,offset_fitted)
        modelZshift=gaussian(xyz, amplitude_fitted,mux_fitted,muy_fitted,muz_fitted+eps,varx_fitted,vary_fitted,varz_fitted,offset_fitted)
        modelAmpshift=gaussian(xyz, amplitude_fitted+10,mux_fitted,muy_fitted,muz_fitted,varx_fitted,vary_fitted,varz_fitted,offset_fitted)
        modelBckgshift=gaussian(xyz, amplitude_fitted,mux_fitted,muy_fitted,muz_fitted,varx_fitted,vary_fitted,varz_fitted,offset_fitted+1)



        serr_partial=[]
        model_partial=[]
        #compute derivative of MSE according to x, y and z:


        model_partial.append((modelXshift-gaussFitted)/eps)
        serr_partial.append(((np.power(np.subtract(modelXshift,gaussFitted),2))-mse)/eps)

        model_partial.append((modelYshift-gaussFitted)/eps)
        serr_partial.append(((np.power(np.subtract(modelYshift,gaussFitted),2))-mse)/eps)

        model_partial.append((modelZshift-gaussFitted)/eps)
        serr_partial.append(((np.power(np.subtract(modelZshift,gaussFitted),2))-mse)/eps)

        model_partial.append((modelAmpshift-gaussFitted)/eps)
        serr_partial.append(((np.power(np.subtract(modelAmpshift,gaussFitted),2))-mse)/10)

        model_partial.append((modelBckgshift-gaussFitted)/eps)
        serr_partial.append(((np.power(np.subtract(modelBckgshift,gaussFitted),2))-mse)/1)

        fisher=np.zeros((5,5))
        #print('model_partial',model_partial)
        #print('serr',serr)
        for i in range(5):
            for j in range(5):
                #fisher[i,j]=np.sum(serr_partial[i]*serr_partial[j])
                #+eps to avoid division 0
                #fisher[i,j]=np.sum(model_partial[i]*model_partial[j]/(mse+eps))
                fisher[i,j]=np.sum(model_partial[i]*model_partial[j]/(gaussFitted))
        #print(fisher)

        thePrecision=None

        try:
            precision=np.linalg.pinv(fisher)
            thePrecision=[np.sqrt(precision[0,0]),np.sqrt(precision[1,1]),np.sqrt(precision[2,2])]

        except Exception:
            thePrecision=None
        return res,mse,thePrecision,[np.sqrt(varx_fitted),np.sqrt(vary_fitted),np.sqrt(varz_fitted)],amplitude_fitted,offset_fitted,distancePixelToInit,distancePixelToEdge,goodnessOfFit


    else:
        return(res,mse,[np.sqrt(varx_fitted),np.sqrt(vary_fitted),np.sqrt(varz_fitted)],amplitude_fitted,offset_fitted,distancePixelToInit,distancePixelToEdge,goodnessOfFit)




#convert X,Y,Z to table to be usable by ChromaticAbberations
def vect2mat(x,y,z,f):
    length=len(x)
    table=[]
    frame=[]
    track=[]
    for i in range(length):
        v=[x[i],y[i],z[i]]
        if i==0:
            table=v
            frame=f[i]
            track=[i]*len(x[i])
        else:
            table=np.concatenate((table,v),axis=1)
            frame=np.concatenate((frame,f[i]))
            track=np.concatenate((track,[i]*len(x[i])))
    table=np.transpose(table)
    return table,frame,track

def mat2vect(table,frame,track):
    x=[]
    y=[]
    z=[]
    f=[]
    for i in np.unique(track):
        frame=np.array(frame)
        index=np.array(track==i)
        xyz=table[index]
        x.append(xyz[:,0])
        y.append(xyz[:,1])
        z.append(xyz[:,2])
        print(len(frame),len(index))
        f.append(frame[index])
    return x,y,z,f



############################# convert sigma in pixel unit


if fixedStdRed is not None:
    fixedStdRed=np.array(fixedStdRed)
    fixedStdRed=fixedStdRed/voxelsize
if fixedStdGreen is not None:
    fixedStdGreen=np.array(fixedStdGreen)
    fixedStdGreen=fixedStdGreen/voxelsize




############################# open image:

data=sk.io.imread(imagepath)




imscale=1
sizeX=np.shape(data)[-1]
sizeY=np.shape(data)[-2]
offsetX=0
offsetY=0





datared=np.array(data[:,:,0,offsetY:offsetY+sizeY,offsetX:offsetX+sizeX],dtype='float32')
datagreen=np.array(data[:,:,1,offsetY:offsetY+sizeY,offsetX:offsetX+sizeX],dtype='float32')
data=None


datared+=photon_intercept
datared*=photon_rescale
datagreen+=photon_intercept
datagreen*=photon_rescale
datared = np.where(datared<0, 0.00000001, datared)
datagreen = np.where(datagreen<0, 0.00000001, datagreen)



############################# open CSV file:

Frame,X1,Y1,Z1,X2,Y2,Z2,Track=readLoc(locpath,voxelsize,imscale=imscale,offsetX=offsetX,offsetY=offsetY,delimiter=delimiter)


'''for i,x1 in enumerate(X1):
    X1[i]=np.add(X1[i],0.5)
for i,y1 in enumerate(Y1):
    Y1[i]=np.add(Y1[i],0.5)
for i,z1 in enumerate(Z1):
    Z1[i]=np.add(Z1[i],0.5)
for i,x2 in enumerate(X2):
    X2[i]=np.add(X2[i],0.5)
for i,y2 in enumerate(Y2):
    Y2[i]=np.add(Y2[i],0.5)
for i,z2 in enumerate(Z2):
    Z2[i]=np.add(Z2[i],0.5)'''




############################# Reverse chromatic correction ####################
sh=np.shape(datagreen)
w=sh[-1]
h=sh[-2]
d=sh[-3]


tableGreen,frameGreen,trackGreen=vect2mat(X2,Y2,Z2,Frame)

ca=ChromaticAberrations(chromagnonpath,driftpath,chromagnonPixSizeXYZ=[-1,-1,-1])

ca.setTable(tableGreen,frameGreen,w,h,d)

ca.trackmate_applyDriftCorrection(-1)

ca.trackmate_ReverseChromaCorrection()

ca.trackmate_applyDriftCorrection(+1)

tableGreen,frameGreen=ca.getTable()

X2,Y2,Z2,Frameb=mat2vect(tableGreen,frameGreen,trackGreen)

##############################################################################





############################## FIT SPOTS with Gaussian model #################


show=False

patchsize=shiftFit*2+1

xyz=make_XYZ_axis(patchsize)

prec_red=[]
prec_green=[]
var_red=[]
var_green=[]



#fixedStd=None

X1new=[]
X2new=[]
Y1new=[]
Y2new=[]
Z1new=[]
Z2new=[]

precX1=[]
precX2=[]
precY1=[]
precY2=[]
precZ1=[]
precZ2=[]

stdX1=[]
stdY1=[]
stdZ1=[]
stdX2=[]
stdY2=[]
stdZ2=[]


photons1=[]
photons2=[]
background1=[]
background2=[]

distancePixelToInit1=[]
distancePixelToEdge1=[]
distancePixelToInit2=[]
distancePixelToEdge2=[]
goodnessOfFit1=[]
goodnessOfFit2=[]

count=0
maxcountthreshold=1





for k in range(len(Frame)):
    #print('track processed=',k,Track[k][0])

    print('track processed=',Track[k][0])

    '''miniX=np.min((np.min(X1[k]),np.min(X2[k])))
    miniY=np.min((np.min(Y1[k]),np.min(Y2[k])))
    maxiX=np.max((np.min(X1[k]),np.max(X2[k])))
    maxiY=np.max((np.min(Y1[k]),np.max(Y2[k])))
    shift=30
    shiftX=int(miniX-shift//2)
    shiftY=int(miniY-shift//2)
    dataredsub=datared[:,:,shiftY:int(maxiY+shift//2),shiftX:int(maxiX+shift//2)]
    datagreensub=datagreen[:,:,shiftY:int(maxiY+shift//2),shiftX:int(maxiX+shift//2)]
    '''

    shred=np.shape(datared)
    shgreen=np.shape(datagreen)
    if (shred!=shgreen):
        print('ERROR image dimension mismatch')
    miniX=0
    miniY=0
    maxiX=shred[-1]
    maxiY=shred[-2]
    shift=0
    shiftX=int(miniX-shift//2)
    shiftY=int(miniY-shift//2)


    dataredsub=datared
    datagreensub=datagreen

    X1newsub=[]
    X2newsub=[]
    Y1newsub=[]
    Y2newsub=[]
    Z1newsub=[]
    Z2newsub=[]

    stdX1sub=[]
    stdY1sub=[]
    stdZ1sub=[]
    stdX2sub=[]
    stdY2sub=[]
    stdZ2sub=[]

    precX1sub=[]
    precX2sub=[]
    precY1sub=[]
    precY2sub=[]
    precZ1sub=[]
    precZ2sub=[]

    photons1sub=[]
    photons2sub=[]
    background1sub=[]
    background2sub=[]

    distancePixelToInit1sub=[]
    distancePixelToEdge1sub=[]
    distancePixelToInit2sub=[]
    distancePixelToEdge2sub=[]

    goodnessOfFit1sub=[]
    goodnessOfFit2sub=[]

    indexFrame=0
    #print(Frame[k][indexFrame],X1[k][indexFrame],Y1[k][indexFrame])
    for indexFrame in range(len(Frame[k])):



        fitRes1,mse,prec1,stand,photon,bckg,distancePixelToInit,distancePixelToEdge,goodnessOfFit=fitLocMLEGaussian(X1[k][indexFrame],Y1[k][indexFrame],Z1[k][indexFrame],xyz,dataredsub[Frame[k][indexFrame]],shiftFit=shiftFit,sigma=fixedStdRed,show=show)



        if fitRes1 is not None:

            X1newsub.append(fitRes1[0]+shiftX)
            Y1newsub.append(fitRes1[1]+shiftY)
            Z1newsub.append(fitRes1[2])
            photons1sub.append(photon)
            background1sub.append(bckg)
            distancePixelToInit1sub.append(distancePixelToInit)
            distancePixelToEdge1sub.append(distancePixelToEdge)
            goodnessOfFit1sub.append(goodnessOfFit)

        else:
            X1newsub.append(X1[k][indexFrame])
            Y1newsub.append(Y1[k][indexFrame])
            Z1newsub.append(Z1[k][indexFrame])
            photons1sub.append(float('nan'))
            background1sub.append(float('nan'))
            distancePixelToInit1sub.append(0)
            distancePixelToEdge1sub.append(float('nan'))
            goodnessOfFit1sub.append(float('nan'))

        if prec1 is not None:
            precX1sub.append(prec1[0])
            precY1sub.append(prec1[1])
            precZ1sub.append(prec1[2])
        else:
            precX1sub.append(float('nan'))
            precY1sub.append(float('nan'))
            precZ1sub.append(float('nan'))

        if stand is not None:
            stdX1sub.append(stand[0])
            stdY1sub.append(stand[1])
            stdZ1sub.append(stand[2])
        else:
            stdX1sub.append(float('nan'))
            stdY1sub.append(float('nan'))
            stdZ1sub.append(float('nan'))


        #print('red ok',show)

        fitRes2,mse,prec2,stand,photon,bckg,distancePixelToInit,distancePixelToEdge,goodnessOfFit=fitLocMLEGaussian(X2[k][indexFrame],Y2[k][indexFrame],Z2[k][indexFrame],xyz,datagreensub[Frame[k][indexFrame]],shiftFit=shiftFit,sigma=fixedStdGreen,show=show)


        if fitRes2 is not None:
            X2newsub.append(fitRes2[0]+shiftX)
            Y2newsub.append(fitRes2[1]+shiftY)
            Z2newsub.append(fitRes2[2])
            photons2sub.append(photon)
            background2sub.append(bckg)
            distancePixelToInit2sub.append(distancePixelToInit)
            distancePixelToEdge2sub.append(distancePixelToEdge)
            goodnessOfFit2sub.append(goodnessOfFit)
        else:
            X2newsub.append(X2[k][indexFrame])
            Y2newsub.append(Y2[k][indexFrame])
            Z2newsub.append(Z2[k][indexFrame])
            photons2sub.append(float('nan'))
            background2sub.append(float('nan'))
            distancePixelToInit2sub.append(0)
            distancePixelToEdge2sub.append(float('nan'))
            goodnessOfFit2sub.append(float('nan'))

        #print(prec2)
        if prec2 is not None:
            precX2sub.append(prec2[0])
            precY2sub.append(prec2[1])
            precZ2sub.append(prec2[2])
        else:
            precX2sub.append(float('nan'))
            precY2sub.append(float('nan'))
            precZ2sub.append(float('nan'))

        if stand is not None:
            stdX2sub.append(stand[0])
            stdY2sub.append(stand[1])
            stdZ2sub.append(stand[2])
        else:
            stdX2sub.append(float('nan'))
            stdY2sub.append(float('nan'))
            stdZ2sub.append(float('nan'))

        #print('green ok ',show)


        count+=1


    if refit_with_averaged_intensity:
        try:
            ########################### START REFIT WITH FIXED INTENSITY #######################
            nonnan=np.invert(np.isnan(photons1sub))

            photons1_tmp=np.array(photons1sub)
            bckg1_tmp=np.array(background1sub)

            photons1_tmp=photons1_tmp[nonnan]
            bckg1_tmp=bckg1_tmp[nonnan]

            f_tmp=Frame[k]

            f_tmp=np.array(f_tmp)
            f_tmp_nonnan=f_tmp[nonnan]

            kr = KernelReg(photons1_tmp,f_tmp_nonnan,'c')
            photons1_fixed, y_std = kr.fit(f_tmp)

            kr = KernelReg(bckg1_tmp,f_tmp_nonnan,'c')
            bckg1_fixed, y_std = kr.fit(f_tmp)



            '''plt.figure()
            plt.plot(f_tmp, photons1sub)
            plt.plot(f_tmp, photons1_fixed)
            plt.show()



            plt.figure()
            plt.plot(f_tmp, background1sub)
            plt.plot(f_tmp, bckg1_fixed)
            plt.show()'''


            nonnan=np.invert(np.isnan(photons2sub))

            photons2_tmp=np.array(photons2sub)
            bckg2_tmp=np.array(background2sub)

            photons2_tmp=photons2_tmp[nonnan]
            bckg2_tmp=bckg2_tmp[nonnan]

            f_tmp=Frame[k]

            f_tmp=np.array(f_tmp)
            f_tmp_nonnan=f_tmp[nonnan]

            kr = KernelReg(photons2_tmp,f_tmp_nonnan,'c')
            photons2_fixed, y_std = kr.fit(f_tmp)

            kr = KernelReg(bckg2_tmp,f_tmp_nonnan,'c')
            bckg2_fixed, y_std = kr.fit(f_tmp)


            '''plt.figure()
            plt.plot(f_tmp, photons2sub)
            plt.plot(f_tmp, photons2_fixed)
            plt.show()



            plt.figure()
            plt.plot(f_tmp, background2sub)
            plt.plot(f_tmp, bckg2_fixed)
            plt.show()'''






            shred=np.shape(datared)
            shgreen=np.shape(datagreen)
            if (shred!=shgreen):
                print('ERROR image dimension mismatch')
            miniX=0
            miniY=0
            maxiX=shred[-1]
            maxiY=shred[-2]
            shift=0
            shiftX=int(miniX-shift//2)
            shiftY=int(miniY-shift//2)


            dataredsub=datared
            datagreensub=datagreen

            X1newsub=[]
            X2newsub=[]
            Y1newsub=[]
            Y2newsub=[]
            Z1newsub=[]
            Z2newsub=[]

            stdX1sub=[]
            stdY1sub=[]
            stdZ1sub=[]
            stdX2sub=[]
            stdY2sub=[]
            stdZ2sub=[]

            precX1sub=[]
            precX2sub=[]
            precY1sub=[]
            precY2sub=[]
            precZ1sub=[]
            precZ2sub=[]

            photons1sub=[]
            photons2sub=[]
            background1sub=[]
            background2sub=[]

            distancePixelToInit1sub=[]
            distancePixelToEdge1sub=[]
            distancePixelToInit2sub=[]
            distancePixelToEdge2sub=[]

            indexFrame=0
            #print(Frame[k][indexFrame],X1[k][indexFrame],Y1[k][indexFrame])
            for indexFrame in range(len(Frame[k])):



                fitRes1,mse,prec1,stand,photon,bckg,distancePixelToInit,distancePixelToEdge,goodnessOfFit=fitLocMLEGaussian(X1[k][indexFrame],Y1[k][indexFrame],Z1[k][indexFrame],xyz,dataredsub[Frame[k][indexFrame]],shiftFit=shiftFit,sigma=fixedStdRed,photonoffset=[photons1_fixed[indexFrame],bckg1_fixed[indexFrame]],show=show)

                #fitRes1,mse,prec1,stand,photon,bckg,distancePixelToInit,distancePixelToEdge=fitLocMLEGaussian(X1[k][indexFrame],Y1[k][indexFrame],Z1[k][indexFrame],xyz,dataredsub[Frame[k][indexFrame]],shiftFit=shiftFit,sigma=fixedStdRed,photonoffset=[photon,bckg],show=show)



                if fitRes1 is not None:

                    X1newsub.append(fitRes1[0]+shiftX)
                    Y1newsub.append(fitRes1[1]+shiftY)
                    Z1newsub.append(fitRes1[2])
                    photons1sub.append(photon)
                    background1sub.append(bckg)
                    distancePixelToInit1sub.append(distancePixelToInit)
                    distancePixelToEdge1sub.append(distancePixelToEdge)
                    goodnessOfFit1sub.append(goodnessOfFit)
                else:
                    X1newsub.append(X1[k][indexFrame])
                    Y1newsub.append(Y1[k][indexFrame])
                    Z1newsub.append(Z1[k][indexFrame])
                    photons1sub.append(float('nan'))
                    background1sub.append(float('nan'))
                    distancePixelToInit1sub.append(0)
                    distancePixelToEdge1sub.append(float('nan'))
                    goodnessOfFit1sub.append(float('nan'))

                if prec1 is not None:
                    precX1sub.append(prec1[0])
                    precY1sub.append(prec1[1])
                    precZ1sub.append(prec1[2])
                else:
                    precX1sub.append(float('nan'))
                    precY1sub.append(float('nan'))
                    precZ1sub.append(float('nan'))

                if stand is not None:
                    stdX1sub.append(stand[0])
                    stdY1sub.append(stand[1])
                    stdZ1sub.append(stand[2])
                else:
                    stdX1sub.append(float('nan'))
                    stdY1sub.append(float('nan'))
                    stdZ1sub.append(float('nan'))






                fitRes2,mse,prec2,stand,photon,bckg,distancePixelToInit,distancePixelToEdge,goodnessOfFit=fitLocMLEGaussian(X2[k][indexFrame],Y2[k][indexFrame],Z2[k][indexFrame],xyz,datagreensub[Frame[k][indexFrame]],shiftFit=shiftFit,sigma=fixedStdGreen,photonoffset=[photons2_fixed[indexFrame],bckg2_fixed[indexFrame]],show=show)

                if fitRes2 is not None:
                    X2newsub.append(fitRes2[0]+shiftX)
                    Y2newsub.append(fitRes2[1]+shiftY)
                    Z2newsub.append(fitRes2[2])
                    photons2sub.append(photon)
                    background2sub.append(bckg)
                    distancePixelToInit2sub.append(distancePixelToInit)
                    distancePixelToEdge2sub.append(distancePixelToEdge)
                    goodnessOfFit2sub.append(goodnessOfFit)
                else:
                    X2newsub.append(X2[k][indexFrame])
                    Y2newsub.append(Y2[k][indexFrame])
                    Z2newsub.append(Z2[k][indexFrame])
                    photons2sub.append(float('nan'))
                    background2sub.append(float('nan'))
                    distancePixelToInit2sub.append(0)
                    distancePixelToEdge2sub.append(float('nan'))
                    goodnessOfFit2sub.append(float('nan'))

                #print(prec2)
                if prec2 is not None:
                    precX2sub.append(prec2[0])
                    precY2sub.append(prec2[1])
                    precZ2sub.append(prec2[2])
                else:
                    precX2sub.append(float('nan'))
                    precY2sub.append(float('nan'))
                    precZ2sub.append(float('nan'))

                if stand is not None:
                    stdX2sub.append(stand[0])
                    stdY2sub.append(stand[1])
                    stdZ2sub.append(stand[2])
                else:
                    stdX2sub.append(float('nan'))
                    stdY2sub.append(float('nan'))
                    stdZ2sub.append(float('nan'))
        except:
            print("error relocalization ",photons1sub)
            print("error relocalization ",photons2sub)

    ########################### END REFIT WITH FIXED INTENSITY #######################


    #print(precX1sub,precY1sub,precZ1sub)
    #print(precX2sub,precY2sub,precZ2sub)

    X1new.append(X1newsub)
    Y1new.append(Y1newsub)
    Z1new.append(Z1newsub)
    X2new.append(X2newsub)
    Y2new.append(Y2newsub)
    Z2new.append(Z2newsub)

    goodnessOfFit1.append(goodnessOfFit1sub)
    goodnessOfFit2.append(goodnessOfFit2sub)

    precX1.append(precX1sub)
    precY1.append(precY1sub)
    precZ1.append(precZ1sub)
    precX2.append(precX2sub)
    precY2.append(precY2sub)
    precZ2.append(precZ2sub)

    stdX1.append(stdX1sub)
    stdY1.append(stdY1sub)
    stdZ1.append(stdZ1sub)
    stdX2.append(stdX2sub)
    stdY2.append(stdY2sub)
    stdZ2.append(stdZ2sub)

    photons1.append(photons1sub)
    photons2.append(photons2sub)
    background1.append(background1sub)
    background2.append(background2sub)

    distancePixelToInit1.append(distancePixelToInit1sub)
    distancePixelToEdge1.append(distancePixelToEdge1sub)
    distancePixelToInit2.append(distancePixelToInit2sub)
    distancePixelToEdge2.append(distancePixelToEdge2sub)


##############################################################################

##############################################################################




'''for i,x1 in enumerate(X1new):
    X1new[i]=np.add(X1new[i],-0.5)
for i,y1 in enumerate(Y1new):
    Y1new[i]=np.add(Y1new[i],-0.5)
for i,z1 in enumerate(Z1new):
    Z1new[i]=np.add(Z1new[i],-0.5)

for i,x2 in enumerate(X2new):
    X2new[i]=np.add(X2new[i],-0.5)
for i,y2 in enumerate(Y2new):
    Y2new[i]=np.add(Y2new[i],-0.5)
for i,z2 in enumerate(Z2new):
    Z2new[i]=np.add(Z2new[i],-0.5)'''

############################# re-apply  chromatic correction #################
tableGreen,frameGreen,trackGreen=vect2mat(X2new,Y2new,Z2new,Frame)

ca=ChromaticAberrations(chromagnonpath,driftpath,chromagnonPixSizeXYZ=[-1,-1,-1])

ca.setTable(tableGreen,frameGreen,w,h,d)

ca.trackmate_applyDriftCorrection(-1)

ca.trackmate_ChromaCorrection()

ca.trackmate_applyDriftCorrection(+1)

tableGreen,frameGreen=ca.getTable()

X2new,Y2new,Z2new,Frameb=mat2vect(tableGreen,frameGreen,trackGreen)

##############################################################################






if refit_with_averaged_intensity:
    locpathnew=locpath[:-4]+'_GaussianFit_refit.csv'

    if fixedStdRed is None and fixedStdGreen is None :
        locpathnew=locpath[:-4]+'_GaussianFit_StdFitted_refit.csv'
    elif fixedStdRed is None:
        locpathnew=locpath[:-4]+'_GaussianFit_StdRedFitted_refit.csv'
    elif fixedStdGreen is None:
        locpathnew=locpath[:-4]+'_GaussianFit_StdGreenFitted_refit.csv'
else:
    locpathnew=locpath[:-4]+'_GaussianFit.csv'

    if fixedStdRed is None and fixedStdGreen is None :
        locpathnew=locpath[:-4]+'_GaussianFit_StdFitted.csv'
    elif fixedStdRed is None:
        locpathnew=locpath[:-4]+'_GaussianFit_StdRedFitted.csv'
    elif fixedStdGreen is None:
        locpathnew=locpath[:-4]+'_GaussianFit_StdGreenFitted.csv'


saveLoc(locpath,locpathnew,voxelsize,X1new,Y1new,Z1new,precX1,precY1,precZ1,stdX1,stdY1,stdZ1,X2new,Y2new,Z2new,precX2,precY2,precZ2,stdX2,stdY2,stdZ2,photons1,photons2,background1,background2,distancePixelToInit1,distancePixelToInit2,distancePixelToEdge1,distancePixelToEdge2,goodnessOfFit1,goodnessOfFit2,imscale=imscale,offsetX=offsetX,offsetY=offsetY,delimiter=delimiter)
