import numpy as np
from scipy.optimize import curve_fit
import random


class Utils_coordinate_difference(object):

    @staticmethod
    def transform_data(distribution0input,distribution1input,shift=0):
        #input: 2 anchor coordinates
        #output: difference of anchor coordinates
        distribution0=np.vstack(distribution0input)
        distribution1=np.vstack(distribution1input)
        data_tmp=np.subtract(distribution0,distribution1)
        #get spherical coordinate to change norm, and go back to cartesian coordinate
        if shift>0:
            p=np.sqrt(np.sum(np.multiply(data_tmp,data_tmp),axis=1))
            theta=np.arccos(np.divide(data_tmp[:,2],p))
            phi=np.arctan2(data_tmp[:,1],data_tmp[:,0])

            p=np.abs(p-shift)

            data_x=np.multiply(p,np.multiply(np.sin(theta),np.cos(phi)))
            data_y=np.multiply(p,np.multiply(np.sin(theta),np.sin(phi)))
            data_z=np.multiply(p,np.cos(theta))
            data_x=np.transpose(data_x[np.newaxis])
            data_y=np.transpose(data_y[np.newaxis])
            data_z=np.transpose(data_z[np.newaxis])
            data_new=np.concatenate((data_x,data_y,data_z),axis=1)
        else:
            data_new=data_tmp
        return(data_new)


    @staticmethod
    def myArrayFlatten(X,Y,Z):
        F=np.transpose(np.array((X[0],Y[0],Z[0])))
        for i in range(1, len(X)):
            F=np.concatenate((F,np.transpose(np.array((X[i],Y[i],Z[i])))),axis=0)
        return F


    @staticmethod
    def readLoc_withoutImport(loc,imscale=1,offsetX=0,offsetY=0):

        x1=np.array(loc['Spot_1_X'].tolist())
        y1=np.array(loc['Spot_1_Y'].tolist())
        z1=np.array(loc['Spot_1_Z'].tolist())
        x2=np.array(loc['Spot_2_X'].tolist())
        y2=np.array(loc['Spot_2_Y'].tolist())
        z2=np.array(loc['Spot_2_Z'].tolist())
        frame=np.array(loc['Frame'].tolist())
        track1=np.array(loc['Track_pair'].tolist())

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
        return Frame,X1,Y1,Z1,X2,Y2,Z2


    @staticmethod
    def makeHistogram(data,bins):
        # convert 3D distribution to 3 concatenated 1D histograms
        tmp=[]
        for dim in [0,1,2]:
            tmp.append(np.histogram(data[:,dim], bins = (bins))[0])
        f=np.asarray(tmp).ravel()
        return(f)

