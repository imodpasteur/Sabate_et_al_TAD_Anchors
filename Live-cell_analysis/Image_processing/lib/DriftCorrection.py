import sys
import skimage as sk
from skimage.util import random_noise
import numpy as np
import imageio
import pandas as pd
import scipy.ndimage as scnd
import scipy.signal as scsig
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import os
import tifffile

class DriftCorrection:

    @staticmethod
    def averageColors(image):
        #input image is 5D: [time, Z,color,  Y, X]
        #output image is 4D: [time, Z, Y, X]
        return np.sum(image,axis=-3)

    @staticmethod
    def averageZ(image):
        #input image is 4D: [time, Z, Y, X]
        #output image is 3D: [time, Y, X]
        sh=np.shape(image)
        imageout=np.zeros((sh[0],sh[-2],sh[-1]))
        for i in range(sh[1]):
            tmp=np.subtract(image[:,i,:,:],np.min(image[:,i,:,:]))
            imageout=np.add(imageout,np.divide(tmp,np.max(tmp)))
        return imageout


    @staticmethod
    def computeDrift(image,refImage=0):
        #image is 3D: [time, Y, X]

        #sizeX=np.shape(image)[-1]
        #sizeY=np.shape(image)[-2]
        #image=image[:,int(sizeY/6):int(sizeY-sizeY/6),int(sizeX/6):int(sizeX-sizeX/6)]
        #viewer.add_image(image,name='image',gamma=0,contrast_limits=[np.min(image),np.max(image)])

        padding=1
        frameNumber=np.shape(image)[0]
        sizeX=np.shape(image)[-1]
        sizeY=np.shape(image)[-2]
        imA=np.zeros((sizeY*padding,sizeX*padding))
        imB=np.zeros((sizeY*padding,sizeX*padding))
        imA[sizeY*padding//2-sizeY//2:sizeY*padding//2-sizeY//2+sizeY,sizeX*padding//2-sizeX//2:sizeX*padding//2-sizeX//2+sizeX]=image[refImage,:,:]



        ftA=np.fft.fftn(imA)
        drift=np.zeros((frameNumber,3))
        for i in range(frameNumber):
            imB=np.zeros((sizeY*padding,sizeX*padding))
            imB[sizeY*padding//2-sizeY//2:sizeY*padding//2-sizeY//2+sizeY,sizeX*padding//2-sizeX//2:sizeX*padding//2-sizeX//2+sizeX]=image[i,:,:]
            imB=np.divide(imB,np.max(imB))


            ftB=np.fft.fftn(imB)

            #out=np.absolute(np.multiply(ftA,ftB))

            #viewer.add_image(out,name='red',gamma=0,contrast_limits=[np.min(out),np.max(out)])
            mul=np.multiply(ftA,np.conjugate(ftB))

            spatial=np.absolute(np.fft.ifft2(mul))
            spatial=np.fft.fftshift(spatial)

            #viewer.add_image(spatial,name='spatial',gamma=0,contrast_limits=[np.min(spatial),np.max(spatial)])
            ind=np.subtract(np.unravel_index(np.argmax(spatial, axis=None), spatial.shape),[padding*sizeY/2,padding*sizeX/2])
            drift[i]=[i,ind[0],ind[1]]
            #print((frameNumber-i),ind)
        return drift

    @staticmethod
    def saveDrift(drift,filepath):
        np.savetxt(filepath, drift, fmt='%s', delimiter=",")

    @staticmethod
    def loadDrift(filepath):
        drift=np.loadtxt(filepath, delimiter=",")
        return drift




    @staticmethod
    def shiftImage(image,drift):
        imshifted=image.copy()

        if np.shape(drift)[0] != np.shape(image)[0]:
            print('ERROR: drift frame number mismatch')
            return None

        for i in range(np.shape(imshifted)[0]):

            dx=int(drift[i,-1])
            dy=int(drift[i,-2])
            #print(i,dx,dy)
            imshifted[i]=np.roll(imshifted[i], dx, axis=-1)
            imshifted[i]=np.roll(imshifted[i], dy, axis=-2)

            if dx>0:
                imshifted[i,...,0:dx]=0
            elif dx<0:
                imshifted[i,...,dx-1:-1]=0
                imshifted[i,...,-1]=0

            if dy>0:
                imshifted[i,...,0:dy,:]=0
            elif dy<0:
                imshifted[i,...,dy-1:-1,:]=0
                imshifted[i,...,-1,:]=0
        return imshifted
        #input image dimension : [time , ... , Y , X]
        #drift is 3-columns: [time , Ydrift , Xdrift]
