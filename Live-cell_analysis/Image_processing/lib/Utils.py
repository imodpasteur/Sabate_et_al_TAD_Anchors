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


class Utils:
    @staticmethod
    def saveImage(filepath,image,imagej_metadata,resolution):
        tifffile.imwrite(filepath, image, imagej=True,metadata=imagej_metadata,resolution=resolution)

    @staticmethod
    def openImage(filepath):
        with tifffile.TiffFile(filepath) as tif:
            imagej_hyperstack = tif.asarray()
            imagej_metadata = tif.imagej_metadata

            axes=tif.series[0].axes
            resolution=tif.pages[0].resolution

            imagej_metadata['axes']=axes
        return [imagej_hyperstack,imagej_metadata,resolution]

    @staticmethod
    def flipYaxis(image):
        return np.flip(image, axis=1)

    @staticmethod
    def image2coordinates(image,flipYaxis=True):
        #convert an image to a 3d matrix (to check registration)
        #by default, Y axis is inverted in chromagnon software
        x=[]
        y=[]
        z=[]
        amp=[]
        sh=np.shape(image)
        for u in range(sh[0]):
            for uu in range(sh[1]):
                for uuu in range(sh[2]):
                    if flipYaxis:
                        z.append(u)
                        y.append(uu)
                        x.append(uuu)
                        amp.append(image[u][sh[1]-uu-1][uuu])
                    else:
                        z.append(u)
                        y.append(uu)
                        x.append(uuu)
                        amp.append(image[u][uu][uuu])
        return np.array([x,y,z,amp]),sh

    @staticmethod
    def coordinates2image(coord,imshape,flipYaxis=True):
        #convert an image to a 3d matrix (to check registration)
        #by default, Y axis is inverted in chromagnon software
        #There is no interpolation in this function. So result might be weird (it is just to check transformation results)
        x=[]
        y=[]
        z=[]
        amp=[]
        image=np.zeros(imshape)
        for u in range(len(coord[0])):
            x=coord[0][u]
            y=coord[1][u]
            z=coord[2][u]
            amp=coord[3][u]
            if flipYaxis:
                if x>=0 and y>=0 and z>=0 and x<imshape[2]-1 and y<imshape[1]-1 and z<imshape[0]-1:
                    image[int(z+.5),imshape[1]-int(y+.5)-1,int(x+.5)]=amp
            else:
                if x>=0 and y>=0 and z>=0 and x<imshape[2]-1 and y<imshape[1]-1 and z<imshape[0]-1:
                    image[int(z+.5),int(y+.5),int(x+.5)]=amp
        image = scnd.maximum_filter(image, size=2)
        return image

    @staticmethod
    def testClass():
        print('test passed')
