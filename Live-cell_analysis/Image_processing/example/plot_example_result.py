

import sys


import numpy as np

import pandas as pd

import os

import matplotlib.pyplot as plt



path_tracks=sys.argv[1]







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
    distance=np.array(loc["Distance"].tolist())
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
    Distance=[]
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
        Distance.append(distance[track1==t].tolist())
    return Frame,X1,Y1,Z1,X2,Y2,Z2,Track1,Distance





######################################### PARAMETERS :
#size of pixel [X,Y,Z] (in µm)
voxelsize=[0.1205,0.1205,0.2897]#um
imscale=1
offsetX=0
offsetY=0
delimiter="," # even ',' or ';'

Frame,X1,Y1,Z1,X2,Y2,Z2,Track,Distance=readLoc(path_tracks,voxelsize,imscale=imscale,offsetX=offsetX,offsetY=offsetY,delimiter=delimiter)




numb=0
plt.figure()
plt.plot(Frame[numb],Distance[numb])
plt.xlabel('frame')
plt.ylabel('distance (µm)')
plt.title('Distance between anchors for the track number '+str(numb))
plt.show()
