

import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import cairo #pip install pycairo
import igraph as ig
import math
import pandas as pd
import os
import shutil
import pickle
import json
import sys
from utils import *
from Chromatin import *
from datetime import datetime
import time as timelib

processivityKb=300
cohesinNumberPerMb=8
extrusionSpeedKb=0.5
avg_CTCF_attach_time=150
proba_CTCF_occupancy=0.5
sizeKb=566
totalAreaMb=2.6
show=True
folderpath='./exampleOutput'
pathCTCF="./CTCF_positions/CTCF_sites_Cell_Lines.pkl"

deltaT=10
totaltime=100 #simulate 100 time points

print('the total simulated time is ',totaltime,' and the step time is ',deltaT)

os.makedirs(folderpath,exist_ok=True)


def readCTCFpositions(pathCTCF,sizeKb):

    fileCTCF = open(pathCTCF,'rb')
    object_ctcf = pickle.load(fileCTCF)

    cell_line=utils.getCellLine(sizeKb)



    ctcf_cl=object_ctcf[cell_line]
    position=ctcf_cl['CTCF_position']
    orientation=ctcf_cl['CTCF_orientation']
    probability=ctcf_cl['CTCF_proba']
    anchors=ctcf_cl['Anchors']
    anchors=np.array(anchors,dtype='int32')
    return position,orientation,probability,anchors



def simulateAnchor2AnchorTrajectory(pathCTCF,totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy,totalAreaMb,folderpath,withAnchor2AnchorDistanceComputation=False,indexSim=0,show=False):
    chr=None
    chr=Chromatin(processivityKb=processivityKb,cohesinNumberPerMb=cohesinNumberPerMb,extrusionSpeedKb=extrusionSpeedKb,avg_CTCF_attach_time=avg_CTCF_attach_time,proba_CTCF_occupancy=proba_CTCF_occupancy,deltaT=deltaT,totalAreaMb=totalAreaMb,folderSaveMovie=folderpath,indexSim=indexSim,show=show)

    beadSizeResolutionKbInCTCFfile=2
    ctcf_position,ctcf_orientation,ctcf_probability,anchors=readCTCFpositions(pathCTCF,sizeKb)
    ctcf_probability[0]=np.inf
    ctcf_probability[-1]=np.inf
    ctcf_position_kb=np.multiply(np.subtract(ctcf_position,1),beadSizeResolutionKbInCTCFfile)
    anchorsKb=np.multiply(np.subtract(anchors,1),beadSizeResolutionKbInCTCFfile)
    chr.setAnchors(anchorsKb)
    for i,pos in enumerate(ctcf_position_kb):
        if ctcf_orientation[i]=='+':
            if i==0 or i==len(ctcf_probability)-1:
                chr.addCTCFsite(ctcf_position_kb[i]*1000,False,ctcf_probability[i],isBound=True)
            else:
                chr.addCTCFsite(ctcf_position_kb[i]*1000,False,ctcf_probability[i],isBound=True)
        else:
            if i==0 or i==len(ctcf_probability)-1:
                chr.addCTCFsite(ctcf_position_kb[i]*1000,True,ctcf_probability[i],isBound=True)
            else:
                chr.addCTCFsite(ctcf_position_kb[i]*1000,True,ctcf_probability[i],isBound=True)



    chr.run(10000)#simulation stabilisation



    k_on=chr.k_on
    k_off=chr.k_off
    anchor2anchor_genomic_dist=[]
    cohesin_Full_List=[]
    CTCF_List=chr.getCTCFList()

    cohResidenceTime=chr.getCohesinResidenceTime()


    time=np.arange(0,totaltime,deltaT)
    for i,t in enumerate(time):
        chr.run(1,withAnchor2AnchorDistanceComputation=withAnchor2AnchorDistanceComputation)
        if withAnchor2AnchorDistanceComputation:
            d=chr.get_anchor2anchor_distanceKb()

        else:
            d=-1


        anchor2anchor_genomic_dist.append(d)

        coh_nb=chr.get_cohesin_number()

        newcohlist=[]
        for c in chr.getCohesinList():
            newcohlist.append(c.copy())
        cohesin_Full_List.append(newcohlist)



        if show :
            #print('plot at t=',t)
            if i%1==0:
                chr.plot()

    return time,k_on,k_off,anchor2anchor_genomic_dist,cohesin_Full_List,CTCF_List,cohResidenceTime,chr.anchorA,chr.anchorB



def saveFullSimulation(folderpath,time,anchor2anchor_genomic_dist,cohesin_Full_List,CTCF_List,k_on,k_off,totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy,cohResidenceTime,totalAreaMb,anchors,simIndex=0,extension_name=''):
    filename='simIndex='+str(simIndex)
    fullpath=os.path.join(folderpath,filename)
    for ctcf in CTCF_List:
        ctcf['position']=int(ctcf['position'])
    for coh in cohesin_Full_List:

        r= range(len(coh))
        for ii in r:
            cc=coh[ii]
            cc['position1']=int(cc['position1'])
            cc['position2']=int(cc['position2'])

    save_file = open(fullpath+extension_name+".json", "w")
    x={
        'totalSize(Mb)':totalAreaMb,
        'loopSize(Kb)':sizeKb,
        'totaltime(s)':totaltime,
        'anchors':(np.array(anchors,dtype="float32")).tolist(),
        'deltaT(s)':deltaT,
        'Kon(1/s)':k_on,
        'Koff(1/s)':k_off,
        'processivity(Kb)':processivityKb,
        'cohesinNumber(1/Mb)':cohesinNumberPerMb,
        'cohesinResidenceTime(s)':cohResidenceTime,
        'extrusionSpeed(Kb/s)':extrusionSpeedKb,
        'avg_CTCF_attach_time(s)':avg_CTCF_attach_time,
        'proba_CTCF_occupancy':proba_CTCF_occupancy,
        'time(s)':(np.array(time,dtype="float32")).tolist(),
        'ctcf':CTCF_List,
        'anchor2anchor_genomic_dist':anchor2anchor_genomic_dist,
        'cohesin_positions':cohesin_Full_List,
    }

    json.dump(x, save_file, indent = 4)
    save_file.close()



indexSim=0

time,k_on,k_off,anchor2anchor_genomic_dist,cohesin_Full_List,CTCF_List,cohResidenceTime,anchorA,anchorB=simulateAnchor2AnchorTrajectory(pathCTCF,totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy,totalAreaMb,folderpath,withAnchor2AnchorDistanceComputation=True,show=show,indexSim=indexSim)

anchors=[anchorA,anchorB]

saveFullSimulation(folderpath,time,anchor2anchor_genomic_dist,cohesin_Full_List,CTCF_List,k_on,k_off,totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy,cohResidenceTime,totalAreaMb,anchors,indexSim,extension_name='_raw')

print('simulation finished. It is saved in ',folderpath)
