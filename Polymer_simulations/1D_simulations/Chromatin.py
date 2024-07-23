# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import cairo #pip install pycairo
import igraph as ig
import math
import pandas as pd
plt.rcParams.update({'font.size': 22})
linewidth=2
import os
import shutil
import pickle
import json
import sys
from utils import *
from datetime import datetime
import time as timelib



class Chromatin:



    #sizeKb: size of simulated loop in Kb
    #processivityKb average extrusion before release of cohesin
    #cohesinNumberPerMb (#/Mb)
    #extrusionSpeedKb (kb/s)
    #residentime_ctcf: time ctcf sites stay on chromatin
    def __init__(self,processivityKb=150,cohesinNumberPerMb=1,extrusionSpeedKb=1,avg_CTCF_attach_time=1,proba_CTCF_occupancy=1,deltaT=1,totalAreaMb=1.8,show=False,folderSaveMovie='.',indexSim=0,isDiffusing=False):

        self.numberOfCohesin=0

        self.cohesinSizeBPequivalent=1

        self.version='ChromatinSimulator_0.0.1'

        self.isDiffusing=isDiffusing

        self.totalArea=totalAreaMb*1000000#2Mb


        self.dirplot=os.path.join(folderSaveMovie,'movie'+str(indexSim))#to save simulated plots

        self.index=0

        self.list_Cohesin=[]

        self.list_CTCF=[]
        self.isCTCFleft=True
        self.isCTCFright=True
        self.proba_CTCF_occupancy=proba_CTCF_occupancy

        self.time=0 #(s)

        self.deltaT=deltaT # simulation step time

        self.avg_CTCF_attach_time=avg_CTCF_attach_time

        if show:
            if os.path.exists(self.dirplot):
                shutil.rmtree(self.dirplot)
            os.mkdir(self.dirplot)

        self.deltaT=deltaT
        self.processivity=processivityKb*1000
        self.cohesinNumberPerBase=cohesinNumberPerMb/1000000
        self.extrusionSpeed=int(extrusionSpeedKb*1000*deltaT)


        self.cohesinresidencetime=self.processivity/(self.extrusionSpeed*2)# *2 because extrusion twice faster if no ctcf site reached


        self.k_off=(1/self.cohesinresidencetime)/self.deltaT
        self.k_on=(self.cohesinNumberPerBase*self.totalArea/self.cohesinresidencetime)/self.deltaT

        numberCohesin=self.cohesinNumberPerBase*self.totalArea
        for i in range(int(np.ceil(numberCohesin))):
            self._makeNewCohesin()
        if show:
            print('cohesinresidencetime=',self.cohesinresidencetime/self.deltaT,'second')
            print('extrusionSpeed=',(self.extrusionSpeed/1000)/self.deltaT,'kb per second')
            print('K_off: probability release of cohesin per second=',self.k_off)
            print('K_on: probability new cohesin per second=',self.k_on,' over ',totalAreaMb,'Mb area')

        self.listEffectiveProcessivity=[]
        self.listAnchor2AnchorDistance=[]

    def setCohesinSizeBPequivalent(self,cohSizeBp):
        self.cohesinSizeBPequivalent=cohSizeBp

    def setAnchors(self,anchorsKb):
        self.anchorA=anchorsKb[0]*1000
        self.anchorB=anchorsKb[1]*1000

    def getCohesinResidenceTime(self):
        return self.cohesinresidencetime*self.deltaT

    def _getProbaCTCF(self,weightCTCF_chipseq):

        p=(1-1/(weightCTCF_chipseq*self.avg_CTCF_attach_time))**self.deltaT

        return p

    def addCTCFsite(self,position,isOrientedLeft,weightCTCF_chipseq,isBound=True,probaOccupancy=-1):
        probaCTCF=self._getProbaCTCF(weightCTCF_chipseq)
        if probaOccupancy < 0:
            #ctcf={'position':position,'isOrientedLeft':isOrientedLeft,'isBound':isBound,'probaCTCF':probaCTCF,'time':weightCTCF_chipseq*self.avg_CTCF_attach_time,'probaOccupancy':self.proba_CTCF_occupancy}
            self.addCTCFsitePrecomputed(position,isOrientedLeft,isBound,probaCTCF,weightCTCF_chipseq*self.avg_CTCF_attach_time,self.proba_CTCF_occupancy)
        else:
            self.addCTCFsitePrecomputed(position,isOrientedLeft,isBound,probaCTCF,weightCTCF_chipseq*self.avg_CTCF_attach_time,probaOccupancy)
            #ctcf={'position':position,'isOrientedLeft':isOrientedLeft,'isBound':isBound,'probaCTCF':probaCTCF,'time':weightCTCF_chipseq*self.avg_CTCF_attach_time,'probaOccupancy':probaOccupancy}
        #self.list_CTCF.append(ctcf)

    def addCTCFsitePrecomputed(self,position,isOrientedLeft,isBound,probaCTCF,time,probaOccupancy):


        ctcf={'position':position,'isOrientedLeft':isOrientedLeft,'isBound':isBound,'probaCTCF':probaCTCF,'time':time,'probaOccupancy':probaOccupancy}
        self.list_CTCF.append(ctcf)

    def resetCohesins(self):
        self.list_Cohesin=[]

    def addCohesin(self,index,position1,position2,isExtrudingLeft,isExtrudingRight,indexCTCFLeft,indexCTCFRight):
        coh={'index':index,'position1':position1,'position2':position2,'isExtrudingLeft':isExtrudingLeft,'isExtrudingRight':isExtrudingRight,'indexCTCFLeft':indexCTCFLeft,'indexCTCFRight':indexCTCFRight}
        self.list_Cohesin.append(coh)
        self.index=index+1
        self.numberOfCohesin+=1

    def _makeNewCohesin(self):
        #isExtrudingLeft =0: in extrusion ; =1: blocked by CTCF ; =-1: passed a ctcf site
        position=random.randint(0,int(self.totalArea-1))
        coh={'index':self.index,'position1':position,'position2':position+1,'isExtrudingLeft':True,'isExtrudingRight':True,'indexCTCFLeft':-1,'indexCTCFRight':-1}
        self.list_Cohesin.append(coh)
        self.index+=1
        self.numberOfCohesin+=1


    def getCohesinList(self):
        return self.list_Cohesin

    def getCTCFList(self):
        return self.list_CTCF


    def _extrudeAllCohesins(self):


        for i,coh in enumerate(self.list_Cohesin):
            if self.isDiffusing:
                sideLeft=int(int(np.random.rand()*2)*2-1)
                leftExtrusion=coh['position1']-self.extrusionSpeed*sideLeft
                sideRight=int(int(np.random.rand()*2)*2-1)
                rightExtrusion=coh['position2']+self.extrusionSpeed*sideRight
            else:
                leftExtrusion=coh['position1']-self.extrusionSpeed
                rightExtrusion=coh['position2']+self.extrusionSpeed

            #check if already reached CTCF site
            if not coh['isExtrudingRight']:
                if random.random()>=self.list_CTCF[coh['indexCTCFRight']]['probaCTCF']:#check if released
                    coh['isExtrudingRight']=True#unblocked
                    coh['indexCTCFRight']=-1


            if not coh['isExtrudingLeft']:
                if random.random()>=self.list_CTCF[coh['indexCTCFLeft']]['probaCTCF']:#check if released
                    coh['isExtrudingLeft']=True#unblocked
                    coh['indexCTCFLeft']=-1

            #check if cohesin reaches a new CTCF site:
            if coh['isExtrudingRight']:
                rightIndexFound=-1
                for j,ctcf in enumerate(self.list_CTCF):
                    if ctcf['isBound']:
                        if not self.isDiffusing:
                            if ctcf['isOrientedLeft']:
                                if rightExtrusion>=ctcf['position'] and coh['position2']<ctcf['position']:
                                    if random.random()<ctcf['probaOccupancy']:#check if new CTCF block
                                        if random.random()<ctcf['probaCTCF']:#check if new CTCF block
                                            if rightIndexFound<0:
                                                rightIndexFound=j
                                            #else if another closed CTCF site at the same time --> select the closest
                                            elif ctcf['position']<self.list_CTCF[rightIndexFound]['position']:
                                                rightIndexFound=j
                        else:#we have to check double side in that case
                            if sideRight>0:#normal side
                                if ctcf['isOrientedLeft']:
                                    if rightExtrusion>=ctcf['position'] and coh['position2']<ctcf['position']:
                                        if random.random()<ctcf['probaOccupancy']:#check if new CTCF block
                                            if random.random()<ctcf['probaCTCF']:#check if new CTCF block
                                                if rightIndexFound<0:
                                                    rightIndexFound=j
                                                #else if another closed CTCF site at the same time --> select the closest
                                                elif ctcf['position']<self.list_CTCF[rightIndexFound]['position']:
                                                    rightIndexFound=j
                            else:
                                if not ctcf['isOrientedLeft']:
                                    if rightExtrusion<ctcf['position'] and coh['position2']>ctcf['position']:
                                        if random.random()<ctcf['probaOccupancy']:#check if new CTCF block
                                            if random.random()<ctcf['probaCTCF']:#check if new CTCF block
                                                if rightIndexFound<0:
                                                    rightIndexFound=j
                                                #else if another closed CTCF site at the same time --> select the closest
                                                elif ctcf['position']<self.list_CTCF[rightIndexFound]['position']:
                                                    rightIndexFound=j
                if rightIndexFound>=0:

                    coh['position2']=self.list_CTCF[rightIndexFound]['position']
                    coh['isExtrudingRight']=False
                    coh['indexCTCFRight']=rightIndexFound


            if coh['isExtrudingLeft']:
                leftIndexFound=-1
                for j,ctcf in enumerate(self.list_CTCF):
                    if ctcf['isBound']:
                        if not self.isDiffusing:
                            if not ctcf['isOrientedLeft']:
                                if leftExtrusion<=ctcf['position']  and coh['position1']>ctcf['position']:
                                    if random.random()<ctcf['probaOccupancy']:#check if new CTCF block
                                        if random.random()<ctcf['probaCTCF']:#check if new CTCF block
                                            if leftIndexFound<0:
                                                leftIndexFound=j
                                            #else if another closed CTCF site at the same time --> select the closest
                                            elif ctcf['position']>self.list_CTCF[leftIndexFound]['position']:
                                                leftIndexFound=j
                        else:#we have to check double side in that case
                            if sideLeft>0:
                                if not ctcf['isOrientedLeft']:
                                    if leftExtrusion<=ctcf['position']  and coh['position1']>ctcf['position']:
                                        if random.random()<ctcf['probaOccupancy']:#check if new CTCF block
                                            if random.random()<ctcf['probaCTCF']:#check if new CTCF block
                                                if leftIndexFound<0:
                                                    leftIndexFound=j
                                                #else if another closed CTCF site at the same time --> select the closest
                                                elif ctcf['position']>self.list_CTCF[leftIndexFound]['position']:
                                                    leftIndexFound=j
                            else:
                                if ctcf['isOrientedLeft']:
                                    if leftExtrusion>=ctcf['position']  and coh['position1']<ctcf['position']:
                                        if random.random()<ctcf['probaOccupancy']:#check if new CTCF block
                                            if random.random()<ctcf['probaCTCF']:#check if new CTCF block
                                                if leftIndexFound<0:
                                                    leftIndexFound=j
                                                #else if another closed CTCF site at the same time --> select the closest
                                                elif ctcf['position']>self.list_CTCF[leftIndexFound]['position']:
                                                    leftIndexFound=j
                if leftIndexFound>=0:
                    coh['position1']=self.list_CTCF[leftIndexFound]['position']
                    coh['isExtrudingLeft']=False
                    coh['indexCTCFLeft']=leftIndexFound

            if coh['isExtrudingRight']:
                coh['position2']=rightExtrusion#can be higher than size-1
            if coh['isExtrudingLeft']:
                coh['position1']=leftExtrusion#can be negative value




        #print(coh['index'],coh['position1'],coh['position2'])
        #break

    def _simulateOneStep(self):

        #####check if atached cohesin should detach
        length=len(self.list_Cohesin)
        for i,rc in enumerate(reversed(self.list_Cohesin)):
            probaCohesinDisappear=1/self.cohesinresidencetime
            isDetached=random.random()<probaCohesinDisappear

            if self.list_Cohesin[length-1-i]['position1']<0 or self.list_Cohesin[length-1-i]['position1']>self.totalArea or self.list_Cohesin[length-1-i]['position2']<0 or self.list_Cohesin[length-1-i]['position2']>self.totalArea:
                isDetached=True
                #detach is position is outside of total simulated area

            if isDetached:
                effectiveProcessivity=np.abs(self.list_Cohesin[length-1-i]['position2']-self.list_Cohesin[length-1-i]['position1'])
                self.listEffectiveProcessivity.append(effectiveProcessivity)
                del self.list_Cohesin[length-1-i]






        #####Extrude existing cohesins
        self._extrudeAllCohesins()


        #####check if new cohesin should atach
        probaNewCohesinAppear=self.cohesinNumberPerBase*(self.totalArea)/self.cohesinresidencetime
        #loop because proba can be >1
        for p in range(int(np.ceil(probaNewCohesinAppear))):
            probatmp=probaNewCohesinAppear-p#the last loop might be 0<probatmp<1
            isNewCohesin=random.random()<probatmp
            if isNewCohesin:
                self._makeNewCohesin()






        self.time+=self.deltaT

    def _test(self):

        coh={'index':0,'position1':10000,'position2':60000}
        self.list_Cohesin.append(coh)

        coh={'index':1,'position1':50000,'position2':200000}
        self.list_Cohesin.append(coh)

        coh={'index':1,'position1':20000,'position2':30000}
        self.list_Cohesin.append(coh)

        coh={'index':1,'position1':210000,'position2':230000}
        self.list_Cohesin.append(coh)

        coh={'index':1,'position1':299999,'position2':250000}
        self.list_Cohesin.append(coh)


    def _computeAnchorDistance(self):

        cohesinSizeBPequivalent=self.cohesinSizeBPequivalent

        #first, check if one cohesin is bound to both CTCF sites:
        for j,coh in enumerate(self.list_Cohesin):
            if coh['position2']==self.anchorB and coh['position1']==self.anchorA:
                self.anchor2anchor_shortest_path=[{'current':self.anchorA,'next':self.anchorB,'isLoop':True}]
                self.anchor_distance=cohesinSizeBPequivalent
                return







        edges=[]
        weight=[]
        isLoop=[]


        nbcoh=len(self.list_Cohesin)
        indexA=0
        indexB=1
        shiftPosit1=2
        shiftPosit2=nbcoh+2

        edges.append((indexA,indexB))
        weight.append((self.anchorB-self.anchorA))
        isLoop.append(False)

        #compute nodes between anchor1 and all cohesin:
        for j,coh in enumerate(self.list_Cohesin):
            edges.append((indexA,j+shiftPosit1))
            weight.append(abs(coh['position1']-(self.anchorA)))
            isLoop.append(False)
            edges.append((indexA,j+shiftPosit2))
            weight.append(abs(coh['position2']-(self.anchorA)))
            isLoop.append(False)

        #compute nodes between anchor2 and all cohesin:
        for j,coh in enumerate(self.list_Cohesin):
            edges.append((indexB,j+shiftPosit2))
            weight.append(abs(coh['position2']-(self.anchorB)))
            isLoop.append(False)
            edges.append((indexB,j+shiftPosit1))
            weight.append(abs(coh['position1']-(self.anchorB)))
            isLoop.append(False)

        #compute nodes between cohesins
        for i,cohA in enumerate(self.list_Cohesin):
            for j in range(i,len(self.list_Cohesin)):
                cohB=self.list_Cohesin[j]
                if i==j:#if same sohesin
                    edges.append((i+shiftPosit1,i+shiftPosit2))
                    weight.append(cohesinSizeBPequivalent)
                    isLoop.append(True)
                else:
                    dist=abs(cohB['position2']-cohA['position2'])
                    edges.append((j+shiftPosit2,i+shiftPosit2))
                    weight.append(dist)
                    isLoop.append(False)

                    dist=abs(cohB['position1']-cohA['position1'])
                    edges.append((j+shiftPosit1,i+shiftPosit1))
                    weight.append(dist)
                    isLoop.append(False)

                    dist=abs(cohB['position1']-cohA['position2'])
                    edges.append((j+shiftPosit1,i+shiftPosit2))
                    weight.append(dist)
                    isLoop.append(False)

                    dist=abs(cohB['position2']-cohA['position1'])
                    edges.append((j+shiftPosit2,i+shiftPosit1))
                    weight.append(dist)
                    isLoop.append(False)





        g = ig.Graph(nbcoh+2,edges,directed=False)

        path = g.get_shortest_paths(indexA, to=indexB, weights=weight,output='epath')[0]


        #for i,p in enumerate(edges):
        #    print('ed',i,p)

        #for i,p in enumerate(path):
        #    print('p',p,edges[p])


        self.anchor_distance=0
        #smallest_path=[]
        self.anchor2anchor_shortest_path=[]
        nextnode=set()
        for i,p in enumerate(path):
            if i==0:#anchor A



                nextnode=set(edges[p])-set([indexA])
                nextEl=list(nextnode)[0]
                #print('i',i,nextEl,shiftPosit1,shiftPosit2)
                nextElposition=self.anchorB
                if nextEl<shiftPosit1:
                    nextElposition=self.anchorB
                elif nextEl<shiftPosit2:
                    nextElposition=self.list_Cohesin[nextEl-shiftPosit1]["position1"]
                else:

                    nextElposition=self.list_Cohesin[nextEl-shiftPosit2]["position2"]

                self.anchor2anchor_shortest_path.append({'current':self.anchorA,'next':nextElposition,'isLoop':isLoop[p]})
                self.anchor_distance+=weight[p]

            elif (i<len(path)-1):
                current=nextnode.copy()
                nextnode=set(edges[p])-current
                currentEl=list(current)[0]
                nextEl=list(nextnode)[0]
                #r=[currentEl,nextEl]
                self.anchor_distance+=weight[p]
                #smallest_path.append(r)


                nextElposition=self.anchorB
                if nextEl<shiftPosit1:
                    print('WARNING: should never be')
                    nextElposition=self.anchorB
                elif nextEl<shiftPosit2:
                    nextElposition=self.list_Cohesin[nextEl-shiftPosit1]["position1"]
                else:
                    nextElposition=self.list_Cohesin[nextEl-shiftPosit2]["position2"]

                currentElPosition=0
                if currentEl<shiftPosit2:
                    currentElposition=self.list_Cohesin[currentEl-shiftPosit1]["position1"]
                else:
                    currentElposition=self.list_Cohesin[currentEl-shiftPosit2]["position2"]

                self.anchor2anchor_shortest_path.append({'current':currentElposition,'next':nextElposition,'isLoop':isLoop[p]})

            else:
                self.anchor_distance+=weight[p]

                self.anchor2anchor_shortest_path.append({'current':nextElposition,'next':self.anchorB,'isLoop':isLoop[p]})

        self.listAnchor2AnchorDistance.append(self.anchor_distance)




    def computeAndGet_HiC_map(self,resolutionBP=10000):

        cohesinSizeBPequivalent=self.cohesinSizeBPequivalent

        distanceMap=np.zeros((int(self.totalArea/resolutionBP),int(self.totalArea/resolutionBP)))

        for indexMapi,anchorA in enumerate(np.arange(0,self.totalArea,resolutionBP)):
            for indexMap_shift,anchorB in enumerate(np.arange(anchorA,self.totalArea,resolutionBP)):
                indexMapii=indexMapi+indexMap_shift


                edges=[]
                weight=[]
                isLoop=[]


                nbcoh=len(self.list_Cohesin)
                indexA=0
                indexB=1
                shiftPosit1=2
                shiftPosit2=nbcoh+2

                edges.append((indexA,indexB))
                weight.append((anchorB-anchorA))
                isLoop.append(False)

                #compute nodes between anchor1 and all cohesin:
                for j,coh in enumerate(self.list_Cohesin):
                    edges.append((indexA,j+shiftPosit1))
                    weight.append(abs(coh['position1']-(anchorA)))
                    isLoop.append(False)
                    edges.append((indexA,j+shiftPosit2))
                    weight.append(abs(coh['position2']-(anchorA)))
                    isLoop.append(False)

                #compute nodes between anchor2 and all cohesin:
                for j,coh in enumerate(self.list_Cohesin):
                    edges.append((indexB,j+shiftPosit2))
                    weight.append(abs(coh['position2']-(anchorB)))
                    isLoop.append(False)
                    edges.append((indexB,j+shiftPosit1))
                    weight.append(abs(coh['position1']-(anchorB)))
                    isLoop.append(False)

                #compute nodes between cohesins
                for i,cohA in enumerate(self.list_Cohesin):
                    for j in range(i,len(self.list_Cohesin)):
                        cohB=self.list_Cohesin[j]
                        if i==j:#if same sohesin
                            edges.append((i+shiftPosit1,i+shiftPosit2))
                            weight.append(cohesinSizeBPequivalent)
                            isLoop.append(True)
                        else:
                            dist=abs(cohB['position2']-cohA['position2'])
                            edges.append((j+shiftPosit2,i+shiftPosit2))
                            weight.append(dist)
                            isLoop.append(False)

                            dist=abs(cohB['position1']-cohA['position1'])
                            edges.append((j+shiftPosit1,i+shiftPosit1))
                            weight.append(dist)
                            isLoop.append(False)

                            dist=abs(cohB['position1']-cohA['position2'])
                            edges.append((j+shiftPosit1,i+shiftPosit2))
                            weight.append(dist)
                            isLoop.append(False)

                            dist=abs(cohB['position2']-cohA['position1'])
                            edges.append((j+shiftPosit2,i+shiftPosit1))
                            weight.append(dist)
                            isLoop.append(False)





                g = ig.Graph(nbcoh+2,edges,directed=False)

                path = g.get_shortest_paths(indexA, to=indexB, weights=weight,output='epath')[0]


                #for i,p in enumerate(edges):
                #    print('ed',i,p)

                #for i,p in enumerate(path):
                #    print('p',p,edges[p])


                anchor_distance=0
                #smallest_path=[]
                anchor2anchor_shortest_path=[]
                nextnode=set()
                for i,p in enumerate(path):
                    if i==0:#anchor A



                        nextnode=set(edges[p])-set([indexA])
                        nextEl=list(nextnode)[0]
                        #print('i',i,nextEl,shiftPosit1,shiftPosit2)
                        nextElposition=anchorB
                        if nextEl<shiftPosit1:
                            nextElposition=anchorB
                        elif nextEl<shiftPosit2:
                            nextElposition=self.list_Cohesin[nextEl-shiftPosit1]["position1"]
                        else:

                            nextElposition=self.list_Cohesin[nextEl-shiftPosit2]["position2"]

                        anchor2anchor_shortest_path.append({'current':anchorA,'next':nextElposition,'isLoop':isLoop[p]})
                        anchor_distance+=weight[p]

                    elif (i<len(path)-1):
                        current=nextnode.copy()
                        nextnode=set(edges[p])-current
                        currentEl=list(current)[0]
                        nextEl=list(nextnode)[0]
                        #r=[currentEl,nextEl]
                        anchor_distance+=weight[p]
                        #smallest_path.append(r)


                        nextElposition=anchorB
                        if nextEl<shiftPosit1:
                            print('WARNING: should never be')
                            nextElposition=anchorB
                        elif nextEl<shiftPosit2:
                            nextElposition=self.list_Cohesin[nextEl-shiftPosit1]["position1"]
                        else:
                            nextElposition=self.list_Cohesin[nextEl-shiftPosit2]["position2"]

                        currentElPosition=0
                        if currentEl<shiftPosit2:
                            currentElposition=self.list_Cohesin[currentEl-shiftPosit1]["position1"]
                        else:
                            currentElposition=self.list_Cohesin[currentEl-shiftPosit2]["position2"]

                        anchor2anchor_shortest_path.append({'current':currentElposition,'next':nextElposition,'isLoop':isLoop[p]})

                    else:
                        anchor_distance+=weight[p]

                        anchor2anchor_shortest_path.append({'current':nextElposition,'next':anchorB,'isLoop':isLoop[p]})

                distanceMap[indexMapi,indexMapii]=anchor_distance
        for i in range(len(distanceMap)):
            distanceMap[i:,i]=distanceMap[i,i:]
        return distanceMap

    def run(self,step_number,withAnchor2AnchorDistanceComputation=False):
        self.withAnchor2AnchorDistanceComputation=withAnchor2AnchorDistanceComputation
        for t in range(step_number):
            self._simulateOneStep()
            if withAnchor2AnchorDistanceComputation:
                self._computeAnchorDistance()


    def computeAnchor2AnchorDistance(self):
        self.withAnchor2AnchorDistanceComputation=True
        self._computeAnchorDistance()


    def get_anchor2anchor_distanceKb(self):
        if self.withAnchor2AnchorDistanceComputation:
            return(self.anchor_distance)/1000
        else:
            print('anchor 2 anchor distance not computed')
            return None

    def get_cohesin_number(self):
        return(len(self.list_Cohesin))

    def plot(self,show=False):
        self._computeAnchorDistance()
        fig=plt.figure()

        maxictcftime=0
        for ctcf in self.list_CTCF:
            if ctcf['time']>maxictcftime:
                if ctcf['time']<np.inf:
                    maxictcftime=ctcf['time']

        #compute macimum arc to get maxi
        x,y=self._getArc(0,self.totalArea)
        maxi=np.max(y)

        for j,coh in enumerate(self.list_Cohesin):
            #plt.plot(coh['position1'],0,'bo')
            #plt.plot(coh['position2'],0,'bo')
            #if coh['position1']>=0 and coh['position1']<=self.totalArea and coh['position2']>=0 and coh['position2']<=self.totalArea :
            x,y=self._getArc(coh['position1'],coh['position2'])
            plt.plot(x/1000,y,'k')
        if self.withAnchor2AnchorDistanceComputation:
            for p in self.anchor2anchor_shortest_path:
                #if p['current']>=0 and p['current']<=self.totalArea and p['next']>=0 and p['next']<=self.totalArea :
                if p['isLoop']:
                    x,y=self._getArc(p['current'],p['next'])
                    if np.max(y)>maxi:
                        maxi=np.max(y)
                    plt.plot(x/1000,y,'r')
                else:
                    plt.plot([p['current']/1000,p['next']/1000],[0,0],'r')


        for j,ctcf in enumerate(self.list_CTCF):
            if ctcf['isBound']:
                if ctcf['isOrientedLeft']:
                    if ctcf['time']>=np.inf:
                        plt.plot(ctcf['position']/1000,0,'g<')
                    else:
                        plt.plot(ctcf['position']/1000,ctcf['time']*maxi/maxictcftime,'g<')
                        plt.plot((ctcf['position']/1000,ctcf['position']/1000),(ctcf['position']/1000,ctcf['time']*maxi/maxictcftime),'g')
                else:
                    if ctcf['time']>=np.inf:
                        plt.plot(ctcf['position']/1000,0,'b>')
                    else:
                        plt.plot(ctcf['position']/1000,ctcf['time']*maxi/maxictcftime,'b>')
                        plt.plot((ctcf['position']/1000,ctcf['position']/1000),(ctcf['position']/1000,ctcf['time']*maxi/maxictcftime),'b')

        plt.plot((self.anchorA/1000,self.anchorA/1000),(0,maxi),'r--')
        plt.plot((self.anchorB/1000,self.anchorB/1000),(0,maxi),'r--')
        plt.xlim((-50, int(self.totalArea/1000)+50))
        plt.gca().axes.get_yaxis().set_visible(False)

        if show:
            plt.show
        else:
            plt.savefig(os.path.join(self.dirplot,'plot_'+str(self.time)+'.png'))
            plt.close(fig)




    def _getArc(self,mini,maxi):
        angles = np.linspace(0 , np.pi, 100 )
        xs =  (np.cos(angles)/2+.5)*(maxi-mini)+mini
        ys =  np.abs(np.sin(angles)*(maxi-mini))
        return xs, ys

    def show(self):
        print('number of cohesin:',len(self.list_Cohesin))
        print('anchor distance:',self.anchor_distance)
        #for i,c in enumerate(self.list_Cohesin):
        #    print('cohesin',i,c)
