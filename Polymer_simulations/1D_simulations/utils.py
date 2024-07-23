#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 15:16:58 2024

@author: benoit
"""
import os
import numpy as np
from sklearn.linear_model import LinearRegression
import pickle
class utils:
    
    
    
    @staticmethod
    def getCellLine(sizeKb):
                
                
        if sizeKb==345:
            cell_line='Lv1'
        elif sizeKb==566:
            cell_line='Lv4'
        elif sizeKb==918:
            cell_line='Tv1'
        elif sizeKb==576:
            cell_line='NoLoop'
        else:
            print('ERROR : polymer size / cell line not found ',sizeKb)
            cell_line='ErrorCellLineNotFound'
        return cell_line
    
    @staticmethod
    def getsizeKb(cell_line):
                
                
        if cell_line=='Lv1':
            sizeKb=345
        elif cell_line=='Lv4':
            sizeKb=566
        elif cell_line=='Tv1':
            sizeKb=918
        elif cell_line=='NoLoop':
            sizeKb=576
        else:
            print('ERROR : polymer size / cell line not found ',sizeKb)
            sizeKb='ErrorCellLineNotFound'
        return sizeKb



    @staticmethod
    def getPathBatchLammps(sizeKb,folderLammpsFiles):
                
                
        if sizeKb==345:
            folderLammpsFiles2=os.path.join(folderLammpsFiles,'runLV1')
        elif sizeKb==566:
            folderLammpsFiles2=os.path.join(folderLammpsFiles,'runLV4')
        elif sizeKb==918:
            folderLammpsFiles2=os.path.join(folderLammpsFiles,'runTV1')
        elif sizeKb==576:
            folderLammpsFiles2=os.path.join(folderLammpsFiles,'runNL')
        else:
            print('ERROR : polymer size / cell line not found ',sizeKb)
            folderLammpsFiles2='ErrorCellLineNotFound'
        return folderLammpsFiles2
    
    
        
        
    @staticmethod
    def lin_fit(mean_tracks_extr,deltaX,nbpts,withIntercept):
    
    
        xaxis=np.arange(-(nbpts-1)*deltaX,deltaX,deltaX)
        yaxis=mean_tracks_extr[len(mean_tracks_extr)-nbpts:]
        reg = LinearRegression(fit_intercept = withIntercept).fit(np.expand_dims(xaxis,-1), yaxis)
        coef=reg.coef_
        intercept=reg.intercept_
    
        yfit=xaxis*reg.coef_+intercept
        return xaxis,yfit,coef,intercept
    
    
    
    @staticmethod
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

    
    