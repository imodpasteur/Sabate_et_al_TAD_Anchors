#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

plt.rcParams.update({'font.size': 22})
linewidth=2


processivityKb=float(sys.argv[1])
cohesinNumberPerMb=float(sys.argv[2])
avg_CTCF_attach_time=float(sys.argv[3])
sizeKb=int(sys.argv[4])    #918  #566
extrusionSpeedKb=float(sys.argv[5])
proba_CTCF_occupancy=float(sys.argv[6])
output_simulation_path=sys.argv[8]

def strtobool (val):
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
    else:
        raise ValueError("invalid truth value %r" % (val,))
showPlot=bool(strtobool(sys.argv[9]))

print("show:",showPlot)

print('path:',output_simulation_path)
#output_simulation_path='./A2Asimulations/240129_A2Asimulations_condition_test'

os.makedirs(output_simulation_path,exist_ok=True)



v=int(np.round(datetime.now().timestamp()*1000000))%(2**32)
random.seed(v)
np.random.seed(v)

# In[ ]:




# In[6]:


folderLammpsFiles='./Files_To_Copy_In_Each_Output_Folder'

folderLammpsFiles1=os.path.join(folderLammpsFiles,'base')

cell_line=utils.getCellLine(sizeKb)
folderLammpsFiles2=utils.getPathBatchLammps(sizeKb,folderLammpsFiles)

output_simulation_path=os.path.join(output_simulation_path,cell_line)
pathCTCF="./CTCF_positions/CTCF_sites_Cell_Lines.pkl"


###### lammpsRatio=3/2000 #3 seconds for 2000 time steps


###### deltaT=3*2705/2000#second

###### totaltime=deltaT*4100#second   2500


beadSizeResolutionKb=2

lammpsRatio=3/2000 #3 seconds for 2000 time steps

expectedExtrusionKbPerTimeStep=beadSizeResolutionKb*2 #2 beads / time step



effectiveExtrusionSpeedKb=extrusionSpeedKb*2 # each cohesin side extrude at extrusionSpeedKb


#deltaT is measured such as there are beadSizeResolutionKb*2 extruded per cohesin at each time step
deltaT=expectedExtrusionKbPerTimeStep/effectiveExtrusionSpeedKb





lammps_deltaT=int(deltaT/lammpsRatio)





lammps_run_1=100#run at the begining
lammps_run_1_ends_at_step=400#run at the begining

lammps_run_2=500#run at the begining
lammps_run_2_ends_at_step=800#run at the begining

lammps_run=lammps_deltaT
print('lammps_run',lammps_run)

lammps_recordingStartsAt=2100


#totaltime is half fixed such as we record 1800 time step when deltaT=4 (happens when extrusionSpeedKb=.5)
totaltime=lammps_recordingStartsAt*deltaT+1800*4




print('total time:',totaltime)

print('deltaT:',deltaT)

time=np.arange(0,totaltime,deltaT)


withAnchor2AnchorDistanceComputation=True
b=1


###### beadSizeResolutionKb=2





totalAreaMb=2.6#Mb






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


def makefolderpath(totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy):
    foldername='A2Asimulation_totaltime='+str(totaltime)+'_deltaT='+str(deltaT)+'_sizeKb='+str(sizeKb)+'_processivityKb='+str(processivityKb)+'_cohesinNumberPerMb='+str(cohesinNumberPerMb)+'_extrusionSpeedKb='+str(extrusionSpeedKb)+'_avg_CTCF_attach_time='+str(avg_CTCF_attach_time)+'_proba_CTCF_occupancy='+str(proba_CTCF_occupancy)

    #if not os.path.exists(os.path.join(output_simulation_path,foldername)):
    os.makedirs(os.path.join(output_simulation_path,foldername),exist_ok=True)
    return os.path.join(output_simulation_path,foldername)




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





def loadFullSimulation(totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy,simIndex=0):
    #filename='precomp_totaltime='+str(totaltime)+'_deltaT='+str(deltaT)+'_sizeKb='+str(sizeKb)+'_processivityKb='+str(processivityKb)+'_cohesinNumberPerMb='+str(cohesinNumberPerMb)+'_extrusionSpeedKb='+str(extrusionSpeedKb)+'_avg_CTCF_attach_time='+str(avg_CTCF_attach_time)+'_simIndex='+str(simIndex)
    foldername='A2Asimulation_totaltime='+str(totaltime)+'_deltaT='+str(deltaT)+'_sizeKb='+str(sizeKb)+'_processivityKb='+str(processivityKb)+'_cohesinNumberPerMb='+str(cohesinNumberPerMb)+'_extrusionSpeedKb='+str(extrusionSpeedKb)+'_avg_CTCF_attach_time='+str(avg_CTCF_attach_time)+'_proba_CTCF_occupancy='+str(proba_CTCF_occupancy)
    filename='simIndex='+str(simIndex)+'.json'
    fullpath=os.path.join(output_simulation_path,foldername,filename)
    with open(fullpath, 'rb') as file:
        # A new file will be created
        #pickle.dump(variable, file)
        variable=pickle.load(file)
    [time,anchor2anchor_genomic_dist,anchor2anchor_squared_dist,cohesin_number,k_on,k_off]=variable
    return time,anchor2anchor_genomic_dist,anchor2anchor_squared_dist,cohesin_number,k_on,k_off







#shift cohesin position to start at 0
#divide by beadSizeResolution
#remove cohesin loaded on itself and out of totalAreaMb range
def filterAndWeightCohesinAndCtcf(cohesinList,ctcfList,sizeKb,totalAreaMb,beadSizeResolutionKb,anchors):
    beadSize=beadSizeResolutionKb*1000
    shift=0#(totalAreaMb*1000000-sizeKb*1000)/2
    totalArea=int(totalAreaMb*1000000/beadSize)
    shiftPositions=1#we add 1 to not start at 0 (for lammps simulations)

    anchors[0]/=beadSize
    anchors[0]=np.round(anchors[0])
    anchors[0]+=(shift//beadSize)+shiftPositions #we add shiftPositions to not start at 0
    anchors[1]/=beadSize
    anchors[1]=np.round(anchors[1])
    anchors[1]+=(shift//beadSize)+shiftPositions #we add shiftPositions to not start at 0

    for ctcf in ctcfList:





        ctcf['position']/=beadSize
        ctcf['position']=np.round(ctcf['position'])
        ctcf['position']+=(shift//beadSize)+shiftPositions #we add shiftPositions to not start at 0


    for i,c in enumerate(cohesinList):
        r= range(len(c))

        for ii in r[::-1]:
            cc=c[ii]
            p1=cc['position1']
            p1/=beadSize
            p1=int(np.round(p1))
            p1+=(shift//beadSize)+shiftPositions #we add shiftPositions to not start at 0
            cc['position1']=int(p1)

            p2=cc['position2']
            p2/=beadSize
            p2=int(np.round(p2))
            p2+=(shift//beadSize)+shiftPositions #we add shiftPositions to not start at 0
            cc['position2']=int(p2)
            if int(np.abs(p1-p2))<=1:
                del c[ii]
                continue

            if p1<0+shiftPositions or p2<0+shiftPositions or p1>=totalArea+shiftPositions or p2>=totalArea+shiftPositions:
                del c[ii]
                continue










#def saveToLammpsFormat(folderpath,time,anchor2anchor_genomic_dist,cohesin_Full_List,CTCF_List,k_on,k_off,totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy,cohResidenceTime,totalAreaMb,lammps_recordingStartsAt,lammps_run_1,lammps_run_1_ends_at_step,lammps_run_2,lammps_run_2_ends_at_step,lammps_run,simIndex=0):
def saveToLammpsFormat(folderpath,cohesin_Full_List,lammps_recordingStartsAt,lammps_run_1,lammps_run_1_ends_at_step,lammps_run_2,lammps_run_2_ends_at_step,lammps_run,simIndex=0):

    filename='simIndex='+str(simIndex)

    fullpath=os.path.join(folderpath,filename)
    save_file = open(fullpath+".in", "w")


    def writeText(myfile,text):
        myfile.write('\n\n############################################################\n# '+str(text)+"\n############################################################\n\n")

    writeText(save_file,'Simulation parameters')
    hyperparam1='lattice fcc 4\nunits\t\tlj\nboundary\tf f f\natom_style\tmolecular\nlog \t\tlog.txt\nread_data\tinitial_conformation.txt\nneighbor 2.0\tmulti\n\ninclude\tinteractions\ninclude param'
    save_file.write(hyperparam1)

    writeText(save_file,'Definition of nucleus and its interaction')
    hyperparam2='region mySphere sphere 0.0 0.0 0.0 ${Radius} side in\n\nfix wall1 all wall/region mySphere lj126 1 0.5 0.5'
    save_file.write(hyperparam2)

    writeText(save_file,'Equilibration steps')
    hyperparam3='velocity \tall create 1.0 $v\nfix             1 all nve/limit 0.01\nfix\t\tlang all langevin 1.0 1.0 1.0 $l\nrun\t\t5000\nunfix 1\n\nfix             1 all nve/limit 0.02\nrun\t\t5000\nunfix 1\n\n\nthermo_style\tcustom step temp\nthermo          1000\nfix\t\t1 all nve/limit 0.05\ntimestep\t0.005\nrun\t\t10000000\nunfix 1\n\nspecial_bonds lj 1.0 0.0 0.0 extra 1000'
    save_file.write(hyperparam3)


    def writeTimePoint(myfile,timepoint,listOfBond,limit,run,timestep,thermo):
        writeText(myfile,'t = '+str(timepoint))
        for bond in listOfBond:
            myfile.write('fix ring'+str(bond[0])+' all bond/create 20 '+str(bond[1])+' '+str(bond[2])+' 4.499999999999999999 2\n')

        myfile.write('\n\nthermo_style\tcustom step temp\n')
        myfile.write('thermo          '+str(thermo)+'\n')
        myfile.write('fix\t\t1 all nve/limit '+str(limit)+'\n')
        myfile.write('timestep\t'+str(timestep)+'\n')
        myfile.write('run\t\t'+str(run)+'\n')

        for bond in listOfBond:
            myfile.write('unfix ring'+str(bond[0])+'\n')
        myfile.write('unfix 1\n\n')

        myfile.write('fix ringb1 all bond/break  1  2  0.0\ntimestep\t0.005\nrun\t\t5\nunfix ringb1\n\n')


    ringIndex=1
    limit=0.05

    timestep=0.005
    thermo=1000

    writeText(save_file,'Equilibriation of cohesin loading cycles without recording')

    for i,cohlist in enumerate(cohesin_Full_List):
        timepoint=i+1

        if (timepoint)==lammps_recordingStartsAt:#startRecording
            writeText(save_file,'Equilibriation of cohesin loading cycles with recording')
            save_file.write('dump  init all dcd 2000 '+filename+'.dcd\n\n')

        listOfBond=[]
        for coh in cohlist:
            p1=np.min((coh["position1"],coh["position2"]))
            p2=np.max((coh["position1"],coh["position2"]))
            #+1 because bonds start at 1
            listOfBond.append((ringIndex,p1,p2))
            ringIndex+=1

        if timepoint<=lammps_run_1_ends_at_step:
            writeTimePoint(save_file,timepoint,listOfBond,limit,lammps_run_1,timestep,thermo)
        elif timepoint<=lammps_run_2_ends_at_step:
            writeTimePoint(save_file,timepoint,listOfBond,limit,lammps_run_2,timestep,thermo)
        else:
            writeTimePoint(save_file,timepoint,listOfBond,limit,lammps_run-5,timestep,thermo)#-5 because there will be 5 more steps of ringb1 after










    save_file.close()

def copyFolderFiles(src_dir,dest_dir):


    # getting all the files in the source directory
    files = os.listdir(src_dir)

    for fname in files:
        shutil.copy2(os.path.join(src_dir,fname), dest_dir)







def simulateAnchor2AnchorTrajectory(pathCTCF,totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy,totalAreaMb,folderpath,withAnchor2AnchorDistanceComputation=False,indexSim=0,show=False):
    chr=None
    chr=Chromatin(processivityKb=processivityKb,cohesinNumberPerMb=cohesinNumberPerMb,extrusionSpeedKb=extrusionSpeedKb,avg_CTCF_attach_time=avg_CTCF_attach_time,proba_CTCF_occupancy=proba_CTCF_occupancy,deltaT=deltaT,totalAreaMb=totalAreaMb,folderSaveMovie=folderpath,indexSim=indexSim,show=show)

    ctcf_position,ctcf_orientation,ctcf_probability,anchors=readCTCFpositions(pathCTCF,sizeKb)
    ctcf_probability[0]=np.inf
    ctcf_probability[-1]=np.inf
    ctcf_position_kb=np.multiply(np.subtract(ctcf_position,1),beadSizeResolutionKb)
    anchorsKb=np.multiply(np.subtract(anchors,1),beadSizeResolutionKb)
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
            if i%10==0:
                chr.plot()

    return time,k_on,k_off,anchor2anchor_genomic_dist,cohesin_Full_List,CTCF_List,cohResidenceTime,chr.anchorA,chr.anchorB

'''

    if show:
        plt.figure(figsize=(10, 6), dpi=80)
        plt.plot(time/3600,anchor2anchor_genomic_dist,label='genomic anchor-anchor distance')
        plt.title('anchor-anchor distance')
        plt.show()



        plt.figure(figsize=(10, 6), dpi=80)
        plt.plot(time/3600,cohesin_number,label='cohesin number')
        plt.title('cohesin number')
        plt.legend()
        plt.show()
    return time,anchor2anchor_genomic_dist,k_on,k_off'''



#np.random.seed(None)
#random.seed(None)





#no loop : batch here

s=int(sys.argv[7])



print('simul',s,avg_CTCF_attach_time)
folderpath=makefolderpath(totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy)
if s==0:
    copyFolderFiles(folderLammpsFiles1,folderpath)
    copyFolderFiles(folderLammpsFiles2,folderpath)
else:
    timelib.sleep(.1)#to be sure the s==0 is already done
if s<2:
    show=showPlot
else:
    show=False
print('run simulation')
time,k_on,k_off,anchor2anchor_genomic_dist,cohesin_Full_List,CTCF_List,cohResidenceTime,anchorA,anchorB=simulateAnchor2AnchorTrajectory(pathCTCF,totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy,totalAreaMb,folderpath,withAnchor2AnchorDistanceComputation=withAnchor2AnchorDistanceComputation,show=show,indexSim=s)

print('simulation ok')

anchors=[anchorA,anchorB]

saveFullSimulation(folderpath,time,anchor2anchor_genomic_dist,cohesin_Full_List,CTCF_List,k_on,k_off,totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy,cohResidenceTime,totalAreaMb,anchors,s,extension_name='_raw')
filterAndWeightCohesinAndCtcf(cohesin_Full_List,CTCF_List,sizeKb,totalAreaMb,beadSizeResolutionKb,anchors)
saveToLammpsFormat(folderpath,cohesin_Full_List,lammps_recordingStartsAt,lammps_run_1,lammps_run_1_ends_at_step,lammps_run_2,lammps_run_2_ends_at_step,lammps_run,s)
saveFullSimulation(folderpath,time,anchor2anchor_genomic_dist,cohesin_Full_List,CTCF_List,k_on,k_off,totaltime,deltaT,sizeKb,processivityKb,cohesinNumberPerMb,extrusionSpeedKb,avg_CTCF_attach_time,proba_CTCF_occupancy,cohResidenceTime,totalAreaMb,anchors,s)

# In[ ]:


print('ok')
