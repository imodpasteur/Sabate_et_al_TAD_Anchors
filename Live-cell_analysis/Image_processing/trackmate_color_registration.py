

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
import xml.etree.ElementTree as ET

from lib.ChromaticAberrations import *
from lib.Utils import *
from lib.DriftCorrection import *



trackRedCh0path=sys.argv[1]

chromagnonfile=sys.argv[2]

folderDrift=sys.argv[3]

folderImage=sys.argv[4]



folder,name=os.path.split(trackRedCh0path)


trackGreenCh1path=trackRedCh0path[:-8]+'GFP.xml'



driftFilePath=os.path.join(folderDrift,name[:-11])+'drift.csv'


imagepath=os.path.join(folderImage,name[:-8]+"CHAB.tif")




if not os.path.isfile(trackRedCh0path):
    print('FILE NOT FOUND',trackRedCh0path)
    exit(0)

if not os.path.isfile(trackGreenCh1path):
    print('FILE NOT FOUND',trackGreenCh1path)
    exit(0)

if not os.path.isfile(driftFilePath):
    print('FILE NOT FOUND',driftFilePath)
    exit(0)

if not os.path.isfile(chromagnonfile):
    print('FILE NOT FOUND',chromagnonfile)
    exit(0)

if not os.path.isfile(imagepath):
    print('FILE NOT FOUND',imagepath)
    exit(0)

print('processing ',name[:-8])




trackRedCh0Outputpath=trackRedCh0path[:-4]+'_CHAB.xml'
trackGreenCh1Outputpath=trackGreenCh1path[:-4]+'_CHAB.xml'



ca=ChromaticAberrations(chromagnonfile,driftFilePath,chromagnonPixSizeXYZ=[-1,-1,-1])
ca.show()



def updateFileName(ca):
    filenameImage=ca.trackmate_getfilename()
    filenameImagenew=filenameImage[:-4]+'_CHAB'+filenameImage[-4:]
    ca.trackmate_setfilename(filenameImagenew)





#GREEN COLOR:

ca.trackmate_FileLoading(trackGreenCh1path)

ca.trackmate_TableLoading()

ca.trackmate_setfoldername(folderImage)
updateFileName(ca)

ca.trackmate_applyDriftCorrection(-1)

ca.trackmate_ChromaCorrection()

ca.trackmate_applyDriftCorrection(+1)

ca.trackmate_updateTable()

ca.trackmate_FileSaving(trackGreenCh1Outputpath)




#RED COLOR (just update image path):

ca.trackmate_FileLoading(trackRedCh0path)

ca.trackmate_setfoldername(folderImage)
updateFileName(ca)

ca.trackmate_FileSaving(trackRedCh0Outputpath)






























print('END')
