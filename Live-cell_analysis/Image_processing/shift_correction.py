

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

from lib.ChromaticAberrations import *
from lib.Utils import *
from lib.DriftCorrection import *







# In[24]:

f=sys.argv[1]
chromagnonfile=sys.argv[2]
print('file ',f)
print('image dimension order has to be: TZCYX')
folder, imagename = os.path.split(f)


driftfilename=imagename[:-4]+'_drift.csv'
outputimagename=imagename[:-4]+'_DC.tif' #drift corrected image
outputimagename2=imagename[:-4]+'_DC_CHAB.tif' #CHromatic ABerration corrected image




imagepath=os.path.join(folder,imagename)
driftpath=os.path.join(folder,driftfilename)
outputimagepath=os.path.join(folder,outputimagename)
outputimagepath2=os.path.join(folder,outputimagename2)






[image5D,imagej_metadata,resolution]=Utils.openImage(imagepath)







frameNumber=np.shape(image5D)[0]






image4D=DriftCorrection.averageColors(image5D)
image3D=DriftCorrection.averageZ(image4D)

drift=DriftCorrection.computeDrift(image3D)

print('drift correction computed')


DriftCorrection.saveDrift(drift,driftpath)




imshifted=DriftCorrection.shiftImage(image5D,drift)
print('hyperStack image shifted')
Utils.saveImage(outputimagepath,imshifted,imagej_metadata,resolution)

print("drift correction terminated successfully.")

image4D=None
image3D=None






ca=ChromaticAberrations(chromagnonfile,driftpath,chromagnonPixSizeXYZ=[-1,-1,-1])
ca.show()


imageA=image5D[:,:,0,:,:]
imageB=image5D[:,:,1,:,:]



imageBreg=ca.transformimage(imageB)

image5Dres=[imageA,imageBreg]

image5Dres=np.moveaxis(image5Dres,0,2)


image5D=None

print("shift",drift)

imshiftednew=DriftCorrection.shiftImage(image5Dres,drift)

Utils.saveImage(outputimagepath2,imshiftednew,imagej_metadata,resolution)








































print('END')
