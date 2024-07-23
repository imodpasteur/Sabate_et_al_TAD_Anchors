#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
#import gc
from PIL import Image
import sys
from numpy import linalg as LA

# In[194]:

pathimage=sys.argv[1]

folder,imagename=os.path.split(pathimage)




colorChannel=0
pathimage=os.path.join(folder,imagename)
thresholdcurv=1.1
thresholdLoG=1


# In[195]:





# In[214]:


class Utils:


    def saveTif(self,filepath,image):

        tifffile.imwrite(filepath, image)


    def saveTifImageJ(self,filepath,image):

        tifffile.imwrite(filepath, image, imagej=True,metadata=self.imagej_metadata)#,resolution=self.resolution)


    def openImage(self,filepath):
        with tifffile.TiffFile(filepath) as tif:
            imagej_hyperstack = tif.asarray()
            self.imagej_metadata = tif.imagej_metadata

            axes=tif.series[0].axes
            self.resolution=tif.pages[0].resolution
            self.imagej_metadata['axes']=axes
        return imagej_hyperstack


def removeAllLayers():
    viewer.layers.select_all()
    viewer.layers.remove_selected()



#performs 3D polynomial fit from a crop assuming center at the middle
def fitLoc(data,threshold=0,show=False):
    sh=np.shape(data)
    x0=sh[-1]//2
    y0=sh[-2]//2
    z0=sh[-3]//2

    x=[]
    #xx=[]
    #xxx=[]
    y=[]
    nbX=0
    maxiValue=data[x0,y0,z0]

    for u in range(sh[-1]):
            nbX+=1
            nbY=0
            for uu in range(sh[-2]):
                    nbY+=1
                    nbZ=0
                    for uuu in range(sh[-3]):
                            nbZ+=1
                            if np.abs(u-x0)<1.5 and np.abs(uu-y0)<1.5 and np.abs(uuu-z0)<1.5:
                                #if 3*3*3 middle window
                                y.append(data[uuu,uu,u])
                                x.append([u,uu,uuu])
                            elif (data[uuu,uu,u]>threshold):
                                y.append(data[uuu,uu,u])
                                x.append([u,uu,uuu])

                            #xx.append(uu)
                            #xxx.append(uuu)
    if show:
        im=data

        rangeplot=[np.min(np.sum(im,axis=0)),np.max(np.sum(im,axis=0))]
        fig=plt.figure()
        imageplt=plt.imshow(np.sum(im,axis=0), vmin=rangeplot[0], vmax=rangeplot[1])
        fig.colorbar(imageplt, location='right')



    poly = PolynomialFeatures(degree=2)
    X_t = poly.fit_transform(x)

    #linearRegression with fit_intercept=False is used to solve linear system
    #fit_intercept has to be false ; otherwise , 2 intercept are fitted (poli + linear) and intercept of polynomial becomes 0
    clf = LinearRegression(fit_intercept=False)
    clf.fit(X_t, y)

    pointC=np.transpose(np.array([[4,4,4]]))
    pointC_t = poly.fit_transform(np.transpose(pointC))
    pointsamp=clf.predict(pointC_t)[0]



    predmodel=clf.predict(X_t)

    serr=np.power(np.subtract(y,predmodel),2)
    mse=np.mean(serr)

    if show:
        model=np.zeros((nbZ,nbY,nbX))
        nb=0
        for u in range(nbX):
            for uu in range(nbY):
                for uuu in range(nbZ):

                    if np.abs(u-x0)<1.5 and np.abs(uu-y0)<1.5 and np.abs(uuu-z0)<1.5:
                        #if 3*3*3 middle window
                        model[uuu,uu,u]=predmodel[nb]
                        nb+=1
                    elif (data[uuu,uu,u]>threshold):
                        model[uuu,uu,u]=predmodel[nb]
                        nb+=1
        fig=plt.figure()
        imageplt=plt.imshow(np.sum(model,axis=0))#, vmin=rangeplot[0], vmax=rangeplot[1])
        fig.colorbar(imageplt, location='right')
        #print('msefit=',np.sqrt(np.mean(np.power(np.subtract(model,im),2))))


    coef=clf.coef_.copy()

    #print('coef',np.array(clf.coef_))
    #print(poly.get_feature_names_out())
    #print(clf.intercept_)
    #print("")

    #['1' 'x0'  'x1'  'x2' 'x0^2' 'x0 x1' 'x0 x2' 'x1^2' 'x1 x2' 'x2^2']
    #[0     1    2     3     4      5        6       7      8      9   ]
    #[a     b    c     d     e      f        g       h      i      j   ]

    #(x-A+B)


    #compute partial derivative to estimate parabola pic in 3d (equivalent to (-b/2a) in 1D):

    partialY=[-coef[1],-coef[2],-coef[3]]
    partialX=[[2*coef[4],coef[5],coef[6]],[coef[5],2*coef[7],coef[8]],[coef[6],coef[8],2*coef[9]]]
    center=np.linalg.solve(partialX, partialY)

    #print('center',center)

    xy=np.abs(coef[5])
    xz=np.abs(coef[6])
    yz=np.abs(coef[8])


    pointC=np.transpose(np.array([center]))
    pointC_t = poly.fit_transform(np.transpose(pointC))
    pointsamp=clf.predict(pointC_t)[0]


    #amplitude spot
    A=pointsamp*.9

    ################################################################################################

    elongationXY=-1
    elongationXZ=-1
    elongationYZ=-1
    elongationXYZ=-1

    for side in range(3):

        if side==0:#XY

            cx=center[0]
            cy=center[1]
            cz=center[2]

            coef=clf.coef_.copy()

        elif side==1:#XZ (switch Y and Z)

            cx=center[0]
            cz=center[1]
            cy=center[2]

            coef=clf.coef_.copy()
            coeftmp=clf.coef_.copy()

            coef[2]=coeftmp[3]
            coef[3]=coeftmp[2]

            coef[5]=coeftmp[6]
            coef[6]=coeftmp[5]

            coef[7]=coeftmp[9]
            coef[9]=coeftmp[7]

        elif side==2:#YZ (switch X and Z)

            cz=center[0]
            cy=center[1]
            cx=center[2]

            coef=clf.coef_.copy()
            coeftmp=clf.coef_.copy()

            coef[1]=coeftmp[3]
            coef[3]=coeftmp[1]

            coef[5]=coeftmp[8]
            coef[8]=coeftmp[5]

            coef[4]=coeftmp[9]
            coef[9]=coeftmp[4]


        allPoints=[]




        #get -45° value:




        #['1' 'x0'  'x1'  'x2' 'x0^2' 'x0 x1' 'x0 x2' 'x1^2' 'x1 x2' 'x2^2']
        #[0     1    2     3     4      5        6       7      8      9   ]
        #[a     b    c     d     e      f        g       h      i      j   ]





        #get -45° x coord for y=-x & z=center[2]  & f=A:
        #a+b*x+c*(-x+A+B)+d*z+e*x*x+f*x*(-x+A+B)+g*x*z+h*(-x+A+B)*(-x+A+B)+i*(-x+A+B)*z+j*z*z-O=0
        ap=coef[4]-coef[5]+coef[7]#y^2
        bp=+coef[1]-coef[2]+cx*coef[5]-2*cx*coef[7] +cy*coef[5]-2*cy*coef[7]+coef[6]*cz-coef[8]*cz
        cp=coef[0] + cx*cx*coef[7] + 2*cx*cy*coef[7] + cx*coef[2] + cx*coef[8]*cz + cy*cy*coef[7] + cy*coef[2] + cy*coef[8]*cz + coef[3]*cz + coef[9]*cz*cz -A
        delta=(bp**2)-(4*ap*cp)
        if delta<0:
            return -1,-1,-1,-1
        yx1=(-bp-np.sqrt(delta))/(2*ap)
        yx2=(-bp+np.sqrt(delta))/(2*ap)
        #print('x(-y)found',yx1,yx2)


        allPoints.append([yx1-cx,(-yx1+cy+cx)-cy])
        allPoints.append([yx2-cx,(-yx2+cy+cx)-cy])



        #get +45° x coord for y=x & z=center[2]  & f=A:
        #a+b*x+c*(x-A+B)+d*z+e*x*x+f*x*(x-A+B)+g*x*z+h*(x-A+B)*(x-A+B)+i*(x-A+B)*z+j*z*z-O=0
        ap=coef[4]+coef[5]+coef[7]#y^2
        bp=coef[1]+coef[2]-cx*coef[5]-2*cx*coef[7]+cy*coef[5]+2*cy*coef[7]+coef[6]*cz+coef[8]*cz
        cp=coef[0] + cx*cx*coef[7] - 2*cx*cy*coef[7] - cx*coef[2] - cx*coef[8]*cz + cy*cy*coef[7] + cy*coef[2] + cy*coef[8]*cz + coef[3]*cz + coef[9]*cz*cz -A
        delta=(bp**2)-(4*ap*cp)
        if delta<0:
            return -1,-1,-1,-1
        xy1=(-bp-np.sqrt(delta))/(2*ap)
        xy2=(-bp+np.sqrt(delta))/(2*ap)
        #print('xyfound',xy1,xy2)

        allPoints.append([xy1-cx,xy1+cy-cx-cy])
        allPoints.append([xy2-cx,xy2+cy-cx-cy])







        #get 0° x coord for y=center[1] & z=center[2]  & f=A:
        ap=coef[4]#y^2
        bp=coef[1]+coef[5]*cy+coef[6]*cz#y
        cp=coef[0]+(coef[2]*cy)+(coef[3]*cz)+(coef[7]*(cy**2))+(coef[8]*cy*cz)+(coef[9]*(cz**2)) -A#y
        delta=(bp**2)-(4*ap*cp)
        if delta<0:
            return -1,-1,-1,-1
        x1=(-bp-np.sqrt(delta))/(2*ap)
        x2=(-bp+np.sqrt(delta))/(2*ap)
        #print('xfound',x1,x2)


        allPoints.append([x1-cx,0])
        allPoints.append([x2-cx,0])






        #get 90° y coord for x=center[0] & z=center[2]  & f=A:
        ap=coef[7]#y^2
        bp=coef[2]+coef[5]*cx+coef[8]*cz#y
        cp=coef[0]+(coef[1]*cx)+(coef[3]*cz)+(coef[4]*(cx**2))+(coef[6]*cx*cz)+(coef[9]*(cz**2)) -A#y
        delta=(bp**2)-(4*ap*cp)
        if delta<0:
            return -1,-1,-1,-1
        y1=(-bp-np.sqrt(delta))/(2*ap)
        y2=(-bp+np.sqrt(delta))/(2*ap)
        #print('yfound',y1,y2)


        allPoints.append([0,y1-cy])
        allPoints.append([0,y2-cy])

        allPoints=np.array(allPoints)





        pts=np.transpose(allPoints)
        c=np.cov(pts)
        if not np.isnan(np.max(c)):
            w, v = LA.eig(c)


            if side==0:#XY

                if (np.min(w)>0):
                    elongationXY=np.max(w)/np.min(w)
                else:
                    elongationXY=100000
                elongationXY-=1

            elif side==1:#XZ (switch Y and Z)

                if (np.min(w)>0):
                    elongationXZ=np.max(w)/np.min(w)
                else:
                    elongationXZ=100000
                elongationXZ-=1

            elif side==2:#YZ (switch X and Z)

                if (np.min(w)>0):
                    elongationYZ=np.max(w)/np.min(w)
                else:
                    elongationYZ=100000
                elongationYZ-=1
        else:
            if side==0:#XY
                elongationXY=-1

            elif side==1:#XZ (switch Y and Z)

                elongationXZ=-1

            elif side==2:#YZ (switch X and Z)

                elongationYZ=-1



        if show:


            plt.figure()
            plt.plot(allPoints[:,0],allPoints[:,1],'+')
            if False:
                print('CHECK COORDINATES:')
                print('')


                #print([yx1,(-yx1+cy+cx),center[2]])
                #print([yx2,(-yx2+cy+cx),center[2]])
                pointN=np.transpose(np.array([[yx1,(-yx1+cy+cx),cz],[yx2,(-yx2+cy+cx),cz]]))
                #print('pointsShift',pointN)
                pointN_t = poly.fit_transform(np.transpose(pointN))
                #print('pointsShift_t',pointN_t)
                pointsampN=clf.predict(pointN_t)
                print('pointsamp YX',pointsampN)

                print('')





                print('')


                print([xy1,xy2,center[2]])
                print([xy2,xy2,center[2]])
                pointN=np.transpose(np.array([[xy1,xy1+cy-cx,cz],[xy2,xy2+cy-cx,cz]]))
                #print('pointsShift',pointN)
                pointN_t = poly.fit_transform(np.transpose(pointN))
                #print('pointsShift_t',pointN_t)
                pointsampN=clf.predict(pointN_t)
                print('pointsamp XY',pointsampN)

                print('')






                print('')


                #print([x1,center[1],center[2]])
                #print([x2,center[1],center[2]])
                pointN=np.transpose(np.array([[x1,cy,cz],[x2,cy,cz]]))
                #print('pointsShift',pointN)
                pointN_t = poly.fit_transform(np.transpose(pointN))
                #print('pointsShift_t',pointN_t)
                pointsampN=clf.predict(pointN_t)
                print('pointsampX',pointsampN)

                print('')





                #print([center[0],y1,center[2]])
                #print([center[0],y2,center[2]])
                pointN=np.transpose(np.array([[cx,y1,cz],[cx,y2,cz]]))
                #print('pointsShift',pointN)
                pointN_t = poly.fit_transform(np.transpose(pointN))
                #print('pointsShift_t',pointN_t)
                pointsampN=clf.predict(pointN_t)
                print('pointsampY',pointsampN)


                print('')

    elongationXYZ=np.max((elongationXY,elongationYZ,elongationXZ))

    return(center,mse,elongationXY,elongationXYZ)




# In[215]:


utils=Utils()


image5D=utils.openImage(pathimage)
image4D=image5D[:,:,colorChannel,:,:]
image5D=None
image4D=np.multiply(np.divide(image4D,np.percentile(image4D,90)),1.)


# In[216]:


shape=np.shape(image4D)
width=shape[-1]
height=shape[-2]


# In[219]:



import colorsys
test_color = colorsys.hsv_to_rgb(359/360.0, 1, 1)
def makecolorLUT():

    rgb=[colorsys.hsv_to_rgb(1,1,0)]
    start=120
    for i in range(1,256):
        angle=i/256;
        H =120.-((120.)*angle);
        if (H<0):
            H+=360;

        #rgb.append(H)
        rgb.append(colorsys.hsv_to_rgb(H/360,1,1))
    rgb=np.array(np.multiply(rgb,255),dtype='int')
    return rgb
def makegrayLUT():

    rgb=np.zeros((256,3))
    for i in range(256):
        rgb[i]=[i,i,i]
    return rgb

rgbcolor=makecolorLUT()
rgbgray=makegrayLUT()
#rgbtmp=[]
#for i in range(50):
#    rgbtmp.append(rgb)
#imLut = Image.fromarray(np.array(rgbtmp,dtype=np.uint8), 'RGB')
#imLut.show()


# In[ ]:





# In[220]:


toRGB=False
resultXY=np.zeros((height,width))
resultXYZ=np.zeros((height,width))
thresholdLoG=1
count=0
toto=[]
for a in image4D:
    tmp=np.array(-scnd.gaussian_laplace(a, 1))
    maxiFilter=scnd.maximum_filter(tmp,size=3)
    localMax=1*((maxiFilter==tmp) & (tmp>thresholdLoG))
    [z_d,y_d,x_d]=np.where(localMax>.5)
    shift=int(4)

    for i in range(len(z_d)):
        sh=np.shape(tmp)
        zz=int(z_d[i])
        yy=int(y_d[i])
        xx=int(x_d[i])
        if zz-shift>=0 and yy-shift>=0 and xx-shift>=0 and zz+shift+1<sh[0] and yy+shift+1<sh[1] and xx+shift+1<sh[2]:
            subImage=tmp[zz-shift:zz+shift+1,yy-shift:yy+shift+1,xx-shift:xx+shift+1]
            fitRes,msetmp,elongationXY,elongationXYZ=fitLoc(subImage,threshold=thresholdLoG,show=False)
            #print('ELONGATION',elongation)
            '''elongationXY/=thresholdcurv
            if elongationXY<0:
                elongationXY=0
            if elongationXY>1:
                elongationXY=1'''

            #result[yy,xx]=rgb[int(elongation)]
            if elongationXY>0:
                resultXY[yy,xx]=elongationXY
            if elongationXYZ>0:
                resultXYZ[yy,xx]=elongationXYZ
            #print(elongationXY)
            count+=1


# In[221]:



proj1=np.max(image4D,axis=1)
proj2=np.max(proj1,axis=0)


# In[222]:


if toRGB:
    im = Image.fromarray(np.array(resultXY,dtype=np.uint8), 'RGB')
    im.save(os.path.join(folder,imagename[:-4]+'_ReplicationMap.png'))
else:
    #result=np.multiply(result,1000)

    met=utils.imagej_metadata

    metnew={'ImageJ': '1.53t', 'images': 2, 'channels': 2}
    newLuts=[np.array(np.transpose(rgbgray),dtype=np.uint8),np.array(np.transpose(rgbcolor),dtype=np.uint8)]
    metnew['LUTs']=newLuts
    metnew['hyperstack']= True
    metnew['mode']='color'
    metnew['Composite mode']='composite'
    metnew['Ranges'] = (0, 0, 0, 1)




    imXY=np.array([proj2,resultXY],dtype=np.float32)
    path=os.path.join(folder,imagename[:-4]+'_ReplicationMapXY.tif')

    tifffile.imwrite(path, imXY, imagej=True,metadata=metnew)#,resolution=self.resolution)



    imXYZ=np.array([proj2,resultXYZ],dtype=np.float32)
    path=os.path.join(folder,imagename[:-4]+'_ReplicationMapXYZ.tif')

    tifffile.imwrite(path, imXYZ, imagej=True,metadata=metnew)#,resolution=self.resolution)



# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:
