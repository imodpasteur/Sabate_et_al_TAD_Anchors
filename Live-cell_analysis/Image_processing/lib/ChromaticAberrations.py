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

class ChromaticAberrations:


    #sometime, output pixel size of chromagnon changes. In that case, set it using chromagnonPixSizeXYZ variable
    def __init__(self,chromagnon_filepath,drift_filepath,chromagnonPixSizeXYZ=[-1,-1,-1]):

        self.chromagnonPixSizeX=chromagnonPixSizeXYZ[0]
        self.chromagnonPixSizeY=chromagnonPixSizeXYZ[1]
        self.chromagnonPixSizeZ=chromagnonPixSizeXYZ[2]


        self.readChromagnonFile(chromagnon_filepath)
        if drift_filepath is not None:
            self.loadDrift(drift_filepath)




    def readChromagnonFile(self,filepath):
        f = open(filepath,"r")
        lines = f.readlines()

        #pix_size
        line_1=lines[2][:-1]
        line_2=line_1.split(',')
        self.pix_size_z=float(line_2[1])
        self.pix_size_y=float(line_2[2])
        self.pix_size_x=float(line_2[3])
        #refwave
        line_1=lines[3][:-1]
        line_2=line_1.split(',')
        self.refwave=int(line_2[1])
        #registration

        line_1=lines[6-self.refwave][:-1]
        line_2=line_1.split(',')
        if self.chromagnonPixSizeZ>0:
            self.tz=float(line_2[2])*self.pix_size_z/self.chromagnonPixSizeZ
        else:
            self.tz=float(line_2[2])

        if self.chromagnonPixSizeY>0:
            self.ty=-float(line_2[3])*self.pix_size_y/self.chromagnonPixSizeY
        else:
            self.ty=-float(line_2[3])

        if self.chromagnonPixSizeX>0:
            self.tx=float(line_2[4])*self.pix_size_x/self.chromagnonPixSizeX
        else:
            self.tx=float(line_2[4])


        self.rz=-float(line_2[5])
        self.mz=float(line_2[6])
        self.my=float(line_2[7])
        self.mx=float(line_2[8])


    def loadDrift(self,filepath):
        # drift dimension is is [frame,Y,X]
        self.drift=np.loadtxt(filepath, delimiter=",")





    def removeDriftCorrection(self):
        #input 2D table: [time, color, Z, Y, X]
        print('todo')

    def show(self):
        print('reference index:',self.refwave)
        print('pix size:',self.pix_size_z,self.pix_size_y,self.pix_size_x)
        print('trZYX:',self.tz,self.ty,self.tx)
        print('rz:',self.rz)
        print('magnZYX:',self.mz,self.my,self.mx)





    def translate(self,matrix):
        matrix[0]=np.add(matrix[0],self.tx)
        matrix[1]=np.add(matrix[1],self.ty)
        matrix[2]=np.add(matrix[2],self.tz)



    def rotate(self,matrix,sizeX,sizeY,sizeZ):



        cy=sizeY//2
        cx=sizeX//2



        #print('cx cy',cx,cy)
        xx=np.subtract(matrix[0],cx)
        yy=np.subtract(matrix[1],cy)

        d=np.sqrt(np.add(np.power(xx,2),np.power(yy,2)))
        angle=np.arctan2(yy,xx)
        angle=np.add(angle,self.rz*np.pi/180)

        matrix[0]=np.add(np.multiply(d,np.cos(angle)),cx)
        matrix[1]=np.add(np.multiply(d,np.sin(angle)),cy)


    def magn(self,matrix,sizeX,sizeY,sizeZ):



        cy=sizeY//2
        cx=sizeX//2
        cz=sizeZ//2

        xx=np.subtract(matrix[0],cx)
        yy=np.subtract(matrix[1],cy)
        zz=np.subtract(matrix[2],cz)

        xx=np.multiply(xx,self.mx)
        yy=np.multiply(yy,self.my)
        zz=np.multiply(zz,self.mz)



        #angle=np.add(angle,self.rz*np.pi/180)
        matrix[0]=np.add(xx,cx)
        matrix[1]=np.add(yy,cy)
        matrix[2]=np.add(zz,cz)

    def transform(self,matrix,sizeX,sizeY,sizeZ):
        #transform a matrix of coordinates
        #matrix contains at least 3 vectors [x, y, z, ...]
        #imshape is the image shape (depth, height, width)





        self.translate(matrix)

        self.magn(matrix,sizeX,sizeY,sizeZ)

        self.rotate(matrix,sizeX,sizeY,sizeZ)
        #
        #self.translate2(matrix)

    def reverseTransform(self,matrix,sizeX,sizeY,sizeZ):

        self.rz=-self.rz
        self.tx=-self.tx
        self.ty=-self.ty
        self.tz=-self.tz
        self.mx=1/self.mx
        self.my=1/self.my
        self.mz=1/self.mz

        self.rotate(matrix,sizeX,sizeY,sizeZ)

        self.magn(matrix,sizeX,sizeY,sizeZ)

        self.translate(matrix)

        self.rz=-self.rz
        self.tx=-self.tx
        self.ty=-self.ty
        self.tz=-self.tz
        self.mx=1/self.mx
        self.my=1/self.my
        self.mz=1/self.mz


    def translateimage(self,image):

        sh=np.shape(image)

        s=np.ones(np.shape(sh))




        if len(s)>=3:
            s[-2]=self.ty
            s[-1]=self.tx
            s[-3]=self.tz
            shifted=scnd.shift(image,s,order=1)
        else:
            shifted=scnd.shift(image,[self.ty,self.tx],order=1)
        return shifted

    def rotateimage(self,image):

        rot=scnd.rotate(image,-self.rz, axes=(-2, -1),reshape=False,order=1)
        return rot

    def magnimage(self,image):#different scale at the end

        sh=np.shape(image)
        pz=np.max((0,np.ceil(sh[-3]-sh[-3]*self.mz)))
        py=np.max((0,np.ceil(sh[-2]-sh[-2]*self.my)))
        px=np.max((0,np.ceil(sh[-1]-sh[-1]*self.mx)))

        s=np.ones(np.shape(sh))

        s[-3]=self.mz
        s[-2]=self.my
        s[-1]=self.mx

        AA=np.zeros((len(s),2))
        AA[-3]=[int(pz),int(pz)]
        AA[-2]=[int(py),int(py)]
        AA[-1]=[int(px),int(px)]
        AA=AA.tolist()
        aa=np.shape(AA)
        for i in range(aa[0]):
            for ii in range(aa[1]):
                AA[i][ii]=int(AA[i][ii])
        #print(AA)
        #for safety, we pad twice more than needed
        imagepad=np.pad(image,AA, 'constant', constant_values=(0))
        result = scnd.zoom(imagepad, s,order=1)
        shr=np.shape(result)
        r=result[...,shr[-3]//2-sh[-3]//2:shr[-3]//2-sh[-3]//2+sh[-3],shr[-2]//2-sh[-2]//2:shr[-2]//2-sh[-2]//2+sh[-2],shr[-1]//2-sh[-1]//2:shr[-1]//2-sh[-1]//2+sh[-1]]
        return r

    def transformimage3D(self,image):
        #transform an image

        image=self.magnimage(image)

        image=self.translateimage(image)

        image=self.rotateimage(image)



        return image



    def transformimage(self,image):
        #transform an image
        sh=np.shape(image)
        if len(sh)<=3:
            return self.transformimage3D(image)

        else:
            result=[]
            for i in range(sh[0]):
                print('registration',(sh[0]-i))
                result.append(self.transformimage3D(image[i]))

            return np.array(result)






    def transformimageUsingCoordinates(self,image,flip=False):
        coordB,imshapeB=ChromaticAbberations.image2coordinates(image,flipYaxis=flip)

        self.transform(coordB,imshapeB)

        return ChromaticAbberations.coordinates2image(coordB,imshapeB,flipYaxis=flip)


    def trackmate_FileLoading(self,filename):
        # importing element tree


        # Pass the path of the xml document
        self.tree = ET.parse(filename)

        # get the parent tag
        self.root = self.tree.getroot()


        imageData=self.root.find('Settings').find('ImageData')
        self.width=int(imageData.get('width'))
        self.height=int(imageData.get('height'))
        self.nslices=int(imageData.get('nslices'))
        self.nframes=int(imageData.get('nframes'))
        self.pixelwidth=float(imageData.get('pixelwidth'))
        self.pixelheight=float(imageData.get('pixelheight'))
        self.voxeldepth=float(imageData.get('voxeldepth'))
        self.timeinterval=float(imageData.get('timeinterval'))
        self.nspots=int(self.root.find('Model').find('AllSpots').get('nspots'))

        self.filename=imageData.get('filename')
        print('filename',self.filename)


        print('error if different: pix size X: ',self.pixelwidth,self.pix_size_x)
        print('error if different: pix size Y: ',self.pixelheight,self.pix_size_y)
        print('error if different: pix size Z: ',self.voxeldepth,self.pix_size_z)



    def trackmate_TableLoading(self):


        self.tableXYZ=[]
        self.tableFrame=[]

        for i,elt in enumerate(self.root.iter('Spot')):

            frame=int(elt.get('FRAME'))
            x=float(elt.get('POSITION_X'))
            y=float(elt.get('POSITION_Y'))
            z=float(elt.get('POSITION_Z'))
            x_pix=x/self.pixelwidth
            y_pix=y/self.pixelheight
            z_pix=z/self.voxeldepth


            self.tableXYZ.append([x_pix,y_pix,z_pix])
            print(np.shape(self.tableXYZ))
            self.tableFrame.append(frame)

    def trackmate_getfilename(self):
        imageData=self.root.find('Settings').find('ImageData')
        self.filename=imageData.get('filename')
        return self.filename

    def trackmate_setfilename(self,newfilename):
        imageData=self.root.find('Settings').find('ImageData')
        imageData.set('filename',newfilename)

    def trackmate_setfoldername(self,newfoldername):
        imageData=self.root.find('Settings').find('ImageData')
        imageData.set('folder',newfoldername)


    def trackmate_updateTable(self):

        for i,elt in enumerate(self.root.iter('Spot')):


            elt.set('POSITION_X',str(self.tableXYZ[i][0]*self.pixelwidth))
            elt.set('POSITION_Y',str(self.tableXYZ[i][1]*self.pixelheight))
            elt.set('POSITION_Z',str(self.tableXYZ[i][2]*self.voxeldepth))



    def trackmate_FileSaving(self,filepath):


        self.tree.write(filepath)


    def trackmate_applyDriftCorrection(self,sign):
        #sign=-1: remove drift correction
        #sign=+1: apply drift correction
        # drift dimension is [frame,Y,X]
        # tableXYZ dimension is [X,Y,Z]
        for i,f in enumerate(self.tableFrame):
            drift=self.drift[f]
            self.tableXYZ[i][0]+=sign*drift[2]
            self.tableXYZ[i][1]+=sign*drift[1]



    def trackmate_ChromaCorrection(self):

        matrix=np.transpose(self.tableXYZ)


        self.transform(matrix,self.width,self.height,self.nslices)

        self.tableXYZ=np.transpose(matrix)


    def trackmate_ReverseChromaCorrection(self):

        matrix=np.transpose(self.tableXYZ)


        self.reverseTransform(matrix,self.width,self.height,self.nslices)

        self.tableXYZ=np.transpose(matrix)


    def setTable(self,tableXYZ,tableFrame,w,h,d):
        self.width=w
        self.height=h
        self.nslices=d
        self.tableXYZ=tableXYZ
        self.tableFrame=tableFrame

    def getTable(self):
        return self.tableXYZ,self.tableFrame
