#===============================================================================
# 
#  License: GPL
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License 2
#  as published by the Free Software Foundation.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# 
#===============================================================================
 
 
import Image
import ImageEnhance
import numpy as np
import threading
import os
import Queue
from CenterRectangle import CenterRectangle
from matplotlib.lines import Line2D
#implicity this relies upon matplotlib.axis matplotlib.AxisImage matplotlib.bar 


#my custom 2d correlation function for numpy 2d matrices.. 
def mycorrelate2d(fixed,moved,skip=1):
    """a 2d correlation function for numpy 2d matrices
    
    arguments
    fixed) is the larger matrix which should stay still 
    moved) is the smaller matrix which should move left/right up/down and sample the correlation
    skip) is the number of positions to skip over when sampling, 
    so if skip =3 it will sample at shift 0,0 skip,0 2*skip,0... skip,0 skip,skip...
    
    returns
    corrmat) the 2d matrix with the corresponding correlation coefficents of the data at that offset
    note the 0,0 entry of corrmat corresponds to moved(0,0) corresponding to fixed(0,0)
    and the 1,1 entry of corrmat corresponds to moved(0,0) corresponding to fixed(skip,skip)
    NOTE) the height of corrmat is given by corrmat.height=ceil((fixed.height-moved.height)/skip)
    and the width in a corresonding manner.
    NOTE)the standard deviation is measured over the entire dataset, so particular c values can be above 1.0
    if the variance in the subsampled region of fixed is lower than the variance of the entire matrix
    
    """
    
    (fh,fw)=fixed.shape
    (mh,mw)=moved.shape
    deltah=(fh-mh)
    deltaw=(fw-mw)
    if (deltah<1 or deltaw<1):
        return
    fixed=fixed-fixed.mean()
    fixed=fixed/fixed.std()
    moved=moved-moved.mean()
    moved=moved/moved.std()
    ch=np.ceil(deltah*1.0/skip)
    cw=np.ceil(deltaw*1.0/skip)
    
    corrmat=np.zeros((ch,cw))
    
    #print (fh,fw,mh,mw,ch,cw,skip,deltah,deltaw)
    for shiftx in range(0,deltaw,skip):
        for shifty in range(0,deltah,skip):
            fixcut=fixed[shifty:shifty+mh,shiftx:shiftx+mw]
            corrmat[shifty/skip,shiftx/skip]=(fixcut*moved).sum()
           
    corrmat=corrmat/(mh*mw)
    
    return corrmat


#thread for making a cropped version of the big image... not very efficent    
class ImageCutThread(threading.Thread):
        def __init__(self, queue):
            threading.Thread.__init__(self)
            self.queue = queue        
        def run(self):
            while True:
                #grabs host from queue
                (filename,rect,i) = self.queue.get()
                image=Image.open(filename)
                image=image.crop(rect)
                (path,file)=os.path.split(filename)
                path=os.path.join(path,"previewstack")
                if not os.path.exists(path):
                    os.path.os.mkdir(path)
                cutfile=os.path.splitext(file)[0]+"stack%3d.tif"%i        
                cutfile=os.path.join(path,cutfile)
                image.save(cutfile)
                #signals to queue job is done
                self.queue.task_done()
                     
class MosaicImage():
    """A class for storing the a large mosaic imagein a matplotlib axis. Also contains functions for finding corresponding points
    in the larger mosaic image, and plotting informative graphs about that process in different axis"""
    def __init__(self,axis,one_axis,two_axis,corr_axis,imagefile,imagematrix,extent=None,flipVert=False,fullRes=False):
        """initialization function which will plot the imagematrix passed in and set the bounds according the bounds specified by extent
        
        keywords)
        axis)the matplotlib axis to plot the image into
        one_axis) the matplotlib axis to plot the cutout of the fixed point when using the corresponding point functionality
        two_axis) the matplotlib axis to plot the cutout of the point that should be moved when using the corresponding point functionality
        corr_axis) the matplotlib axis to plot out the matrix of correlation values found when using the corresponding point functionality
        imagefile) a string with the path of the file which contains the full resolution image that should be used when calculating the corresponding point funcationality
         currently the reading of the image is using PIL so the path specified must be an image which is PIL readable
        imagematrix) a numpy 2d matrix containing a low resolution varsion of the full resolution image, for the purposes of faster plotting/memory management
        extent) a list [minx,maxx,miny,maxy] of the corners of the image.  This will specify the scale of the image, and allow the corresponding point functionality
        to specify how much the movable point should be shifted in the units given by this extent.  If omitted the units will be in pixels and extent will default to
        [0,width,height,0].
       
        """
        #define the attributes of this class
        self.axis=axis
        self.one_axis=one_axis
        self.two_axis=two_axis
        self.corr_axis=corr_axis
        self.imagefile=imagefile
        self.flipVert=flipVert
        self.imagematrix=imagematrix
        
        #read in the full resolution height/width using PIL
        image = Image.open(imagefile)
        (self.originalwidth,self.originalheight)=image.size
        
        #if extent was not specified default to units of pixels with 0,0 in the upper left
        if extent==None:
            if flipVert:
                self.extent=[0,self.originalwidth,0,self.originalheight]
            else:
                self.extent=[0,self.originalwidth,self.originalheight,0]
        else:
            self.extent=extent
 
        #calculate the width of the image (calling it _um assuming its in units of microns)
        #from now on I will assume the units are in microns, though if they were in some other unit it would just carry through
        width_um=abs(self.extent[1]-self.extent[0])
        #height_um=self.extent[2]-self.extent[3]
        
        #calculate the pixels/micron of full resolution picture
        self.orig_um_per_pix=width_um/self.originalwidth    
        #calculate the pixels/micron of the downsampled matrix 
        (matrix_height,matrix_width)=imagematrix.shape
        self.matrix_scale=matrix_width/width_um
        
        #plot the image using paintImage
        self.paintImage()
        
        #initialize the images for the various subplots as None
        self.oneImage=None
        self.twoImage=None
        self.corrImage=None
        self.set_maxval(self.imagematrix.max(axis=None))
        self.fullRes=fullRes
        self.axis.set_title('Mosaic Image')
        
     
    def set_extent(self,extent):
        self.extent=extent
        width_um=abs(self.extent[1]-self.extent[0])
        #height_um=self.extent[2]-self.extent[3]
        
        #calculate the pixels/micron of full resolution picture
        self.orig_um_per_pix=width_um/self.originalwidth    
        #calculate the pixels/micron of the downsampled matrix 
        (matrix_height,matrix_width)=self.imagematrix.shape
        self.matrix_scale=matrix_width/width_um  
        self.Image.set_extent(self.extent)
        self.axis.set_xlim(self.extent[0],self.extent[1])
        self.axis.set_ylim(self.extent[2],self.extent[3])
        self.axis.set_xlabel('X Position (microns)')
        self.axis.set_ylabel('Y Position (microns)')
                  
    def paintImage(self):
        """plots self.imagematrix in self.axis using self.extent to define the boundaries"""
        self.Image=self.axis.imshow(self.imagematrix,cmap='gray',extent=self.extent)
        (minval,maxval)=self.Image.get_clim()
        self.maxvalue=maxval
        #self.axis.canvas.get_toolbar().slider.SetSelection(minval,self.maxvalue)
        self.axis.autoscale(False)
        self.axis.set_xlabel('X Position (pixels)')
        self.axis.set_ylabel('Y Position (pixels)')
        self.Image.set_clim(0,25000)
    
    def set_maxval(self,maxvalue):
        """set the maximum value in the image colormap"""
        self.maxvalue=maxvalue;
        self.repaint()
        
    def repaint(self):
        """sets the new clim for the Image using self.maxvalue as the new maximum value"""
        (minval,maxval)=self.Image.get_clim()
        self.Image.set_clim(minval,self.maxvalue)
        if self.oneImage!=None:
            self.oneImage.set_clim(minval,self.maxvalue)
        if self.twoImage!=None:
            self.twoImage.set_clim(minval,self.maxvalue)
    
    def paintImageCenter(self,cut,theaxis,xc=0,yc=0,skip=1,cmap='gray',scale=1):
        """paints an image and redefines the coordinates such that 0,0 is at the center
        
        keywords
        cut)the 2d numpy matrix with the image data
        the axis)the matplotlib axis to plot it in
        skip)the factor to rescale the axis by so that 1 entry in the cut, is equal to skip units on the axis (default=1)
        cmap)the colormap designation to use for the plot (default 'gray')
        
        """
        theaxis.cla()
        (h,w)=cut.shape
        dh=skip*1.0*(h-1)/2       
        dw=skip*1.0*(w-1)/2
        dh=dh*scale;
        dw=dw*scale;
        
        if self.extent[0]<self.extent[1]:
            left=xc-dw
            right=xc+dw
        else:
            left=xc+dw
            right=xc-dw
        if self.extent[2]<self.extent[3]:
            top=yc-dh
            bot=yc+dh
        else:
            top=yc+dh
            bot=yc-dh
            
        ext=[left,right,top,bot]

        image=theaxis.imshow(cut,cmap=cmap,extent=ext)
        
        if self.extent[0]<self.extent[1]:
            theaxis.set_xlim(xc-dw,xc+dw)
        else:
            theaxis.set_xlim(xc+dw,xc-dw)
        if self.extent[2]<self.extent[3]:
            theaxis.set_ylim(yc-dh,yc+dh)
        else:
            theaxis.set_ylim(yc+dh,yc-dh)
            
        theaxis.hold(True)

        return image 
    
    def updateImageCenter(self,cut,theimage,theaxis,xc=0,yc=0,skip=1,scale=1):
        """updates an image with a new image
        
        keywords
        cut) the 2d numpy matrix with the image data 
        theimage) the image to update
        theaxis) the axis that the image is in
        skip)the factor to rescale the axis by so that 1 entry in the cut, is equal to skip units on the axis (default=1)

        """
        (h,w)=cut.shape
        dh=skip*1.0*(h-1)/2       
        dw=skip*1.0*(w-1)/2
        dh=dh*scale;
        dw=dw*scale;
        theimage.set_array(cut)
        if self.extent[0]<self.extent[1]:
            left=xc-dw
            right=xc+dw
            theaxis.set_xlim(xc-dw,xc+dw)
        else:
            left=xc+dw
            right=xc-dw
            theaxis.set_xlim(xc+dw,xc-dw)
        if self.extent[2]<self.extent[3]:
            top=yc-dh
            bot=yc+dh
            theaxis.set_ylim(yc-dh,yc+dh)
        else:
            top=yc+dh
            bot=yc-dh
            theaxis.set_ylim(yc+dh,yc-dh)
        ext=[left,right,top,bot]
        theimage.set_extent(ext)
         
         
    def paintImageOne(self,cut,xy=(0,0),dxy_pix=(0,0),window=0):
        """paints an image in the self.one_axis axis, plotting a box of size 2*window+1 around that point
        
        keywords
        cut) the 2d numpy matrix with the image data
        dxy_pix) the center of the box to be drawn given as an (x,y) tuple
        window)the size of the box, where the height is 2*window+1
        
        """ 
        (xc,yc)=xy  
        (dx,dy)=dxy_pix
        dx=dx*self.orig_um_per_pix;
        dy=dy*self.orig_um_per_pix;
        #the size of the cutout box in microns
        boxsize_um=(2*window+1)*self.orig_um_per_pix;
        
        #if there is no image yet, create one and a box
        if self.oneImage==None:
            self.oneImage=self.paintImageCenter(cut, self.one_axis,xc=xc,yc=yc,scale=self.orig_um_per_pix)
            self.oneBox=CenterRectangle((xc+dx,yc+dy),width=50,height=50,edgecolor='r',linewidth=1.5,fill=False)
            self.one_axis.add_patch(self.oneBox)
            self.one_axis_center=Line2D([xc],[yc],marker='+',markersize=7,markeredgewidth=1.5,markeredgecolor='r')
            self.one_axis.add_line(self.one_axis_center) 
            self.one_axis.set_title('Point 1')
            self.one_axis.set_ylabel('Microns')
            self.one_axis.autoscale(False)
            self.oneImage.set_clim(0,self.maxvalue)     
        #if there is an image update it and the self.oneBox
        else:
            self.updateImageCenter(cut, self.oneImage, self.one_axis,xc=xc,yc=yc,scale=self.orig_um_per_pix)
            self.oneBox.set_center((dx+xc,dy+yc))
            self.oneBox.set_height(boxsize_um)
            self.oneBox.set_width(boxsize_um)
            self.one_axis_center.set_xdata([xc])
            self.one_axis_center.set_ydata([yc])
    
        
    def paintImageTwo(self,cut,xy=(0,0)):
        """paints an image in the self.two_axis, with 0,0 at the center cut=the 2d numpy"""
        #create or update appropriately
        (xc,yc)=xy
        if self.twoImage==None:
            self.twoImage=self.paintImageCenter(cut, self.two_axis,xc=xc,yc=yc,scale=self.orig_um_per_pix)
            self.two_axis_center=Line2D([xc],[yc],marker='+',markersize=7,markeredgewidth=1.5,markeredgecolor='r')
            self.two_axis.add_line(self.two_axis_center) 
            self.two_axis.set_title('Point 2')
            self.two_axis.set_ylabel('Pixels from point 2')
            self.two_axis.autoscale(False)
            self.twoImage.set_clim(0,self.maxvalue)
 
        else:
            self.updateImageCenter(cut, self.twoImage, self.two_axis,xc=xc,yc=yc,scale=self.orig_um_per_pix)
            self.two_axis_center.set_xdata([xc])
            self.two_axis_center.set_ydata([yc])
    
    def paintCorrImage(self,corrmat,dxy_pix,skip):
        """paints an image in the self.corr_axis, with 0,0 at the center and rescaled by skip, plotting a point at dxy_pix
        
        keywords)
        corrmat) the 2d numpy matrix with the image data
        dxy_pix) the offset in pixels from the center of the image to plot the point
        skip) the factor to rescale the axis by, so that when corrmat was produced by mycorrelate2d with a certain skip value, 
        the axis will be in units of pixels
        
        """
        #unpack the values
        (dx,dy)=dxy_pix
        #update or create new
        if self.corrImage==None:
            self.corrImage=self.paintImageCenter(corrmat, self.corr_axis,skip=skip,cmap='jet')             
            self.maxcorrPoint,=self.corr_axis.plot(dx,dy,'ro')

            self.colorbar=self.corr_axis.figure.colorbar(self.corrImage,shrink=.9)
            self.corr_axis.set_title('Cross Correlation')
            self.corr_axis.set_ylabel('Pixels shifted')
          
        else:
            self.updateImageCenter(corrmat, self.corrImage, self.corr_axis,skip=skip)
            self.maxcorrPoint.set_data(dx,dy)   
        #hard code the correlation maximum at .5
        self.corrImage.set_clim(0,.5)

    def convert_pos_to_orig_ind(self,x,y):
        """converts a position in original units (usually microns) to indices in the original full resolution image
        
        keywords)
        x)x position in microns
        y)y position in microns
        
        returns) (x_pix,y_pix) the indices in pixels of that location
        
        """
        print "(x,y) before"
        print (x,y)
        #calculate distance from left hand edge
        
        x=abs(x-self.extent[0])
        print self.flipVert
        print self.fullRes
        if self.flipVert and not self.fullRes:
            y=abs(self.extent[2]-y)
            print "alternative distance"
        else:
            print "regular distance"
            y=abs(self.extent[3]-y)
        #calculate distance from top of image
        print "x,y relative"
        print (x,y)
        x_pix=int(round(x/self.orig_um_per_pix))
        y_pix=int(round(y/self.orig_um_per_pix))
        print (x_pix,y_pix)
        return (x_pix,y_pix)
    
    def cutout_window(self,x,y,window):
        """returns a cutout of the original image at a certain location and size
        
        keywords)
        x)x position in microns
        y)y position in microns
        window) size of the patch to cutout, will cutout +/- window in both vertical and horizontal dimensions
        note.. behavior not well specified at edges, may crash
        
        function uses PIL to read in image and crop it appropriately
        returns) cut: a 2d numpy matrix containing the removed patch
        
        """
        print (x,y)
        (xpx,ypx)=self.convert_pos_to_orig_ind(x,y)
        print "centered at %d %d"%(x,y)
        if self.fullRes==True:
            cut=self.imagematrix[ypx-window:ypx+window,xpx-window:xpx+window]
            return cut
        else:
            image=Image.open(self.imagefile)
            if image.mode != 'P':
                image=image.convert('P')
            #print "image mode is %s"%image.mode
            image=image.crop([xpx-window,ypx-window,xpx+window,ypx+window])
            #enh = ImageEnhance.Contrast(image)
            #image=enh.enhance(1.5)   
            (width,height)=image.size
            #print "about to get data"
            thedata=image.getdata()
            #print "this is the data"
            #print thedata
            cut=np.reshape(np.array(thedata,np.dtype('uint8')),(height,width))
            if self.flipVert:
                cut=np.flipud(cut)
            return cut                      
    
    def cross_correlate_two_to_one(self,xy1,xy2,window=100,delta=75,skip=3):
        """take two points in the image, and calculate the 2d cross correlation function of the image around those two points
        
        keywords)
        xy1) a (x,y) tuple specifying point 1, the point that should be fixed
        xy2) a (x,y) tuple specifiying point 2, the point that should be moved
        window) the size of the patch to cutout (+/- window around the points) for calculating the correlation (default = 100 pixels)
        delta) the size of the maximal shift +/- delta from no shift to calculate
        skip) the number of integer pixels to skip over when calculating the correlation
        
        returns (one_cut,two_cut,corrmat)
        one_cut) the patch cutout around point 1
        two_cut) the patch cutout around point 2
        corrmat) the matrix of correlation values measured with 0,0 being a shift of -delta,-delta
        
        """
        (x1,y1)=xy1
        (x2,y2)=xy2
        one_cut=self.cutout_window(x1,y1,window+delta)
        two_cut=self.cutout_window(x2,y2,window)
        #return (target_cut,source_cut,mycorrelate2d(target_cut,source_cut,mode='valid'))
        return (one_cut,two_cut,mycorrelate2d(one_cut,two_cut,skip))
       
    def align_by_correlation(self,xy1,xy2,window=100,delta=75,skip=3):
        """take two points in the image, and calculate the 2d cross correlation function of the image around those two points
        plots the results in the appropriate axis, and returns the shift which aligns the two points given in microns
        
        keywords)
        xy1) a (x,y) tuple specifying point 1, the point that should be fixed
        xy2) a (x,y) tuple specifiying point 2, the point that should be moved
        window) the size of the patch to cutout (+/- window around the points) for calculating the correlation (default = 100 pixels)
        delta) the size of the maximal shift +/- delta from no shift to calculate
        skip) the number of integer pixels to skip over when calculating the correlation
        
        returns) (maxC,dxy_um)
        maxC)the maximal correlation measured
        dxy_um) the (x,y) tuple which contains the shift in microns necessary to align point xy2 with point xy1
        
        """
        #calculate the cutout patches and the correlation matrix
        (one_cut,two_cut,corrmat)=self.cross_correlate_two_to_one(xy1,xy2,window,delta,skip)
        #find the peak of the matrix
        maxind=corrmat.argmax()
        #determine the indices of that peak
        (max_i,max_j)=np.unravel_index(maxind,corrmat.shape)
        #calculate the shift for that index in pixels
        dy_pix=max_i*skip-delta
        dx_pix=max_j*skip-delta
        #convert those indices into microns
        dy_um=dy_pix*self.orig_um_per_pix
        dx_um=dx_pix*self.orig_um_per_pix
        #pack up the shifts into tuples
        dxy_pix=(dx_pix,dy_pix)
        dxy_um=(dx_um,dy_um)
        #calculate what the maximal correlation was
        corrval=corrmat.max()
        
        print "(correlation,(dx,dy))="
        print (corrval,dxy_pix)
        #paint the patch around the first point in its axis, with a box of size of the two_cut centered around where we found it
        self.paintImageOne(one_cut,xy=xy1,dxy_pix=dxy_pix, window=window)
        #paint the patch around the second point in its axis
        self.paintImageTwo(two_cut,xy=xy2)
        #paint the correlation matrix in its axis
        self.paintCorrImage(corrmat, dxy_pix,skip)
        return (corrmat.max(),dxy_um)
    
    def paintPointsOneTwo(self,xy1,xy2,window):
        (x1,y1)=xy1
        (x2,y2)=xy2
        one_cut=self.cutout_window(x1,y1,window)
        two_cut=self.cutout_window(x2,y2,window)
        self.paintImageOne(one_cut,xy1)
        #paint the patch around the second point in its axis
        self.paintImageTwo(two_cut,xy2)
        
    def make_preview_stack(self,xpos,ypos,width,height,directory):
        print "make a preview stack"
   
        hw_pix=int(round(width*.5/self.orig_um_per_pix))
        hh_pix=int(round(height*.5/self.orig_um_per_pix))
        queue = Queue.Queue()
          
        #spawn a pool of threads, and pass them queue instance 
        for i in range(4):
            t = ImageCutThread(queue)
            t.setDaemon(True)
            t.start()
              
        for i in range(len(self.mosaicArray.xpos)):
            (cx_pix,cy_pix)=self.convert_pos_to_ind(xpos[i],ypos[i])
            rect=[cx_pix-hw_pix,cy_pix-hh_pix,cx_pix+hw_pix,cy_pix+hh_pix]
            queue.put((self.imagefile,rect,i))
        queue.join()