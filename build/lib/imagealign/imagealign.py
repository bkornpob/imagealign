# Kornpob Bhirombhakdi
# kbhirombhakdi@stsci.edu

import copy
import pandas as pd
from recentroid.recentroid import Recentroid
import scipy.ndimage
from astropy.io import fits
import numpy as np

class ImageAlign:
    """
    ImageAlign is a class to facilitate aligning two images by matching centroids of sources in pixel (x,y) coordinate.
    + How to use:
      - Instantiate ImageAlign object.
      - Call public methods to compute outputs.
      - Call save() to save outputs.
    + Attributes:
      - self.image1, self.image2 = 2D array of images
      - self.source1, self.source2 = a dict with each key as source name and value as a dict with keys = {'X','Y'} and values as pixel coordinate (x,y) for each corresdponding key. 
        > source1 and source2 must run parallelly.
        > source1 and source2 are recorded as they are regardless of index scheme.
        > They will be changed to zero-indexing scheme after any calling method such as compute_to_zero_index or compute_recentroid.
      - self.source_zero_indexing = boolean. True if the source1 and source2 as inputs are in zero-indexing scheme (i.e., python standard). False if they are one-indexing scheme (e.g., DS9 standard).
        > If any call changes source1 and source2 to zero-indexing scheme, source_zero_indexing would also be updated to True.
      - self.source_shift = simply source1 - source2 for the shift in (x,y) coordinate available after running self.compute_shift()
      - self.shift = shift values in pixel unit available after running self.compute_shift() 
    + Methods (public):
      - self.compute_to_zero_index() = simply convert source1 and source from one-indexing to zero-indexing by (x,y) -= (1.,1.).
        > This also updates self.source_zero_indexing to True.
        > Calling compute_recentroid() also implements compute_to_zero_index().
      - self.compute_recentroid(box_size,centroid_func) = use the initial source1 and source2, and re-centroid. New centroids will replace source1 and source2, and source_zero_indexing sets to True.
        > Package Recentroid (pip install recentroid) is used in this step.
        > box_size = a search area for a new centroid specified by a square of size box_size centered on the initial centroid.
        > centroid_func = centroid function. See the package Recentroid (pip install recentroid) for the support.
          + If None, centroid_func = photutils.centroid_2dg
      - self.compute_shift() = use source1 and source2 to compute shifts in pixel unit.
        > self.source_shift = source1 - source2
        > self.shift = mean and std from source_shift
      - self.make_shift_image(shift_more, order, fill_value) = use self.shift and shift_more to shift image2. Output as self.image2_shifted.
        > Package scipy.ndimage.shift is used in this process.
        > shift_more = a tuple (x,y) specifying arbitrary shift in addition to self.shift: image2_shifted = image2 + self.shift + shift_more.
        > order = integer for spline interpolation
        > fill_value = float to be assigned to any none finite value in image2 before shift.
        > self.image2_shifted is the output.
      - self.save(container) = save all outputs as name convention ./savefolder/saveprefix_savesuffix.extension where savefolder and saveprefix are specified by Container class.
        > container = imagealign.container.Container class
        > Other switches include save_zero_indexing, save_recentroid, save_shift, save_shifted_image, save_shifted_image_overwrite. 
        > Outputs include 
          + recentroid1.reg, recentroid2.reg = new centroids from source1 and source2 respectively. Each file has two columns of (x,y) coordinate with space separation.
          + shift.csv = self.shift in csv format
          + shifted.fits = fits file with self.image2_shifted in EXT1 with EXT0 as an empty primary HDU.
    + Methods (private):
      - self._to_zero_index() = converting one-indexing to zero-indexing schemes, and also updating source1, source2, and source_zero_indexing = True.
        > simply performing x_zero_index = x_one_index - 1., and vice versa for y.
    """
    def __init__(self,image1,image2,source1,source2,source_zero_indexing=True):
        self.image1 = image1
        self.image2 = image2
        self.source1 = source1
        self.source2 = source2
        self.source_zero_indexing = source_zero_indexing
    def compute_recentroid(self,box_size=10,centroid_func=None):
        # update to zero indexing scheme
        if not self.source_zero_indexing:
            self._to_zero_index()
        #####
        t1 = Recentroid(self.source1,self.image1,box_size,centroid_func)
        t1.compute()
        t2 = Recentroid(self.source2,self.image2,box_size,centroid_func)
        t2.compute()
        #####
        self.source1 = t1.source_table_new.T.to_dict()
        self.source2 = t2.source_table_new.T.to_dict()
        print('Recentroid source1 and source2.')
    def compute_shift(self):
        # update to zero indexing scheme
        if not self.source_zero_indexing:
            self._to_zero_index()
        ##### self.source_shift
        source1 = pd.DataFrame(self.source1).T
        source2 = pd.DataFrame(self.source2).T
        t = source1 - source2
        self.source_shift = copy.deepcopy(t)
        print('Register self.source_shift')
        ##### self.shift
        shift = {}
        shift['X'] = {'MEAN':self.source_shift['X'].values.mean(), 'STD':self.source_shift['X'].values.std()}
        shift['Y'] = {'MEAN':self.source_shift['Y'].values.mean(), 'STD':self.source_shift['Y'].values.std()}
        self.shift = shift
        print('Register self.shift')
    def make_shifted_image(self,shift_more=(0.,0.),order=3,fill_value=0.):
        # register shift_more
        self.shift['X']['SHIFT_MORE'] = shift_more[0]
        self.shift['Y']['SHIFT_MORE'] = shift_more[1]
        print('Register SHIFT_MORE to self.shift')
        #####
        shift_x = self.shift['X']['MEAN'] + shift_more[0]
        shift_y = self.shift['Y']['MEAN'] + shift_more[1]
        shift_total = np.array([shift_x,shift_y])
        #####
        image2 = self.image2.copy()
        m = ~np.isfinite(image2)
        image2[m] = fill_value
        image2_shifted = scipy.ndimage.shift(input=image2,shift=shift_total,order=order)
        self.image2_shifted = image2_shifted.copy()
        print('Register self.image2_shifted')
    def save(self,container=None,save_zero_indexing=True,save_recentroid=True,save_shift=True,save_shifted_image=True,save_shifted_image_overwrite=True):
        if container is None:
            raise ValueError('container must be specified. See imagealign.container.Container class.')
        #####
        if save_recentroid:
            string1 = './{0}/{1}_recentroid1.reg'.format(container.data['savefolder'],container.data['saveprefix'])
            string2 = './{0}/{1}_recentroid2.reg'.format(container.data['savefolder'],container.data['saveprefix'])
            t1 = pd.DataFrame(self.source1).T
            t2 = pd.DataFrame(self.source2).T
            if not save_zero_indexing:
                t1['X'] += 1.
                t1['Y'] += 1.
                t2['X'] += 1.
                t2['Y'] += 1.
            t1.to_csv(string1,sep=' ',index=False,header=False)
            t2.to_csv(string2,sep=' ',index=False,header=False)
            print('Save {0}'.format(string1))
            print('Save {0}'.format(string2))
        #####
        if save_shift:
            string = './{0}/{1}_shift.csv'.format(container.data['savefolder'],container.data['saveprefix'])
            t = pd.DataFrame(self.shift).T
            t.to_csv(string)
            print('Save {0}'.format(string))
        #####
        if save_shifted_image:
            string = './{0}/{1}_shifted.fits'.format(container.data['savefolder'],container.data['saveprefix'])
            phdu = fits.PrimaryHDU()
            ihdu = fits.ImageHDU()
            hdul = fits.HDUList([phdu,ihdu])
            hdul[1].data = self.image2_shifted.copy()
            hdul.writeto(string,overwrite=save_shifted_image_overwrite)
            print('Save {0}'.format(string))
    def _to_zero_index(self):
        if self.source_zero_indexing:
            print('source_zero_indexing == True. Terminate compute_to_zero_index.')
            return
        #####
        source1 = copy.deepcopy(self.source1)
        source1 = pd.DataFrame(source1).T
        source1['X'] -= 1.
        source1['Y'] -= 1.
        self.source1 = source1.T.to_dict()
        print('Update source1 to zero-indexing scheme.')
        #####
        source2 = copy.deepcopy(self.source2)
        source2 = pd.DataFrame(source2).T
        source2['X'] -= 1.
        source2['Y'] -= 1.
        self.source2 = source2.T.to_dict()
        print('Update source2 to zero-indexing scheme.')
        #####
        self.source_zero_indexing = True
        print('Update source_zero_indexing to True.')
    