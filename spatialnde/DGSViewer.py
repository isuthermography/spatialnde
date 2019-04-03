###!!!! BUG: The name of this module is inappropriate (has nothing specific to do with DGS) and should be renamed. 

import numpy as np
import cv2
import sys
from matplotlib import pyplot as plt
import matplotlib

import copy
import itertools
import os

#sys.path.insert(0, 'spatialnde/')
#import calibration_image
#import intrinsic_calibration

from spatialnde import calibration_image
from spatialnde import intrinsic_calibration

from limatix import dc_value

matplotlib.rcParams['toolbar'] = 'None'

COVERAGE_ON = 0
COVERAGE_CURRENT = 1
COVERAGE_OFF = 2

class VIEWER():


    def __init__(self, params, calib_doc, pattern_max):
        
        self.calib_doc = calib_doc # Calibration source
        self.pattern_max = pattern_max # Max pattern to search
        self.current_source = calibration_image.calibration_image(calibration_image.calibration_image_params(**params))
        self.current_searched = None # Store calibration image of searched pattern. So frame number, offset, etc. can be changed (updating self.current_source), but if intrinsic calibration is written, it uses the last searched pattern stored in self.current_pattern
        
        self.calibration_images = [] # Intrinsic calibration views written to file

        self.intrinsic_calibration = None 
        self.undistortion_processor = None

        self.figure = None
        self.axes = None
        self.image = None

        self.fig_text_dict = {} # Text displayed on viewer
        self.imout_files = [] # Saved snapshots

        self.vmin = 0.0
        self.vmax = 1.0

        self.pattern_num_x = None
        self.pattern_num_y = None
        self.reopen = False

        self.picked_points = []
        self.scatter_dict = {} # Pattern points found
        self.coverage_dict = {} # Grid point coverage

        self.toggle_coverage = COVERAGE_ON

        self.im_size = self.current_source.numpy_image.shape[:2]
        self.is_undistorted = False # If True, use undistorted images and pattern points
        
        self.associate_calib_source() # Bring in camera calibration from calibsource (or create new document) and create intrinsic calibration and undistortion processor


    def associate_calib_source(self):

        assert(isinstance(self.calib_doc, dc_value.hrefvalue))
        
        try:
            self.intrinsic_calibration = intrinsic_calibration.intrinsic_calibration.fromfile(self.calib_doc.absurl())
            self.undistortion_processor = intrinsic_calibration.undistortion_processor(self.intrinsic_calibration,self.intrinsic_calibration.sizex,self.intrinsic_calibration.sizey,self.intrinsic_calibration.sizex,self.intrinsic_calibration.sizey,interpolation=cv2.INTER_LINEAR,alpha=1.0)
            print("Pulling calibration information from {f}").format(f=self.calib_doc.absurl())
            

            for image in self.intrinsic_calibration.calibration_images:
                print("Found calibration data generated from {f} on {t}").format(f=image.params.input_href.absurl(), t=image.timestamp)
                self.calibration_images.append(image)
                pass

        except IOError:
            print("No calibration source... Creating calibration file {f}").format(f=self.calib_doc.absurl())
            pass


    def initialize_plot(self):
        self.figure = plt.figure(figsize=(14,10), dpi=80)
        self.axes = self.figure.add_axes([.25, .05, .8*self.im_size[0]/self.im_size[1], .9])
        self.image = self.axes.imshow(self.yield_image(), vmin=self.vmin, vmax=self.vmax, cmap='gray')
        self.axes.autoscale(False)

        self.axes.tick_params(labelbottom='off', labelleft='off')
       
        self.fig_text_dict.update({
            'current_frame': self.figure.text(.01,.8,"Frame number: {f}".format(f=self.current_source.params.frame_number+1), fontsize=10), 
            'units_per_intensity': self.figure.text(.01,.7,"Units per Intensity: {g} [log(deg. Kelvin)]".format(g=self.current_source.params.units_per_intensity), fontsize=10), 
            'offset': self.figure.text(.01,.65,"Offset: {b} [log(deg. Kelvin)]".format(b=self.current_source.params.offset), fontsize=10)})

        if not self.intrinsic_calibration is None:
            self.fig_text_dict.update({'reprojection_error': self.figure.text(.01,.60,"Reprojection error: {r} [pixels]".format(r=self.intrinsic_calibration.reprojection_error), fontsize=10)})
        else:
            self.fig_text_dict.update({'reprojection_error': self.figure.text(.01,.60,"Reprojection error: ? [pixels]")})


        self.fig_text_dict.update({
            'command_<>': self.figure.text(.01,.55,"</>: Frame down/up", fontsize=10), 
            'command_lr': self.figure.text(.01,.525,"left/right arrow: Offset down/up", fontsize=10), 
            'command_ud': self.figure.text(.01,.5,"up/down arrow: Units per intensity down/up", fontsize=10), 
            'command_t': self.figure.text(.01,.45,"t: Toggle coverage plotting", fontsize=10), 
            'command_c': self.figure.text(.01,.425,"c: Search calibration patterns", fontsize=10), 
            'command_u': self.figure.text(.01,.4,"u: Undistort using available calibration", fontsize=10),
            'command_r': self.figure.text(.01, .375, 'r: Reset found pattern points', fontsize=10),
            'command_s': self.figure.text(.01,.35,"s: Snapshot of current screen", fontsize=10), 
            'command_w': self.figure.text(.01,.325,"w: Append calibation to file", fontsize=10), 
            'command_esc': self.figure.text(.01,.3,"<Esc>: Leave interactive mode", fontsize=10),
            'command_m': self.figure.text(.01,.25,"m: Calibrate based on points chosen with mouse", fontsize=10),
            'desc_m': self.figure.text(.01,.225,"TAKE NOTE OF GRID SIZE (ROWS and COLS) BEFORE PRESSING 'm'.", fontsize=10),
            'command_lclick': self.figure.text(.01,.2,"left click: Manually chose grid point.", fontsize=10),
            'desc_lclick1': self.figure.text(.01,.175,"CHOOSE POINTS STARTING AT TOP LEFT AND", fontsize=10),
            'desc_lclick2': self.figure.text(.01,.15,"MOVING HORIZONTALLY IN RASTER SCAN.", fontsize=10),
            'command_rclick': self.figure.text(.01,.125,"right click: Remove last grid point chosen.", fontsize=10),
            })

        self.clear_coverage()
        self.plot_coverage()
        self.plot_picked(self.picked_points)


    def enter_interactive_mode(self):
        cid = self.figure.canvas.mpl_connect('key_press_event', self.on_key)
        cid2 = self.figure.canvas.mpl_connect('button_press_event', self.on_click)
        plt.show()
        
        return self.reopen

    
    def update_axes(self, im_data=None):

        if im_data is None:
            im_data = self.yield_image()

        self.image.set_data(im_data)
        self.image.set_clim(vmin=self.vmin, vmax=self.vmax)

        self.figure.texts.remove(self.fig_text_dict.get('current_frame'))
        self.figure.texts.remove(self.fig_text_dict.get('units_per_intensity'))
        self.figure.texts.remove(self.fig_text_dict.get('offset'))
        self.figure.texts.remove(self.fig_text_dict.get('reprojection_error'))

        self.fig_text_dict.update({
            'current_frame': self.figure.text(.01,.8,"Frame number: {f}".format(f=self.current_source.params.frame_number+1), size='small'), 
            'units_per_intensity': self.figure.text(.01,.7,"Units per Intensity: {g} [log(deg. Kelvin)]".format(g=self.current_source.params.units_per_intensity), size='small'), 
            'offset': self.figure.text(.01,.65,"Offset: {b} [log(deg. Kelvin)]".format(b=self.current_source.params.offset), size='small')})

        if not self.intrinsic_calibration is None:
            self.fig_text_dict.update({'reprojection_error': self.figure.text(.01,.60,"Reprojection error: {r} [pixels]".format(r=self.intrinsic_calibration.reprojection_error), size='small')})
        else:
            self.fig_text_dict.update({'reprojection_error': self.figure.text(.01,.60,"Reprojection error: ? [pixels]")})

        plt.draw()


    def plot_found(self):
        
        if not self.current_searched is None:

            if self.is_undistorted:
                image_points = cv2.undistortPoints(self.current_searched.image_points, self.intrinsic_calibration.cameramatrix, self.intrinsic_calibration.distcoeffs, P=self.undistortion_processor.newcameramatrix)
            else:
                image_points = self.current_searched.image_points
                pass

            for pt in image_points.reshape(-1,2):
                self.scatter_dict[tuple(pt)] = self.axes.scatter(x=pt[0], y=pt[1], c='w', s=40)
            
            plt.draw()

        else:
            return

    def plot_picked(self, points):

        for pt in points:
            self.scatter_dict.update({tuple(pt): self.axes.scatter(x=pt[0], y=pt[1], c='w', s=40)})
            pass
            
        plt.draw()
        return
    
    def clear_picked(self, points):

        for pt in points:
            val = self.scatter_dict.pop(tuple(pt))
            val.remove()
            pass

        plt.draw()
        return

    def plot_coverage(self):
        img_points = []
        
        if self.toggle_coverage == COVERAGE_ON:
            for calibration_image in self.calibration_images:
                for pt in calibration_image.image_points.reshape(-1,2):
                    img_points.append((pt[0],pt[1]))
                    pass
                pass
            pass
        
        if self.toggle_coverage == COVERAGE_ON or self.toggle_coverage == COVERAGE_CURRENT:
            if not self.current_searched is None:
                for pt in self.current_searched.image_points.reshape(-1,2):
                    img_points.append((pt[0],pt[1]))
            else:
                pass
            pass
        
        if len(img_points) > 0:
            img_points = np.swapaxes(np.array([img_points]), 0, 1)
            if self.is_undistorted:
                undistorted_image_points = cv2.undistortPoints(img_points, self.intrinsic_calibration.cameramatrix, self.intrinsic_calibration.distcoeffs, P=self.undistortion_processor.newcameramatrix)
            else:
                undistorted_image_points = img_points
                pass

            for pt in undistorted_image_points.reshape(-1,2):
                self.coverage_dict[tuple(pt)] = self.axes.scatter(x=pt[0], y=pt[1], edgecolor='r', s=200, facecolor='none', lw=2)
                pass
            pass
            
        plt.draw()
        
        return


    def clear_found(self):

        for key, val in self.scatter_dict.items():
            val.remove()

        self.scatter_dict = {}

        plt.draw()


    def clear_coverage(self):

        for key, val in self.coverage_dict.items():
            val.remove()

        self.coverage_dict = {}

        plt.draw()

    def yield_image(self):
        im_data = self.current_source.numpy_image
        scaled = calibration_image.normalize_float(im_data, self.current_source.params.units_per_intensity, self.current_source.params.offset)
        
        if self.is_undistorted:
            return self.undistortion_processor.undistort_numpy_image(scaled)
        else:
            return scaled


    def search_pattern(self):

        print("Searching for ideal pattern...")

        maxx = self.pattern_max[0] # Based on the maximum pattern size visible
        maxy = self.pattern_max[1]
        
        if self.current_source.params.pattern_type == 'chessboard':
        	minx = 3
        	best_pattern = (3,3) # Store the best pattern size
       		pass
        else:
        	minx = 2
        	best_pattern = (1,1) # Store the best pattern size
        	pass

        for (pattern_num_x, pattern_num_y) in ((p,p) for p in range(minx,min(maxx+1, maxy+1))):
            # Start by checking square grids

            new_params = self.current_source.params.copy_and_update(pattern_num_x=pattern_num_x, pattern_num_y=pattern_num_y)
            self.current_source = calibration_image.calibration_image(new_params)

            if not self.current_source.image_points is None:
                print("pattern size ({p},{q}) found".format(p=pattern_num_x,q=pattern_num_y))
                best_pattern = (pattern_num_x, pattern_num_y)
            else:
                print("Pattern size ({p},{q}) not found".format(p=pattern_num_x,q=pattern_num_y))
                if pattern_num_x == pattern_num_y:
                    # Termination criteria, the best pattern has already been found
                    # If the square pattern (x,x) cannot be found, then it will not find (x+1,x+1) or (x,x+1)
                    break
                pass
        
        for (pattern_num_x, pattern_num_y) in itertools.chain(itertools.product([best_pattern[0]],range(best_pattern[0]+1,maxy+1)),itertools.product(range(best_pattern[0]+1,maxx+1),[best_pattern[0]])): 
            # Check rectangular grids            
 
            new_params = self.current_source.params.copy_and_update(pattern_num_x=pattern_num_x, pattern_num_y=pattern_num_y)
            self.current_source = calibration_image.calibration_image(new_params)

            if not self.current_source.image_points is None:
                print("pattern size ({p},{q}) found".format(p=pattern_num_x,q=pattern_num_y))
                if pattern_num_x*pattern_num_y > best_pattern[0]*best_pattern[1]:
                    best_pattern = (pattern_num_x, pattern_num_y)
                    pass
            else:
                print("Pattern size ({p},{q}) not found".format(p=pattern_num_x,q=pattern_num_y))

        if best_pattern == (1,1):
            return False

        (pattern_num_x, pattern_num_y) = best_pattern

        new_params = self.current_source.params.copy_and_update(pattern_num_x=pattern_num_x, pattern_num_y=pattern_num_y)
        self.current_source = calibration_image.calibration_image(new_params)

        self.current_searched = copy.copy(self.current_source)

        return True


    def set_offset(self, offset):

        """Set image offset, update limits"""
        new_params = self.current_source.params.copy_and_update(offset=offset)
        self.current_source = calibration_image.calibration_image(new_params)
        pass


    def set_units_per_intensity(self, units_per_intensity):

        """Set units per intensity, update limites"""

        new_params = self.current_source.params.copy_and_update(units_per_intensity=units_per_intensity)
        self.current_source = calibration_image.calibration_image(new_params)
        pass


    def offset_step(self):

        """Offset increment"""

        return self.current_source.params.units_per_intensity*0.125


    def units_per_intensity_step(self):
    
        """Units per intensity increment"""    
        # Scope gives units per intensity .1, .2, .5, 1, 2, 5, 10, 20, 50, 100...
        
        return 2.0


    def save_fig(self, imout_filename=None):

        """Save .png image and append filename to imout_files"""
        
        # Filename strings
        frame_str = 'F' + ("%d" % (self.current_source.params.frame_number+1))
        units_per_intensity_str = 'G' + str(self.current_source.params.units_per_intensity).replace('.','p')
        offset_str = 'B' + str(self.current_source.params.offset).replace('.','p').replace('-','n')

        # Save full snapshot
        if imout_filename is not None:
            imout_href=dc_value.hrefvalue(quote("%s%s%s%s.png" % (posixpath.split(self.current_source.params.input_href.get_bare_unquoted_filename())[0],frame_str,units_per_intensity_str,offset_str)),contexthref=self.current_source.params.input_href)
            imout_filename=imout_href.getpath()
            pass
        print("Saving snapshot {s}".format(s=imout_filename))
        
        self.figure.savefig(imout_filename, dpi=600)
        self.imout_files.append(imout_filename)
        
        # Save cropped snapshot
        imout_filename_cropped = imout_filename[:-4] + 'cropped.png'
        print("Saving cropped snapshot {s}".format(s=imout_filename_cropped))
        
        cropped_extent = self.axes.get_window_extent().transformed(self.figure.dpi_scale_trans.inverted())
        self.figure.savefig(imout_filename_cropped, bbox_inches=cropped_extent)
        self.imout_files.append(imout_filename_cropped)


    def on_click(self, event):

        if event.button == 1:
            x = event.xdata
            y = event.ydata

            if x is not None and y is not None:
                self.picked_points.append((x, y))
                self.plot_picked([(x,y)])
                pass

            pass

        if event.button == 3:
            
            # This can be used to remove points by mouse location
            #dist = [np.sqrt((event.xdata-pt[0])**2 + (event.ydata-pt[1])**2) for pt in self.picked_points]
            #x = event.xdata
            #y = event.ydata
            #pt = self.picked_points.pop(dist.index(min(dist)))

			# Here we are removing the most recent chosen point
            pt = self.picked_points.pop(-1)
            self.clear_picked([pt])
            pass

        return 


    def on_key(self, event):

        # Event handling

        if event.key.lower() == '.':
            # Increment frame number
            new_frame = self.current_source.params.frame_number + 1
            new_params = self.current_source.params.copy_and_update(frame_number=new_frame)
            self.current_source = calibration_image.calibration_image(new_params)

        if event.key.lower() == ',':
            # Decrement frame number
            if self.current_source.params.frame_number -1 < 0:
                pass
            else:
                new_frame = self.current_source.params.frame_number - 1
                new_params = self.current_source.params.copy_and_update(frame_number=new_frame)
                self.current_source = calibration_image.calibration_image(new_params)

        if event.key.lower() == 'up':
            # Increment units per intensity
            newval = self.current_source.params.units_per_intensity/self.units_per_intensity_step()
            self.set_units_per_intensity(newval)

        if event.key.lower() == 'down':
            # Decrement units per intensity
            newval = self.current_source.params.units_per_intensity*self.units_per_intensity_step()
            self.set_units_per_intensity(newval)

        if event.key.lower() == 'right':
            # Increment offset
            self.set_offset(self.current_source.params.offset + self.offset_step())

        if event.key.lower() == 'left':
            # Decrement offset
            self.set_offset(self.current_source.params.offset - self.offset_step())

        if event.key.lower() == 'u':
            # Undistort using current intrinsic calibration
            if not self.undistortion_processor:
                print("No calibration data from source...")
                pass
            else:
                if self.is_undistorted:
                    self.is_undistorted = False
                    pass
                else:
                    self.is_undistorted = True
                    pass

                self.clear_found()
                self.plot_found()
                self.clear_coverage()
                self.plot_coverage()

            pass

        if event.key.lower() == 'c':
            # Search pattern and update intrinsic calibration
            nx = self.im_size[1]
            ny = self.im_size[0]

            self.clear_found()
            self.clear_coverage()
            found = self.search_pattern()
            
            if found:
                self.intrinsic_calibration = intrinsic_calibration.intrinsic_calibration(self.calibration_images + [self.current_searched],nx,ny)
                self.undistortion_processor = intrinsic_calibration.undistortion_processor(self.intrinsic_calibration,nx,ny,nx,ny,interpolation=cv2.INTER_LINEAR,alpha=1.0)
                self.plot_coverage()
                self.plot_found()

        if event.key.lower() == 'r':
            # Clear intrinsic calibration not written to file
            self.clear_picked(self.picked_points)
            self.picked_points = []

            self.clear_coverage()
            self.clear_found()
            self.current_searched = None
           
            if not len(self.calibration_images) > 0:
                self.intrinsic_calibration=None
                self.undistortion_processor=None
                self.is_undistorted=False
            else:
                nx = self.im_size[1]
                ny = self.im_size[0]

                self.intrinsic_calibration = intrinsic_calibration.intrinsic_calibration(self.calibration_images,nx,ny)
                self.undistortion_processor = intrinsic_calibration.undistortion_processor(self.intrinsic_calibration,nx,ny,nx,ny,interpolation=cv2.INTER_LINEAR,alpha=1.0)

            self.plot_coverage()

            pass

        if event.key.lower() == 'w':
            # Write current view to file
            if self.intrinsic_calibration is None:
                print("No calibration data available for writing...")
                pass
            else:                
                if not self.current_source.image_points is None:
                    print("Appending calibration data to {f}".format(f=self.calib_doc.absurl()))
                    self.calibration_images.append(self.current_searched)
                    self.intrinsic_calibration.savetofile(self.calib_doc)
                    pass
                else:
                    print("No calibration data available for writing!")


        if event.key.lower() == 's':
            # Save snapshot of view
            self.save_fig()
        
        if event.key.lower() == 'm':

            self.reopen = True

            plt.close()

            self.pattern_num_x=int(raw_input("%s: " % ('Enter number of grid columns')).strip())
            self.pattern_num_y=int(raw_input("%s: " % ('Enter number of grid rows')).strip())
            pass
        
        if event.key.lower() == 't':

            if self.toggle_coverage + 1 <= COVERAGE_OFF:
                self.toggle_coverage += 1
                pass
            else:
                self.toggle_coverage = COVERAGE_ON
                pass

            self.clear_coverage()
            self.plot_coverage()

        if event.key.lower() == 'escape':
            # Exit viewer
            plt.close()

        self.update_axes()

        pass






