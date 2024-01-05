4# -*- coding: utf-8 -*-
"""
#     #    #    ####### ######                           ###        #####  
##   ##   # #   #     # #     # #    #  ####  #####     #   #      #     # 
# # # #  #   #  #     # #     # #    # #    #   #      #     #           # 
#  #  # #     # #     # ######  ###### #    #   #      #     #      #####  
#     # ####### #     # #       #    # #    #   #      #     # ### #       
#     # #     # #     # #       #    # #    #   #       #   #  ### #       
#     # #     # ####### #       #    #  ####    #        ###   ### ####### 

Welcome to MAOPhot 0.2, a PSF Photometry tool using Astropy and Photutils.psf

    0.2 Revision
    - Added Non Iterative PSF Photometry option which uses class
       PSFPhotometry 
    - check_star boolean in DataFrame is bool type throughout, not str
    - Upgraded to Photutils 1.10.0 and Panda 2.1.4, etc.
    - changed "Close" button in Settings window to "OK"
    - report JD is DATE-OBS + EXPOSURE/2 (like AAVSO report)
    - minor changes to Notes in report to reflect current AAVSO report
    - added AMASS to AAVSO report

MAOPhot calculates stellar magnitudes from 2 dimensional digital photographs. 
It produces an extended AAVSO (American Association of Variable Star Observers)
 report which can be submitted to the AAVSO using the online tool WebObs 
 (http://www.aavso.org/webobs).

There are many photometry measuring programs available such as VPhot 
(http://www.aavso.org/vphot) and AstroImageJ (University of Louisville). VPhot
uses the aperture photometry method. 

MAOPhot uses the PSF Photometry method exclusively. PSF (point spread function)
modeling is well suited for measuring stellar magnitudes in crowded  fields,
or the magnitude of a star that has a close companion, e.g., Z Tau.
(See https://www.aavso.org/lpv-double-trouble-campaign-0)

MAOPhot is written in Python. It uses many Python 'astropy' 
(https://www.astropy.org/) libraries. The astropy package contains key 
functionality and common tools for performing astronomy and astrophysics
 with Python. Included in the package is Photutils.psf.  See "PSF Photometry" 
 (https://photutils.readthedocs.io/en/stable/psf.html) which describes many of 
 the classes and methods used in MAOPhot 

This program was derived from MetroPSF by Maxym Usatov.  
It has been redesigned for AAVSO reporting only and includes, but not limited 
to the following enhancements:

- Generation of Effective PSF model, and ability to create a ‘rejection list’
- option to use an Integrated Gaussian PRF (Pixel Response Function) as model 
- PSF Photometry using an iterative algorithm to perform point spread function 
photometry in crowded fields
- Photometry using an ensemble of comparison stars or a single comp star
- Generation of Two-Color Photometry (B, V), (V, R) or (V, I), and Single Image
 Photometry reports in AAVSO extended format 
- Use of telescope Transformation Coefficients (needed for Two Color Photometry)
- Image display shows comp star AAVSO label number and name of any found VSX 
objects in image field
- Intermediate results are saved as .csv files 
- User can optionally enter a AAVSO Chart ID when retrieving comparison star
 data
- User can specify check star and list of comp stars to use



    More about Single Image Photometry

        Single Image Photometry does not utilize the Transformation
        coefficients. Simple differential photometry is used.
        Only a single comp star is used (which must be the case if the AAVSO 
        VPhot tool, 'Transform Applier' is to be used).
        
        Var mag = Var IM - Comp IM + Comp (known) mag
        Check star mag = Check star IM - Comp IM + Comp (known) mag

        where IM is the instrumental magnitude 
            -2.5 * np.log10(self.results_tab_df["flux_fit"] / exptime)
        self.results_tab_df["flux_fit"] represents the fitted flux for the 
        star (in that row)

        

    More about Two Color Photometry

        MAOPhot mimics VPhot's "Two Color Photometry" (Usually B and V)
            
        See spreadsheet: ProcessingMaoImages_202281V1117Her.xlsx 

        It includes formulas to generate "two color photometry". 
        The method generate_aavso_report_2color() uses the same formula which
		uses the telescope Transformation Coefficients.

        The following coefficients are used and are added to the list of "pink"
        fields. (Pink fields required for report generation)
            Tbv
            Tv_bv
            Tb_bv
            
        These are associated with a particular telescope, eg. MAO
        
        Required data (Pink fields):
            Tbv
            Tv_bv
            Tb_bv
            "Select Comp Stars (AAVSO Label)"
            "Use Check star (AAVSO Label)"
            ...
            etc.
                       
            
        The following formulas are calculated for each comp star and 
		var = check star;
        then calculated again for each comp star and var = "Object Name". 
		E.g. V1117 Her
        
        Δv = vvar - vcomp
        Δ(B-V) = Tbv * Δ(b-v)
        Vvar = Δv + Tv_bv * Δ(B-V) + Vcomp
        
        where:
        Δv  = IM of variable  - IM comp
        Vcomp = published V-magnitude of comp
        Δ(B-V) = Tbv * difference between standard color of var and standard 
		color of comp
        
        To calculate these formulas and then generate a report, two (2) sets
        of results_tab_df in csv format must exist, one for B and one for V and
        must have been derived from the B and V images of the Var under 
		analysis.
        
        When generate_aavso_report_2color is called, by the menu item,
        'Two Color Photometry->Two Color Photometry (B,V)', the user will 
        be asked to specify the 2 aformentioned csv files.
        
        From these files/Panda databases, the formulas are calculated, and 
		results are displayed. 


        
        Error Estimation :
        MAOPhot mimics VPhot when calculating error estimation. 
        From VPhot documentation:
                In an ensemble solution with more than two comp stars, 
            the magnitude is estimated as the average of the individual 
            comp stars estimate [of the check star], and the error is taken as 
            the standard deviation of this sample. 
        
            If one or two comp stars are used, the error estimate is
            based on the SNR of each measurement (the target measurement 
            and the comp stars measurements). The standard error of a 
            measurement is defined as 2.5 * np.log10(1 + 1 / SNR)
            [The errors are added in quadruture.]

        --
        get_comparison_stars
        
        Given the specified Check Star in field "Use Check Star (AAVSO label)"
        and the list of comp stars in "Select Comp Stars (AAVSO Label)", these
        stars are gotten from the AAVSO and each has an additional attribute,
        is_check_star, added to the results_tab_df DB. 
        


    MAOPhot derived from original MetroPSF (C) Copyright 2021, Maxym Usatov 
    <maxim.usatov@bcsatellite.net> Refer to metropsf.pdf for license information.

    This research made use of Photutils, an Astropy package for
    detection and photometry of astronomical sources (Bradley et al. 20XX).

"""
from ast import Assert
from astropy.stats import SigmaClip
from astropy.stats import sigma_clipped_stats
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.table import Table
from astropy.modeling.fitting import LevMarLSQFitter, SLSQPLSQFitter, SimplexLSQFitter
from astropy.nddata import NDData
from astropy.nddata import Cutout2D
from astropy.visualization import SqrtStretch, LogStretch, AsinhStretch, simple_norm
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.io import fits
from photutils.background import Background2D
from photutils.background import MedianBackground
from photutils.background import MADStdBackgroundRMS
from photutils.background import LocalBackground
from photutils.detection import find_peaks
from photutils.psf import IterativePSFPhotometry, PSFPhotometry
from photutils.psf import IntegratedGaussianPRF, SourceGrouper
from photutils.detection import IRAFStarFinder
from photutils.psf import extract_stars
from photutils.psf import EPSFBuilder
from tkinter import filedialog as fd
from tkinter import ttk
import tkinter as tk
from mpl_toolkits.mplot3d import Axes3D
import math
from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import matplotlib
from astroquery.vizier import Vizier
from astroquery.astrometry_net import AstrometryNet
from tqdm import tqdm
from PIL import Image, ImageTk, ImageMath
import pandas as pd
import logging
import sys
import requests
import csv
import os.path
import numpy as np
import warnings
import datetime

warnings.filterwarnings("ignore")
matplotlib.use("TkAgg")


# Photometry

# No imports beyond this line.


def save_background_image(stretch_min, stretch_max, zoom_level, image_data):
    global FITS_minimum
    global FITS_maximum
    global background_image
    background_image = Image.fromarray(image_data)
    width, height = background_image.size
    new_size = (int(width * zoom_level), int(height * zoom_level))
    background_image = background_image.resize(new_size, Image.LANCZOS)
    background_image = ImageMath.eval("(a + " + str(stretch_min / 100 * FITS_maximum) +
                                      ") * 255 / " + str(stretch_max / 100 * FITS_maximum), a=background_image)
    background_image = ImageMath.eval("convert(a, 'L')", a=background_image)
    background_image.save('background.jpg')


def save_image(stretch_min, stretch_max, zoom_level, image_data, filename):
    global FITS_minimum
    global FITS_maximum
    _image = Image.fromarray(image_data)
    width, height = _image.size
    new_size = (int(width * zoom_level), int(height * zoom_level))
    _image = _image.resize(new_size, Image.LANCZOS)
    _image = ImageMath.eval("(a + " + str(stretch_min / 100 * FITS_maximum) +
                            ") * 255 / " +
                            str(stretch_max / 100 * FITS_maximum),
                            a=background_image)
    _image = ImageMath.eval("convert(a, 'L')", a=background_image)
    _image.save(filename)


def generate_FITS_thumbnail(stretch_min, stretch_max, zoom_level,
                            stretching_stringvar):
    global generated_image
    global image_data
    global FITS_minimum
    global FITS_maximum
    converted_data = image_data.astype(float)
    if stretching_stringvar == "Square Root":
        stretch = SqrtStretch()
        converted_data = (converted_data - np.min(converted_data)
                          ) / np.ptp(converted_data)
        converted_data = stretch(converted_data)

    if stretching_stringvar == "Log":
        stretch = LogStretch()
        converted_data = (converted_data - np.min(converted_data)
                          ) / np.ptp(converted_data)
        converted_data = stretch(converted_data)

    if stretching_stringvar == "Asinh":
        stretch = AsinhStretch()
        converted_data = (converted_data - np.min(converted_data)
                          ) / np.ptp(converted_data)
        converted_data = stretch(converted_data)

    generated_image = Image.fromarray(converted_data)
    FITS_maximum = np.max(converted_data)
    width, height = generated_image.size
    new_size = (int(width * zoom_level), int(height * zoom_level))
    generated_image = generated_image.resize(new_size, Image.LANCZOS)
    generated_image = ImageMath.eval("(a + " +
                                     str(stretch_min / 100 * FITS_maximum) +
                                     ") * 255 / " +
                                     str(stretch_max / 100 * FITS_maximum),
                                     a=generated_image)

image_file = ""
image_data = []
       


class MyGUI:

    zoom_level = 1
    linreg_error = 0
    zoom_step = 0.5
    photometry_results_plotted = False
    ePSF_results_plotted = False
    results_tab_df = pd.DataFrame()
    image_bkg_value = 0
    bkg2D = None #if fetched, the Background2D object
    fit_shape = 21
    error_raised = False
    histogram_slider_low = 0
    histogram_slider_high = 5
    last_clicked_x = 0
    last_clicked_y = 0
    ensemble_size = 0
    jd = 0
    image_file = ""
    photometry_circles = {}
    valid_parameter_list = {}
    ePSF_rejection_list = pd.DataFrame({'x':[],'y':[]})
    epsf_model = None
    stars_tbl = None 

    # Parameter declaration  and init 
    photometry_aperture_entry = None
    max_ensemble_magnitude_entry = None
    fwhm_entry = None
    star_detection_threshold_entry = None
    photometry_iterations_entry = None
    sharplo_entry = None
    bkg_filter_size_entry = None
    matching_radius_entry = None
    aavso_obscode_entry = None
    telescope_entry = None
    tbv_entry = None
    tv_bv_entry = None
    tb_bv_entry = None
    tvr_entry = None
    tv_vr_entry = None
    tr_vr_entry = None
    tvi_entry = None
    tv_vi_entry = None
    ti_vi_entry = None
    catalog_stringvar = None
    vizier_catalog_entry = None
    fitter_stringvar = None
    astrometrynet_entry = None
    astrometrynet_key_entry = None
    object_kref_entry = None
    object_sel_comp_entry = None
    object_name_entry = None
    object_notes_entry = None

    #
    # The TopLoevel window containing the settings
    es_top = None

    def console_msg(self, MAOPhot_message, level=logging.INFO):
        # add a time stamp
        message = datetime.datetime.now().strftime("%d %b %Y %H:%M:%S")
        message += "      " + MAOPhot_message
        self.console.insert(tk.END, message+"\n")
        self.console.see(tk.END)
        self.window.update_idletasks()
        
        self.our_logger.log(level=level, msg=MAOPhot_message)
        
        """
        file = open("MAOPhotpsf.log", "a+", encoding="utf-8")
        file.write(str(message) + "\n")
        file.close()
        """

    def display_image(self):
        if len(image_data) > 0:
            self.canvas.delete("all")
            global generated_image
            generate_FITS_thumbnail(self.histogram_slider_low,
                                    self.histogram_slider_high,
                                    self.zoom_level,
                                    self.stretching_stringvar.get())
            self.image = ImageTk.PhotoImage(generated_image)
            self.canvas.create_image(0, 0, anchor=tk.NW, image=self.image)
            self.canvas.config(scrollregion=self.canvas.bbox(tk.ALL))
            self.canvas.bind("<Button-1>", self.mouse_click)
            if self.ePSF_results_plotted:
                self.plot_ePSF()
            self.plot_photometry()

    def load_FITS(self, image_file):
        global image_figure
        global image_data
        global image
        global image_width
        global image_height
        global header
        global FITS_minimum
        global FITS_maximum
        global generated_image
        try:
            self.console_msg("Loading FITS: " + image_file)
            with fits.open(image_file) as image:
                self.image_file = image_file
                self.filename_label['text'] = "FITS: " + image_file
                self.canvas.delete("all")
                self.zoom_level = 1
                self.photometry_results_plotted = False
                self.ePSF_results_plotted = False
                #Load previous work if it exists
                if os.path.isfile(self.image_file+".csv"):
                    self.results_tab_df = pd.read_csv(self.image_file + ".csv")
                else:
                    self.results_tab_df = pd.DataFrame()
                header = image[0].header
                image_data = fits.getdata(image_file)
                image_width = image_data.shape[1]
                image_height = image_data.shape[0]
                self.wcs_header = WCS(image[0].header)

                FITS_minimum = np.min(image_data)
                FITS_maximum = np.max(image_data)
                self.console_msg("Width: " + str(image_width) +
                                 " Height: " + str(image_height))
                self.console_msg(
                    "FITS Minimum: " + str(FITS_minimum) + " Maximum: " +
                    str(FITS_maximum))
                if 'filter' in header:
                    self.filter = str(header['filter'])
                    self.set_entry_text(self.filter_entry, self.filter)
                    self.console_msg("Filter: " + self.filter)
                else:
                    self.console_msg(
                        "Filter name not in FITS header. Set filter manually.")
                if 'airmass' in header:
                    self.airmass = str(header['airmass'])
                    self.set_entry_text(self.airmass_entry, self.airmass)
                    self.console_msg("Airmass: " + self.airmass)
                else:
                    self.console_msg(
                        "Airmass not in FITS header. Airmass may be required for AAVSO report. Set Airmass manually.")
                if 'exptime' in header:
                    exptime = header['exptime']
                    self.set_entry_text(self.exposure_entry, exptime)
                    self.console_msg("Exposure: " + str(exptime))
                else:
                    self.console_msg(
                        "Exposure (EXPTIME) not in FITS header. Set exposure manually.")

                self.jd = 0

                if 'date-obs' in header:
                    try:
                        date_obs = Time(header['date-obs'])
                        self.jd = Time(date_obs, format='jd')
                        self.console_msg("DATE-OBS: " + str(self.jd.to_value('iso')) + " UTC; JD: " + str(self.jd))
                        self.set_entry_text(self.date_obs_entry, str(self.jd))
                        
                    except Exception as e:
                        self.error_raised = True
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)
                        pass

                if 'jd' in header:
                    jd = header['jd']
                    self.console_msg(
                        "Julian date at the start of exposure (from JD): " + str(jd))
                    self.jd = jd
                    self.date_obs_entry.delete(0, tk.END)
                    self.date_obs_entry.insert(0, str(self.jd))

                self.image_bkg_value = np.median(image_data)
                self.console_msg("Median background level, ADU: " 
                    + str(self.image_bkg_value))
                self.console_msg("Ready")

                
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) 
                +" "+str(e), level=logging.ERROR)
            pass

    def open_FITS_file(self):
        global header

        options = {}
        options['defaultextension'] = '.fit'
        options['filetypes'] = [('FIT', '.fit'),('FITS', '.fits'), ('FTS', '.fts')]
        options['title'] = 'Open FIT image file...'

        image_file = fd.askopenfilename(**options)
        
        if len(image_file) > 0:
            try:
                
                self.load_FITS(image_file)
                self.display_image()
                self.clear_psf_label()
                self.clear_epsf_plot()
                
            except Exception as e:
                self.error_raised = True
                exc_type, exc_obj, exc_tb = sys.exc_info()
                self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
            

    def save_FITS_file_as(self):
        global image
        global image_data
        global header

        options = {}
        options['defaultextension'] = '.fit'
        options['filetypes'] = [('FIT', '.fit'),('FITS', '.fits'), ('FTS', '.fts')]
        options['title'] = 'Save FIT image file as...'

        file_name = fd.asksaveasfile(**options)
        
        try:
            if len(str(file_name)) > 0:
                self.console_msg("Saving FITS as " + str(file_name.name))
                fits.writeto(file_name.name, image_data,
                             header, overwrite=True)
                self.console_msg("Saved.")
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  +" "+str(e), level=logging.ERROR)

    def save_FITS_file(self):
        global image
        global image_data
        global header
        file_name = self.image_file
        try:
            if len(str(file_name)) > 0:
                self.console_msg("Saving FITS as " + str(file_name))
                fits.writeto(file_name, image_data, header, overwrite=True)
                self.console_msg("Saved.")
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)   +" "+str(e), level=logging.ERROR)

    def create_ePSF(self):
        global header
        self.console_msg("Initiating Effective PSF building...")

        #make sure an image is loaded
        if len(image_data) == 0:
            self.console_msg("Cannot proceed; an image must be loaded first; use File->Open...")
            return

        try:
            # Find local peaks in an image that are above a specified 
            # threshold value. Threshold determined 
            # 
            #determine the background, it will be subtrated from image_data, later
            bkg_filter_size = int(self.bkg_filter_size_entry.get())
            #make sure this is not an even number 
            if bkg_filter_size % 2 == 0:
                bkg_filter_size += 1
            sigma_clip = SigmaClip(sigma=3.0)
            bkg_estimator = MedianBackground()
            self.console_msg("Estimating background...")
            self.bkg2D = Background2D(
                                image_data, 
                                box_size=(self.fit_shape * 10, self.fit_shape * 10),
                                filter_size=(bkg_filter_size, bkg_filter_size), 
                                edge_method='crop',
                                sigma_clip=sigma_clip,
                                bkg_estimator=bkg_estimator)

            self.console_msg("Median Background2D level: " 
                + str(self.bkg2D.background_median))

            peaks_tbl = find_peaks(image_data, threshold=self.bkg2D.background_median * 3)

            peaks_tbl_len = len(peaks_tbl)
            self.console_msg("ePSF: found " + str(peaks_tbl_len) + " peaks.")

            #mask out peaks near the boundary; use twice the aperature entry
            #
            ##Calculate size of cutouts for EPSFBuilder
            ## make it 10 times the aperture entry
            ##same size used in plot_ePSF
            self.fit_shape = int(self.photometry_aperture_entry.get())
            size = 10*self.fit_shape + 1
            hsize = (size - 1)/2
            x = peaks_tbl['x_peak']  
            y = peaks_tbl['y_peak']  
            _image = Image.fromarray(image_data)
            width, height = _image.size
            mask = ((x > hsize) & (x < (width -1 - hsize)) &
                    (y > hsize) & (y < (height -1 - hsize)))  

            #prelim_stars_tbl are inbound stars (not to close to the edge)
            prelim_stars_tbl = Table()
            prelim_stars_tbl['x'] = x[mask]  
            prelim_stars_tbl['y'] = y[mask]  
            prelim_stars_tbl['rejected'] = False #init

            prelim_stars_tbl_len = len(prelim_stars_tbl)
            self.console_msg("ePSF: found and removed " + str(peaks_tbl_len - prelim_stars_tbl_len) + " peaks on edge.")
            self.console_msg("ePSF: " + str(prelim_stars_tbl_len) + " peaks remain.")

            # now set 'rejected' to True for any stars that are proximate to 
            # another in the same list
            for i in range(len(prelim_stars_tbl)):
                for ii in range(len(prelim_stars_tbl)):
                    if ii == i:
                        continue
                    i_x = prelim_stars_tbl[i]['x']
                    i_y = prelim_stars_tbl[i]['y']
                    ii_x = prelim_stars_tbl[ii]['x']
                    ii_y = prelim_stars_tbl[ii]['y']
                    if math.dist([i_x, i_y], [ii_x, ii_y]) <= hsize:
                        #reject this because it is too close to that companion
                        prelim_stars_tbl[i]['rejected'] = True
                        prelim_stars_tbl[ii]['rejected'] = True
    
            x = prelim_stars_tbl['x']  
            y = prelim_stars_tbl['y']  
            reject_this = prelim_stars_tbl['rejected']

            mask = reject_this == False  # only keep ones we don't reject

            isolated_stars_tbl = Table()
            isolated_stars_tbl['x'] = x[mask]  
            isolated_stars_tbl['y'] = y[mask]  
            isolated_stars_tbl['rejected'] = False #init

            isolated_stars_tbl_len = len(isolated_stars_tbl)
            self.console_msg("ePSF: found and removed " + str(prelim_stars_tbl_len - isolated_stars_tbl_len) + " close companions.")
            self.console_msg("ePSF: " + str(isolated_stars_tbl_len) + " peaks remain.")

            # now set 'rejected' to True for any stars that are proximate to a 
            # coordinate in ePSF_rejection_list
            for psf_x, psf_y, psf_reject in isolated_stars_tbl:
                isolated_stars_tbl.add_index('x')
                for index, row in self.ePSF_rejection_list.iterrows():
                    reject_x = row['x']
                    reject_y = row['y']
                    if abs(reject_x - psf_x) <= hsize and abs(reject_y - psf_y) <= hsize:
                        #user does not want this one
                        isolated_stars_tbl.loc[psf_x]['rejected'] = True
                        break
    
            x = isolated_stars_tbl['x']  
            y = isolated_stars_tbl['y']  
            reject_this = isolated_stars_tbl['rejected']

            mask = reject_this == False  # only keep ones we don't reject

            self.stars_tbl = Table()
            self.stars_tbl['x'] = x[mask]  
            self.stars_tbl['y'] = y[mask]  

            stars_tbl_len = len(self.stars_tbl)
            self.console_msg("ePSF: found and removed " + str(isolated_stars_tbl_len - stars_tbl_len) + " peaks rejected by user.")
            self.console_msg("ePSF: " + str(stars_tbl_len) + " peaks remain for EPSF Builder.")

            # Subtract background from image_data
            # bkg.background is a 2D ndarray of background image
            clean_image = image_data-self.bkg2D.background

            working_image = NDData(data=clean_image)

            candidate_stars = extract_stars(working_image, self.stars_tbl, size=size)  
            epsf_builder = EPSFBuilder(oversampling=4, maxiters=10, progress_bar=True)
            self.epsf_model, fitted_stars = epsf_builder(candidate_stars)  

            model_length = len(self.epsf_model.data)
            x = np.arange(0, model_length, 1)
            y = np.arange(0, model_length, 1)
            x, y = np.meshgrid(x, y)

            #plot it in the psf
            self.psf_plot.clear()
            self.psf_plot.plot_surface(x, y, self.epsf_model.data, cmap=cm.jet)
            self.psf_plot_canvas.draw()

            #plot it in the ePSF
            self.clear_epsf_plot()
            self.plot_ePSF()

            norm = simple_norm(self.epsf_model.data, 'log', percent=99.)
            im = self.ePSF_plot.imshow(self.epsf_model.data, norm=norm, origin='lower', cmap='viridis')
            self.ePSF_plot.set_title("Effective PSF")
            self.fig_ePSF.colorbar(im, ax=self.ePSF_plot)

            self.ePSF_plot_canvas.draw()

            self.ePSF_results_plotted = True

            self.console_msg("Ready")

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  +" "+str(e), level=logging.ERROR)

    def clear_ePSF(self):
        global header
        self.console_msg("Clearing ePSF model, Rejection List, plot...")
        #drop all the rows but keep the 'x' and 'y' column
        self.ePSF_rejection_list.drop(self.ePSF_rejection_list.index, inplace=True)
        self.epsf_model = None #reset
        self.clear_psf_label()
        self.clear_epsf_plot()
        self.ePSF_results_plotted = False
        self.display_image()
        self.console_msg("Ready")
        return

    def load_ePSF_rejection_list(self):
        global header

        try:
            options = {}
            options['defaultextension'] = '.csv'
            options['filetypes'] = [('CSV', '.csv')]
            options['title'] = 'Choose the ' + self.image_file + '-rejection.csv'

            file_name = fd.askopenfilename(**options)

            if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
                self.console_msg("Loading Rejection list from: " + str(file_name))
                self.ePSF_rejection_list = pd.read_csv(str(file_name))
                self.display_image()
            else:
                return
            self.console_msg("Ready")

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)

        return

    def save_as_ePSF_rejection_list(self):
        global header

        try:
            options = {}
            options['defaultextension'] = '.csv'
            options['filetypes'] = [('CSV', '.csv')]
            options['title'] = 'Save Rejection List As... ' + self.image_file + '-rejection.csv'
            options['initialfile'] = self.object_name_entry.get() + '-rejection.csv'

            file_name = fd.asksaveasfile(**options)

            if len(str(file_name)) > 0:
                self.console_msg("Saving Rejection List as " + file_name.name)
                dir_path = os.path.dirname(os.path.realpath(file_name.name)) + "\\"
                self.ePSF_rejection_list.to_csv(file_name.name, index=False)

            self.console_msg("Ready")

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)

        return

    def execute_noniterative_psf_photometry(self):
        global header
        self.console_msg(
            "Initiating Non-iterative PSF photometry...")
        if len(image_data) == 0:
            self.console_msg("Cannot proceed; an image must be loaded first; use File->Open...")
            return
        try:
            #determine the background, it will be subtrated from image_data, later
            bkg_filter_size = int(self.bkg_filter_size_entry.get())
            #make sure this is not an even number 
            if bkg_filter_size % 2 == 0:
                bkg_filter_size += 1
            self.fit_shape = int(self.photometry_aperture_entry.get())
            sigma_clip = SigmaClip(sigma=3.0)
            bkg_estimator = MedianBackground()
            self.console_msg("Estimating background...")
            #Backgound2D used because background may vary over FOV
            self.bkg2D = Background2D(
                                image_data, 
                                box_size=(self.fit_shape * 10, self.fit_shape * 10),
                                filter_size=(bkg_filter_size, bkg_filter_size), 
                                edge_method='crop',
                                sigma_clip=sigma_clip,
                                bkg_estimator=bkg_estimator)

            self.console_msg("Median Background2D level: "
                             + str(self.bkg2D.background_median))

            # Subtract background from image_data
            # bkg.background is a 2D ndarray of background image
            clean_image = image_data-self.bkg2D.background

            working_image = NDData(data=clean_image)

            fwhm = float(self.fwhm_entry.get())
            star_detection_threshold = float(
                self.star_detection_threshold_entry.get())

            sharplo = float(self.sharplo_entry.get())

            self.console_msg("Finding stars...")
            star_find = IRAFStarFinder(threshold = star_detection_threshold,
                                        fwhm = fwhm,
                                        exclude_border=True,
                                        roundhi = 3.0,
                                        roundlo = -5.0,
                                        sharplo = sharplo,
                                        sharphi = 2.0)

           
            # the 2.5 in the following is from 
            # https://photutils.readthedocs.io/en/stable/grouping.html#getting-started
            #grouper = SourceGrouper(2.5 * fwhm)

            local_bkg = LocalBackground(inner_radius=fwhm*4, outer_radius=fwhm*8)

            if self.fitter_stringvar.get() == "Sequential LS Programming":
                self.console_msg(
                    "Setting fitter to Sequential Least Squares Programming")
                selected_fitter = SLSQPLSQFitter()

            elif self.fitter_stringvar.get() == "Simplex LS":
                self.console_msg(
                    "Setting fitter to Simplex and Least Squares Statistic")
                selected_fitter = SimplexLSQFitter()

            else: # self.fitter_stringvar.get() == "Levenberg-Marquardt":
                self.console_msg("Setting fitter to Levenberg-Marquardt")
                selected_fitter = LevMarLSQFitter()

            if self.epsf_model != None:
                self.console_msg("Using derived Effective PSF Model")
                psf_model = self.epsf_model
            else:
                #Use Gausian
                sigma = 2.00 * fwhm / gaussian_sigma_to_fwhm 
                self.console_msg("Using Gaussian PRF for model; sigma = "+format(sigma, '.3f')+\
                                 "; fwhm = "+format(fwhm, '.1f')+\
                                 "; fitting width/height = "+str(self.fit_shape))
                psf_model = IntegratedGaussianPRF(sigma)


            photometry = PSFPhotometry(
                                        psf_model = psf_model,
                                        fit_shape = self.fit_shape,
                                        finder = star_find,
                                        grouper = None,
                                        fitter = selected_fitter,
                                        fitter_maxiters = 100,
                                        localbkg_estimator = local_bkg,
                                        aperture_radius=1.5*fwhm,
                                        progress_bar=True
                                        )
            


            self.console_msg("Performing photometry...")
            result_tab = photometry(data=working_image)

            if 'message' in selected_fitter.fit_info:
                self.console_msg("Done. PSF fitter message(s): " + str(selected_fitter.fit_info['message']))
            else:
                self.console_msg("Done. PSF fitter; no message available")

            self.results_tab_df = result_tab[result_tab.colnames[0:14]].to_pandas()  # only need first 15 columns
            self.results_tab_df["removed_from_ensemble"] = False
            self.results_tab_df["date-obs"] = float(self.date_obs_entry.get())
            if self.airmass_entry.get().isnumeric():
                self.results_tab_df["AMASS"] = float(self.airmass_entry.get())
            else:
                self.results_tab_df["AMASS"] = "na"
                


            # Calculate instrumental magnitudes
            # Following added for "True" inst mag used in AAVSO report
            image_exposure_time = float(self.exposure_entry.get())
            self.results_tab_df["inst_mag"] = -2.5 * np.log10(self.results_tab_df["flux_fit"] / image_exposure_time)

            #record for later 
            self.results_tab_df["exposure"] = image_exposure_time

            self.results_tab_df.to_csv(self.image_file + ".csv", index=False)
            self.console_msg("Photometry saved to " + str(self.image_file + ".csv") + "; len = " + str(len(self.results_tab_df)))
            self.plot_photometry()

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  +" "+str(e), level=logging.ERROR)

    def execute_iterative_psf_photometry(self):
        global header
        self.console_msg(
            "Initiating iterative PSF photometry...")
        if len(image_data) == 0:
            self.console_msg("Cannot proceed; an image must be loaded first; use File->Open...")
            return
        try:
            #determine the background, it will be subtrated from image_data, later
            bkg_filter_size = int(self.bkg_filter_size_entry.get())
            #make sure this is not an even number 
            if bkg_filter_size % 2 == 0:
                bkg_filter_size += 1

            self.fit_shape = int(self.photometry_aperture_entry.get())
            sigma_clip = SigmaClip(sigma=3.0)
            bkg_estimator = MedianBackground()
            self.console_msg("Estimating background...")
            #Backgound2D used because background may vary over FOV
            self.bkg2D = Background2D(
                                image_data, 
                                box_size=(self.fit_shape * 10, self.fit_shape * 10),
                                filter_size=(bkg_filter_size, bkg_filter_size), 
                                edge_method='crop',
                                sigma_clip=sigma_clip,
                                bkg_estimator=bkg_estimator)

            self.console_msg("Median Background2D level: "
                             + str(self.bkg2D.background_median))

            # Subtract background from image_data
            # bkg.background is a 2D ndarray of background image
            clean_image = image_data-self.bkg2D.background

            working_image = NDData(data=clean_image)

            fwhm = float(self.fwhm_entry.get())
            star_detection_threshold = float(
                self.star_detection_threshold_entry.get())

            iterations = int(self.photometry_iterations_entry.get()) if self.photometry_iterations_entry.get().isnumeric() else None

            sharplo = float(self.sharplo_entry.get())

            #bkg_filter_size = int(self.bkg_filter_size_entry.get())
            self.console_msg("Finding stars...")
            star_find = IRAFStarFinder(threshold = star_detection_threshold,
                                        fwhm = fwhm,
                                        exclude_border=True,
                                        roundhi = 3.0,
                                        roundlo = -5.0,
                                        sharplo = sharplo,
                                        sharphi = 2.0)

           
            local_bkg = LocalBackground(inner_radius=fwhm*4, outer_radius=fwhm*8)

            if self.fitter_stringvar.get() == "Sequential LS Programming":
                self.console_msg(
                    "Setting fitter to Sequential Least Squares Programming")
                selected_fitter = SLSQPLSQFitter()

            elif self.fitter_stringvar.get() == "Simplex LS":
                self.console_msg(
                    "Setting fitter to Simplex and Least Squares Statistic")
                selected_fitter = SimplexLSQFitter()

            else: # self.fitter_stringvar.get() == "Levenberg-Marquardt":
                self.console_msg("Setting fitter to Levenberg-Marquardt")
                selected_fitter = LevMarLSQFitter(calc_uncertainties=True)

            if self.epsf_model != None:
                self.console_msg("Using derived Effective PSF Model")
                psf_model = self.epsf_model
            else:
                #Use Gausian
                sigma = 2.00 * fwhm / gaussian_sigma_to_fwhm 
                self.console_msg("Using Gaussian PRF for model; sigma = "+format(sigma, '.3f')+\
                                 "; fwhm = "+format(fwhm, '.1f')+\
                                 "; fitting width/height = "+str(self.fit_shape))
                psf_model = IntegratedGaussianPRF(sigma)

            photometry = IterativePSFPhotometry(
                                                psf_model = psf_model,
                                                fit_shape = self.fit_shape,
                                                finder = star_find,
                                                grouper = None,
                                                fitter = selected_fitter,
                                                fitter_maxiters = 100,
                                                maxiters = iterations,
                                                localbkg_estimator = local_bkg,
                                                aperture_radius=1.5*fwhm,
                                                sub_shape=None, #defaults to fit_shape
                                                progress_bar=True
                                                )

            result_tab = photometry(data=working_image)

            if 'message' in selected_fitter.fit_info:
                self.console_msg("Done. PSF fitter message(s): " + str(selected_fitter.fit_info['message']))
            else:
                self.console_msg("Done. PSF fitter; no message available")

            #self.results_tab_df = result_tab[result_tab.colnames[0:14]].to_pandas()  # only need first 15 columns
            self.results_tab_df = result_tab.to_pandas()
            self.results_tab_df["removed_from_ensemble"] = False
            self.results_tab_df["date-obs"] = float(self.date_obs_entry.get())
            if self.airmass_entry.get().isnumeric():
                self.results_tab_df["AMASS"] = float(self.airmass_entry.get())
            else:
                self.results_tab_df["AMASS"] = "na"

            # Calculate instrumental magnitudes
            # Following added for "True" inst mag used in AAVSO report
            image_exposure_time = float(self.exposure_entry.get())
            self.results_tab_df["inst_mag"] = -2.5 * np.log10(self.results_tab_df["flux_fit"] / image_exposure_time)

            #record for later 
            self.results_tab_df["exposure"] = image_exposure_time

            self.results_tab_df.to_csv(self.image_file + ".csv", index=False)
            self.console_msg("Photometry saved to " + str(self.image_file + ".csv") + "; len = " + str(len(self.results_tab_df)))
            self.plot_photometry()

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  +" "+str(e), level=logging.ERROR)



    # center coordinates, radius
    def create_circle(self, x, y, r, canvas_name, outline="grey50"):
        x0 = x - r
        y0 = y - r
        x1 = x + r
        y1 = y + r
        self.photometry_circles[str(x)+str(y)] = canvas_name.create_oval(x0, y0, x1, y1, outline=outline)

    def create_text(self, x, y, r, canvas_name, anchor, text='', fill='white'):
        x1 = x 
        y1 = y - 1.25*r
        canvas_name.create_text(x1, y1, fill=fill, anchor=anchor, text=text)

    def plot_ePSF(self):
        try:
            #display the non-rejected stars as white, and rejected as red circles

            ## make it twice the aperture entry
            ##same size used in create_ePSF
            self.fit_shape = int(self.photometry_aperture_entry.get())
            size = 2*self.fit_shape + 1
            hsize = (size - 1)/2
            for psf_x, psf_y in self.stars_tbl.iterrows('x', 'y'):
                color = 'white' # it is a white circle until a reject match is found
                for index, row in self.ePSF_rejection_list.iterrows():
                    reject_x = row['x']
                    reject_y = row['y']
                    if abs(reject_x - psf_x) <= hsize and abs(reject_y - psf_y) <= hsize:
                        color = 'red'  #paint rejects red
                        break;
                self.create_circle(x=psf_x * self.zoom_level, y=psf_y * self.zoom_level,
                                      r=hsize * self.zoom_level, canvas_name=self.canvas, outline=color)

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) + " " + str(e), level=logging.ERROR)

    def plot_photometry(self):
        try:
            self.console_msg("Plotting Photometry...")

            labels_in_photometry_table = "label" in self.results_tab_df
            vsx_ids_in_photometry_table = "vsx_id" in self.results_tab_df

            if os.path.isfile(self.image_file+".csv"):
                self.fit_shape = int(self.photometry_aperture_entry.get())
                self.results_tab_df = pd.read_csv(self.image_file + ".csv")
                if "removed_from_ensemble" not in self.results_tab_df:
                    # This prefilling is required for backwards compatibility to read .phot
		            #(now called .csv) files from older versions.
                    self.results_tab_df["removed_from_ensemble"] = False

                self.ePSF_results_plotted = False
                self.photometry_results_plotted = True

                for index, row in self.results_tab_df.iterrows():
                    outline = "grey50"
                    if labels_in_photometry_table:
                        if len(str(row["label"])) > 0 and str(row["label"]) != "nan" and str(row["label"]).isnumeric():
                            outline = "green"
                            self.create_circle(x=row["x_fit"] * self.zoom_level,
                                y=row["y_fit"] * self.zoom_level,
                                r=self.fit_shape / 2 * self.zoom_level,
                                canvas_name=self.canvas, outline=outline)
                            self.create_text(  x=row["x_fit"] * self.zoom_level,
                                y=row["y_fit"] * self.zoom_level, 
                                r=self.fit_shape / 2 * self.zoom_level,
                                canvas_name=self.canvas,
                                anchor=tk.CENTER,
                                text=str(int(row["label"])))
                            continue

                    if row["removed_from_ensemble"]:
                        assert False, "Found an entry 'removed from ensemble???!'"

                    if vsx_ids_in_photometry_table:
                        if len(str(row["vsx_id"])) > 0 and str(row["vsx_id"]) != "nan":
                            outline = "yellow"
                            self.create_circle(x=row["x_fit"] * self.zoom_level,
                                y=row["y_fit"] * self.zoom_level,
                                r=self.fit_shape / 2 * self.zoom_level,
                                canvas_name=self.canvas,
                                outline=outline)
                            self.create_text(  x=row["x_fit"] * self.zoom_level,
                                y=row["y_fit"] * self.zoom_level, 
                                r=self.fit_shape / 2 * self.zoom_level,
                                canvas_name=self.canvas,
                                anchor=tk.CENTER,
                                text=str(row["vsx_id"]).strip(),
                                fill='yellow')

            self.console_msg("Plotting Photometry...complete")
            self.console_msg("Ready")
            
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) + " " + str(e), level=logging.ERROR)

    def match_photometry_table(self, x, y, r=5):
        x_criterion = self.results_tab_df['x_fit'] < (x + r)
        matched_objects = self.results_tab_df[x_criterion]
        x_criterion = self.results_tab_df['x_fit'] > (x - r)
        matched_objects = matched_objects[x_criterion]
        y_criterion = self.results_tab_df['y_fit'] < (y + r)
        matched_objects = matched_objects[y_criterion]
        y_criterion = self.results_tab_df['y_fit'] > (y - r)
        matched_objects = matched_objects[y_criterion]
        if len(matched_objects) > 0:
            return(matched_objects.iloc[0]["x_fit"],
                   matched_objects.iloc[0]["y_fit"],
                   matched_objects.iloc[0]["flux_fit"],
                   matched_objects.iloc[0]["inst_mag"]
                   )
        else:
            return(0, 0, 0, 0, 0, 0) #, 0, 0, 0)

    def mouse_click(self, event):
        global image_data
        x = int(self.canvas.canvasx(event.x) / self.zoom_level)
        y = int(self.canvas.canvasy(event.y) / self.zoom_level)
        self.last_clicked_x = x
        self.last_clicked_y = y
        ADU = image_data[y-1, x-1]
        sky = self.wcs_header.pixel_to_world(x, y)
        sky_coordinate_string = ""


        #clear plot label
        self.plotname_label['text'] = "Plot: "

        if hasattr(sky, 'ra'):
            c = SkyCoord(ra=sky.ra, dec=sky.dec)
            sky_coordinate_string = "α δ: "+c.to_string("hmsdms")
        self.console_msg("Position X: "+str(x)+"\t Y: "+str(y) +
                         "\t ADU: "+str(ADU) + "\t\t\t" + sky_coordinate_string)
        psf_canvas_x = x
        psf_canvas_y = y

        if self.ePSF_results_plotted:
            #add the selected coordinate into the ePSF_rejection_list
            self.ePSF_rejection_list.loc[len(self.ePSF_rejection_list.index)] = [x, y]
            #indicate the rejected ones
            self.display_image()

        elif self.photometry_results_plotted:
            vsx_ids_in_photometry_table = "vsx_id" in self.results_tab_df
            self.display_image()
            self.console_msg("")

            x_fit, y_fit, flux_fit, inst_mag = self.match_photometry_table(x, y)
            sky = self.wcs_header.pixel_to_world(x_fit, y_fit)
            sky_coordinate_string = ""
            if hasattr(sky, 'ra'):
                c = SkyCoord(ra=sky.ra, dec=sky.dec)
                sky_coordinate_string = " α δ: " + c.to_string("hmsdms")
            if x_fit != 0 and y_fit != 0:
                psf_canvas_x = x_fit
                psf_canvas_y = y_fit
            if str(x_fit)+str(y_fit) in self.photometry_circles:
                self.canvas.delete(
                    self.photometry_circles[str(x_fit)+str(y_fit)])
            self.canvas.create_line(x_fit*self.zoom_level, y_fit*self.zoom_level - 35*self.zoom_level, x_fit *
                                    self.zoom_level, y_fit*self.zoom_level - 10*self.zoom_level, fill="white")  # Draw "target" lines
            self.canvas.create_line(x_fit*self.zoom_level+35*self.zoom_level, y_fit*self.zoom_level,
                                    x_fit*self.zoom_level + 10*self.zoom_level, y_fit*self.zoom_level, fill="white")
            self.console_msg("Photometry fits, X: " + str(round(x_fit, 2)) + " Y: " + str(round(y_fit, 2)) + " Flux (ADU): " + str(
                round(flux_fit, 2)) + " Instrumental magnitude: " + str(round(inst_mag, 3)) + " " + sky_coordinate_string)

            # Reset object name field in the left panel to avoid user mistakes
            self.set_entry_text(self.object_name_entry, "")
            if "match_id" in self.results_tab_df:
                matching_star_criterion = (self.results_tab_df["x_fit"] == x_fit) & (
                    self.results_tab_df["y_fit"] == y_fit)
                if len(self.results_tab_df[matching_star_criterion]) > 0:
                    matching_star = self.results_tab_df[matching_star_criterion].iloc[0]
                    if type(matching_star["match_id"]) in (str, int, np.float64):
                        self.console_msg(
                            "Matching catalog source ID: " + str(matching_star["match_id"]) + 
                                "; label: " + str(int(matching_star["label"])) +
                                " magnitude: " + str(matching_star["match_mag"]))
                        self.set_entry_text(self.object_name_entry, str(matching_star["match_id"]))
                        
                        #update plot label
                        self.plotname_label['text'] = "Plot: " + str(matching_star["match_id"]) + \
                            "; " + str(int(matching_star["label"]))
                        
                    if vsx_ids_in_photometry_table:
                        if len(str(matching_star["vsx_id"])) > 0 and str(matching_star["vsx_id"]) != "nan":
                            self.console_msg(
                                "Matching VSX Source: " + str(matching_star["vsx_id"]))
                            self.set_entry_text(
                                self.object_name_entry, str(matching_star["vsx_id"]))
                            #update plot label
                            self.plotname_label['text'] = "Plot: " + str(matching_star["vsx_id"])

            self.update_PSF_canvas(psf_canvas_x, psf_canvas_y)

    def update_PSF_canvas(self, x, y):
        global image_data
        global FITS_minimum
        global FITS_maximum
        try:
            if len(image_data) > 0:
                self.fit_shape = int(self.photometry_aperture_entry.get())
                x0 = int(x - (self.fit_shape - 1) / 2)
                y0 = int(y - (self.fit_shape - 1) / 2)
                x1 = int(x + (self.fit_shape - 1) / 2)
                y1 = int(y + (self.fit_shape - 1) / 2)
                position = (x, y)
                size = (self.fit_shape - 1, self.fit_shape - 1)
                data = Cutout2D(image_data, position, size).data
                x = np.arange(x0, x1, 1)
                y = np.arange(y0, y1, 1)
                x, y = np.meshgrid(x, y)
                self.psf_plot.clear()
                self.psf_plot.plot_surface(x, y, data, cmap=cm.jet)
                self.psf_plot_canvas.draw()
        except Exception as e:
            self.error_raised = True
            pass

    def clear_epsf_plot(self):
        self.ePSF_plot.clear()
        self.fig_ePSF.clear()
        self.ePSF_plot_canvas.draw()
        self.fig_ePSF, self.ePSF_plot = plt.subplots()
        self.ePSF_plot_canvas = FigureCanvasTkAgg(self.fig_ePSF, self.right_frame)
        self.ePSF_plot_canvas.draw()
        self.ePSF_canvas = self.ePSF_plot_canvas.get_tk_widget()
        self.ePSF_canvas.config(width=int(self.screen_width/8.5), height=int(self.screen_width/8.5))
        # Allocate small PSF canvas to a new grid inside the right_frame
        self.ePSF_canvas.grid(row=3, column=0)   #was row0

    def clear_psf_label(self):
        #clear plot label
        self.plotname_label['text'] = "Plot: "
        self.psf_plot.clear()
        self.psf_plot_canvas.draw()

    def update_PSF_canvas_2d(self, x, y):
        global image_data
        global FITS_minimum
        global FITS_maximum
        image_crop = Image.fromarray(image_data)
        self.fit_shape = int(self.photometry_aperture_entry.get())
        x0 = int(x - (self.fit_shape-1)/2)
        y0 = int(y - (self.fit_shape-1)/2)
        x1 = int(x + (self.fit_shape-1)/2)
        y1 = int(y + (self.fit_shape-1)/2)
        image_crop = image_crop.crop((x0, y0, x1, y1))
        image_crop = image_crop.resize((300, 300), resample=0)
        image_crop = ImageMath.eval(
            "a * 255 / " + str(self.histogram_slider_high / 100 * FITS_maximum), a=image_crop)
        self.image_crop = ImageTk.PhotoImage(image_crop)
        self.psf_canvas.create_image(0, 0, anchor=tk.NW, image=self.image_crop)

    def zoom_in(self):
        self.canvas.scale("all", 0, 0, 1+self.zoom_step, 1+self.zoom_step)
        self.zoom_level = self.zoom_level * (1+self.zoom_step)
        self.console_msg("Zoom: "+str(self.zoom_level))
        self.display_image()

    def zoom_out(self):
        self.canvas.scale("all", 0, 0, 1-self.zoom_step, 1-self.zoom_step)
        self.zoom_level = self.zoom_level * (1-self.zoom_step)
        self.console_msg("Zoom: " + str(self.zoom_level))
        self.display_image()

    def zoom_100(self):
        self.canvas.scale("all", 0, 0, 1, 1)
        self.zoom_level = 1
        self.console_msg("Zoom: " + str(self.zoom_level))
        self.display_image()

    def solve_image(self):
        global generated_image
        global header

        if "flux_fit" not in self.results_tab_df:
            self.console_msg("Cannot plate solve before executing PSF Photometry; execute 'Photometry->Iterative PSF Photometry' first.")
            self.console_msg("Ready")
            return

        self.console_msg("Solving via Astrometry.Net...")
        try:
            # First check if a phtometry table exists.

            ast = AstrometryNet()
            ast.api_key = self.astrometrynet_key_entry.get()
            ast.URL = "http://" + self.astrometrynet_entry.get()
            ast.API_URL = "http://" + self.astrometrynet_entry.get() + "/api"

            sources_df = self.results_tab_df.sort_values("flux_fit", ascending=False)
            width, height = generated_image.size

            self.wcs_header = ast.solve_from_source_list(sources_df['x_fit'], sources_df['y_fit'],
                                                         width, height,
                                                         solve_timeout=360)
            self.console_msg("Astrometry.Net solution reference point RA: " + 
                                str(self.wcs_header["CRVAL1"]) + " Dec: " + 
                                str(self.wcs_header["CRVAL2"]))
            self.console_msg("Ready")

            header = header + self.wcs_header
            self.wcs_header = WCS(header)
            
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
            

    #
    #   BV_two_color_photometry: called from Menu selection  Two Color Photometry->(B-V)
    #
    def BV_two_color_photometry(self):
        self.two_color_photometry('B-V')

    #
    #   VR_two_color_photometry: called from Menu selection  Two Color Photometry->(V-R)
    #
    def VR_two_color_photometry(self):
        self.two_color_photometry('V-R')

    #
    #   VI_two_color_photometry: called from Menu selection  Two Color Photometry->(V-I)
    #
    def VI_two_color_photometry(self):
        self.two_color_photometry('V-I')

    #
    # two_color_phtometry
    # 
    # This calculates the two color photometry process as executed in AAVSO VPhot
    # TwoColorPhotometry. 
    # The parameter input_color is either 'B-V','V-R', or 'V-I' The formulae are the same/
    #
    def two_color_photometry(self, input_color):
        try:
            """
            #
            # NOTE! Since B-V was implemented first, the B-V filternames are
            # still used even though imput_color may be V-R; 
            # Only when it counts does the B-V change to the real V-R or V-I
            #

            Two Color Photometry requires 'Object Name' to be filled; eg. 'V1117 Her'

            These comments illustrate the case when input_color is B-V
            
            The following formulas are calculated for each comp star and var = check star;
            then calculated again for each comp star and var = "Object Name". Eg. V1117 Her
            
            Δv = vvar - vcomp
            Δ(B-V) = Tbv * Δ(b-v)
            Vvar = Δv + Tv_bv * Δ(B-V) + Vcomp
            
            where:
            Δv  = IM of variable  - IM comp
            Vcomp = published V-magnitude of comp
            Δ(B-V) = Tbv * difference between standard color of var and standard color of comp


            Build 2 tables, "check" and "var", with the following columns (see 
            E:\Astronomy\AAVSO\Reports\AAVSO Reports\MAO\2022 6 4 V1117 Her/
            TwoColor V1117_Her 2022 6 4.xlsx)
            
            type                     "check" or "var"
            name                     <check star label> or <var name>
            label                    label of comp star
            IMB                      instrumental mag of comp (B)
            IMV                      instrumental mag of comp (V)
            B                        published B mag
            V                        published V mag
            delta_b_minus_v          var(b-v) - comp_b_minus_v
            delta_B_minus_V          Tbv *  delta_b_minus_v
            delta_v                  var_IMV - IMV
            delta_b                  var_IMB - IMB
            comp_b_minus_v           IMB - IMV
            B_star                   Bvar = Δb + Tv_bv * Δ(B-V) + Bcomp
            V_star                   Vvar = Δv + Tv_bv * Δ(B-V) + Vcomp

            """
            # use first_filter, amd second_filter dict to index into appropriate filter
            first_filter = {'B-V': 'B', 'V-R': 'V', 'V-I': 'V'}
            second_filter = {'B-V': 'V', 'V-R': 'R', 'V-I': 'I'}

            #
            # NOTE! Since B-V was implemented first, the B-V filternames are
            # still used even though imput_color may be V-R; 
            # Only when it counts does the B-V change to the real V-R or V-I
            #

            variable_star = self.object_name_entry.get().strip()
            if len(variable_star) == 0:
                self.console_msg(
                    "Two Color Photometry requires 'Object Name' to be filled; eg. 'V1117 Her'")
                return
            
            # Ask for the B and V CSV files
            options = {}
            options['defaultextension'] = '.csv'
            options['filetypes'] = [('CSV', '.csv')]
            options['title'] = 'Choose a file for filter ' + first_filter[input_color]

            file_name = fd.askopenfilename(**options)

            if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
                self.console_msg("Loading filter " + first_filter[input_color] + " from " + str(file_name))
                self.results_tab_df_colorB = pd.read_csv(str(file_name), dtype={'check_star': bool})
            else:
                return

            #Test to make sure csv file is ready                
            if "label" not in self.results_tab_df_colorB:
                self.console_msg("Cannot proceed; run 'Photometry->Get Comparison Stars' first.")
                return

            options['title'] = 'Choose a file for filter ' + second_filter[input_color]
            file_name = fd.askopenfilename(**options)
            if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
                self.console_msg("Loading filter " + second_filter[input_color] + " from " + str(file_name))
                self.results_tab_df_colorV = pd.read_csv(str(file_name), dtype={'check_star': bool})
            else:
                return
            
            #Test to make sure csv file is ready                
            if "label" not in self.results_tab_df_colorV:
                self.console_msg("Cannot proceed; run 'Photometry->Get Comparison Stars' first.")
                return

            self.console_msg("Performing Two Color Photometry...")
            
            #get the transformation coefficients
            #
            # NOTE! Since B-V was implemented first, the B-V filternames are
            # still used even though imput_color may be V-R; 
            # Only when it counts does the B-V change to the real V-R or V-I
            #

            if input_color == 'B-V':
                tbv_coefficient = float(self.tbv_entry.get())
                tb_bv_coefficient = float(self.tb_bv_entry.get())
                tv_bv_coefficient = float(self.tv_bv_entry.get())
            elif input_color == 'V-R':
                tbv_coefficient = float(self.tvr_entry.get())
                tb_bv_coefficient = float(self.tv_vr_entry.get())
                tv_bv_coefficient = float(self.tr_vr_entry.get())
            elif input_color == 'V-I':
                tbv_coefficient = float(self.tvi_entry.get())
                tb_bv_coefficient = float(self.tv_vi_entry.get())
                tv_bv_coefficient = float(self.ti_vi_entry.get())
            else:
                raise Exception("two_color_photometry: unknown imput_color entered")

         
            #
            #   CHECK STAR Calculations
            #
            # Find the check star; 
            check_star_B = self.results_tab_df_colorB[self.results_tab_df_colorB["check_star"] == True].iloc[0]

            check_star_label = int(check_star_B["label"])

            self.console_msg("Using check star " + str(check_star_label))

            check_IMB = check_star_B["inst_mag"]
            check_B = check_star_B["match_mag"]

            check_star_V = self.results_tab_df_colorV[self.results_tab_df_colorV["check_star"] == True].iloc[0]
            check_IMV = check_star_V["inst_mag"]
            check_V = check_star_V["match_mag"]
            
            # Find the variable_star; 
            var_star_B = self.results_tab_df_colorB[self.results_tab_df_colorB["vsx_id"] == variable_star].iloc[0]
            var_IMB = var_star_B["inst_mag"]

            var_star_label = var_star_B["vsx_id"]

            var_star_V = self.results_tab_df_colorV[self.results_tab_df_colorV["vsx_id"] == variable_star].iloc[0]
            var_IMV = var_star_V["inst_mag"]


            # AAVSO adds EXPOSURE/2 to this time and sets it in the report
            # E.g., Z Tau,2460300.57931,13.167,0.018,V,YES,STD,ENSEMBLE,na,...
            date_obs_B = self.results_tab_df_colorB[self.results_tab_df_colorB["check_star"] == True].iloc[0]["date-obs"]
            date_obs_V = self.results_tab_df_colorV[self.results_tab_df_colorV["check_star"] == True].iloc[0]["date-obs"]
            
            # add EXPOSURE/2 
            half_exposure_B = (self.results_tab_df_colorB[self.results_tab_df_colorB["check_star"] == True].iloc[0]["exposure"])/2
            half_exposure_V = (self.results_tab_df_colorV[self.results_tab_df_colorV["check_star"] == True].iloc[0]["exposure"])/2

            #use Julian Datw
            _obs_B = Time(date_obs_B, format='jd') + TimeDelta(half_exposure_B, format='sec')
            _obs_V = Time(date_obs_V, format='jd') + TimeDelta(half_exposure_V, format='sec')

            date_obs_B = _obs_B.jd
            date_obs_V = _obs_V.jd

            #AIRMASS
            amass_B = self.results_tab_df_colorB[self.results_tab_df_colorB["check_star"] == True].iloc[0]["AMASS"]
            amass_V = self.results_tab_df_colorV[self.results_tab_df_colorV["check_star"] == True].iloc[0]["AMASS"]

            """
            Build the result_check_star and result_var_table which are similar to the aforementioed 
            spreadsheet ProcessingMaoImages_202281V1117Her.xlsx
            """
            
            if input_color == 'B-V':
                result_check_star = pd.DataFrame(columns=["type", "name", "comp", "IMB", "IMV", "B", "V", "delta_b_minus_v", "delta_B_minus_V",
                                                      "delta_b", "delta_v", "comp_b_minus_v", "B_star", "V_star", "outlier"])
            
                result_var_star = pd.DataFrame(columns=["type", "name", "comp", "IMB", "IMV", "B", "V", "delta_b_minus_v", "delta_B_minus_V",
                                                    "delta_b", "delta_v", "comp_b_minus_v", "B_star", "V_star"])
            elif input_color == 'V-R':
                result_check_star = pd.DataFrame(columns=["type", "name", "comp", "IMV", "IMR", "V", "R", "delta_v_minus_r", "delta_V_minus_R",
                                                      "delta_v", "delta_r", "comp_v_minus_r", "V_star", "R_star", "outlier"])
            
                result_var_star = pd.DataFrame(columns=["type", "name", "comp", "IMV", "IMR", "V", "R", "delta_v_minus_r", "delta_V_minus_R",
                                                      "delta_v", "delta_r", "comp_v_minus_r", "V_star", "R_star"])
            elif input_color == 'V-I':
                result_check_star = pd.DataFrame(columns=["type", "name", "comp", "IMV", "IMI", "V", "I", "delta_v_minus_i", "delta_V_minus_I",
                                                      "delta_v", "delta_i", "comp_v_minus_i", "V_star", "I_star", "outlier"])
            
                result_var_star = pd.DataFrame(columns=["type", "name", "comp", "IMV", "IMI", "V", "I", "delta_v_minus_i", "delta_V_minus_I",
                                                      "delta_v", "delta_i", "comp_v_minus_i", "V_star", "I_star"])
            else:
                raise Exception("two_color_photometry: unknown imput_color entered")

            

            """
            loop through all the selected comp stars and fill the * table
            """
            sel_comps = [] #init
            sel_comps_to_use = self.object_sel_comp_entry.get().strip().split(',')

            #make array of int csalled sel_comps            
            for comp in sel_comps_to_use:
                sel_comps.append(comp.strip())
            
            for comp in sel_comps:

                #Dont use the check star 
                if comp == check_star_label:
                    continue

                #selected comp must be in both tables
                if str(comp) in self.results_tab_df_colorB["label"].values:
                    comp_star_B = self.results_tab_df_colorB[self.results_tab_df_colorB["label"] == str(comp)].iloc[0]
                elif int(comp) in self.results_tab_df_colorB["label"].values:
                    comp_star_B = self.results_tab_df_colorB[self.results_tab_df_colorB["label"] == int(comp)].iloc[0]
                else:
                    self.console_msg("Comp star: "+ str(int(comp)) + " not in " + first_filter[input_color] + " table")
                    continue

                if str(comp) in self.results_tab_df_colorV["label"].values:
                    comp_star_V = self.results_tab_df_colorV[self.results_tab_df_colorV["label"] == str(comp)].iloc[0]
                elif int(comp) in self.results_tab_df_colorV["label"].values:
                    comp_star_V = self.results_tab_df_colorV[self.results_tab_df_colorV["label"] == int(comp)].iloc[0]
                else:
                    self.console_msg("Comp star: "+ str(int(comp)) + " not in " + second_filter[input_color] + " table")
                    continue
                

                comp_b_minus_v = comp_star_B["inst_mag"] - comp_star_V["inst_mag"]
                delta_b_minus_v = (check_IMB - check_IMV) - comp_b_minus_v
                delta_B_minus_V = tbv_coefficient*delta_b_minus_v
                delta_b = check_IMB - comp_star_B["inst_mag"]
                delta_v = check_IMV - comp_star_V["inst_mag"]
                B_star = delta_b + (tb_bv_coefficient*delta_B_minus_V) + float(comp_star_B["match_mag"])
                V_star = delta_v + (tv_bv_coefficient*delta_B_minus_V) + float(comp_star_V["match_mag"])
                
                
                if input_color == 'B-V':
                    result_check_star.loc[len(result_check_star)] =\
                        {
                        "type": "check",
                        "name": check_star_label,
                        "comp": int(comp),
                        "IMB": comp_star_B["inst_mag"],
                        "IMV": comp_star_V["inst_mag"],
                        "B": comp_star_B["match_mag"],
                        "V": comp_star_V["match_mag"],
                        "delta_b_minus_v": delta_b_minus_v,
                        "delta_B_minus_V": delta_B_minus_V,
                        "delta_b": delta_b,
                        "delta_v": delta_v,
                        "comp_b_minus_v": comp_b_minus_v,
                        "B_star": B_star,
                        "V_star": V_star,
                        "outlier": ''
                        }
                elif input_color == 'V-R':
                    result_check_star.loc[len(result_check_star)] =\
                        {
                        "type": "check",
                        "name": check_star_label,
                        "comp": int(comp),
                        "IMV": comp_star_B["inst_mag"],
                        "IMR": comp_star_V["inst_mag"],
                        "V": comp_star_B["match_mag"],
                        "R": comp_star_V["match_mag"],
                        "delta_v_minus_r": delta_b_minus_v,
                        "delta_V_minus_R": delta_B_minus_V,
                        "delta_v": delta_b,
                        "delta_r": delta_v,
                        "comp_v_minus_r": comp_b_minus_v,
                        "V_star": B_star,
                        "R_star": V_star,
                        "outlier": ''
                        }
                elif input_color == 'V-I':
                    result_check_star.loc[len(result_check_star)] =\
                        {
                        "type": "check",
                        "name": check_star_label,
                        "comp": int(comp),
                        "IMV": comp_star_B["inst_mag"],
                        "IMI": comp_star_V["inst_mag"],
                        "V": comp_star_B["match_mag"],
                        "I": comp_star_V["match_mag"],
                        "delta_v_minus_i": delta_b_minus_v,
                        "delta_V_minus_I": delta_B_minus_V,
                        "delta_v": delta_b,
                        "delta_i": delta_v,
                        "comp_v_minus_i": comp_b_minus_v,
                        "V_star": B_star,
                        "I_star": V_star,
                        "outlier": ''
                        }
                else:
                    raise Exception("two_color_photometry: unknown imput_color entered")

                #
                # NOTE! Since B-V was implemented first, the B-V filternames are
                # still used even though imput_color may be V-R; 
                # Only when it counts does the B-V change to the real V-R or V-I
                #

                # comp_b_minus_v , calculated above
                delta_b_minus_v = (var_IMB - var_IMV) - comp_b_minus_v
                delta_B_minus_V = tbv_coefficient*delta_b_minus_v
                delta_b = var_IMB - comp_star_B["inst_mag"]
                delta_v = var_IMV - comp_star_V["inst_mag"]
                B_star = delta_b + (tb_bv_coefficient*delta_B_minus_V) + float(comp_star_B["match_mag"])
                V_star = delta_v + (tv_bv_coefficient*delta_B_minus_V) + float(comp_star_V["match_mag"])
                
                
                if input_color == 'B-V':
                    result_var_star.loc[len(result_var_star)] =\
                        {
                        "type": "var",
                        "name": var_star_label,
                        "comp": int(comp),
                        "IMB": comp_star_B["inst_mag"],
                        "IMV": comp_star_V["inst_mag"],
                        "B": comp_star_B["match_mag"],
                        "V": comp_star_V["match_mag"],
                        "delta_b_minus_v": delta_b_minus_v,
                        "delta_B_minus_V": delta_B_minus_V,
                        "delta_b": delta_b,
                        "delta_v": delta_v,
                        "comp_b_minus_v": comp_b_minus_v,
                        "B_star": B_star,
                        "V_star": V_star
                        }
                elif input_color == 'V-R':
                    result_var_star.loc[len(result_var_star)] =\
                        {
                        "type": "var",
                        "name": var_star_label,
                        "comp": int(comp),
                        "IMV": comp_star_B["inst_mag"],
                        "IMR": comp_star_V["inst_mag"],
                        "V": comp_star_B["match_mag"],
                        "R": comp_star_V["match_mag"],
                        "delta_v_minus_r": delta_b_minus_v,
                        "delta_V_minus_R": delta_B_minus_V,
                        "delta_v": delta_b,
                        "delta_r": delta_v,
                        "comp_v_minus_r": comp_b_minus_v,
                        "V_star": B_star,
                        "R_star": V_star
                        }
                elif input_color == 'V-I':
                    result_var_star.loc[len(result_var_star)] =\
                        {
                        "type": "var",
                        "name": var_star_label,
                        "comp": int(comp),
                        "IMV": comp_star_B["inst_mag"],
                        "IMI": comp_star_V["inst_mag"],
                        "V": comp_star_B["match_mag"],
                        "I": comp_star_V["match_mag"],
                        "delta_v_minus_i": delta_b_minus_v,
                        "delta_V_minus_I": delta_B_minus_V,
                        "delta_v": delta_b,
                        "delta_i": delta_v,
                        "comp_v_minus_i": comp_b_minus_v,
                        "V_star": B_star,
                        "I_star": V_star
                        }
                else:
                    raise Exception("two_color_photometry: unknown imput_color entered")

            #
            # NOTE! Since B-V was implemented first, the B-V filternames are
            # still used even though imput_color may be V-R; 
            # Only when it counts does the B-V change to the real V-R or V-I
            #
          
            B_mean_check = result_check_star[first_filter[input_color] + "_star"].mean()
            V_mean_check = result_check_star[second_filter[input_color] + "_star"].mean()
            B_mean_var = result_var_star[first_filter[input_color] + "_star"].mean()
            V_mean_var = result_var_star[second_filter[input_color] + "_star"].mean()
            
            B_std_check = result_check_star[first_filter[input_color] + "_star"].std()
            V_std_check = result_check_star[second_filter[input_color] + "_star"].std()
            B_std_var = result_var_star[first_filter[input_color] + "_star"].std()
            V_std_var = result_var_star[second_filter[input_color] + "_star"].std()

            #Find IQR and iqr 
            q3, q1 = np.percentile(result_check_star[first_filter[input_color] + "_star"], [75 ,25])
            b_iqr = q3 - q1
            b_upper_limit = q3 + (b_iqr * 1.5)
            b_lower_limit = q1 - (b_iqr * 1.5)

            q3, q1 = np.percentile(result_check_star[second_filter[input_color] + "_star"], [75 ,25])
            v_iqr = q3 - q1
            v_upper_limit = q3 + (v_iqr * 1.5)
            v_lower_limit = q1 - (v_iqr * 1.5)

            #record any outliers
            
            if (result_check_star[first_filter[input_color] + "_star"] < b_lower_limit).any():
                result_check_star.loc[result_check_star.loc[result_check_star[first_filter[input_color] + "_star"] < b_lower_limit].index, "outlier"] = "<--OUTLIER"
            
            if (result_check_star[first_filter[input_color] + "_star"] > b_upper_limit).any():
                result_check_star.loc[result_check_star.loc[result_check_star[first_filter[input_color] + "_star"] > b_upper_limit].index, "outlier"] = "<--OUTLIER"
            
            if (result_check_star[second_filter[input_color] + "_star"] < v_lower_limit).any():
                result_check_star.loc[result_check_star.loc[result_check_star[second_filter[input_color] + "_star"] < v_lower_limit].index, "outlier"] = "<--OUTLIER"
            
            if (result_check_star[second_filter[input_color] + "_star"] > v_upper_limit).any():
                result_check_star.loc[result_check_star.loc[result_check_star[second_filter[input_color] + "_star"] > v_upper_limit].index, "outlier"] = "<--OUTLIER"
            
            self.console_msg("\n")
            #
            # NOTE! Since B-V was implemented first, the B-V filternames are
            # still used even though imput_color may be V-R; 
            # Only when it counts does the B-V change to the real V-R or V-I, 
            # like in the following..

            if input_color == 'B-V':
                #create an aux table containing misc data needed for AAVSO report
                #this data is appended to the notes section 
                #(See E:\Astronomy\AAVSO\Reports\AAVSO Reports\MAO\2022 8 1 V1117 Her\AAVSOReport_V1117-Her_B_20220802.txt)
                result_aux_report = pd.DataFrame(columns=["color", "JD", "KMAGS", "KMAGINS", "KREFMAG", "T_bv", "Tv_bv", "VMAGINS", "Date-Obs", "KNAME", "AMASS"])
                
                result_aux_report.loc[len(result_aux_report)] =\
                    {
                    "color" : "B",
                    "JD" : comp_star_B["date-obs"], #orig
                    "KMAGS" : B_mean_check,
                    "KMAGINS" : check_IMB,
                    "KREFMAG" : check_B,
                    "T_bv" : tbv_coefficient,
                    "Tv_bv" : tv_bv_coefficient,
                    "VMAGINS" : var_IMB,
                    "Date-Obs" : date_obs_B, #orig plus EXPOSURE/2
                    "KNAME" : check_star_label,
                    "AMASS" : amass_B
                    }
                
                result_aux_report.loc[len(result_aux_report)] =\
                    {
                    "color" : "V",
                    "JD" : comp_star_V["date-obs"], #orig 
                    "KMAGS" : V_mean_check,
                    "KMAGINS" : check_IMV,
                    "KREFMAG" : check_V,
                    "T_bv" : tbv_coefficient,
                    "Tv_bv" : tv_bv_coefficient,
                    "VMAGINS" : var_IMV,
                    "Date-Obs" : date_obs_V,  #orig plus EXPOSURE/2
                    "KNAME" : check_star_label,
                    "AMASS" : amass_V
                    }
            
                self.console_msg("Check Star Estimates using check star: " + str(check_star_label) + " (B: " + str(check_B) +")" + " (V: " + str(check_V) +")" "\n" +
                                result_check_star.sort_values(by="name").to_string() +
                                '\n' +
                                ("B* Ave: " + format(B_mean_check, ' >6.3f') +
                                "  V* Ave: " + format(V_mean_check, ' >6.3f')).rjust(137) +
                                '\n' +
                                ("B* Std: " + format(B_std_check, ' >6.3f') +
                                "  V* Std: " + format(V_std_check, ' >6.3f')).rjust(137))

                self.console_msg(("Check Star IQR limit for B*: " + format(b_lower_limit, ' >6.3f') + ';' + format(b_upper_limit, ' >6.3f')).rjust(123))
                self.console_msg(("Check Star IQR limit for V*: " + format(v_lower_limit, ' >6.3f') + ';' + format(v_upper_limit, ' >6.3f')).rjust(123))
                self.console_msg('\n')

                self.console_msg("Variable Star Estimates of Var: " + var_star_B["vsx_id"] + "\n" +
                                result_var_star.sort_values(by="name").to_string() +
                                '\n' + 
                                ("B* Ave: " + format(B_mean_var, ' >6.3f') +
                                "  V* Ave: " + format(V_mean_var, ' >6.3f')).rjust(137) +
                                '\n' +
                                ("B* Std: " + format(B_std_var, ' >6.3f') +  
                                "  V* Std: " + format(V_std_var, ' >6.3f')).rjust(137))
            elif input_color == 'V-R':
                #create an aux table containing misc data needed for AAVSO report
                #this data is appended to the notes section 
                #(See E:\Astronomy\AAVSO\Reports\AAVSO Reports\MAO\2022 8 1 V1117 Her\AAVSOReport_V1117-Her_B_20220802.txt)
                result_aux_report = pd.DataFrame(columns=["color", "JD", "KMAGS", "KMAGINS", "KREFMAG", "T_vr", "Tr_vr", "VMAGINS", "Date-Obs", "KNAME", "AMASS"])
                
                result_aux_report.loc[len(result_aux_report)] =\
                    {
                    "color" : "V",
                    "JD" : comp_star_B["date-obs"], #orig 
                    "KMAGS" : B_mean_check,
                    "KMAGINS" : check_IMB,
                    "KREFMAG" : check_B,
                    "T_vr" : tbv_coefficient,
                    "Tr_vr" : tv_bv_coefficient,
                    "VMAGINS" : var_IMB,
                    "Date-Obs" : date_obs_B, #orig plus EXPOSURE/2
                    "KNAME" : check_star_label,
                    "AMASS" : amass_B
                    }
                
                result_aux_report.loc[len(result_aux_report)] =\
                    {
                    "color" : "R",
                    "JD" : comp_star_V["date-obs"], #orig 
                    "KMAGS" : V_mean_check,
                    "KMAGINS" : check_IMV,
                    "KREFMAG" : check_V,
                    "T_vr" : tbv_coefficient,
                    "Tr_vr" : tv_bv_coefficient,
                    "VMAGINS" : var_IMV,
                    "Date-Obs" : date_obs_V, #orig plus EXPOSURE/2
                    "KNAME" : check_star_label,
                    "AMASS" : amass_V
                    }
            
                self.console_msg("Check Star Estimates using check star: " + str(check_star_label) + " (V: " + str(check_B) +")" + " (R: " + str(check_V) +")" "\n" +
                                result_check_star.sort_values(by="name").to_string() +
                                '\n' +
                                ("V* Ave: " + format(B_mean_check, ' >6.3f') +
                                "  R* Ave: " + format(V_mean_check, ' >6.3f')).rjust(137) +
                                '\n' +
                                ("V* Std: " + format(B_std_check, ' >6.3f') +
                                "  R* Std: " + format(V_std_check, ' >6.3f')).rjust(137))

                self.console_msg(("Check Star IQR limit for V*: " + format(b_lower_limit, ' >6.3f') + ';' + format(b_upper_limit, ' >6.3f')).rjust(123))
                self.console_msg(("Check Star IQR limit for R*: " + format(v_lower_limit, ' >6.3f') + ';' + format(v_upper_limit, ' >6.3f')).rjust(123))
                self.console_msg('\n')

                self.console_msg("Variable Star Estimates of Var: " + var_star_B["vsx_id"] + "\n" +
                                result_var_star.sort_values(by="name").to_string() +
                                '\n' + 
                                ("V* Ave: " + format(B_mean_var, ' >6.3f') +
                                "  R* Ave: " + format(V_mean_var, ' >6.3f')).rjust(137) +
                                '\n' +
                                ("V* Std: " + format(B_std_var, ' >6.3f') +  
                                "  R* Std: " + format(V_std_var, ' >6.3f')).rjust(137))
            elif input_color == 'V-I':
                #create an aux table containing misc data needed for AAVSO report
                #this data is appended to the notes section 
                #(See E:\Astronomy\AAVSO\Reports\AAVSO Reports\MAO\2022 8 1 V1117 Her\AAVSOReport_V1117-Her_B_20220802.txt)
                result_aux_report = pd.DataFrame(columns=["color", "JD", "KMAGS", "KMAGINS", "KREFMAG", "T_vi", "Ti_vi", "VMAGINS", "Date-Obs", "KNAME", "AMASS"])
                
                result_aux_report.loc[len(result_aux_report)] =\
                    {
                    "color" : "V",
                    "JD" : comp_star_B["date-obs"], #orig 
                    "KMAGS" : B_mean_check,
                    "KMAGINS" : check_IMB,
                    "KREFMAG" : check_B,
                    "T_vi" : tbv_coefficient,
                    "Ti_vi" : tv_bv_coefficient,
                    "VMAGINS" : var_IMB,
                    "Date-Obs" : date_obs_B, #orig plus EXPOSURE/2
                    "KNAME" : check_star_label,
                    "AMASS" : amass_B
                    }
                
                result_aux_report.loc[len(result_aux_report)] =\
                    {
                    "color" : "I",
                    "JD" : comp_star_V["date-obs"], #orig 
                    "KMAGS" : V_mean_check,
                    "KMAGINS" : check_IMV,
                    "KREFMAG" : check_V,
                    "T_vi" : tbv_coefficient,
                    "Ti_vi" : tv_bv_coefficient,
                    "VMAGINS" : var_IMV,
                    "Date-Obs" : date_obs_V, #orig plus EXPOSURE/2
                    "KNAME" : check_star_label,
                    "AMASS" : amass_V
                    }
            
                self.console_msg("Check Star Estimates using check star: " + str(check_star_label) + " (V: " + str(check_B) +")" + " (I: " + str(check_V) +")" "\n" +
                                result_check_star.sort_values(by="name").to_string() +
                                '\n' +
                                ("V* Ave: " + format(B_mean_check, ' >6.3f') +
                                "  I* Ave: " + format(V_mean_check, ' >6.3f')).rjust(137) +
                                '\n' +
                                ("V* Std: " + format(B_std_check, ' >6.3f') +
                                "  I* Std: " + format(V_std_check, ' >6.3f')).rjust(137))

                self.console_msg(("Check Star IQR limit for V*: " + format(b_lower_limit, ' >6.3f') + ';' + format(b_upper_limit, ' >6.3f')).rjust(123))
                self.console_msg(("Check Star IQR limit for I*: " + format(v_lower_limit, ' >6.3f') + ';' + format(v_upper_limit, ' >6.3f')).rjust(123))
                self.console_msg('\n')

                self.console_msg("Variable Star Estimates of Var: " + var_star_B["vsx_id"] + "\n" +
                                result_var_star.sort_values(by="name").to_string() +
                                '\n' + 
                                ("V* Ave: " + format(B_mean_var, ' >6.3f') +
                                "  I* Ave: " + format(V_mean_var, ' >6.3f')).rjust(137) +
                                '\n' +
                                ("V* Std: " + format(B_std_var, ' >6.3f') +  
                                "  I* Std: " + format(V_std_var, ' >6.3f')).rjust(137))
            else:
                raise Exception("two_color_photometry: unknown imput_color entered")
                             
            #concat result_aux_report, result_check_star, and result_var_star
            master_frames = [result_check_star, result_var_star, result_aux_report]

            master_report = pd.concat(master_frames, keys=['check', 'var', 'aux'])
            
            dir_path = os.path.dirname(os.path.realpath(file_name)) + "\\"
            
            master_report.to_csv(dir_path + self.object_name_entry.get() + "-" + input_color + "-Master-Report.csv", index=False)
            self.console_msg("Master Report saved to " + str(dir_path + self.object_name_entry.get() + "-" + input_color + "-Master-Report.csv"))

            pass
            
            
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
     


    def get_comparison_stars(self):
        global image_width, image_height
        comp_stars_used = [] #init
        try:
            self.filter = self.filter_entry.get()

            frame_center = self.wcs_header.pixel_to_world(
                int(image_width / 2), (image_height / 2))

            if not hasattr(frame_center, "ra"):
                #file needs to be plate solved first
                self.console_msg("Image needs to be plate solved first! execute 'Photometry->Solve Image'")
                self.console_msg("Ready")
                return

            frame_center_coordinates = SkyCoord(
                ra=frame_center.ra, dec=frame_center.dec)
            frame_edge = self.wcs_header.pixel_to_world(
                int(image_width), (image_height / 2))
            frame_edge_coordinates = SkyCoord(
                ra=frame_edge.ra, dec=frame_edge.dec)
            frame_radius = frame_edge.separation(frame_center)
            
            frame_top_left = self.wcs_header.pixel_to_world(0, 0)
            frame_top_right = self.wcs_header.pixel_to_world(image_width, 0)

            fov_horizontal = frame_top_right.separation(frame_top_left).arcminute

            self.console_msg("Updating photometry table with Sky coordinates...")
            for index, row in self.results_tab_df.iterrows():
                sky = self.wcs_header.pixel_to_world(row["x_fit"], row["y_fit"])
                c = SkyCoord(ra=sky.ra, dec=sky.dec)
                self.results_tab_df.loc[index, "ra_fit"] = c.ra / u.deg
                self.results_tab_df.loc[index, "dec_fit"] = c.dec / u.deg
            self.results_tab_df.to_csv(self.image_file + ".csv", index=False)
            self.console_msg("Photometry table saved to " + str(self.image_file + ".csv"))


            catalog_selection = self.catalog_stringvar.get()


            """
                Most if not all the time we will be using AAVSO catalog for comp and check stars
            """
    
            if catalog_selection == "AAVSO":
                using_aavso_catalog = True
                
                #Here because we want to use the comparison stars that the 
                #AAVSO has picked out; not all APASS DR9 'comp' stars
                #have a AUID
                self.console_msg("Getting AAVSO Comparison Stars...")
                comparison_stars = \
                    self.aavso_get_comparison_stars(frame_center, filter_band=str(
                    self.filter_entry.get()),
                    field_of_view=fov_horizontal,
                    maglimit=self.max_ensemble_magnitude_entry.get())
                
                ra_column_name = "RA"
                dec_column_name = "Dec"
                mag_column_name = "Mag"
                
 
            
            else:
                using_aavso_catalog = False
                
                #Use VizierR catalog
                if catalog_selection == "APASS DR9":
                    catalog = "II/336"
                    #catalog = "J/other/JAD/20.4"
                    ra_column_name = "RAJ2000"
                    dec_column_name = "DEJ2000"
    
                elif catalog_selection == "URAT1":
                    catalog = "I/329"
                    ra_column_name = "RAJ2000"
                    dec_column_name = "DEJ2000"
    
                elif catalog_selection == "USNO-B1.0":
                    catalog = "I/284"
                    ra_column_name = "RAJ2000"
                    dec_column_name = "DEJ2000"
    
                elif catalog_selection == "VizieR Catalog":
                    catalog = self.vizier_catalog_entry.get()
                    ra_column_name = "RAJ2000"
                    dec_column_name = "DEJ2000"
    
                elif catalog_selection == "Gaia DR2":
                    catalog = "I/345"
                    ra_column_name = "RA_ICRS"
                    dec_column_name = "DE_ICRS"
    
                else:
                    self.console_msg("UNEXPECTED ERROR UNKOWN CATALOG SELECTION!!!!!", level=logging.ERROR)
                    return

                self.console_msg(
                    "Inquiring VizieR Catalog " + catalog_selection + ", center α δ " +
                        frame_center.to_string("hmsdms") +
                        ", radius " + str(frame_radius))
    
                
                mag_column_name = self.filter + "mag"
    
                comparison_stars = Vizier(
                    catalog=catalog, row_limit=-1).query_region(frame_center, frame_radius)[0]


                # print(comparison_stars)
                if mag_column_name not in comparison_stars.colnames:
                    self.console_msg(
                        "Catalog " + self.catalog_stringvar.get() + " does not list " + self.filter + " magnitudes.")
                    return

            if len(comparison_stars) == 0:
                self.console_msg(
                    "NO Comparison stars found in the field; make sure filter and chart Id is correct.")
                return
            else:                
                self.console_msg(
                    "Found " + str(len(comparison_stars)) + " objects in the field.")


            self.console_msg("Matching image to catalog...")
            matching_radius = float(
                self.matching_radius_entry.get()) * 0.000277778  # arcsec to degrees
            
            if using_aavso_catalog:
                catalog_comparison = SkyCoord(
                    comparison_stars[ra_column_name],
                    comparison_stars[dec_column_name],
                    unit=(u.hourangle, u.deg))
            else:
                catalog_comparison = SkyCoord(
                    comparison_stars[ra_column_name],
                    comparison_stars[dec_column_name])

                
            for index, row in self.results_tab_df.iterrows():
                photometry_star_coordinates = SkyCoord(
                    ra=row["ra_fit"] * u.deg, dec=row["dec_fit"] * u.deg)
                match_index, d2d_match, d3d_match = photometry_star_coordinates.match_to_catalog_sky(
                    catalog_comparison)
                # print(str(photometry_star_coordinates))
                # print(match_index)
                # Name of the catalog
                if using_aavso_catalog:
                    match_id = comparison_stars.iloc[match_index]["AUID"]
                    match_label = comparison_stars.iloc[match_index]["Label"]
                    match_ra = catalog_comparison[match_index].ra.degree
                    match_dec= catalog_comparison[match_index].dec.degree
                    match_mag = comparison_stars.iloc[match_index][mag_column_name]
                    match_is_check = comparison_stars.iloc[match_index]["Check Star"]
                else:
                    match_id = comparison_stars[match_index][0]
                    match_label = ""
                    match_ra = comparison_stars[match_index][ra_column_name]
                    match_dec = comparison_stars[match_index][dec_column_name]
                    match_mag = comparison_stars[match_index][mag_column_name]
                    
                    
                match_coordinates = SkyCoord(
                    ra=match_ra * u.deg, dec=match_dec * u.deg)
                separation = photometry_star_coordinates.separation(
                    match_coordinates)
                

                if separation < matching_radius * u.deg:
                    # Sometimes the flux_fit is negative.
                    # That causes a blank inst_mag (can't take a log of neg number) 
                    # Check for this and if so, ignore
                    # 
                    if self.results_tab_df.loc[index, "flux_fit"] < 0:
                        self.results_tab_df.loc[index, "check_star"] = False # all values in this column must be of boolean
                        continue

                    # Make sure the earliest iteration is used
                    # The first iteration is not necessarily iter_detected = 1
                    # The first or earliest detection could be any int up to max iterations
                    
                    # "iter_detected" may not exisit becase Non-iterative PSF Photometry was used
                    if "label" in self.results_tab_df and "iter_detected" in self.results_tab_df:
                        already_gotten = self.results_tab_df.loc[self.results_tab_df["label"] == int(match_label)]    
                        if not already_gotten.empty:
                            # Confirm that this already gotten one has an EARLIER iteration,
                            # and if so, then we can continue, else
                            # we have to erase already gotten and use the new one
                            if already_gotten.iloc[0]["iter_detected"] < self.results_tab_df.loc[index, "iter_detected"]:
                                self.results_tab_df.loc[index, "label"] = ""
                                self.results_tab_df.loc[index, "check_star"] = False # all values in this column must be of boolean
                                continue
                            else:
                                #erase already_gotten label 
                                already_gotten.iloc[0]['label'] = ""
                                # fall thru and label the new index 

                    #Here if not already_gotten or the already_gotten has 
                    # a "later/stale" iteration
                    self.results_tab_df.loc[index, "match_id"] = \
                        str(self.catalog_stringvar.get()) + \
                            " " + str(match_id)
                    self.results_tab_df.loc[index, "label"] = int(match_label)
                    self.results_tab_df.loc[index, "match_ra"] = match_ra
                    self.results_tab_df.loc[index, "match_dec"] = match_dec
                    self.results_tab_df.loc[index, "match_mag"] = match_mag
                    self.results_tab_df.loc[index, "check_star"] = match_is_check
                    
                    #record comp stars used for console if AAVSO comp stars
                    comp_stars_used.append((match_label, match_is_check))
                    
                else:
                    #Here if separation >= matching_radius
                    self.results_tab_df.loc[index, "match_id"] = ""
                    self.results_tab_df.loc[index, "label"] = ""
                    self.results_tab_df.loc[index, "match_ra"] = ""
                    self.results_tab_df.loc[index, "match_dec"] = ""
                    self.results_tab_df.loc[index, "match_mag"] = ""
                    self.results_tab_df.loc[index, "check_star"] = False # all values in this column must be of boolean 

            self.console_msg(
                "Inquiring VizieR (B/vsx/vsx) for VSX variables in the field...")
            vsx_result = Vizier(catalog="B/vsx/vsx", row_limit=-1).query_region(frame_center, frame_radius)
            # print(vsx_result)

            """
            Look for any and all VSX stars
            """
            if len(vsx_result) > 0:
                vsx_stars = vsx_result[0]
                self.console_msg(
                    "Found " + str(len(vsx_stars)) + " VSX sources in the field. Matching...")
                catalog_vsx = SkyCoord(
                    vsx_stars["RAJ2000"], vsx_stars["DEJ2000"])
                for index, row in self.results_tab_df.iterrows():
                    photometry_star_coordinates = SkyCoord(
                        ra=row["ra_fit"] * u.deg, dec=row["dec_fit"] * u.deg)
                    match_index, d2d_match, d3d_match = photometry_star_coordinates.match_to_catalog_sky(
                        catalog_vsx)
                    
                    match_id = vsx_stars[match_index]["Name"]
                    match_ra = vsx_stars[match_index]["RAJ2000"]
                    match_dec = vsx_stars[match_index]["DEJ2000"]
                    match_coordinates = SkyCoord(
                        ra=match_ra * u.deg, dec=match_dec * u.deg)
                    separation = photometry_star_coordinates.separation(
                        match_coordinates)
                    if separation < matching_radius * u.deg:
                        # Sometimes the flux_fit is negative out of the IterativelySubtractedPSFPhotometry.
                        # That causes a blank inst_mag (can't take a log of neg number) 
                        # Check for this and if so, ignore
                        # (IterativelySubtractedPSFPhotometry is now deprecated and replaced by 
                        # IterativePSFPhotometry; not known if following is still needed)
                        # 
                        if self.results_tab_df.loc[index, "flux_fit"] < 0:
                            continue

                        # Make sure the earliest iteration is used
                        # The first iteration is not necessarily iter_detected = 1
                        # The first or earliest detection could be any int up to max iterations
                        if "vsx_id" in self.results_tab_df:
                            already_gotten = self.results_tab_df.loc[self.results_tab_df["vsx_id"] == str(match_id)]    
                            if not already_gotten.empty:
                                # Confirm that this already gotten one has an EARLIER iteration,
                                # and if so, then we can continue, else
                                # we have to erase already gotten and use the new one
                                if already_gotten.iloc[0]["iter_detected"] < self.results_tab_df.loc[index, "iter_detected"]:
                                    self.results_tab_df.loc[index, "vsx_id"] = ""
                                    continue
                                else:
                                    #erase already_gotten label 
                                    already_gotten.iloc[0]['vsx_id'] = ""
                                    # fall thru and label the new index 

                        #Found a match...should be earliest iteration
                        self.results_tab_df.loc[index, "vsx_id"] = str(match_id)
                        self.console_msg("Match VSX source: " + str(match_id))
                    else:
                        self.results_tab_df.loc[index, "vsx_id"] = ""
            else:
                self.console_msg("Found no VSX sources in the field.")
                
            self.results_tab_df.to_csv(self.image_file + ".csv", index=False)
            self.console_msg("Photometry table saved (with following comp stars) to " + str(self.image_file + ".csv"))

            if using_aavso_catalog:
                comp_list = ''
                #output comp_stars_used
                found_check = False #init
                for comp in comp_stars_used:
                    (label, ischeck) = comp
                    found_check |= ischeck == True
                    comp_list += str(label) + ', ' 
                self.console_msg("AAVSO comp stars used: " + comp_list)
                    
                check_star = self.object_kref_entry.get().strip()
                if not found_check and check_star != '':
                    self.console_msg("WARNING!! The requested check_star: " + check_star + " was NOT FOUND in image! Choose another.", level=logging.WARNING)
                    
            self.display_image()
            self.console_msg("Ready")
        
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)

    def set_entry_text(self, entry, text):
        entry.delete(0, tk.END)
        entry.insert(0, text)

    def update_display(self):
        self.display_image()
        self.update_PSF_canvas(self.last_clicked_x, self.last_clicked_y)

    def safe_float_convert(self, x):
        try:
            z = float(x)
            if np.isnan(z):
                return False    # Nan!
            return True  # numeric, success!
        except ValueError:
            return False  # not numeric
        except TypeError:
            return False  # null type

    def update_histogram_low(self, value):
        self.histogram_slider_low = int(value)
        self.update_display()

    def update_histogram_high(self, value):
        self.histogram_slider_high = int(value)
        self.update_display()


    def save_settings_as(self):
        options = {}
        options['defaultextension'] = '.txt'
        options['filetypes'] = [('TXT', '.txt')]
        #options['initialfile'] = ''
        options['title'] = 'Save MAOPhot settings as...'

        file_name = fd.asksaveasfile(**options)

        try:
            if len(str(file_name)) > 0:
                self.console_msg("Saving settings as " + str(file_name.name))
                settings = {}
                '''
                Use the valid_parameter list which contains the official list
                of user parameters
                '''
                mao_parameters = list(self.valid_parameter_list)
                for param in mao_parameters:
                    settings.update({param : self.valid_parameter_list[param].get()})

                with open(str(file_name.name), 'w') as f:
                    w = csv.DictWriter(f, settings.keys())
                    w.writeheader()
                    w.writerow(settings)
                self.console_msg("Saved.")
        
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
        finally:
            # bring the settings window back to the front                
            self.es_top.lift()



    def load_settings(self):
        try:
            options = {}
            options['defaultextension'] = '.txt'
            options['filetypes'] = [('TXT', '.txt')]
            #options['initialfile'] = ''
            options['title'] = 'Load MAOPhot settings...'

            file_name = fd.askopenfilename(**options)
            
            if len(str(file_name)) > 0:
                self.console_msg("Loading settings from " + str(file_name))
                settings = {}
                with open(str(file_name)) as f:
                    r = csv.DictReader(f)
                    for row in r:
                        # dict from OrderedDict, required for Python < 3.8 as DictReader behavior changed
                        row = dict(row)
                        # append settings dictionary with the read row
                        settings.update(row)
                    for key in settings:
                        if key not in self.valid_parameter_list:
                            continue
                        if type(getattr(self, key)) == tk.Entry:
                            self.set_entry_text(
                                getattr(self, key), settings[key])
                        if type(getattr(self, key)) == tk.StringVar:
                            getattr(self, key).set(settings[key])
                        if type(getattr(self, key)) == tk.BooleanVar:
                            getattr(self, key).set(settings[key])

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
        finally:
            # bring the settings window back to the front                
            self.es_top.lift()


    #
    # show_settings; From File->Edit Settings.... 
    #
    def show_settings(self):
        self.es_top.deiconify()

    #
    # edit_setting; command for "File->Edit Settings..." menu item which creates popup window
    #
    def launch_settings(self):
        try:
            es_ = tk.Toplevel(self.window, padx=15, pady=15, takefocus=True)

            es_.geometry(str(int(self.screen_width*.2)) + "x" + str(int(self.screen_height*.55)))
            es_.title("Settings")
            self.es_top = es_

            self.es_top.protocol("WM_DELETE_WINDOW", self.es_top.withdraw)

            self.es_top.minsize(width=int(self.screen_width*.2), height=int(self.screen_height*.55))

            tk.Grid.columnconfigure(self.es_top, 0, weight=1)

            settings_entry_width = 6
            settings_entry_pad = 0
            extended_settings_entry_width = 30
            extended_settings_entry_pad = 0

            row = 0

            photometry_aperture_label = tk.Label(
                self.es_top, text="Fitting Width/Height, px (odd only):")
            photometry_aperture_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.photometry_aperture_entry = tk.Entry(
                self.es_top, width=settings_entry_width)
            self.photometry_aperture_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            max_ensemble_magnitude_label = tk.Label(
                self.es_top, text="Maximum Ensemble Magnitude:")
            max_ensemble_magnitude_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.max_ensemble_magnitude_entry = tk.Entry(
                self.es_top, width=settings_entry_width)
            self.max_ensemble_magnitude_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            fwhm_label = tk.Label(self.es_top, text="FWHM, px:")
            fwhm_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.fwhm_entry = tk.Entry(self.es_top, width=settings_entry_width)
            self.fwhm_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1


            star_detection_threshold_label = tk.Label(
                self.es_top, text="IRAFStarFinder Threshold:")
            star_detection_threshold_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.star_detection_threshold_entry = tk.Entry(
                self.es_top, width=settings_entry_width)
            self.star_detection_threshold_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            photometry_iterations_label = tk.Label(
                self.es_top, text="Photometry Iterations:")
            photometry_iterations_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.photometry_iterations_entry = tk.Entry(
                self.es_top, width=settings_entry_width)
            self.photometry_iterations_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            sharplo_label = tk.Label(
                self.es_top, text="Lower Bound for Sharpness:")
            sharplo_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.sharplo_entry = tk.Entry(
                self.es_top, width=settings_entry_width)
            self.sharplo_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            bkg_filter_size_label = tk.Label(
                self.es_top, text="Background Median Filter, px:")
            bkg_filter_size_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.bkg_filter_size_entry = tk.Entry(
                self.es_top, width=settings_entry_width)
            self.bkg_filter_size_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            matching_radius_label = tk.Label(
                self.es_top, text="Matching Radius, arcsec:")
            matching_radius_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.matching_radius_entry = tk.Entry(
                self.es_top, width=settings_entry_width)
            self.matching_radius_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            astrometrynet_label = tk.Label(
                self.es_top, text="Astrometry.net Server:")
            astrometrynet_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.astrometrynet_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width)
            self.astrometrynet_entry.grid(row=row, column=2, ipadx=extended_settings_entry_pad)
            row += 1

            astrometrynet_key_label = tk.Label(
                self.es_top, text="Astrometry.net API Key:")
            astrometrynet_key_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.astrometrynet_key_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width)
            self.astrometrynet_key_entry.grid(row=row, column=2, ipadx=extended_settings_entry_pad)
            self.astrometrynet_key_entry.config(show="*")
            row += 1

            separator_telescope = ttk.Separator(self.es_top, orient='horizontal')
            separator_telescope.grid(row=row, columnspan=3, pady=5, sticky=tk.EW)
            row += 1

            telescope_label = tk.Label(self.es_top, text="Telescope:")
            telescope_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.telescope_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width)
            self.telescope_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            #
            # Tbv, Tb_bv, Tv_bv
            #

            tbv_label = tk.Label(self.es_top, text="Tbv:")
            tbv_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.tbv_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.tbv_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            tb_bv_label = tk.Label(self.es_top, text="Tb_bv:")
            tb_bv_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.tb_bv_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.tb_bv_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            tv_bv_label = tk.Label(self.es_top, text="Tv_bv:")
            tv_bv_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.tv_bv_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.tv_bv_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            #
            # Tvr, Tr_vr, Tv_vr
            #
            tvr_label = tk.Label(self.es_top, text="Tvr:")
            tvr_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.tvr_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.tvr_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            tv_vr_label = tk.Label(self.es_top, text="Tv_vr:")
            tv_vr_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.tv_vr_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.tv_vr_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            tr_vr_label = tk.Label(self.es_top, text="Tr_vr:")
            tr_vr_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.tr_vr_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.tr_vr_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            #
            # Tvi, Ti_vi, Ti_vi
            #
            tvi_label = tk.Label(self.es_top, text="Tvi:")
            tvi_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.tvi_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.tvi_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            tv_vi_label = tk.Label(self.es_top, text="Tv_vi:")
            tv_vi_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.tv_vi_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.tv_vi_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            ti_vi_label = tk.Label(self.es_top, text="Ti_vi:")
            ti_vi_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.ti_vi_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.ti_vi_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            separator_telescope_ = ttk.Separator(self.es_top, orient='horizontal')
            separator_telescope_.grid(row=row, columnspan=3, pady=5, sticky=tk.EW)
            row += 1

            aavso_obscode_label = tk.Label(
                self.es_top, text="AAVSO Observer Code:")
            aavso_obscode_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.aavso_obscode_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.aavso_obscode_entry.grid(row=row, column=2, ipadx=extended_settings_entry_pad)
            row += 1

            exposure_label = tk.Label(
                self.es_top, text="Exposure Time:")
            exposure_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.exposure_entry = tk.Entry(
                self.es_top, width=settings_entry_width, background='pink')
            self.exposure_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            filter_label = tk.Label(self.es_top, text="CCD Filter:")
            filter_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.filter_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.filter_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            airmass_label = tk.Label(self.es_top, text="Airmass:")
            airmass_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.airmass_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.airmass_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            date_obs_label = tk.Label(
                self.es_top, text="Date-Obs (JD):")
            date_obs_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.date_obs_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.date_obs_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            object_name_label = tk.Label(
                self.es_top, text="Object Name:")
            object_name_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.object_name_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.object_name_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            object_notes_label = tk.Label(self.es_top, text="Notes:")
            object_notes_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.object_notes_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.object_notes_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            object_kref_label = tk.Label(self.es_top, text="Use Check Star (AAVSO Label):")
            object_kref_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.object_kref_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.object_kref_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            object_sel_comp_label = tk.Label(self.es_top, text="Select Comp Stars (AAVSO Label):")
            object_sel_comp_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.object_sel_comp_entry = tk.Entry(
                self.es_top, width=extended_settings_entry_width, background='pink')
            self.object_sel_comp_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            catalog_label = tk.Label(self.es_top, text="Comparison Catalog:")
            catalog_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            catalog_dropdown = tk.OptionMenu(
                self.es_top, self.catalog_stringvar, "AAVSO")
                # , "APASS DR9", "URAT1", "USNO-B1.0", "Gaia DR2", "VizieR Catalog") <--not supporting now
            catalog_dropdown.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            vizier_catalog_label = tk.Label(self.es_top, text="AAVSO ChartID:")
            vizier_catalog_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.vizier_catalog_entry = tk.Entry(self.es_top,width=extended_settings_entry_width)
            self.vizier_catalog_entry.grid(row=row, column=2, sticky=tk.E)
            row += 1

            fitter_label = tk.Label(self.es_top, text="PSF Fitter:")
            fitter_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            fitter_dropdown = tk.OptionMenu(self.es_top, self.fitter_stringvar,
                                                "Levenberg-Marquardt", "Sequential LS Programming", "Simplex LS")
            fitter_dropdown.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            separator_ = ttk.Separator(self.es_top, orient='horizontal')
            separator_.grid(row=row, columnspan=3, pady=5, sticky=tk.EW)
            row += 1

            #
            # Buttons 
            #
            load_settings_button = tk.Button(self.es_top, text="Load...", command=self.load_settings)
            load_settings_button.grid(row=row, column=0, padx=20, sticky=tk.W)
            save_settings_button = tk.Button(self.es_top, text="Save As...", command=self.save_settings_as)
            save_settings_button.grid(row=row, column=1, padx=20) #, sticky=tk.W)
            close_settings_button = tk.Button(self.es_top, text="  OK  ", command=self.es_top.withdraw) #hide
            close_settings_button.grid(row=row, column=2, padx=20, sticky=tk.E)
            row += 1

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)

    def aavso_get_comparison_stars(self, frame_center, filter_band='V', field_of_view=18.5, maglimit=20):
        try:
            #Some telescopes use 'I' instead of 'Ic', but AAVSO charts use Ic
            if filter_band == 'I':
                filter_band = 'Ic'
            if filter_band == 'R':
                filter_band = 'Rc'
            
            ra = frame_center.to_string("hmsdms").split()[0].replace(
                "h", " ").replace("m", " ").replace("s", "")
            dec = frame_center.to_string("hmsdms").split()[1].replace(
                "d", " ").replace("m", " ").replace("s", "").replace("+", "")

            aavso_chartId_to_use = self.vizier_catalog_entry.get().strip()
            
            if aavso_chartId_to_use == "":
                r = requests.get(
                    'https://www.aavso.org/apps/vsp/api/chart/?format=json&fov=' +
                    str(field_of_view) + '&ra=' + str(ra) + '&dec=' + str(dec) +
                    '&maglimit=' + str(maglimit))
            else:
                #Use the following when just specifying known chartID
                r = requests.get('https://app.aavso.org/vsp/api/chart/' + 
                                 aavso_chartId_to_use + '/?format=json')


            if not 'chartid' in r.json():
                self.console_msg("Invalid chartId: " + aavso_chartId_to_use)
                return None


            #See if there are any specific comps to be used only
            #if a check star is specified, then add it to list, it will marked a check star later
            sel_comps = [] #init
            sel_comps_to_use = self.object_sel_comp_entry.get().strip().split(',')
                
            for comp in sel_comps_to_use:
                sel_comps.append(int(comp.strip()))
                
            check_star_to_use = self.object_kref_entry.get().strip()
            if check_star_to_use.isnumeric():
                #add it to the list of comps
                check_star = int(check_star_to_use)
                sel_comps.append(check_star)
            else:
                check_star = -1 #this value never matches a label

            chart_id = r.json()['chartid']
            self.console_msg(
                'Downloaded AAVSO Comparison Star Chart ID ' + str(chart_id))
            
            result = pd.DataFrame(columns=["AUID", "RA", "Dec", "Mag", "Label", "Chart ID", "Check Star"])

            for star in r.json()['photometry']:
                auid = star['auid']
                ra = star['ra']
                dec = star['dec']
                label = int(star['label'])

                if label not in sel_comps:
                    continue #skip
                
                # if this label is the selected check star, mark it
                is_check_star = (label == check_star)
                
                # init mag here because if 
                # chart doesn't have the mag for 
                # this comp or check star, then it is not used
                mag = None
                for band in star['bands']:
                    if band['band'] == filter_band:
                        mag = band['mag']

                if mag == None:
                    self.console_msg("label: " + str(label) + " has no mag for " + filter_band + "..skipping")
                    continue #skip this one
                        
                result.loc[len(result)] = {"AUID": auid,
                    "RA": ra,
                    "Dec": dec,
                    "Mag": mag,
                    "Label": label,
                    "Chart ID": chart_id,
                    "Check Star": is_check_star}

            return result
        
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)

    def generate_aavso_report_1image(self):
        global image_width, image_height

        """
        Typical Single Image Report is shown below (note only 1 check and 1 comp star)

        Requirements to run:
         - Fits file loaded, plate solved and comp and vsx star shown and entered in settings
         If there are more than 1 comp stars in the list 'Select Comp Stars', then the first one in the list
         is used.
        """

        self.console_msg("Beginning Generate AAVSO Single Image Photometry Report...")

        var_star_name = self.object_name_entry.get().strip()
        if len(var_star_name) == 0:
            self.console_msg(
                "'Object Name' must be specified; eg. 'W Her'")
            return

        check_star_name = self.object_kref_entry.get().strip()
        if len(check_star_name) == 0:
            self.console_msg(
                "'Use Check Star (AAVSO Label)' must be specified; eg. '144'")
            return

        comp_star_list = self.object_sel_comp_entry.get().strip()
        if len(comp_star_list) == 0:
            self.console_msg(
                "'Select Comp Stars (AAVSO Label)' must be specified; eg. '144'")
            return

        #if there is a list of comp stars, use the first one in the list
        comp_star_name = comp_star_list.split(",")[0]
        if len(comp_star_name) == 0:
            self.console_msg(
                "Cannot read first comp star in the list as a label; eg. '144'")
            return

        report_dir = "aavso_reports"
        if os.path.isdir(os.path.dirname(self.image_file)):
            os.chdir(os.path.dirname(self.image_file))
            if not os.path.exists(report_dir):
                os.mkdir(report_dir)
        else:
            self.console_msg("Dir path: \"" + str(self.image_file) + "\" does not exist; check if file loaded")
            return

        image_basename = os.path.basename(self.image_file)
        report_filename = os.path.join(report_dir, "AAVSO " + os.path.splitext(
            image_basename)[0] + " " + str(self.object_name_entry.get()) + "_single.txt")
        
        """
        Report is generated from 'Saved' <image_file>.csv 

        """
        #Ask user for <image-file>.csv

        options = {}
        options['defaultextension'] = '.csv'
        options['filetypes'] = [('CSV', '.csv')]
        options['title'] = 'Choose the ' + self.image_file + '.csv'

        file_name = fd.askopenfilename(**options)

        if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
            self.console_msg("Loading Single Image data from " + str(file_name))
            self.results_tab_df_color = pd.read_csv(str(file_name))
        else:
            return

        #Test to make sure csv file is ready                
        if "vsx_id" not in self.results_tab_df_color:
            self.console_msg("Cannot proceed; run 'Photometry->Get Comparison Stars' first.")
            return

        try:
            """
            Typical Single Image Report (note only 1 check and 1 comp star)

            #TYPE=EXTENDED
            #OBSCODE=FPIA
            #SOFTWARE=VPhot 4.0.27
            #DELIM=,
            #DATE=JD
            #OBSTYPE=CCD
            #NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
            W Her,2459735.64795,14.599,0.010,B,NO,STD,144,-6.423,139,-6.659,na,na,X27989RA,
            Mittelman ATMoB Observatory|CMAGINS=-6.423|CREFERR=0.061|CREFMAG=14.933|KMAG=14.697
            |KMAGINS=-6.659|KREFERR=0.037|KREFMAG=14.709|VMAGINS=-6.757
            """

            #Check if the Var to report on has been measured
            if not var_star_name in self.results_tab_df_color["vsx_id"].values:
                self.console_msg("Var not found in table; "+ str(var_star_name) + " not found!", level=logging.ERROR)
                return

            #Check if comp star has been measured
            if not int(comp_star_name) in self.results_tab_df_color["label"].values:
                self.console_msg("Comp star not found in table; "+ str(comp_star_name) + " not found!", level=logging.ERROR)
                return


            with open(report_filename, mode='w') as f:
    
                decimal_places = 3 #report is usually 3
    
                f.write("#TYPE=Extended\n")
                f.write("#OBSCODE="+self.aavso_obscode_entry.get()+"\n")
                f.write("#SOFTWARE=Self-developed using photutils.psf; IRAFStarFinder\n") 
                f.write("#DELIM=,\n")
                f.write("#DATE=JD\n")
                f.write("#OBSTYPE=CCD\n")
                f.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")

                #var_star_name, check_star_name, comp_star_name were determined above
                var_star = self.results_tab_df_color[self.results_tab_df_color["vsx_id"] == var_star_name].iloc[0]
                check_star = self.results_tab_df_color[self.results_tab_df_color["label"] == int(check_star_name)].iloc[0]
                comp_star = self.results_tab_df_color[self.results_tab_df_color["label"] == int(comp_star_name)].iloc[0]

                comp_IM = comp_star["inst_mag"]
                comp_star_mag = comp_star["match_mag"]

                check_IM = check_star["inst_mag"]
                check_star_mag = check_IM - comp_IM + comp_star_mag
                check_star_mag_ref = check_star["match_mag"]

                var_IM = var_star["inst_mag"]

                #differential photometry
                var_star_mag = var_IM - comp_IM + comp_star_mag
                check_star_mag = check_IM - comp_IM + comp_star_mag

                #error estimates
                # MERR
                # use the  Background2D object's median level if it was fetched
                if not self.bkg2D == None:
                    snr_var_star = var_star['flux_fit']/self.bkg2D.background_median
                    snr_check_star = check_star['flux_fit']/self.bkg2D.background_median
                    snr_comp_star = comp_star['flux_fit']/self.bkg2D.background_median
                else:
                    snr_var_star = var_star['flux_fit']/self.image_bkg_value
                    snr_check_star = check_star['flux_fit']/self.image_bkg_value
                    snr_comp_star = comp_star['flux_fit']/self.image_bkg_value

                var_star_err = 2.5*np.log10(1 + 1/snr_var_star)
                check_star_err = 2.5*np.log10(1 + 1/snr_check_star)
                comp_star_err = 2.5*np.log10(1 + 1/snr_comp_star)

                err_in_quadrature = math.sqrt(var_star_err**2 + check_star_err**2)

                starid = var_star_name
                date = format(float(self.date_obs_entry.get()), '.5f') 
                mag = str(round(var_star_mag, decimal_places))
                merr = str(round(err_in_quadrature, decimal_places))
                filt = self.filter_entry.get().strip()
                trans = "NO"
                mtype = "STD"
                cname = comp_star_name
                cmag = str(round(comp_IM, decimal_places))
                kname = check_star_name
                kmag = str(round(check_IM, decimal_places)) #not same as KMAG in notes
                amass = self.airmass_entry.get().strip() if self.airmass_entry.get().isnumeric() else "na"
                group = "na"
                chart = self.vizier_catalog_entry.get().strip()
                notes = self.object_notes_entry.get().strip()
                notes += "|CMAGINS=" + cmag + \
                         "|CREFERR=" + str(round(comp_star_err, decimal_places)) +\
                         "|CREFMAG=" + str(round(comp_star_mag, decimal_places)) +\
                         "|KMAG=" + str(round(check_star_mag, decimal_places)) +\
                         "|KMAGINS=" + str(round(check_IM, decimal_places)) +\
                         "|KREFERR=" + str(round(check_star_err, decimal_places)) +\
                         "|KREFMAG=" + str(round(check_star_mag_ref, decimal_places)) +\
                         "|VMAGINS=" + str(round(var_IM, decimal_places))
                
                # Add " " after notes, because TA clobbers last char
                f.write(starid+","+date+","+mag+","+merr+","+filt+","+trans+","+mtype+"," +
                        cname+","+cmag+","+kname+","+kmag+","+amass+","+group+","+chart+","+notes+" \n")

            self.console_msg(
                "AAVSO Photometry report saved to " + str(report_filename))

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)

        return
        
    #
    #   BV_generate_aavso_report_2color: called from Menu selection 
    #   Generate AAVSO Report->Two Color Photometry->(B-V)
    #
    def BV_generate_aavso_report_2color(self):
        self.generate_aavso_report_2color('B-V')

    #
    #   VR_generate_aavso_report_2color: called from Menu selection 
    #   Generate AAVSO Report->Two Color Photometry->(V-R)
    #
    def VR_generate_aavso_report_2color(self):
        self.generate_aavso_report_2color('V-R')

    #
    #   VI_generate_aavso_report_2color: called from Menu selection 
    #   Generate AAVSO Report->Two Color Photometry->(V-I)
    #
    def VI_generate_aavso_report_2color(self):
        self.generate_aavso_report_2color('V-I')

    #
    # generate_aavso_report_2color; this generates a AAVSO report in extended format
    #

    def generate_aavso_report_2color(self, input_color):
        global image_width, image_height
        
        """
        Typical AAVSO Two Color Report:
            
        #TYPE=EXTENDED
        #OBSCODE=FPIA
        #SOFTWARE=VPhot 4.0.44
        #DELIM=,
        #DATE=JD
        #OBSTYPE=CCD
        #NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
        Z Tau,2460300.57931,13.167,0.018,V,YES,STD,ENSEMBLE,na,136,13.592,1.362,na,X29320NP,Mittelman ATMoB Observatory|KMAGINS=-8.083|KMAGSTD=13.592|KREFMAG=13.585|Tvi=1.194|VMAGINS=-8.606
        Z Tau,2460300.58965,8.268,0.023,I,YES,STD,ENSEMBLE,na,136,12.722,1.312,na,X29320NP,Mittelman ATMoB Observatory|KMAGINS=-7.134|KMAGSTD=12.722|KREFMAG=12.731|Ti_vi=-0.138|VMAGINS=-11.032

        Note:
          DATE (in JD) is DATE-OBS + EXPOSURE/2

        """

        # use first_filter, amd second_filter dict to index into appropriate filter
        first_filter = {'B-V': 'B', 'V-R': 'V', 'V-I': 'V'}
        second_filter = {'B-V': 'V', 'V-R': 'R', 'V-I': 'I'}
        t_coefficient_tbv = {'B-V': 'T_bv', 'V-R': 'T_vr', 'V-I': 'T_vi'}
        t_coefficient_tv_bv = {'B-V': 'Tv_bv', 'V-R': 'Tr_vr', 'V-I': 'Ti_vi'}

        t_coefficient_tbv_for_report_only = {'B-V': 'Tbv', 'V-R': 'Tvr', 'V-I': 'Tvi'}

        self.console_msg("Beginning Generate AAVSO Two Color Ensemble Report...")
        var_star_name = self.object_name_entry.get().strip()
        if len(var_star_name) == 0:
            self.console_msg(
                "Object Name must be specified; eg. 'V1117 Her'")
            return

        """
        Report is generated from 'Saved' master_result DataFrame
        that was saved as <object_name_entry>-X-Y-Master-Result.csv
        
        """
        #Ask user for <ObjectID>-Master-Report.csv
        # Ask for the B and V CSV files
        options = {}
        options['defaultextension'] = '.csv'
        options['filetypes'] = [('CSV', '.csv')]
        options['title'] = 'Choose the ' + var_star_name + '-' + input_color + '-Master-Report.csv'

        file_name = fd.askopenfilename(**options)

        if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
            self.console_msg("Loading Master-Report from " + str(file_name))
            master_report = pd.read_csv(str(file_name))
        else:
            return
            
        report_dir = "aavso_reports"
        try:        
            if os.path.isdir(os.path.dirname(str(file_name))):
                os.chdir(os.path.dirname(str(file_name)))
                if not os.path.exists(report_dir):
                    os.mkdir(report_dir)
            else:
                raise Exception("Can't create report_dir: aavso_reports")
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)
            

        image_basename = os.path.basename(self.image_file)
        report_filename = os.path.join(report_dir, "AAVSO " + os.path.splitext(
            image_basename)[0] + " " + str(self.object_name_entry.get()) + ".txt")


        decimal_places = 3 #report is usually 3
        
        #
        # NOTE! Since B-V was implemented first, the B-V filternames are
        # still used even though imput_color may be V-R; 
        # Only when it counts does the B-V change to the real V-R or V-I, 
        # like in the following..

        #extract the estimates
        result_check_star = master_report[master_report["type"] == "check"]
        result_var_star = master_report[master_report["type"] == "var"]

        B_mean_check = result_check_star[first_filter[input_color] + "_star"].mean()
        V_mean_check = result_check_star[second_filter[input_color] + "_star"].mean()
        B_mean_var = result_var_star[first_filter[input_color] + "_star"].mean()
        V_mean_var = result_var_star[second_filter[input_color] + "_star"].mean()
        
        B_std_check = result_check_star[first_filter[input_color] + "_star"].std()
        V_std_check = result_check_star[second_filter[input_color] + "_star"].std()
        B_std_var = result_var_star[first_filter[input_color] + "_star"].std()
        V_std_var = result_var_star[second_filter[input_color] + "_star"].std()
            

        #
        # NOTE! Since B-V was implemented first, the B-V filternames are
        # still used even though imput_color may be V-R; 
        # Only when it counts does the B-V change to the real V-R or V-I, 
        # like in the following..

        # need the Check data label now
        aux_result_first_color = master_report[master_report["color"] == first_filter[input_color]]
        aux_result_second_color = master_report[master_report["color"] == second_filter[input_color]]



        if input_color == 'B-V':
            self.console_msg("Check Star Estimates using check star: " + str(int(aux_result_first_color["KNAME"])) + 
                            " (B: " + str(round(float(aux_result_first_color['KREFMAG']), decimal_places)) +")" + 
                            " (V: " + str(round(float(aux_result_second_color['KREFMAG']), decimal_places)) +")" "\n" + 
                            ("B* Ave: " + format(B_mean_check, ' >6.3f') +
                            "  V* Ave: " + format(V_mean_check, ' >6.3f') +
                            "  (B-V)*: " + format(B_mean_check-V_mean_check, ' >6.3f')).rjust(72) +
                            '\n' +
                            ("B* Std: " + format(B_std_check, ' >6.3f') +
                            "  V* Std: " + format(V_std_check, ' >6.3f')).rjust(56))

            self.console_msg("Variable Star Estimates\n" +
                            ("B* Ave: " + format(B_mean_var, ' >6.3f') +
                            "  V* Ave: " + format(V_mean_var, ' >6.3f') +
                            "  (B-V)*: " + format(B_mean_var-V_mean_var, ' >6.3f')).rjust(72) +
                            '\n' +
                            ("B* Std: " + format(B_std_var, ' >6.3f') +  
                            "  V* Std: " + format(V_std_var, ' >6.3f')).rjust(56))
        elif input_color == 'V-R':
            self.console_msg("Check Star Estimates using check star: " + str(int(aux_result_first_color["KNAME"])) + 
                            " (V: " + str(round(float(aux_result_first_color['KREFMAG']), decimal_places)) +")" + 
                            " (R: " + str(round(float(aux_result_second_color['KREFMAG']), decimal_places)) +")" "\n" + 
                            ("V* Ave: " + format(B_mean_check, ' >6.3f') +
                            "  R* Ave: " + format(V_mean_check, ' >6.3f') +
                            "  (V-R)*: " + format(B_mean_check-V_mean_check, ' >6.3f')).rjust(72) +
                            '\n' +
                            ("V* Std: " + format(B_std_check, ' >6.3f') +
                            "  R* Std: " + format(V_std_check, ' >6.3f')).rjust(56))

            self.console_msg("Variable Star Estimates\n" +
                            ("V* Ave: " + format(B_mean_var, ' >6.3f') +
                            "  R* Ave: " + format(V_mean_var, ' >6.3f') +
                            "  (V-R)*: " + format(B_mean_var-V_mean_var, ' >6.3f')).rjust(72) +
                            '\n' +
                            ("V* Std: " + format(B_std_var, ' >6.3f') +  
                            "  R* Std: " + format(V_std_var, ' >6.3f')).rjust(56))
        elif input_color == 'V-I':
            self.console_msg("Check Star Estimates using check star: " + str(int(aux_result_first_color["KNAME"])) + 
                            " (V: " + str(round(float(aux_result_first_color['KREFMAG']), decimal_places)) +")" + 
                            " (I: " + str(round(float(aux_result_second_color['KREFMAG']), decimal_places)) +")" "\n" + 
                            ("V* Ave: " + format(B_mean_check, ' >6.3f') +
                            "  I* Ave: " + format(V_mean_check, ' >6.3f') +
                            "  (V-I)*: " + format(B_mean_check-V_mean_check, ' >6.3f')).rjust(72) +
                            '\n' +
                            ("V* Std: " + format(B_std_check, ' >6.3f') +
                            "  I* Std: " + format(V_std_check, ' >6.3f')).rjust(56))

            self.console_msg("Variable Star Estimates\n" +
                            ("V* Ave: " + format(B_mean_var, ' >6.3f') +
                            "  I* Ave: " + format(V_mean_var, ' >6.3f') +
                            "  (V-I)*: " + format(B_mean_var-V_mean_var, ' >6.3f')).rjust(72) +
                            '\n' +
                            ("V* Std: " + format(B_std_var, ' >6.3f') +  
                            "  I* Std: " + format(V_std_var, ' >6.3f')).rjust(56))
        else:
            raise Exception("generate_aavso_report_2color: unknown imput_color entered")

        try:
            with open(report_filename, mode='w') as f:
                f.write("#TYPE=Extended\n")
                f.write("#OBSCODE="+self.aavso_obscode_entry.get()+"\n")
                f.write("#SOFTWARE=Self-developed using photutils.psf; IRAFStarFinder\n") 
                f.write("#DELIM=,\n")
                f.write("#DATE=JD\n")
                f.write("#OBSTYPE=CCD\n")
                f.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")

                #B filter
                aux_result = aux_result_first_color
                
                starid = str(self.object_name_entry.get())
                date = format(float(aux_result['Date-Obs']), '.5f') 
                mag = str(round(B_mean_var, decimal_places))
                merr = str(round(B_std_check,decimal_places))
                filt = first_filter[input_color]
                trans = "YES"
                mtype = "STD"
                cname = "ENSEMBLE"
                cmag = "na"
                kname = str(int(aux_result_first_color["KNAME"]))
                kmag = str(round(B_mean_check, decimal_places))
                amass = format(float(aux_result_first_color["AMASS"]), '.3f') if aux_result_first_color["AMASS"].iloc[0].isnumeric() else "na"
                group = "na"
                chart = self.vizier_catalog_entry.get().strip()
                notes = self.object_notes_entry.get().strip()
                notes += "|KMAGINS=" + str(round(float(aux_result['KMAGINS']), decimal_places)) + \
                         "|KMAGSTD=" + str(round(B_mean_check, decimal_places)) + \
                         "|KREFMAG=" + str(round(float(aux_result['KREFMAG']), decimal_places)) + \
                         "|" + t_coefficient_tbv_for_report_only[input_color] + "="  + str(round(float(aux_result[t_coefficient_tbv[input_color]]), decimal_places)) + \
                         "|VMAGINS=" + str(round(float(aux_result['VMAGINS']), decimal_places))
                
                
                # Add " " after notes, because TA clobbers last char
                f.write(starid+","+date+","+mag+","+merr+","+filt+","+trans+","+mtype+"," +
                        cname+","+cmag+","+kname+","+kmag+","+amass+","+group+","+chart+","+notes+" \n")

                #V filter (easier just to repeat most of this)
                
                aux_result = aux_result_second_color
                
                starid = str(self.object_name_entry.get())
                date = format(float(aux_result['Date-Obs']), '.5f') 
                mag = str(round(V_mean_var, decimal_places))
                merr = str(round(V_std_check,decimal_places))
                filt = second_filter[input_color]
                trans = "YES"
                mtype = "STD"
                cname = "ENSEMBLE"
                cmag = "na"
                kname = self.object_kref_entry.get().strip()
                kmag = str(round(V_mean_check, decimal_places))
                amass = format(float(aux_result_second_color["AMASS"]), '.3f') if aux_result_first_color["AMASS"].iloc[0].isnumeric() else "na"
                group = "na"
                chart = self.vizier_catalog_entry.get().strip()
                notes = self.object_notes_entry.get().strip()
                notes += "|KMAGINS=" + str(round(float(aux_result['KMAGINS']), decimal_places)) + \
                         "|KMAGSTD=" + str(round(V_mean_check, decimal_places)) + \
                         "|KREFMAG=" + str(round(float(aux_result['KREFMAG']), decimal_places)) + \
                         "|" + t_coefficient_tv_bv[input_color] + "="  + str(round(float(aux_result[t_coefficient_tv_bv[input_color]]), decimal_places)) + \
                         "|VMAGINS=" + str(round(float(aux_result['VMAGINS']), decimal_places))
                
                # Add " " after notes, because TA clobbers last char
                f.write(starid+","+date+","+mag+","+merr+","+filt+","+trans+","+mtype+"," +
                        cname+","+cmag+","+kname+","+kmag+","+amass+","+group+","+chart+","+notes+" \n")

            self.console_msg(
                "AAVSO Photometry report saved to " + str(report_filename))

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)


    def exit_app(self):
        os._exit(0)

    def __init__(self):
        #Wis heisen Sie?
        self.program_name = "MAOPhot"
        self.program_version = "0.2"
        self.program_name_note = ""
        self.program_full_name = self.program_name + " " + self.program_version + " " + self.program_name_note

        #set the logger up
        self.our_logger = logging.getLogger(self.program_name + self.program_version + ".log")
        self.our_fh = logging.FileHandler(self.program_name + self.program_version + ".log")
        self.our_logger.setLevel(logging.INFO)

        # create formatter and add it to the handlers
        self.our_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        self.our_fh.setFormatter(self.our_formatter)
        self.our_logger.addHandler(self.our_fh)
        
        
        #set up GUI

        self.window = tk.Tk()
        self.screen_width = self.window.winfo_screenwidth()
        self.screen_height = self.window.winfo_screenheight()


        # Matplotlib settings
        matplotlib.rc('xtick', labelsize=7)
        matplotlib.rc('ytick', labelsize=7)

        # Maximize that works everywhere
        m = self.window.maxsize()
        self.window.geometry('{}x{}+0+0'.format(*m))

        self.window.title(self.program_full_name)

        self.menubar = tk.Menu(self.window)
        self.filemenu = tk.Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="Open...", command=self.open_FITS_file)
        self.filemenu.add_command(label="Save", command=self.save_FITS_file)
        self.filemenu.add_command(
            label="Save As...", command=self.save_FITS_file_as)
        self.filemenu.add_separator()
        self.filemenu.add_command(
            label="Edit Settings...", command=self.show_settings)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=self.exit_app)
        self.menubar.add_cascade(label="File", menu=self.filemenu)

        self.viewmenu = tk.Menu(self.menubar, tearoff=0)
        self.viewmenu.add_command(label="Zoom In", command=self.zoom_in)
        self.viewmenu.add_command(label="Zoom Out", command=self.zoom_out)
        self.viewmenu.add_command(label="100% Zoom", command=self.zoom_100)
        self.menubar.add_cascade(label="View", menu=self.viewmenu)

        self.photometrymenu = tk.Menu(self.menubar, tearoff=0)
        self.photometrymenu.add_command(label="Create Effective PSF", command=self.create_ePSF)
        self.photometrymenu.add_separator()
        self.photometrymenu.add_command(label="Load Rejection List...", command=self.load_ePSF_rejection_list)
        self.photometrymenu.add_command(label="Save Rejection List As...", command=self.save_as_ePSF_rejection_list)
        self.photometrymenu.add_command(label="Clear ePSF Data", command=self.clear_ePSF)
        self.photometrymenu.add_separator()

        self.photometrymenu.add_command(
            label="Non-iterative PSF Photometry", command=self.execute_noniterative_psf_photometry)
        self.photometrymenu.add_command(
            label="Iterative PSF Photometry", command=self.execute_iterative_psf_photometry)
        self.photometrymenu.add_separator()

        self.photometrymenu.add_command(
            label="Solve Image", command=self.solve_image)
        self.photometrymenu.add_separator()

        self.photometrymenu.add_command(
            label="Get Comparison Stars", command=self.get_comparison_stars)
        self.menubar.add_cascade(label="Photometry", menu=self.photometrymenu)
        
        self.two_color_photo_menu = tk.Menu(self.menubar, tearoff=0)
        self.two_color_photo_menu.add_command(
            label="(B,V)", command=self.BV_two_color_photometry)
        self.two_color_photo_menu.add_command(
            label="(V,R)", command=self.VR_two_color_photometry)
        self.two_color_photo_menu.add_command(
            label="(V,I)", command=self.VI_two_color_photometry)
        self.menubar.add_cascade(label="Two Color Photometry", menu=self.two_color_photo_menu)

        self.reportmenu = tk.Menu(self.menubar, tearoff=0)
        self.reportmenu.add_command(label="Single Image Photometry", command=self.generate_aavso_report_1image)
        self.two_color_sub_menu = tk.Menu(self.reportmenu, tearoff=False)
        self.two_color_sub_menu.add_command(label = "(B-V)", command=self.BV_generate_aavso_report_2color)
        self.two_color_sub_menu.add_command(label = "(V-R)", command=self.VR_generate_aavso_report_2color)
        self.two_color_sub_menu.add_command(label = "(V-I)", command=self.VI_generate_aavso_report_2color)
        self.reportmenu.add_cascade(label="Two Color Photometry", menu=self.two_color_sub_menu)

        self.menubar.add_cascade(label="Generate AAVSO Report", menu=self.reportmenu)

        self.window.config(menu=self.menubar)

        self.left_half = tk.Frame(self.window)  # Left half of the window
        self.left_half.grid(row=0, column=0, sticky=tk.NSEW)
        self.center = tk.Frame(self.window)  # Center of the window
        self.center.grid(row=0, column=1, sticky=tk.NSEW)
        self.right_half = tk.Frame(self.window)  # Right half of the window
        self.right_half.grid(row=0, column=2, sticky=tk.NSEW)
        # Expand center horizontally
        tk.Grid.columnconfigure(self.window, 1, weight=1)
        # Expand everything vertically
        tk.Grid.rowconfigure(self.window, 0, weight=1)

        self.filename_label = tk.Label(self.center, text="FITS:" + image_file)
        self.filename_label.grid(row=0, column=0)  # Place label

        self.canvas = tk.Canvas(self.center, bg='black')  # Main canvas
        # Place main canvas, sticky to occupy entire
        self.canvas.grid(row=1, column=0, sticky=tk.N+tk.S+tk.E+tk.W)
        # cell dimensions
        # Expand main canvas column to fit whole window
        tk.Grid.columnconfigure(self.center, 0, weight=1)
        # Expand main canvas row to fit whole window
        tk.Grid.rowconfigure(self.center, 1, weight=1)
        self.canvas_scrollbar_V = tk.Scrollbar(
            self.center, orient=tk.VERTICAL)  # Main canvas scrollbars
        self.canvas_scrollbar_V.grid(row=1, column=1)
        self.canvas_scrollbar_V.grid(
            sticky=tk.N+tk.S+tk.E+tk.W, column=1, row=1)
        self.canvas_scrollbar_H = tk.Scrollbar(
            self.center, orient=tk.HORIZONTAL)
        self.canvas_scrollbar_H.grid(row=2, column=0)
        self.canvas_scrollbar_H.grid(
            sticky=tk.N + tk.S + tk.E + tk.W, column=0, row=2)
        self.canvas_scrollbar_H.config(command=self.canvas.xview)
        self.canvas_scrollbar_V.config(command=self.canvas.yview)
        self.canvas.config(xscrollcommand=self.canvas_scrollbar_H.set)
        self.canvas.config(yscrollcommand=self.canvas_scrollbar_V.set)

        # Console below
        self.console = tk.Text(self.center, #height=40,
                               bg='black', fg='white', width=200)
        self.console.grid(sticky=tk.N+tk.S+tk.E+tk.W, column=0, row=3)
        self.console_scrollbar = tk.Scrollbar(self.center)
        self.console_scrollbar.grid(
            sticky=tk.N + tk.S + tk.E + tk.W, column=1, row=3)

        self.console.config(yscrollcommand=self.console_scrollbar.set)
        self.console_scrollbar.config(command=self.console.yview)

        self.console_msg(self.program_full_name)
        self.console_msg("Ready")




        # We will lay out interface things into the new right_frame grid
        self.right_frame = tk.Frame(self.right_half)
        # Place right_frame into the top of the main canvas row, right next to it
        self.right_frame.grid(row=1, column=2, sticky=tk.N)
        
        self.plotname_label = tk.Label(self.right_frame, text="Plot:")
        self.plotname_label.grid(row=0, column=0)  # Place label
        self.fig_psf = Figure()
        self.psf_plot = self.fig_psf.add_subplot(111, projection='3d')
        # PSF 3D plot canvas - Matplotlib wrapper for Tk
        self.psf_plot_canvas = FigureCanvasTkAgg(self.fig_psf, self.right_frame)
        self.psf_plot_canvas.draw()
        self.psf_canvas = self.psf_plot_canvas.get_tk_widget()
        self.psf_canvas.config(width=int(self.screen_width/8.5), height=int(self.screen_width/8.5))
        # Allocate small PSF canvas to a new grid inside the right_frame
        self.psf_canvas.grid(row=1, column=0)   #was row0
        
        #make another canvas for 2D plot of effectivePSF
        self.ePSF_plotname_label = tk.Label(self.right_frame, text="Effective PSF:")
        self.ePSF_plotname_label.grid(row=2, column=0)  # Place label


        self.fig_ePSF, self.ePSF_plot = plt.subplots()
        #self.ePSF_plot = self.fig_ePSF.add_subplot(111)
        #
        self.ePSF_plot_canvas = FigureCanvasTkAgg(self.fig_ePSF, self.right_frame)
        self.ePSF_plot_canvas.draw()
        self.ePSF_canvas = self.ePSF_plot_canvas.get_tk_widget()
        self.ePSF_canvas.config(width=int(self.screen_width/8.5), height=int(self.screen_width/8.5))
        # Allocate small PSF canvas to a new grid inside the right_frame
        self.ePSF_canvas.grid(row=3, column=0)   #was row0

        #
        # Left Frame
        #

        # We will lay out interface things into the new left_frame grid
        self.left_frame = tk.Frame(self.left_half)
        # Place right_frame into the top of the main canvas row, right next to it
        self.left_frame.grid(row=1, column=0, sticky=tk.NSEW)
        

        # Frame to hold "SOME" settings in grid pattern
        self.settings_frame = tk.Frame(self.left_frame)
        self.settings_frame.pack()

        # Expand settings_frame column that holds labels
        tk.Grid.columnconfigure(self.settings_frame, 0, weight=1)

        settings_entry_width = 6
        settings_entry_pad = 0
        extended_settings_entry_width = 30
        extended_settings_entry_pad = 0

        #
        # Some parameter settings need initial values
        #
        self.catalog_stringvar = tk.StringVar()
        self.catalog_stringvar.set("AAVSO")
        self.fitter_stringvar = tk.StringVar()
        self.fitter_stringvar.set("Levenberg-Marquardt")

        row = 0

        self.stretching_label = tk.Label(
            self.settings_frame, text="Image Stretching:")
        self.stretching_label.grid(row=row, column=0, sticky=tk.E)
        self.stretching_stringvar = tk.StringVar()
        self.stretching_stringvar.set("Asinh")
        self.stretching_dropdown = tk.OptionMenu(
            self.settings_frame, self.stretching_stringvar, "None", "Square Root", "Log", "Asinh")
        self.stretching_dropdown.grid(row=row, column=1, sticky=tk.EW)
        row += 1
        self.stretching_stringvar.trace(
            "w", lambda name, index, mode, sv=self.stretching_stringvar: self.update_display())

        # Histogram stretch sliders
        tk.Grid.columnconfigure(self.settings_frame, 1, weight=0)
        self.stretch_label = tk.Label(
            self.settings_frame, text="Histogram Stretch Low/High:")
        self.stretch_label.grid(row=row, column=0, sticky=tk.W)
        row += 1
        self.stretch_low = tk.Scale(
            self.settings_frame, from_=0, to=100, orient=tk.HORIZONTAL, command=self.update_histogram_low)
        self.stretch_low.grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        row += 1
        self.stretch_high = tk.Scale(
            self.settings_frame, from_=0, to=100, orient=tk.HORIZONTAL, command=self.update_histogram_high)
        self.stretch_high.set(5)
        self.stretch_high.grid(row=row, column=0, columnspan=2, sticky=tk.NSEW)

        #
        # laumch_settings; pops up settings window and initializes all settings
        #
        self.launch_settings()

        """
        valid_parameter_list facilitates loading from and saving to a settings file
        buttons in popup window: Settings
        """
        self.valid_parameter_list = {
            'photometry_aperture_entry': self.photometry_aperture_entry,
            'max_ensemble_magnitude_entry': self.max_ensemble_magnitude_entry,
            'fwhm_entry': self.fwhm_entry,
            'star_detection_threshold_entry': self.star_detection_threshold_entry,
            'photometry_iterations_entry': self.photometry_iterations_entry,
            'sharplo_entry': self.sharplo_entry,
            'bkg_filter_size_entry': self.bkg_filter_size_entry,
            'matching_radius_entry': self.matching_radius_entry,
            'aavso_obscode_entry': self.aavso_obscode_entry,
            'telescope_entry': self.telescope_entry,
            'tbv_entry': self.tbv_entry,
            'tv_bv_entry': self.tv_bv_entry,
            'tb_bv_entry': self.tb_bv_entry,
            'tvr_entry': self.tvr_entry,
            'tv_vr_entry': self.tv_vr_entry,
            'tr_vr_entry': self.tr_vr_entry,
            'tvi_entry': self.tvi_entry,
            'tv_vi_entry': self.tv_vi_entry,
            'ti_vi_entry': self.ti_vi_entry,
            'catalog_stringvar': self.catalog_stringvar,
            'vizier_catalog_entry': self.vizier_catalog_entry,
            'fitter_stringvar': self.fitter_stringvar,
            'astrometrynet_entry': self.astrometrynet_entry,
            'astrometrynet_key_entry': self.astrometrynet_key_entry,
            'object_kref_entry': self.object_kref_entry,
            'object_sel_comp_entry': self.object_sel_comp_entry,
            'object_name_entry': self.object_name_entry,
            'object_notes_entry': self.object_notes_entry
            }

        tk.mainloop()


myGUI = MyGUI()

