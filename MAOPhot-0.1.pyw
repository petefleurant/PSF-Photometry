# -*- coding: utf-8 -*-
"""
    MAOPhot 1.0 
    
    MAOPhot calculates stellar magnitudes from 2 dimensional digital astrophotographs.
    It has the option of producing an extended AAVSO (American Association of Variable Star Observers)
    report which can then be submitted to the AAVSO using an online tool WebObs (http://www.aavso.org/webobs).

    There are many photometry measureing programs available such as VPhot (http://www.aavso.org/vphot) 
    and AstroImageJ (University of Louisville).

    MAOPhot uses PSF (point spread function) modeling which is well suited for measuring stellar magnitudes
    in crowded fields, or a particular star that has a close companion, e.g., Z Tau.

    See "PSF Photometry" (https://photutils.readthedocs.io/en/stable/psf.html) whcih describes many of the classes and methods used in MAOPhot
    
    Original MetroPSF (C) Copyright 2021, Maxym Usatov 
    <maxim.usatov@bcsatellite.net> Refer to metropsf.pdf for license information.
 


    General Workflow for Two Color Phtometry and AAVSO report generation:
        - load Settings
        - Open B color fits file
        - Solve Image if not done already
        - Photometry->Interatively Subtracted PSF Phtometry
        - Photometry->Get Comparison Stars
        - Repeat last 4 steps for V color
        - Two Color Photometry->Two Color Photometry (B,V)
        - Report->AAVSO: Generate Report; Two Color Photometry
        


    Modifications by petefleurant@gmail.com, as follows

    2022 6 16

    Slight modificationsa added
    - Add comp star label to .phot file
    - Use log10 in place of plain log (which is natural logarithm) in
        plot_photometry
    - Change extention ".phot" to ".csv"
    - replace chart with AAVSO Chart ID
    - Added pink fields "Notes" and "Check Star"
    - Fixed bug with resizing image before AAVSO report generates
    - Added title "Plot" (plotname_label) above 3d canvas with vsx_id 
        if available
    - Clear canvas if no star selected
    - Add choice for comparison catalog: Add AAVSO; this returns 
        comparison stars goven frame dimnensions. 
    - Make sure comparison stars are unique, uniquify by using largest flux_0
    - plots get cleared when opening a FITS file

    Added
        python logging into console_msg()


        generate_aavso_report_2color()


        2022 9 13 
        
        MAJOR CHANGE
        
        Modifications made to basically have MAOPhot mimic VPhot's "Two
        Color Phtometry" (Usually B and V)
            
        See spreadsheet: 
            E:\Astronomy\Processing MAO images\ProcessingMaoImages_202281V1117Her.xlsx

        This includes formulas to generate "two color photometry". 
        generate_aavso_report_2color() uses the same formula which includes
        using transformation coefficients.

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
            
            
        The following formulas are calculated for each comp star and var = check star;
        then calculated again for each comp star and var = "Object Name". Eg. V1117 Her
        
        Δv = vvar - vcomp
        Δ(B-V) = Tbv * Δ(b-v)
        Vvar = Δv + Tv_bv * Δ(B-V) + Vcomp
        
        where:
        Δv  = IM of variable  - IM comp
        Vcomp = published V-magnitude of comp
        Δ(B-V) = Tbv * difference between standard color of var and standard color of comp
        
        
        In order to calculate these formulas and then generate a report, two (2) sets
        of results_tab_df in csv format must exist, one for B and one for V and
        must have been derived from the B and V images of the Var under analysis.
        
        When generate_aavso_report_2color is called, by the menu item: 
            Report->AAVSO: Generate Report Two Color Photometry
        then the user will be asked to specify the 2 aformentioned csv files.
        
        From these files/Panda databases, the formulas are calculated, and results 
        are displayed. 
        
        
        
        --
        get_comparison_stars
        
        Given the specified Check Star in field "Use Check Star (AAVSO label)"
        and the list of comp stars in "Select Comp Stars (AAVSO Label)", these
        stars are gotten from the AAVSO and each has an additional attribute,
        is_check_star, added to the results_tab_df DB. 
        

    --------------------------------------------------------------------------
        
        Removed use of flux_unc which is not always returned; use poisson_noise_error = 1.0857 / snr

        Select Comp Stars entry list; comp stars only used that are in this list




    Things needed
        - Add label number to red circle in plot; so that one knows which are
        outliers
      
      
 
 
 
 
 
"""
from ast import Assert
from photutils import Background2D, MedianBackground
from astropy.stats import SigmaClip
from astropy.modeling.fitting import LevMarLSQFitter, LinearLSQFitter, SLSQPLSQFitter, SimplexLSQFitter
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils.background import MMMBackground, MADStdBackgroundRMS
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.detection import IRAFStarFinder
from tkinter import filedialog as fd
import tkinter as tk
from astropy.nddata import Cutout2D
from astropy.stats import gaussian_sigma_to_fwhm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import matplotlib
from astroquery.vizier import Vizier
from astroquery.astrometry_net import AstrometryNet
from PIL import Image, ImageTk, ImageMath
from astropy.visualization import SqrtStretch, LogStretch, AsinhStretch
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.io import fits
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
    background_image = background_image.resize(new_size, Image.ANTIALIAS)
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
    _image = _image.resize(new_size, Image.ANTIALIAS)
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
    generated_image = generated_image.resize(new_size, Image.ANTIALIAS)
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
    results_tab_df = pd.DataFrame()
    bkg_value = 0
    fit_shape = 21
    error_raised = False
    histogram_slider_low = 0
    histogram_slider_high = 5
    last_clicked_x = 0
    last_clicked_y = 0
    last_clicked_differential_magnitude = 0
    last_clicked_differential_uncertainty = 0
    last_clicked_true_inst_mag = 0
    ensemble_size = 0
    a = 0
    b = 0
    jd = 0
    image_file = ""
    photometry_circles = {}

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
            if self.photometry_results_plotted:
                self.plot_photometry()


    def display_background(self):
        self.canvas.delete("all")
        self.console_msg("Displaying extracted background.")
        global background_image
        self.image = ImageTk.PhotoImage(background_image)
        self.canvas.create_image(0, 0, anchor=tk.NW, image=self.image)
        self.canvas.config(scrollregion=self.canvas.bbox(tk.ALL))

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
                self.results_tab_df = pd.DataFrame()
                header = image[0].header
                image_data = fits.getdata(image_file)
                image_width = image_data.shape[1]
                image_height = image_data.shape[0]
                self.wcs_header = WCS(image[0].header)

                if self.crop_fits_entry.get() != "100":
                    factor = 100 / int(self.crop_fits_entry.get())
                    new_width = int(image_width / factor)
                    new_height = int(image_height / factor)
                    x0 = int((image_width - new_width) / 2)
                    y0 = int((image_height - new_height) / 2)
                    x1 = x0 + new_width
                    y1 = y0 + new_width
                    image_data = image_data[y0:y1, x0:x1]
                    image_width = new_width
                    image_height = new_height

#               self.console_msg("WCS Header: " + str(self.wcs_header))

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
                if 'exptime' in header:
                    exptime = header['exptime']
                    self.set_entry_text(self.exposure_entry, exptime)
                    self.console_msg("Exposure: " + str(exptime))
                else:
                    self.console_msg(
                        "Exposure (EXPTIME) not in FITS header. Set exposure manually.")

                if 'gain' in header:
                    gain = header['gain']
                    self.set_entry_text(self.ccd_gain, gain)
                    self.console_msg("Gain (e-/ADU): " + str(gain))
                else:
                    self.console_msg(
                        "Gain not in FITS header. Set gain manually for aperture photometry.")

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

                self.bkg_value = np.median(image_data)
                self.console_msg(
                    "Median background level, ADU: " + str(self.bkg_value))
                
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
            pass

    def initialize_debug(self):
        self.set_entry_text(self.crop_fits_entry, "30")
        image_file = "calibrated-T11-blackhaz-CG Dra-20210516-010321-V-BIN1-E-300-012.fit"
        if len(image_file) > 0:
            self.load_FITS(image_file)
            self.display_image()

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
                self.clear_psf_linreg_plot()

                if self.plate_solve_on_open.get():
                    self.console_msg("Solving via Astrometry.Net...")
                    
                    ast = AstrometryNet()
                    ast.api_key = self.astrometrynet_key_entry.get()
                    ast.URL = "http://" + self.astrometrynet_entry.get()
                    ast.API_URL = "http://" + self.astrometrynet_entry.get() + "/api"
        
                    self.wcs_header = ast.solve_from_image(image_file, force_image_upload=True)
                    
                    self.console_msg(
                        "Astrometry.Net solution reference point RA: " + str(self.wcs_header["CRVAL1"]) + " Dec: " + str(
                            self.wcs_header["CRVAL2"]))
                    
                    header = header + self.wcs_header
                    self.wcs_header = WCS(header)
                
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

    def aperture_photometry(self):
        from photutils.aperture import CircularAperture
        from photutils.aperture import aperture_photometry
        from photutils.utils import calc_total_error
        global header
        self.console_msg("Initiating aperture photometry...")
        try:
            self.fit_shape = int(self.photometry_aperture_entry.get())
            fwhm = int(self.fwhm_entry.get())
            star_detection_threshold = float(
                self.star_detection_threshold_entry.get())
            iterations = int(self.photometry_iterations_entry.get())
            bkgrms = MADStdBackgroundRMS()
            sharplo = float(self.sharplo_entry.get())
            bkg_filter_size = int(self.bkg_filter_size_entry.get())
            std = bkgrms(image_data)
            sigma_clip = SigmaClip(sigma=3.0)
            bkg_estimator = MedianBackground()
            self.console_msg("Estimating background...")
            bkg = Background2D(image_data, (self.fit_shape * 1, self.fit_shape * 1),
                               filter_size=(
                                   bkg_filter_size, bkg_filter_size), sigma_clip=sigma_clip,
                               bkg_estimator=bkg_estimator)
            clean_image = image_data-bkg.background
            save_background_image(
                self.histogram_slider_low, self.histogram_slider_high, self.zoom_level, bkg.background)
            self.console_msg("Estimating total error...")
            #effective_gain = int(float(self.exposure_entry.get()))
            effective_gain = float(self.ccd_gain_entry.get())
            error = calc_total_error(
                image_data, bkg.background_rms, effective_gain)

            self.console_msg("Finding stars...")
            iraffind = IRAFStarFinder(threshold=star_detection_threshold * std,
                                      fwhm=fwhm, roundhi=3.0, roundlo=-5.0,
                                      sharplo=sharplo, sharphi=2.0)
            positions_found = iraffind(clean_image)
            positions = np.array(positions_found[['xcentroid', 'ycentroid']])
            # Convert a list of tuples to 2D array
            positions = [list(item) for item in positions]
            aperture = CircularAperture(positions, r=self.fit_shape)
            self.console_msg("Performing photometry...")
            phot_table = aperture_photometry(
                clean_image, aperture, error=error)
            self.results_tab_df = phot_table.to_pandas()
            # Let's rename and add columns so the resulting dataframe will be consistent with that produced by the PSF photometry
            self.results_tab_df.rename(columns={"xcenter": "x_0",
                                                "ycenter": "y_0",
                                                "aperture_sum": "flux_fit",
                                                "aperture_sum_err": "flux_unc"}, inplace=True)
            self.results_tab_df["x_fit"] = self.results_tab_df["x_0"]
            self.results_tab_df["y_fit"] = self.results_tab_df["y_0"]
            self.results_tab_df["flux_0"] = self.results_tab_df["flux_fit"]
            #self.results_tab_df["sigma_fit"] = 0

            self.results_tab_df["removed_from_ensemble"] = False
            self.results_tab_df["date-obs"] = float(self.date_obs_entry.get())
            self.results_tab_df.to_csv(self.image_file + ".csv", index=False)
            self.console_msg("Photometry saved to " + str(self.image_file + ".csv"))
            self.plot_photometry()
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  +" "+str(e), level=logging.ERROR)

    def perform_psf_photometry(self):
        global header
        self.console_msg(
            "Initiating interatively subtracted PSF photometry...")
        try:
            bkgrms = MADStdBackgroundRMS()
            std = bkgrms(image_data)
            self.fit_shape = int(self.photometry_aperture_entry.get())
            fwhm = float(self.fwhm_entry.get())
            star_detection_threshold = float(
                self.star_detection_threshold_entry.get())
            iterations = int(self.photometry_iterations_entry.get())
            sharplo = float(self.sharplo_entry.get())
            bkg_filter_size = int(self.bkg_filter_size_entry.get())
            self.console_msg("Finding stars...")
            iraffind = IRAFStarFinder(threshold = star_detection_threshold * std,
                                      fwhm = fwhm,
                                      roundhi = 3.0,
                                      roundlo = -5.0,
                                      sharplo = sharplo,
                                      sharphi = 2.0)

         
            daogroup = DAOGroup(2.0 * fwhm)
            sigma_clip = SigmaClip(sigma=3.0)
            mmm_bkg = MMMBackground(sigma_clip=sigma_clip)


            """
            
            
              Let IterativelySubtractedPSFPhotometry do the subtraction
            
            sigma_clip = SigmaClip(sigma=3.0)
            bkg_estimator = MedianBackground()
            self.console_msg("Estimating background...")
            bkg = Background2D(image_data, (self.fit_shape * 1, self.fit_shape * 1),
                               filter_size=(
                                   bkg_filter_size, bkg_filter_size), sigma_clip=sigma_clip,
                               bkg_estimator=bkg_estimator)
            clean_image = image_data-bkg.background
            save_background_image(
                self.histogram_slider_low, self.histogram_slider_high, self.zoom_level, bkg.background)
            
            """

            if self.fitter_stringvar.get() == "Linear Least Square":
                self.console_msg("Setting fitter to Linear Least Square")
                selected_fitter = LinearLSQFitter(calc_uncertainties=True)

            elif self.fitter_stringvar.get() == "Sequential LS Programming":
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
                
            # sigma=2 here is the initial guess
            psf_model = IntegratedGaussianPRF(fwhm / gaussian_sigma_to_fwhm)
            #Dont do the next line, see https://photutils.readthedocs.io/en/stable/psf.html)
            #This results in mopre complicated outcomes
            #psf_model.sigma.fixed = False   # This allows to fit Gaussian PRF sigma as well 
            photometry = IterativelySubtractedPSFPhotometry(finder=iraffind,
                                                            group_maker = daogroup,
                                                            psf_model = psf_model,
                                                            bkg_estimator = mmm_bkg,
                                                            fitter = selected_fitter,
                                                            aperture_radius=2*fwhm,
                                                            niters = iterations,
                                                            fitshape = (self.fit_shape, self.fit_shape))
            self.console_msg("Performing photometry...")
            result_tab = photometry(image = image_data)
            residual_image = photometry.get_residual_image()
            fits.writeto("residuals.fits", residual_image, header, overwrite=True)
            self.console_msg("Done. PSF fitter message(s): " + str(selected_fitter.fit_info['message']))

            self.results_tab_df = result_tab.to_pandas()
            self.results_tab_df["removed_from_ensemble"] = False
            self.results_tab_df["date-obs"] = float(self.date_obs_entry.get())
            
            self.results_tab_df.to_csv(self.image_file + ".csv", index=False)
            self.console_msg("Photometry saved to " + str(self.image_file + ".csv"))
            self.plot_photometry()
            self.results_tab_df.to_csv(self.image_file + ".csv", index=False)
            self.console_msg("Photometry saved to " + str(self.image_file + ".csv"))

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  +" "+str(e), level=logging.ERROR)




    # center coordinates, radius
    def create_circle(self, x, y, r, canvas_name, outline_color="grey50"):
        x0 = x - r
        y0 = y - r
        x1 = x + r
        y1 = y + r
        self.photometry_circles[str(x)+str(y)] = canvas_name.create_oval(x0, y0, x1, y1, outline=outline_color)

    def create_text(self, x, y, r, canvas_name, anchor, text=''):
        x1 = x 
        y1 = y - 1.25*r
        canvas_name.create_text(x1, y1, fill='white', anchor=anchor, text=text)

    def plot_photometry(self):
        try:
            self.console_msg("Plotting Photometry...")
            labels_in_photometry_table = "label" in self.results_tab_df
            vsx_ids_in_photometry_table = False
            if "vsx_id" in self.results_tab_df:
                vsx_ids_in_photometry_table = True
            exptime = float(self.exposure_entry.get())
            if os.path.isfile(self.image_file+".csv"):
                self.fit_shape = int(self.photometry_aperture_entry.get())
                self.results_tab_df = pd.read_csv(self.image_file + ".csv")
                if "removed_from_ensemble" not in self.results_tab_df:
                    # This prefilling is required for backwards compatibility to read .phot
		            #(now called .csv) files from older versions.
                    self.results_tab_df["removed_from_ensemble"] = False
                # Calculate instrumental magnitudes
		        # Following added for "True" inst mag used in AAVSO report
                self.results_tab_df["true_inst_mag"] = -2.5 * np.log10(self.results_tab_df["flux_fit"] / exptime)
                self.results_tab_df["inst_mag"] = self.results_tab_df["true_inst_mag"] + 25
                
                #self.results_tab_df["inst_mag_min"] = -2.5 * np.log10((self.results_tab_df["flux_fit"] - self.results_tab_df["flux_unc"]) / exptime) + 25
                #self.results_tab_df["inst_mag_max"] = -2.5 * np.log10((self.results_tab_df["flux_fit"] + self.results_tab_df["flux_unc"]) / exptime) + 25

                self.photometry_results_plotted = True

                for index, row in self.results_tab_df.iterrows():
                    outline_color = "grey50"
                    if labels_in_photometry_table:
                        if len(str(row["label"])) > 0 and str(row["label"]) != "nan":
                            outline_color = "green"
                            self.create_circle(x=row["x_fit"] * self.zoom_level, y=row["y_fit"] * self.zoom_level,
                                r=self.fit_shape / 2 * self.zoom_level, canvas_name=self.canvas, outline_color=outline_color)
                            self.create_text(  x=row["x_fit"] * self.zoom_level, y=row["y_fit"] * self.zoom_level, 
                                r=self.fit_shape / 2 * self.zoom_level, canvas_name=self.canvas, anchor=tk.CENTER, text=str(int(row["label"])))
                            continue
                    if row["removed_from_ensemble"]:
                        assert False, "Found an entry 'removed from ensemble???!'"
                        outline_color = "red"
                    if vsx_ids_in_photometry_table:
                        if len(str(row["vsx_id"])) > 0 and str(row["vsx_id"]) != "nan":
                            outline_color = "yellow"
                            self.create_text(  x=row["x_fit"] * self.zoom_level, y=row["y_fit"] * self.zoom_level, 
                                r=self.fit_shape / 2 * self.zoom_level, canvas_name=self.canvas, anchor=tk.CENTER, text=str(row["vsx_id"]).strip())
                            self.create_circle(x=row["x_fit"] * self.zoom_level, y=row["y_fit"] * self.zoom_level,
                                       r=self.fit_shape / 2 * self.zoom_level, canvas_name=self.canvas, outline_color=outline_color)
            self.console_msg("Plotting Photometry...complete")
            self.console_msg("Ready")
            
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) + " " + str(e), level=logging.ERROR)

    def hide_photometry(self):
        self.photometry_results_plotted = False
        self.display_image()
    
    
    """
    def plot_sigma_heatmap(self):
        global image_width, image_height
        try:
            heatmap = np.empty((image_width, image_height))
            fig, ax = plt.subplots()
            if os.path.isfile(self.image_file+".csv"):
                self.results_tab_df = pd.read_csv(self.image_file + ".csv")
                for index, row in self.results_tab_df.iterrows():
                    heatmap[int(row["x_fit"]), int(
                        row["y_fit"])] = row["sigma_fit"]
                im = ax.imshow(heatmap)
                plt.show()

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
    """
    
    
    
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
                   #matched_objects.iloc[0]["sigma_fit"],
                   matched_objects.iloc[0]["true_inst_mag"],
                   matched_objects.iloc[0]["inst_mag"] #,
                   #matched_objects.iloc[0]["flux_unc"],
                   #matched_objects.iloc[0]["inst_mag_min"],
                   #matched_objects.iloc[0]["inst_mag_max"],
                   )
        else:
            return(0, 0, 0, 0, 0, 0) #, 0, 0, 0)

    def mouse_click(self, event):
        global image_data
        self.last_clicked_differential_magnitude = 0
        self.last_clicked_differential_uncertainty = 0
        self.last_clicked_true_inst_mag = 0
        decimal_places = int(self.decimal_places_entry.get())
        vsx_ids_in_photometry_table = False
        if "vsx_id" in self.results_tab_df:
            vsx_ids_in_photometry_table = True
        self.display_image()
        self.console_msg("")
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
        if self.photometry_results_plotted:
            # OLD x_fit, y_fit, flux_fit, sigma_fit, true_inst_mag, inst_mag, flux_unc, inst_mag_min, inst_mag_max = self.match_photometry_table(x, y)
            x_fit, y_fit, flux_fit, true_inst_mag, inst_mag = self.match_photometry_table(x, y)
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
                round(flux_fit, 2)) + " Instrumental magnitude: " + str(round(true_inst_mag, 3)) + " " + sky_coordinate_string)
                # OLD: round(flux_fit, 2)) + " Instrumental magnitude: " + str(round(true_inst_mag, 3)) + " PSF σ: " + str(round(sigma_fit, 2)) + sky_coordinate_string)
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

                        """    
                        if len(str(matching_star["nearby_vsx_id"])) > 0 and str(
                                matching_star["nearby_vsx_id"]) != "nan" and self.nearby_vsx_var.get():
                            self.console_msg(
                                "Nearby VSX Source: " + str(matching_star["nearby_vsx_id"]) + " separation: " + str(
                                    matching_star["nearby_vsx_separation"]))
                        """
                        
            # Has regression fit been done?
            if self.a != 0 and self.b != 0 and x_fit != 0 and y_fit != 0: 
                #regression fit has been done
                differential_magnitude = round(
                    inst_mag * self.a + self.b, decimal_places)

                total_error = np.nan
                
                if "mag_error" in self.results_tab_df:
                    """
                    Always use poisson noise error because flux_unc is not always available
                    """
                    # Aperture photometry mode - use Poisson noise error
                    snr = flux_fit/self.bkg_value
                    poisson_noise_error = 1.0857 / snr
                    total_error = round(np.sqrt(self.linreg_error ** 2 + poisson_noise_error ** 2), decimal_places)
                    self.console_msg("Linear Regression Error: " + str(round(self.linreg_error, decimal_places)) +
                                     " Poisson Noise Error: " + str(round(poisson_noise_error, decimal_places)))

                obj_name_string = ""
                if type(matching_star["match_id"]) in (str, int, np.float64) and str(matching_star["match_id"]) != "nan":
                    obj_name_string = matching_star["match_id"] + ", "
                if 'label' in self.results_tab_df.columns:
                    if type(matching_star["label"]) in (str, int, np.float64) and str(matching_star["label"]) != "nan":
                        obj_name_string += "label: " + str(round(matching_star["label"])) + ", "
                if vsx_ids_in_photometry_table and len(str(matching_star["vsx_id"])) > 0 and str(matching_star["vsx_id"]) != "nan":
                    obj_name_string = str(matching_star["vsx_id"]) + ", "
                self.console_msg("Differential ensemble photometry magnitude:")
                jd_time = Time(self.jd, format='jd')
                self.console_msg(obj_name_string + str(jd_time.to_value('iso')) + " UTC, " + self.filter_entry.get() + " = " + str(round(differential_magnitude, decimal_places)) + " ± " + str(
                    total_error))

                self.last_clicked_differential_magnitude = round(
                    differential_magnitude, decimal_places)
                self.last_clicked_differential_uncertainty = total_error
                self.last_clicked_true_inst_mag = round(
                    true_inst_mag, decimal_places)

                # Checking PSF sigma
                """
                median_psf_sigma = self.results_tab_df["sigma_fit"].median()
                if sigma_fit > median_psf_sigma * 1.2:
                    self.console_msg(
                        #
                        WARNING: This source has PSF σ exceeding median value by more than 20%.
                        Check PSF shape - likely non-linear sensor regime or improperly
                        subtracted sky background.
                        )
                """
            
        #if self.photometry_results_plotted:
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

    def clear_psf_linreg_plot(self):
        #clear plot label
        self.plotname_label['text'] = "Plot: "
        self.psf_plot.clear()
        self.psf_plot_canvas.draw()
        self.linreg_plot.clear()
        self.plot_canvas.draw()


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
        self.console_msg("Solving via Astrometry.Net...")
        try:
            ast = AstrometryNet()
            ast.api_key = self.astrometrynet_key_entry.get()
            ast.URL = "http://" + self.astrometrynet_entry.get()
            ast.API_URL = "http://" + self.astrometrynet_entry.get() + "/api"

            sources_df = self.results_tab_df.sort_values(
                "flux_fit", ascending=False)
            width, height = generated_image.size

            self.wcs_header = ast.solve_from_source_list(sources_df['x_fit'], sources_df['y_fit'],
                                                         width, height,
                                                         solve_timeout=360)
            self.console_msg(
                "Astrometry.Net solution reference point RA: " + str(self.wcs_header["CRVAL1"]) + " Dec: " + str(
                    self.wcs_header["CRVAL2"]))
            header = header + self.wcs_header
            self.wcs_header = WCS(header)
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
            
            
    def two_color_photometry_BV(self):
        try:
            """
            Two Color Photometry requires 'Object Name' to be filled; eg. 'V1117 Her'
            
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
            
            star                     "check" or "var"
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
            variable_star = self.object_name_entry.get().strip()
            if len(variable_star) == 0:
                self.console_msg(
                    "Two Color Photometry requires 'Object Name' to be filled; eg. 'V1117 Her'")
                return
            
            # Ask for the B and V CSV files
            options = {}
            options['defaultextension'] = '.csv'
            options['filetypes'] = [('CSV', '.csv')]
            #options['initialdir'] = "E:\Astronomy\AstronomyData\MAO\2022\2022 6 4 V1117 Her\PIData\master\MAOPhot"
            #options['initialfile'] = ''
            options['title'] = 'Choose a file for color B'

            file_name = fd.askopenfilename(**options)

            if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
                self.console_msg("Loading color B from " + str(file_name))
                self.results_tab_df_colorB = pd.read_csv(str(file_name))
            else:
                return
                
            options['title'] = 'Choose a file for color V'
            file_name = fd.askopenfilename(**options)
            if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
                self.console_msg("Loading color V from " + str(file_name))
                self.results_tab_df_colorV = pd.read_csv(str(file_name))
            else:
                return
            self.console_msg("Performing Two Color Photometry...")
            
            #get the transformation coefficients
            tbv_coefficient = float(self.tbv_entry.get())
            tb_bv_coefficient = float(self.tb_bv_entry.get())
            tv_bv_coefficient = float(self.tv_bv_entry.get())

         
            #
            #   CHECK STAR Calculations
            #
            # Find the check star; 
            check_star_B = self.results_tab_df_colorB[self.results_tab_df_colorB["check_star"] == True].iloc[0]
            
            self.console_msg("Using check star " + str(int(check_star_B["label"])))

            check_IMB = check_star_B["true_inst_mag"]
            check_B = check_star_B["match_mag"]

            check_star_V = self.results_tab_df_colorV[self.results_tab_df_colorV["check_star"] == True].iloc[0]
            check_IMV = check_star_V["true_inst_mag"]
            check_V = check_star_V["match_mag"]
            
            # Find the variable_star; 
            var_star_B = self.results_tab_df_colorB[self.results_tab_df_colorB["vsx_id"] == variable_star].iloc[0]
            var_IMB = var_star_B["true_inst_mag"]

            var_star_V = self.results_tab_df_colorV[self.results_tab_df_colorV["vsx_id"] == variable_star].iloc[0]
            var_IMV = var_star_V["true_inst_mag"]


            # The proper dae-obs should be set in post-processing, where the number of subs and exposure time is 
            #taken into consideration
            date_obs_B = self.results_tab_df_colorB[self.results_tab_df_colorB["check_star"] == True].iloc[0]["date-obs"]
            date_obs_V = self.results_tab_df_colorV[self.results_tab_df_colorV["check_star"] == True].iloc[0]["date-obs"]
            
            
            """
            Build the result_check_star and result_var_table which are similar to the aforementioed 
            spreadsheet ProcessingMaoImages_202281V1117Her.xlsx
            """
            
            result_check_star = pd.DataFrame(columns=["star", "label", "IMB", "IMV", "B", "V", "delta_b_minus_v", "delta_B_minus_V",
                                                      "delta_b", "delta_v", "comp_b_minus_v", "B_star", "V_star", "outlier"])
            
            result_var_star = pd.DataFrame(columns=["star", "label", "IMB", "IMV", "B", "V", "delta_b_minus_v", "delta_B_minus_V",
                                                    "delta_b", "delta_v", "comp_b_minus_v", "B_star", "V_star"])
            

            """
            loop through all the selected comp stars and fill the result_check_star table
            """
            sel_comps = [] #init
            sel_comps_to_use = self.object_sel_comp_entry.get().strip().split(',')

            #make array of int                
            for comp in sel_comps_to_use:
                sel_comps.append(int(comp.strip()))
            
            for comp in sel_comps:
                #selected comp must be in both tables
                if comp not in self.results_tab_df_colorB["label"].values:
                    self.console_msg("Comp star: "+ str(int(comp)) + " not in B table")
                    continue

                if comp not in self.results_tab_df_colorV["label"].values:
                    self.console_msg("Comp star: "+ str(int(comp)) + " not in V table")
                    continue
                
                #Dont use the check star 
                if comp == int(check_star_B["label"]):
                    continue
                
                comp_star_B = self.results_tab_df_colorB[self.results_tab_df_colorB["label"] == comp].iloc[0]
                comp_star_V = self.results_tab_df_colorV[self.results_tab_df_colorV["label"] == comp].iloc[0]

                comp_b_minus_v = comp_star_B["true_inst_mag"] - comp_star_V["true_inst_mag"]
                delta_b_minus_v = (check_IMB - check_IMV) - comp_b_minus_v
                delta_B_minus_V = tbv_coefficient*delta_b_minus_v
                delta_b = check_IMB - comp_star_B["true_inst_mag"]
                delta_v = check_IMV - comp_star_V["true_inst_mag"]
                B_star = delta_b + (tb_bv_coefficient*delta_B_minus_V) + comp_star_B["match_mag"]
                V_star = delta_v + (tv_bv_coefficient*delta_B_minus_V) + comp_star_V["match_mag"]
                
                
                result_check_star = result_check_star.append({
                    "star": "check",
                    "label": int(comp),
                    "IMB": comp_star_B["true_inst_mag"],
                    "IMV": comp_star_V["true_inst_mag"],
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
                    }, ignore_index=True)
          

                # comp_b_minus_v , calculated above
                delta_b_minus_v = (var_IMB - var_IMV) - comp_b_minus_v
                delta_B_minus_V = tbv_coefficient*delta_b_minus_v
                delta_b = var_IMB - comp_star_B["true_inst_mag"]
                delta_v = var_IMV - comp_star_V["true_inst_mag"]
                B_star = delta_b + (tb_bv_coefficient*delta_B_minus_V) + comp_star_B["match_mag"]
                V_star = delta_v + (tv_bv_coefficient*delta_B_minus_V) + comp_star_V["match_mag"]
                
                
                result_var_star = result_var_star.append({
                    "star": "var",
                    "label": int(comp),
                    "IMB": comp_star_B["true_inst_mag"],
                    "IMV": comp_star_V["true_inst_mag"],
                    "B": comp_star_B["match_mag"],
                    "V": comp_star_V["match_mag"],
                    "delta_b_minus_v": delta_b_minus_v,
                    "delta_B_minus_V": delta_B_minus_V,
                    "delta_b": delta_b,
                    "delta_v": delta_v,
                    "comp_b_minus_v": comp_b_minus_v,
                    "B_star": B_star,
                    "V_star": V_star
                    }, ignore_index=True)
          
            B_mean_check = result_check_star["B_star"].mean()
            V_mean_check = result_check_star["V_star"].mean()
            B_mean_var = result_var_star["B_star"].mean()
            V_mean_var = result_var_star["V_star"].mean()
            
            B_std_check = result_check_star["B_star"].std()
            V_std_check = result_check_star["V_star"].std()
            B_std_var = result_var_star["B_star"].  std()
            V_std_var = result_var_star["V_star"].std()

            #Find IQR and iqr 
            q3, q1 = np.percentile(result_check_star["B_star"], [75 ,25])
            b_iqr = q3 - q1
            b_upper_limit = q3 + (b_iqr * 1.5)
            b_lower_limit = q1 - (b_iqr * 1.5)

            q3, q1 = np.percentile(result_check_star["V_star"], [75 ,25])
            v_iqr = q3 - q1
            v_upper_limit = q3 + (v_iqr * 1.5)
            v_lower_limit = q1 - (v_iqr * 1.5)

            #record any outliers
            
            if (result_check_star["B_star"] < b_lower_limit).any():
                result_check_star.loc[result_check_star.loc[result_check_star["B_star"] < b_lower_limit].index, "outlier"] = "<--OUTLIER"
            
            if (result_check_star["B_star"] > b_upper_limit).any():
                result_check_star.loc[result_check_star.loc[result_check_star["B_star"] > b_upper_limit].index, "outlier"] = "<--OUTLIER"
            
            if (result_check_star["V_star"] < v_lower_limit).any():
                result_check_star.loc[result_check_star.loc[result_check_star["V_star"] < v_lower_limit].index, "outlier"] = "<--OUTLIER"
            
            if (result_check_star["V_star"] > v_upper_limit).any():
                result_check_star.loc[result_check_star.loc[result_check_star["V_star"] > v_upper_limit].index, "outlier"] = "<--OUTLIER"
            
            #create an aux table containing misc data needed for AAVSO report
            #this data is appended to the notes section 
            #(See E:\Astronomy\AAVSO\Reports\AAVSO Reports\MAO\2022 8 1 V1117 Her\AAVSOReport_V1117-Her_B_20220802.txt)
            result_aux_report = pd.DataFrame(columns=["color", "JD", "KMAGS", "KMAGINS", "KREFMAG", "Tb_bv", "T_bv", "VMAGINS", "Date-Obs"])
            
            result_aux_report = result_aux_report.append({
                "color" : "B",
                "JD" : comp_star_B["date-obs"],
                "KMAGS" : B_mean_check,
                "KMAGINS" : check_IMB,
                "KREFMAG" : check_B,
                "T_bv" : tbv_coefficient,
                "Tv_bv" : tv_bv_coefficient,
                "VMAGINS" : var_IMB,
                "Date-Obs" : date_obs_B
                }, ignore_index=True)                     
            
            result_aux_report = result_aux_report.append({
                "color" : "V",
                "JD" : comp_star_V["date-obs"],
                "KMAGS" : V_mean_check,
                "KMAGINS" : check_IMV,
                "KREFMAG" : check_V,
                "T_bv" : tbv_coefficient,
                "Tv_bv" : tv_bv_coefficient,
                "VMAGINS" : var_IMV,
                "Date-Obs" : date_obs_V
                }, ignore_index=True)                     
            
            self.console_msg("\n")
            self.console_msg("Check Star Estimates using check star: " + str(int(check_star_B["label"])) + "\n" +
                             result_check_star.sort_values(by="label").to_string() +
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
                             result_var_star.sort_values(by="label").to_string() +
                             '\n' + 
                             ("B* Ave: " + format(B_mean_var, ' >6.3f') +
                             "  V* Ave: " + format(V_mean_var, ' >6.3f')).rjust(137) +
                             '\n' +
                             ("B* Std: " + format(B_std_var, ' >6.3f') +  
                             "  V* Std: " + format(V_std_var, ' >6.3f')).rjust(137))
                             
            #concat result_aux_report, result_check_star, and result_var_star
            master_frames = [result_check_star, result_var_star, result_aux_report]

            master_report = pd.concat(master_frames, keys=['check', 'var', 'aux'])
            
            dir_path = os.path.dirname(os.path.realpath(file_name)) + "\\"
            
            master_report.to_csv(dir_path + self.object_name_entry.get() + "-Master-Report.csv", index=False)
            self.console_msg("Master Report saved to " + str(dir_path + self.object_name_entry.get() + "-Master-Report.csv"))



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


            if True:
                self.console_msg("Updating photometry table with Sky coordinates...")
                for index, row in self.results_tab_df.iterrows():
                    sky = self.wcs_header.pixel_to_world(row["x_fit"], row["y_fit"])
                    c = SkyCoord(ra=sky.ra, dec=sky.dec)
                    self.results_tab_df.loc[index, "ra_fit"] = c.ra / u.deg
                    self.results_tab_df.loc[index, "dec_fit"] = c.dec / u.deg
                self.results_tab_df.to_csv(self.image_file + ".csv", index=False)


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



            if True:
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
                        #uniq
                        already_gotten = self.results_tab_df.loc[self.results_tab_df["label"] == str(match_label)]    
                        if not already_gotten.empty:
                            #no need to add 'ghost'
                            continue

                        #Here if not already_gotten                        
                        self.results_tab_df.loc[index, "match_id"] = \
                            str(self.catalog_stringvar.get()) + \
                                " " + str(match_id)
                        self.results_tab_df.loc[index, "label"] = str(match_label)
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
                        self.results_tab_df.loc[index, "check_star"] = ""

                if self.remove_vsx_var.get():
                    self.console_msg(
                        "Inquiring VizieR (B/vsx/vsx) for VSX variables in the field...")
                    vsx_result = Vizier(
                        catalog="B/vsx/vsx", row_limit=-1).query_region(frame_center, frame_radius)
                    # print(vsx_result)
                    matched_match_id = 0 #init
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
                                if match_id == matched_match_id:
                                    self.console_msg("WARNING: VSX source: " + str(match_id) + " already matched! XXghostXX")
                                    self.results_tab_df.loc[index, "vsx_id"] = str(match_id) + "XXghostXX"
                                else:
                                    self.results_tab_df.loc[index, "vsx_id"] = str(match_id)
                                    self.console_msg("Match VSX source: " + str(match_id))
                                matched_match_id = match_id
                            else:
                                self.results_tab_df.loc[index, "vsx_id"] = ""

                            """
                            Removing the following because of confusion
                            
                            # Check 30 arcsec vicinity for nearby VSX sources:
                            if separation < 0.00833333 * u.deg:  <----- always true if separation < matching_radius * u.deg
                                self.results_tab_df.loc[index, "nearby_vsx_id"] = str(
                                    match_id)
                                self.results_tab_df.loc[index,
                                                        "nearby_vsx_separation"] = separation
                            else:
                                self.results_tab_df.loc[index,
                                                        "nearby_vsx_id"] = ""
                                self.results_tab_df.loc[index,
                                                        "nearby_vsx_separation"] = ""
                            """
                    else:
                        self.console_msg("Found no VSX sources in the field.")
                    
                self.results_tab_df.to_csv(
                    self.image_file + ".csv", index=False)
                self.console_msg("Photometry table saved (with following comp stars) to " +
                                 str(self.image_file + ".csv"))

                if using_aavso_catalog:
                    comp_list = ''
                    #output comp_stars_used
                    found_check = False #init
                    for comp in comp_stars_used:
                        (label, ischeck) = comp
                        found_check |= ischeck == "True"
                        comp_list += str(label) + ', ' 
                    self.console_msg("AAVSO comp stars used: " + comp_list)
                        
                    check_star = self.object_kref_entry.get().strip()
                    if not found_check and check_star != '':
                        self.console_msg("WARNING!! The requested check_star: " + check_star + " was NOT FOUND in image! Choose another.", level=logging.WARNING)
                        
                self.display_image()
                
                        
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

    def plot_photometry_menu_action(self):
        self.plot_photometry()
        self.update_display()

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

    def find_linear_regression_model(self):
        # Select stars with catalog matches with known magnitudes
        self.console_msg(
            "Finding linear regression model with ensemble stars...")
        vsx_ids_in_photometry_table = False
        if "vsx_id" in self.results_tab_df:
            vsx_ids_in_photometry_table = True
        try:
            min_comparison_magnitude = int(
                self.min_ensemble_magnitude_entry.get())
            max_comparison_magnitude = int(
                self.max_ensemble_magnitude_entry.get())

            # Only stars with known magnitude from catalogs
            mask = self.results_tab_df['match_mag'].map(
                self.safe_float_convert)
            ensemble_stars = self.results_tab_df.loc[mask]

            # Only stars without VSX id
            if self.remove_vsx_var.get() and vsx_ids_in_photometry_table:
                self.console_msg(
                    "Ensemble size before removal of VSX sources: "+str(len(ensemble_stars)))
                mask = ensemble_stars['vsx_id'].isnull()
                ensemble_stars = ensemble_stars.loc[mask]

            # Convert magnitudes to floats
            ensemble_stars["inst_mag"] = ensemble_stars["inst_mag"].astype(
                float)
            ensemble_stars["match_mag"] = ensemble_stars["match_mag"].astype(
                float)

            # Filter by magnitude
            mask = ensemble_stars['match_mag'] < max_comparison_magnitude
            ensemble_stars = ensemble_stars.loc[mask]
            mask = ensemble_stars['match_mag'] > min_comparison_magnitude
            ensemble_stars = ensemble_stars.loc[mask]

            # Remove removed outliers
            mask = ensemble_stars['removed_from_ensemble'] == False
            ensemble_stars = ensemble_stars.loc[mask]

            weights = None

            if self.weighting_stringvar.get() == "Raw Flux":
                ensemble_stars = ensemble_stars.sort_values(
                    by=['flux_0'], ascending=False)
                ensemble_stars = ensemble_stars.head(
                    n=int(self.ensemble_limit_entry.get()))
                weights = ensemble_stars["flux_0"]

            if self.weighting_stringvar.get() == "Instrumental Magnitude":
                ensemble_stars = ensemble_stars.sort_values(
                    by=['inst_mag'], ascending=True)
                ensemble_stars = ensemble_stars.head(
                    n=int(self.ensemble_limit_entry.get()))
                weights = 1/ensemble_stars["inst_mag"]

            if self.weighting_stringvar.get() == "PSF Sigma":
                ensemble_stars = ensemble_stars.sort_values(
                    by=['sigma_fit'], ascending=True)
                ensemble_stars = ensemble_stars.head(
                    n=int(self.ensemble_limit_entry.get()))
                weights = 1/ensemble_stars["sigma_fit"]

            self.console_msg(
                "Using "+str(len(ensemble_stars))+" ensemble stars.")

            x = ensemble_stars["inst_mag"]
            y = ensemble_stars["match_mag"]
            lb = ensemble_stars['label']

            # Fit a 1st degree polynomial
            self.a, self.b = np.polyfit(x, y, deg=1, w=weights)
            fit_fn = np.poly1d((self.a, self.b))
            yhat = fit_fn(x)
            ybar = np.sum(y) / len(y)
            ssreg = np.sum((yhat - ybar) ** 2)
            sstot = np.sum((y - ybar) ** 2)
            r_squared = ssreg / sstot
            self.console_msg(
                "Linear regression fit: a = " + str(round(self.a,4)) +
                " b = " + str(round(self.b,4)) +
                " r^2 = " + str(round(r_squared,4)))

            self.linreg_plot.clear()
            self.linreg_plot.plot(x, y, 'ro', ms=5) #ms is size of red circle
            self.linreg_plot.plot(x, fit_fn(x), '--k', ms=1)
            self.linreg_plot.text(
                0.1, 0.85, "n = "+str(len(ensemble_stars)), transform=self.linreg_plot.transAxes)
            self.plot_canvas.draw()
            self.ensemble_size = len(ensemble_stars)

            for index, row in self.results_tab_df.iterrows():
                self.results_tab_df.loc[index, "differential_mag"] = round(
                    float(row["inst_mag"]) * self.a + self.b, 3)

            for index, row in ensemble_stars.iterrows():
                ensemble_stars.loc[index, "differential_mag"] = round(
                    float(row["inst_mag"]) * self.a + self.b, 3)

            for index, row in self.results_tab_df.iterrows():
                try:
                    self.results_tab_df.loc[index, "mag_error"] = abs(
                        float(row["differential_mag"]) - float(row["match_mag"]))
                except:
                    continue

            for index, row in ensemble_stars.iterrows():
                try:
                    # ensemble_stars.loc[index, "mag_error"] = abs(
                    #    float(row["differential_mag"]) - float(row["match_mag"]))
                    ensemble_stars.loc[index, "mag_error_squared"] = (
                        float(row["differential_mag"]) - float(row["match_mag"]))**2
                except:
                    continue
            self.linreg_error = np.sqrt(ensemble_stars["mag_error_squared"].sum() / (len(ensemble_stars)-1))        # Standard deviation
            decimal_places = int(self.decimal_places_entry.get())
            self.console_msg("Linear Regression Error: " + str(round(self.linreg_error, decimal_places)) +
                             " Number of ensemble stars: " + str(len(ensemble_stars)))

            self.results_tab_df.to_csv(self.image_file + ".csv", index=False)
            self.console_msg("Photometry table with differential magnitudes saved to " + str(self.image_file + ".csv"))
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)

    def remove_fit_outlier(self):
        # Find not yet removed outlier within the ensemble magnitude range
        self.console_msg("Processing...")
        error = 0
        outlier_index = 0
        min_comparison_magnitude = int(self.min_ensemble_magnitude_entry.get())
        max_comparison_magnitude = int(self.max_ensemble_magnitude_entry.get())
        for index, row in self.results_tab_df.iterrows():
            if abs(row["mag_error"]) > error and row["removed_from_ensemble"] == False and float(row[
                    "match_mag"]) > min_comparison_magnitude and float(row["match_mag"]) < max_comparison_magnitude:
                outlier_index = index
                error = abs(row["mag_error"])
        self.results_tab_df.loc[outlier_index, "removed_from_ensemble"] = True
        self.find_linear_regression_model()
        self.update_display()

    def remove_distant_fit_outliers(self):
        # Removing outliers further than X arcseconds from the target
        self.console_msg("Processing...")
        max_separation = int(self.max_outliers_separation_entry.get())
        target_coordinates = self.wcs_header.pixel_to_world(
            self.last_clicked_x, self.last_clicked_y)
        n = 0
        for index, row in self.results_tab_df.iterrows():
            comparison_coordinates = self.wcs_header.pixel_to_world(
                row["x_fit"], row["y_fit"])
            separation = comparison_coordinates.separation(target_coordinates)
            if separation > max_separation * 0.000277778 * u.deg:
                self.results_tab_df.loc[index, "removed_from_ensemble"] = True
                n = n + 1
        self.console_msg("Removed "+str(n) +
                         " sources from the ensemble as distant outliers.")
        self.find_linear_regression_model()
        self.update_display()

    def remove_fit_outliers_until_ensemble_limit(self):
        while self.ensemble_size > int(self.ensemble_limit_entry.get()):
            self.remove_fit_outlier()

    def reset_fit_outliers(self):
        self.console_msg("Processing...")
        self.results_tab_df["removed_from_ensemble"] = False
        self.find_linear_regression_model()
        self.update_display()

    def delete_photometry_file(self):
        if os.path.isfile(self.image_file + ".csv"):
            os.remove(self.image_file + ".csv")
            self.console_msg("Photometry file deleted.")
            self.a = 0
            self.b = 0
            self.update_display()

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
                settings.update(
                    {'photometry_aperture_entry': self.photometry_aperture_entry.get()})
                settings.update(
                    {'min_ensemble_magnitude_entry': self.min_ensemble_magnitude_entry.get()})
                settings.update(
                    {'max_ensemble_magnitude_entry': self.max_ensemble_magnitude_entry.get()})
                settings.update({'fwhm_entry': self.fwhm_entry.get()})
                settings.update(
                    {'star_detection_threshold_entry': self.star_detection_threshold_entry.get()})
                settings.update(
                    {'photometry_iterations_entry': self.photometry_iterations_entry.get()})
                settings.update({'sharplo_entry': self.sharplo_entry.get()})
                settings.update(
                    {'bkg_filter_size_entry': self.bkg_filter_size_entry.get()})
                settings.update(
                    {'matching_radius_entry': self.matching_radius_entry.get()})
                settings.update(
                    {'ensemble_limit_entry': self.ensemble_limit_entry.get()})
                settings.update(
                    {'decimal_places_entry': self.decimal_places_entry.get()})
                settings.update({'obscode_entry': self.obscode_entry.get()})
                settings.update(
                    {'aavso_obscode_entry': self.aavso_obscode_entry.get()})
                settings.update({'latitude_entry': self.latitude_entry.get()})
                settings.update(
                    {'longitude_entry': self.longitude_entry.get()})
                settings.update({'height_entry': self.height_entry.get()})
                settings.update(
                    {'telescope_entry': self.telescope_entry.get()})
                settings.update(
                    {'tbv_entry': self.tbv_entry.get()})
                settings.update(
                    {'tv_bv_entry': self.tv_bv_entry.get()})
                settings.update(
                    {'tb_bv_entry': self.tb_bv_entry.get()})
                settings.update(
                    {'accessory_entry': self.accessory_entry.get()})
                settings.update({'ccd_entry': self.ccd_entry.get()})
                settings.update({'ccd_gain_entry': self.ccd_gain_entry.get()})
                settings.update(
                    {'weighting_stringvar': self.weighting_stringvar.get()})
                settings.update(
                    {'catalog_stringvar': self.catalog_stringvar.get()})
                settings.update(
                    {'vizier_catalog_entry': self.vizier_catalog_entry.get()})
                settings.update({'remove_vsx_var': self.remove_vsx_var.get()})
                settings.update({'nearby_vsx_var': self.nearby_vsx_var.get()})
                settings.update(
                    {'max_outliers_separation_entry': self.max_outliers_separation_entry.get()})
                settings.update(
                    {'fitter_stringvar': self.fitter_stringvar.get()})
                settings.update({'batch_psf_var': self.batch_psf_var.get()})
                settings.update({'plate_solve_on_open': self.plate_solve_on_open.get()})
                settings.update(
                    {'crop_fits_entry': self.crop_fits_entry.get()})
                settings.update(
                    {'astrometrynet_entry': self.astrometrynet_entry.get()})
                settings.update(
                    {'astrometrynet_key_entry': self.astrometrynet_key_entry.get()})
                settings.update(
                    {'object_kref_entry': self.object_kref_entry.get()})
                settings.update(
                    {'object_sel_comp_entry': self.object_sel_comp_entry.get()})
                settings.update(
                    {'object_name_entry': self.object_name_entry.get()})
                settings.update(
                    {'object_notes_entry': self.object_notes_entry.get()})
                with open(str(file_name.name), 'w') as f:
                    w = csv.DictWriter(f, settings.keys())
                    w.writeheader()
                    w.writerow(settings)
                self.console_msg("Saved.")
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)

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
                        if type(getattr(self, key)) == tk.Entry:
                            self.set_entry_text(
                                getattr(self, key), settings[key])
                        if type(getattr(self, key)) == tk.StringVar:
                            getattr(self, key).set(settings[key])
                        if type(getattr(self, key)) == tk.BooleanVar:
                            getattr(self, key).set(settings[key])
                
                pass
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)

    def generate_baa_report(self):
        global image_width, image_height
        report_dir = "baa_reports"
        self.console_msg("Changing folder to " +
                         os.path.dirname(self.image_file))
        os.chdir(os.path.dirname(self.image_file))
        if not os.path.exists(report_dir):
            os.mkdir(report_dir)
        image_basename = os.path.basename(self.image_file)
        report_filename = os.path.join(report_dir, "BAA " + os.path.splitext(
            image_basename)[0] + " " + str(self.object_name_entry.get()) + ".txt")
        try:
            with open(report_filename, mode='w') as f:
                f.write("File Format\tCCD/DSLR v2.01\n")
                f.write("Observation Method\tCCD\n")
                f.write("Variable\t" + str(self.object_name_entry.get()) + "\n")
                # Generate Chart ID
                chart_id = self.catalog_stringvar.get()
                if self.catalog_stringvar.get() == "VizieR Catalog":
                    chart_id = "VizieR " + self.vizier_catalog_entry.get()
                frame_center = self.wcs_header.pixel_to_world(
                    int(image_width / 2), (image_height / 2))
                frame_top_left = self.wcs_header.pixel_to_world(0, 0)
                frame_bottom_left = self.wcs_header.pixel_to_world(
                    0, image_height)
                frame_top_right = self.wcs_header.pixel_to_world(
                    image_width, 0)
                frame_center_string = frame_center.to_string(
                    "hmsdms", precision=0, sep="")
                fov_horizontal = frame_top_right.separation(
                    frame_top_left).arcminute
                fov_vertical = frame_top_left.separation(
                    frame_bottom_left).arcminute
                chart_id = chart_id + " RA Dec:" + frame_center_string + " FOV:" + str(
                    int(fov_horizontal)) + "'x" + str(int(fov_vertical)) + "'"
                if len(chart_id) > 50:
                    self.console_msg(
                        "Chart ID length is longer than 50 bytes - correct manually!")
                f.write("Chart ID\t" + chart_id + "\n")
                f.write("Observer code\t" +
                        str(self.obscode_entry.get()) + "\n")
                f.write(
                    "Location\t" + str(self.latitude_entry.get()) + " " + str(self.longitude_entry.get()) + " H" + str(
                        self.height_entry.get()) + "m\n")
                f.write("Telescope\t" + str(self.telescope_entry.get()
                                            ) + " " + str(self.accessory_entry.get()) + "\n")
                f.write("Tbv\t" + str(self.tbv_entry.get()
                                            ) + " " + str(self.accessory_entry.get()) + "\n")
                f.write("Tv_bv\t" + str(self.tv_bv_entry.get()
                                            ) + " " + str(self.accessory_entry.get()) + "\n")
                f.write("Tb_bv\t" + str(self.tb_bv_entry.get()
                                            ) + " " + str(self.accessory_entry.get()) + "\n")
                f.write("Camera\t" + str(self.ccd_entry.get()) + "\n")
                f.write("Magnitude type\tPreCalculated\n")
                f.write("Photometry software\t" + self.program_full_name + "\n")
                f.write("\n")
                f.write(
                    "JulianDate\tFilter\tVarAbsMag\tVarAbsErr\tExpLen\tCmpStar\n")
                f.write(
                    str(round(float(self.date_obs_entry.get()), 5)) + "\t" + str(self.filter_entry.get()) + "\t" + str(
                        self.last_clicked_differential_magnitude) + "\t" + str(
                        self.last_clicked_differential_uncertainty) +
                    "\t" + str(self.exposure_entry.get()) + "\t" + "Ensemble\n")
            self.console_msg(
                "BAA Photometry Database report saved to " + str(report_filename))
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)

    def aavso_get_comparison_stars(self, frame_center, filter_band='V', field_of_view=18.5, maglimit=20):
        
        try:
            #Some telescopes use 'I' instead of 'Ic', but AAVSO charts use Ic
            if filter_band == 'I':
                filter_band = 'Ic'
            
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
            result = pd.DataFrame(columns=["AUID", "RA", "Dec", "Mag", "Label", "Chart ID"])


            for star in r.json()['photometry']:
                auid = star['auid']
                ra = star['ra']
                dec = star['dec']
                label = star['label']
                
                if label not in sel_comps:
                    continue #skip
                
                # if this label is the selected check star, mark it
                if label == check_star:
                    #mark it a True check star
                    is_check_star = "True"
                else:
                    is_check_star = "False"
                            
                
                for band in star['bands']:
                    if band['band'] == filter_band:
                        mag = band['mag']
                        
                result = result.append({
                    "AUID": auid,
                    "RA": ra,
                    "Dec": dec,
                    "Mag": mag,
                    "Label": label,
                    "Chart ID": chart_id,
                    "Check Star": is_check_star
                }, ignore_index=True)
                
            return result
        
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)

    def generate_aavso_report(self):
        global image_width, image_height
        self.console_msg("Beginning Generate AAVSO Report wCREF...")
        report_dir = "aavso_reports"
        self.console_msg("Changing folder to " +
                         os.path.dirname(self.image_file))
        os.chdir(os.path.dirname(self.image_file))
        if not os.path.exists(report_dir):
            os.mkdir(report_dir)
        image_basename = os.path.basename(self.image_file)
        report_filename = os.path.join(report_dir, "AAVSO " + os.path.splitext(
            image_basename)[0] + " " + str(self.object_name_entry.get()) + ".txt")
        frame_center = self.wcs_header.pixel_to_world(
            int(image_width / 2), (image_height / 2))
        frame_top_left = self.wcs_header.pixel_to_world(0, 0)
        frame_bottom_left = self.wcs_header.pixel_to_world(0, image_height)
        frame_top_right = self.wcs_header.pixel_to_world(image_width, 0)
        fov_horizontal = frame_top_right.separation(frame_top_left).arcminute
        fov_vertical = frame_top_left.separation(frame_bottom_left).arcminute

        self.console_msg("Getting AAVSO Comparison Stars...")
        comparison_stars = self.aavso_get_comparison_stars(frame_center, filter_band=str(
            self.filter_entry.get()), field_of_view=fov_horizontal, maglimit=self.max_ensemble_magnitude_entry.get())

        self.console_msg("Matching comparison stars to the photometry table...")
        matching_radius = float(
            self.matching_radius_entry.get()) * 0.000277778  # arcsec to degrees
        catalog_comparison = SkyCoord(
            comparison_stars["RA"], comparison_stars["Dec"], unit=(u.hourangle, u.deg))
        for index, row in self.results_tab_df.iterrows():
            photometry_star_coordinates = SkyCoord(
                ra=row["ra_fit"] * u.deg, dec=row["dec_fit"] * u.deg)
            match_index, d2d_match, d3d_match = photometry_star_coordinates.match_to_catalog_sky(
                catalog_comparison)
            match_id = comparison_stars.iloc[match_index]["AUID"]
            match_ra = comparison_stars.iloc[match_index]["RA"]
            match_dec = comparison_stars.iloc[match_index]["Dec"]
            match_mag = comparison_stars.iloc[match_index]["Mag"]
            match_label = comparison_stars.iloc[match_index]["Label"]
            match_chartId = comparison_stars.iloc[match_index]["Chart ID"]
            match_coordinates = SkyCoord(
                ra=match_ra, dec=match_dec, unit=(u.hourangle, u.deg))
            separation = photometry_star_coordinates.separation(
                match_coordinates)
            if separation < matching_radius * u.deg:
                self.results_tab_df.loc[index, "auid"] = str(match_id)
                self.results_tab_df.loc[index, "label"] = str(match_label)
                self.results_tab_df.loc[index, "auid_mag"] = float(match_mag)
                self.results_tab_df.loc[index, "chart ID"] = str(match_chartId)
            else:
                self.results_tab_df.loc[index, "auid"] = ""
                self.results_tab_df.loc[index, "label"] = ""
                self.results_tab_df.loc[index, "auid_mag"] = np.nan
                self.results_tab_df.loc[index, "chart ID"] = ""

        self.console_msg("Searching for the best check star for magnitude " +
                         str(self.last_clicked_differential_magnitude))

        self.results_tab_df["difference_between_differential_mag_and_auid_mag"] = abs(
            self.results_tab_df["differential_mag"] - self.results_tab_df["auid_mag"])
        self.results_tab_df["difference_between_target_mag_and_auid_mag"] = abs(float(
            self.last_clicked_differential_magnitude) - self.results_tab_df["auid_mag"])

        # Find comparison star that has the least magntiude difference to target star
        check_star = self.results_tab_df[self.results_tab_df["difference_between_target_mag_and_auid_mag"].eq(
            self.results_tab_df["difference_between_target_mag_and_auid_mag"].min())].iloc[0]
        self.console_msg("Found check star "+str(check_star["label"])+" of magnitude " + str(
            check_star["auid_mag"]) + ", difference with current photometry: " +
                str(round(check_star["difference_between_differential_mag_and_auid_mag"], 4)))
        self.console_msg("Removing this comparison star from ensemble...")
        self.results_tab_df.loc[self.results_tab_df["auid"] ==
                                check_star["auid"], "removed_from_ensemble"] = True
        self.console_msg("Finding new fit and differential magntiude...")
        self.find_linear_regression_model()
        x = self.last_clicked_x
        y = self.last_clicked_y
        self.canvas.xview(tk.MOVETO, 0)     # Reset canvas position
        self.canvas.yview(tk.MOVETO, 0)
        # Simulate a mouse click to redo photometry
        self.canvas.event_generate('<Button-1>', x=x, y=y)
        width, height = generated_image.size
        canvas_width = self.canvas.winfo_width()
        canvas_height = self.canvas.winfo_height()
        # Center canvas on the object
        self.canvas.xview(tk.MOVETO, (x-canvas_width/2)/width)
        self.canvas.yview(tk.MOVETO, (y-canvas_height/2)/height)
        # Updating comparison star with new photometry
        check_star = self.results_tab_df[self.results_tab_df["difference_between_target_mag_and_auid_mag"].eq(
            self.results_tab_df["difference_between_target_mag_and_auid_mag"].min())].iloc[0]

        try:
            with open(report_filename, mode='w') as f:
                f.write("#TYPE=Extended\n")
                f.write("#OBSCODE="+self.aavso_obscode_entry.get()+"\n")
                f.write("#SOFTWARE="+self.program_full_name+"\n")
                f.write("#DELIM=,\n")
                f.write("#DATE=JD\n")
                f.write("#OBSTYPE=CCD\n")
                f.write(
                    "#STARID,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")

                starid = str(self.object_name_entry.get())
                date = str(round(float(self.date_obs_entry.get()), 5))
                mag = str(self.last_clicked_differential_magnitude)
                merr = str(self.last_clicked_differential_uncertainty)
                filt = str(self.filter_entry.get())
                trans = "NO"
                mtype = "STD"
                cname = "ENSEMBLE"
                cmag = "na"
                kname = str(check_star["label"])
                kmag = str(check_star["differential_mag"])
                amass = "na"
                group = "na"
                chart = str(check_star["chart ID"])
                if self.catalog_stringvar.get() == "VizieR Catalog":
                    chart = "VizieR " + self.vizier_catalog_entry.get()
                notes = str(self.object_notes_entry.get())
                
                if notes == "":
                    notes = "na"
                    
                # Add " " after notes, because TA clobbers last char
                f.write(starid+","+date+","+mag+","+merr+","+filt+","+trans+","+mtype+"," +
                        cname+","+cmag+","+kname+","+kmag+","+amass+","+group+","+chart+","+notes+" ")

            self.console_msg(
                "AAVSO Photometry Database report saved to " + str(report_filename))
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)

    def generate_aavso_report_2color(self):
        global image_width, image_height
        
        """
        Typical AAVSO Two Color Report:
            
        #TYPE=EXTENDED
        #OBSCODE=FPIA
        #SOFTWARE=VPhot 4.0.27
        #DELIM=,
        #DATE=JD
        #OBSTYPE=CCD
        #NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
        V1117 Her,2459735.67981,12.891,0.013,B,YES,STD,ENSEMBLE,na,120,13.019,na,na,X27876MZ,Mittelman ATMoB Observatory|KMAG=13.019|KMAGINS=-8.706|KREFMAG=13.019|Tbv=1.182|VMAGINS=-8.809
        V1117 Her,2459735.66950,12.549,0.007,V,YES,STD,ENSEMBLE,na,120,12.029,na,na,X27876MZ,Mittelman ATMoB Observatory|KMAG=12.029|KMAGINS=-9.497|KREFMAG=12.033|Tv_bv=-0.115|VMAGINS=-9.051

        """


        self.console_msg("Beginning Generate AAVSO Two Color Ensemble Report...")
        var_star_name = self.object_name_entry.get().strip()
        if len(var_star_name) == 0:
            self.console_msg(
                "Object Name must be specified; eg. 'V1117 Her'")
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
            image_basename)[0] + " " + str(self.object_name_entry.get()) + ".txt")

        """
        Report is generated from 'Saved' master_result DataFrame
        that was saved as <object_name_entry>-Master-Result.csv
        
        """
        #Ask user for <ObjectID>-Master-Report.csv
        # Ask for the B and V CSV files
        options = {}
        options['defaultextension'] = '.csv'
        options['filetypes'] = [('CSV', '.csv')]
        options['title'] = 'Choose the ' + var_star_name + '-Master-Report.csv'

        file_name = fd.askopenfilename(**options)

        if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
            self.console_msg("Loading Master-Report from " + str(file_name))
            master_report = pd.read_csv(str(file_name))
        else:
            return
            

        decimal_places = 3 #report is usually 3
        
        #extract the estimates
        result_check_star = master_report[master_report["star"] == "check"]
        result_var_star = master_report[master_report["star"] == "var"]

        B_mean_check = result_check_star["B_star"].mean()
        V_mean_check = result_check_star["V_star"].mean()
        B_mean_var = result_var_star["B_star"].mean()
        V_mean_var = result_var_star["V_star"].mean()
        
        B_std_check = result_check_star["B_star"].std()
        V_std_check = result_check_star["V_star"].std()
        B_std_var = result_var_star["B_star"].std()
        V_std_var = result_var_star["V_star"].std()
            
        self.console_msg("Check Star Estimates\n" + 
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
        
        
        pass
                         

        try:
            with open(report_filename, mode='w') as f:
                f.write("#TYPE=Extended\n")
                f.write("#OBSCODE="+self.aavso_obscode_entry.get()+"\n")
                f.write("#SOFTWARE=Self-developed using photoutils.psf; DAOPHOT\n") 
                f.write("#DELIM=,\n")
                f.write("#DATE=JD\n")
                f.write("#OBSTYPE=CCD\n")
                f.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")

                #B filter
                aux_result = master_report[master_report["color"] == "B"]
                
                starid = str(self.object_name_entry.get())
                date = format(float(aux_result['Date-Obs']), '.5f') 
                mag = str(round(B_mean_var, decimal_places))
                merr = str(round(B_std_check,decimal_places))
                filt = "B"
                trans = "YES"
                mtype = "STD"
                cname = "ENSEMBLE"
                cmag = "na"
                kname = self.object_kref_entry.get().strip()
                kmag = str(round(B_mean_check, decimal_places))
                amass = "na"
                group = "na"
                chart = self.vizier_catalog_entry.get().strip()
                notes = self.object_notes_entry.get().strip()
                notes += "|KMAG=" + str(round(B_mean_check, decimal_places)) + \
                         "|KMAGINS=" + str(round(float(aux_result['KMAGINS']), decimal_places)) + \
                         "|KREFMAG=" + str(round(float(aux_result['KREFMAG']), decimal_places)) + \
                         "|Tbv="  + str(round(float(aux_result['T_bv']), decimal_places)) + \
                         "|VMAGINS=" + str(round(float(aux_result['VMAGINS']), decimal_places))
                
                
                # Add " " after notes, because TA clobbers last char
                f.write(starid+","+date+","+mag+","+merr+","+filt+","+trans+","+mtype+"," +
                        cname+","+cmag+","+kname+","+kmag+","+amass+","+group+","+chart+","+notes+" \n")

                #V filter (easier just to repeat most of this)
                
                aux_result = master_report[master_report["color"] == "V"]
                
                starid = str(self.object_name_entry.get())
                date = format(float(aux_result['Date-Obs']), '.5f') 
                mag = str(round(V_mean_var, decimal_places))
                merr = str(round(V_std_check,decimal_places))
                filt = "V"
                trans = "YES"
                mtype = "STD"
                cname = "ENSEMBLE"
                cmag = "na"
                kname = self.object_kref_entry.get().strip()
                kmag = str(round(V_mean_check, decimal_places))
                amass = "na"
                group = "na"
                chart = self.vizier_catalog_entry.get().strip()
                notes = self.object_notes_entry.get().strip()
                notes += "|KMAG=" + str(round(V_mean_check, decimal_places)) + \
                         "|KMAGINS=" + str(round(float(aux_result['KMAGINS']), decimal_places)) + \
                         "|KREFMAG=" + str(round(float(aux_result['KREFMAG']), decimal_places)) + \
                         "|Tv_bv="  + str(round(float(aux_result['Tv_bv']), decimal_places)) + \
                         "|VMAGINS=" + str(round(float(aux_result['VMAGINS']), decimal_places))
                
                # Add " " after notes, because TA clobbers last char
                f.write(starid+","+date+","+mag+","+merr+","+filt+","+trans+","+mtype+"," +
                        cname+","+cmag+","+kname+","+kmag+","+amass+","+group+","+chart+","+notes+" \n")

            self.console_msg(
                "AAVSO Photometry Database report saved to " + str(report_filename))

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)


    def next_vsx_source(self):
        mask = self.results_tab_df['vsx_id'].notnull()
        vsx_sources = self.results_tab_df.loc[mask]
        next_vsx_source_id = ""
        next_row = False
        if self.object_name_entry.get() == "":
            # Reset to first VSX source
            next_vsx_source_id = vsx_sources.iloc[0]['vsx_id']
        if len(vsx_sources) > 0 and self.object_name_entry.get() != "":
            for index, row in vsx_sources.iterrows():
                if next_row:
                    next_vsx_source_id = row['vsx_id']
                    break
                if row['vsx_id'] == self.object_name_entry.get():
                    next_row = True
            if next_vsx_source_id == "":
                # Reset to first VSX source
                next_vsx_source_id = vsx_sources.iloc[0]['vsx_id']
        self.console_msg("Next VSX Source: "+next_vsx_source_id)
        x = int(self.results_tab_df.loc[self.results_tab_df['vsx_id']
                == next_vsx_source_id]['x_0'] * self.zoom_level)
        y = int(self.results_tab_df.loc[self.results_tab_df['vsx_id']
                == next_vsx_source_id]['y_0'] * self.zoom_level)
        self.canvas.xview(tk.MOVETO, 0)     # Reset canvas position
        self.canvas.yview(tk.MOVETO, 0)
        # Simulate a mouse click
        self.canvas.event_generate('<Button-1>', x=x, y=y)
        width, height = generated_image.size
        canvas_width = self.canvas.winfo_width()
        canvas_height = self.canvas.winfo_height()
        # Center canvas on the object
        self.canvas.xview(tk.MOVETO, (x-canvas_width/2)/width)
        self.canvas.yview(tk.MOVETO, (y-canvas_height/2)/height)

    def report_on_all_vsx_sources(self):
        mask = self.results_tab_df['vsx_id'].notnull()
        vsx_sources = self.results_tab_df.loc[mask]
        n = len(vsx_sources)
        for i in range(0, n):
            self.next_vsx_source()
            self.reset_fit_outliers()
            self.remove_distant_fit_outliers()
            self.remove_fit_outliers_until_ensemble_limit()
            x = self.last_clicked_x
            y = self.last_clicked_y
            self.canvas.xview(tk.MOVETO, 0)  # Reset canvas position
            self.canvas.yview(tk.MOVETO, 0)
            # Simulate a mouse click to redo photometry
            self.canvas.event_generate('<Button-1>', x=x, y=y)
            self.console_msg("Generating AAVSO report...")
            self.generate_aavso_report()
            self.console_msg("Generating BAA report...")
            self.generate_baa_report()
        self.console_msg("Report on all VSX sources finished for this frame.")

    def batch_report_on_all_vsx_sources(self):
        self.error_raised = False
        batch_dir = fd.askdirectory()
        if len(batch_dir) > 0:
            self.console_msg("Changing working directory to "+batch_dir)
            os.chdir(batch_dir)
            dir_list = os.listdir(batch_dir)
            for filename in dir_list:
                if os.path.splitext(filename)[1].lower() == ".fits" or os.path.splitext(filename)[1].lower() == ".fit":
                    try:
                        if self.error_raised:
                            self.console_msg(
                                "Aborting batch processing due to an error.")
                            break
                        full_filename = os.path.join(batch_dir, filename)
                        self.console_msg("Processing " + full_filename)
                        self.load_FITS(full_filename)
                        self.display_image()
                        if self.batch_psf_var.get():
                            self.perform_psf_photometry()
                        else:
                            self.aperture_photometry()
                        self.solve_image()
                        self.get_comparison_stars()
                        self.console_msg("Finding linear fit...")
                        self.find_linear_regression_model()
                        self.console_msg("Reporting on all VSX sources...")
                        self.report_on_all_vsx_sources()
                    except Exception as e:
                        self.error_raised = True
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)
                        continue
            self.console_msg("Batch processing finished.")

    def display_curve_from_baa_reports(self):
        if len(self.object_name_entry.get()) == 0:
            self.console_msg(
                "Enter object name or click on a source with name.")
            return
        dir = fd.askdirectory()
        object_name = self.object_name_entry.get()
        filter = self.filter_entry.get()
        light_curve = pd.DataFrame()
        if len(dir) > 0:
            dir_list = os.listdir(dir)
            for filename in dir_list:
                if object_name + ".txt" in filename:
                    self.console_msg("Adding point from "+filename)
                    full_filename = os.path.join(dir, filename)
                    # Skip to the 12th line with data
                    point = pd.read_csv(full_filename, header=10, sep='\t')
                    try:
                        light_curve = light_curve.append(
                            point, ignore_index=True)
                    except Exception as e:
                        self.error_raised = True
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)
                        continue
            fig, ax = plt.subplots(figsize=(12, 7))
            plt.gca().invert_yaxis()
            ax.xaxis.set_tick_params(
                which='major', size=6, width=1, direction='in', top='on')
            ax.xaxis.set_tick_params(
                which='minor', size=6, width=1, direction='in', top='on')
            ax.yaxis.set_tick_params(
                which='major', size=6, width=1, direction='in', right='on')
            ax.yaxis.set_tick_params(
                which='minor', size=6, width=1, direction='in', right='on')
            for axis in ['top', 'bottom', 'left', 'right']:
                ax.spines[axis].set_linewidth(1)
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize('medium')
                item.set_family('serif')
            ax.xaxis.label.set_fontsize('x-large')
            ax.yaxis.label.set_fontsize('x-large')
            formatter = ScalarFormatter(useOffset=False)
            ax.xaxis.set_major_formatter(formatter)
            title = object_name+"\n"
            if len(self.obscode_entry.get()) > 0:
                title = title + "Observer code "+self.obscode_entry.get()+", "
            if len(self.telescope_entry.get()) > 0:
                title = title + self.telescope_entry.get()+", "
            if len(self.tbv_entry.get()) > 0:
                title = title + self.tbv_entry.get()+", "
            if len(self.tv_bv_entry.get()) > 0:
                title = title + self.tv_bv_entry.get()+", "
            if len(self.tb_bv_entry.get()) > 0:
                title = title + self.tb_bv_entry.get()+", "
            if len(self.accessory_entry.get()) > 0:
                title = title + self.accessory_entry.get()+", "
            if len(self.ccd_entry.get()) > 0:
                title = title + self.ccd_entry.get()
            ax.set_title(title)
            ax.title.set_fontsize('x-large')
            ax.errorbar(light_curve["JulianDate"]-2400000.5, light_curve["VarAbsMag"],
                        yerr=light_curve["VarAbsErr"], fmt="o", color="black")
            ax.set_xlabel('MJD', labelpad=10)
            ax.set_ylabel(filter, labelpad=10)
            plt.show()

    def exit_app(self):
        os._exit(0)

    def __init__(self):
        #Wis heisen Sie?
        self.program_name = "MAOPhot"
        self.program_version = "1.0"
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
            label="Load Settings...", command=self.load_settings)
        self.filemenu.add_command(
            label="Save Settings As...", command=self.save_settings_as)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=self.exit_app)
        self.menubar.add_cascade(label="File", menu=self.filemenu)

        self.viewmenu = tk.Menu(self.menubar, tearoff=0)
        self.viewmenu.add_command(label="Update", command=self.update_display)
        self.viewmenu.add_separator()
        self.viewmenu.add_command(label="Zoom In", command=self.zoom_in)
        self.viewmenu.add_command(label="Zoom Out", command=self.zoom_out)
        self.viewmenu.add_command(label="100% Zoom", command=self.zoom_100)
        self.viewmenu.add_separator()
        self.viewmenu.add_command(
            label="Next VSX Source", command=self.next_vsx_source)
        self.menubar.add_cascade(label="View", menu=self.viewmenu)

        self.photometrymenu = tk.Menu(self.menubar, tearoff=0)
        self.photometrymenu.add_command(
            label="Iteratively Subtracted PSF Photometry", command=self.perform_psf_photometry)
        self.photometrymenu.add_command(
            label="Aperture Photometry", command=self.aperture_photometry)
        self.photometrymenu.add_separator()
        self.photometrymenu.add_command(
            label="Plot", command=self.plot_photometry_menu_action)
        self.photometrymenu.add_command(
            label="Hide", command=self.hide_photometry)
        self.photometrymenu.add_separator()
        self.photometrymenu.add_command(
            label="Solve Image", command=self.solve_image)
        self.photometrymenu.add_command(
            label="Get Comparison Stars", command=self.get_comparison_stars)
        self.photometrymenu.add_command(
            label="Find Regression Model", command=self.find_linear_regression_model)
        self.photometrymenu.add_separator()
        self.photometrymenu.add_command(
            label="Remove Fit Outlier", command=self.remove_fit_outlier)
        self.photometrymenu.add_command(
            label="Remove Fit Outliers Until Ensemble Limit", command=self.remove_fit_outliers_until_ensemble_limit)
        self.photometrymenu.add_command(
            label="Remove Fit Outliers Beyond Separation Limit", command=self.remove_distant_fit_outliers)
        self.photometrymenu.add_command(
            label="Reset Fit Outliers", command=self.reset_fit_outliers)
        self.photometrymenu.add_separator()
        self.photometrymenu.add_command(
            label="Delete Photometry File", command=self.delete_photometry_file)
        self.photometrymenu.add_command(
            label="Display Background", command=self.display_background)
        #self.photometrymenu.add_command(label="Plot Sigma Heatmap", command=self.plot_sigma_heatmap)
        self.menubar.add_cascade(label="Photometry", menu=self.photometrymenu)
        
        self.two_color_photometry = tk.Menu(self.menubar, tearoff=0)
        self.two_color_photometry.add_command(
            label="Two Color Photometry (B,V)", command=self.two_color_photometry_BV)
        self.menubar.add_cascade(label="Two Color Photometry", menu=self.two_color_photometry)

        self.reportmenu = tk.Menu(self.menubar, tearoff=0)
        self.reportmenu.add_command(
            label="BAA: Generate Report", command=self.generate_baa_report)
        self.reportmenu.add_command(
            label="BAA: Light Curve from Reports...", command=self.display_curve_from_baa_reports)
        self.reportmenu.add_separator()
        #self.reportmenu.add_command(
        #    label="AAVSO: Generate Report; ensemble", command=self.generate_aavso_report)
        self.reportmenu.add_command(
            label="AAVSO: Generate Report; Two Color Photometry", command=self.generate_aavso_report_2color)
        self.reportmenu.add_separator()
        self.reportmenu.add_command(
            label="BAA/AAVSO Reports on All VSX Sources", command=self.report_on_all_vsx_sources)
        self.reportmenu.add_command(
            label="BAA/AAVSO Batch Reports on All VSX Sources...", command=self.batch_report_on_all_vsx_sources)

        self.menubar.add_cascade(label="Report", menu=self.reportmenu)

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

        # We will lay out interface things into the new right_frame grid
        self.right_frame = tk.Frame(self.right_half)
        # Place right_frame into the top of the main canvas row, right next to it
        self.right_frame.grid(row=1, column=2, sticky=tk.N)
        self.plotname_label = tk.Label(self.right_frame, text="Plot:")
        self.plotname_label.grid(row=0, column=0)  # Place label
        #self.psf_canvas = tk.Canvas(self.right_frame, bg='grey', width=300, height=300) # Small PSF canvas
        self.fig_psf = Figure()
        self.psf_plot = self.fig_psf.add_subplot(111, projection='3d')
        # PSF 3D plot canvas - Matplotlib wrapper for Tk
        self.psf_plot_canvas = FigureCanvasTkAgg(
            self.fig_psf, self.right_frame)
        self.psf_plot_canvas.draw()
        self.psf_canvas = self.psf_plot_canvas.get_tk_widget()
        self.psf_canvas.config(
            width=int(self.screen_width/8.5), height=int(self.screen_width/8.5))
        # Allocate small PSF canvas to a new grid inside the right_frame
        self.psf_canvas.grid(row=1, column=0)   #was row0

        self.fig = Figure()
        self.linreg_plot = self.fig.add_subplot(111)
        # Linear regression canvas - Matplotlib wrapper for Tk
        self.plot_canvas = FigureCanvasTkAgg(self.fig, self.right_frame)
        self.plot_canvas.draw()
        self.linreg_canvas = self.plot_canvas.get_tk_widget()
        self.linreg_canvas.config(
            width=int(self.screen_width/8.5), height=int(self.screen_width/12))
        # Allocate small PSF canvas to a new grid inside the right_frame
        self.linreg_canvas.grid(row=2, column=0) #was row1

        # We will lay out interface things into the new right_frame grid
        self.left_frame = tk.Frame(self.left_half)
        # Place right_frame into the top of the main canvas row, right next to it
        self.left_frame.grid(row=0, column=0, sticky=tk.N)

        # Frame to hold settings grid
        self.settings_frame = tk.Frame(self.left_frame)
        # Settings_frame under the canvas in the right_frame
        self.settings_frame.grid(row=2, column=0, sticky=tk.NSEW)
        # Expand settings_frame column that holds labels
        tk.Grid.columnconfigure(self.settings_frame, 0, weight=1)

        settings_entry_width = 6
        extended_settings_entry_width = 30

        row = 0

        self.photometry_aperture_label = tk.Label(
            self.settings_frame, text="Fitting Width/Height, px:")
        self.photometry_aperture_label.grid(row=row, column=0, sticky=tk.W)
        self.photometry_aperture_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.photometry_aperture_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.photometry_aperture_entry, self.fit_shape)
        row = row + 1

        self.min_ensemble_magnitude_label = tk.Label(
            self.settings_frame, text="Minimum Ensemble Magnitude:")
        self.min_ensemble_magnitude_label.grid(row=row, column=0, sticky=tk.W)
        self.min_ensemble_magnitude_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.min_ensemble_magnitude_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.min_ensemble_magnitude_entry, "7")
        row = row + 1

        self.max_ensemble_magnitude_label = tk.Label(
            self.settings_frame, text="Maximum Ensemble Magnitude:")
        self.max_ensemble_magnitude_label.grid(row=row, column=0, sticky=tk.W)
        self.max_ensemble_magnitude_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.max_ensemble_magnitude_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.max_ensemble_magnitude_entry, "20")
        row = row + 1

        self.fwhm_label = tk.Label(self.settings_frame, text="FWHM, px:")
        self.fwhm_label.grid(row=row, column=0, sticky=tk.W)
        self.fwhm_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.fwhm_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.fwhm_entry, "4")
        row = row + 1

        self.star_detection_threshold_label = tk.Label(
            self.settings_frame, text="StarFinder Threshold k in (k * std):")
        self.star_detection_threshold_label.grid(
            row=row, column=0, sticky=tk.W)
        self.star_detection_threshold_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.star_detection_threshold_entry.grid(
            row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.star_detection_threshold_entry, "3.5")
        row = row + 1

        self.photometry_iterations_label = tk.Label(
            self.settings_frame, text="Photometry Iterations:")
        self.photometry_iterations_label.grid(row=row, column=0, sticky=tk.W)
        self.photometry_iterations_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.photometry_iterations_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.photometry_iterations_entry, "1")
        row = row + 1

        self.sharplo_label = tk.Label(
            self.settings_frame, text="Lower Bound for Sharpness:")
        self.sharplo_label.grid(row=row, column=0, sticky=tk.W)
        self.sharplo_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.sharplo_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.sharplo_entry, "0.5")
        row = row + 1

        self.bkg_filter_size_label = tk.Label(
            self.settings_frame, text="Background Median Filter, px:")
        self.bkg_filter_size_label.grid(row=row, column=0, sticky=tk.W)
        self.bkg_filter_size_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.bkg_filter_size_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.bkg_filter_size_entry, "1")
        row = row + 1

        self.filter_label = tk.Label(self.settings_frame, text="CCD Filter:")
        self.filter_label.grid(row=row, column=0, sticky=tk.W)
        self.filter_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.filter_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.filter_entry, "")
        row = row + 1

        self.exposure_label = tk.Label(
            self.settings_frame, text="Exposure Time:")
        self.exposure_label.grid(row=row, column=0, sticky=tk.W)
        self.exposure_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.exposure_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.exposure_entry, "0")
        row = row + 1

        self.matching_radius_label = tk.Label(
            self.settings_frame, text="Matching Radius, arcsec:")
        self.matching_radius_label.grid(row=row, column=0, sticky=tk.W)
        self.matching_radius_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.matching_radius_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.matching_radius_entry, "2")
        row = row + 1

        self.ensemble_limit_label = tk.Label(
            self.settings_frame, text="Limit Ensemble to:")
        self.ensemble_limit_label.grid(row=row, column=0, sticky=tk.W)
        self.ensemble_limit_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.ensemble_limit_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.ensemble_limit_entry, "1000")
        row = row + 1

        self.decimal_places_label = tk.Label(
            self.settings_frame, text="Decimal Places to Report:")
        self.decimal_places_label.grid(row=row, column=0, sticky=tk.W)
        self.decimal_places_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.decimal_places_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.decimal_places_entry, "3")
        row = row + 1

        self.max_outliers_separation_label = tk.Label(
            self.settings_frame, text="Ensemble Outliers Separation Limit, arcsec:")
        self.max_outliers_separation_label.grid(row=row, column=0, sticky=tk.W)
        self.max_outliers_separation_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.max_outliers_separation_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.max_outliers_separation_entry, "300")
        row = row + 1

        self.crop_fits_label = tk.Label(
            self.settings_frame, text="FITS Crop, %:")
        self.crop_fits_label.grid(row=row, column=0, sticky=tk.W)
        self.crop_fits_entry = tk.Entry(
            self.settings_frame, width=settings_entry_width)
        self.crop_fits_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.crop_fits_entry, "100")
        row = row + 1

        self.astrometrynet_label = tk.Label(
            self.settings_frame, text="Astrometry.net Server:")
        self.astrometrynet_label.grid(row=row, column=0, stick=tk.W)
        self.astrometrynet_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width)
        self.astrometrynet_entry.grid(row=row, column=1, sticky=tk.E)
        self.set_entry_text(self.astrometrynet_entry, "nova.astrometry.net")
        row = row + 1

        self.astrometrynet_key_label = tk.Label(
            self.settings_frame, text="Astrometry.net API Key:")
        self.astrometrynet_key_label.grid(row=row, column=0, stick=tk.W)
        self.astrometrynet_key_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width)
        self.astrometrynet_key_entry.grid(row=row, column=1, sticky=tk.E)
        self.astrometrynet_key_entry.config(show="*")
        self.set_entry_text(self.astrometrynet_key_entry, "pwjgdcpwaugkhkln")
        row = row + 1

        self.obscode_label = tk.Label(
            self.settings_frame, text="BAA Observer Code:")
        self.obscode_label.grid(row=row, column=0, stick=tk.W)
        self.obscode_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.obscode_entry.grid(row=row, column=1, sticky=tk.E)
        row = row + 1

        self.aavso_obscode_label = tk.Label(
            self.settings_frame, text="AAVSO Observer Code:")
        self.aavso_obscode_label.grid(row=row, column=0, stick=tk.W)
        self.aavso_obscode_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.aavso_obscode_entry.grid(row=row, column=1, sticky=tk.E)
        row = row + 1

        self.latitude_label = tk.Label(
            self.settings_frame, text="Observatory Latitude:")
        self.latitude_label.grid(row=row, column=0, stick=tk.W)
        self.latitude_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.latitude_entry.grid(row=row, column=1, sticky=tk.E)
        row = row + 1

        self.longitude_label = tk.Label(
            self.settings_frame, text="Observatory Longitude:")
        self.longitude_label.grid(row=row, column=0, stick=tk.W)
        self.longitude_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.longitude_entry.grid(row=row, column=1, sticky=tk.E)
        row = row + 1

        self.height_label = tk.Label(
            self.settings_frame, text="Observatory Height, m ASL:")
        self.height_label.grid(row=row, column=0, stick=tk.W)
        self.height_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.height_entry.grid(row=row, column=1, sticky=tk.E)
        row = row + 1

        self.telescope_label = tk.Label(self.settings_frame, text="Telescope:")
        self.telescope_label.grid(row=row, column=0, stick=tk.W)
        self.telescope_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.telescope_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1

        self.tbv_label = tk.Label(self.settings_frame, text="Tbv:")
        self.tbv_label.grid(row=row, column=0, stick=tk.W)
        self.tbv_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.tbv_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1

        self.tv_bv_label = tk.Label(self.settings_frame, text="Tv_bv:")
        self.tv_bv_label.grid(row=row, column=0, stick=tk.W)
        self.tv_bv_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.tv_bv_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1

        self.tb_bv_label = tk.Label(self.settings_frame, text="Tb_bv:")
        self.tb_bv_label.grid(row=row, column=0, stick=tk.W)
        self.tb_bv_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.tb_bv_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1

        self.accessory_label = tk.Label(self.settings_frame, text="Accessory:")
        self.accessory_label.grid(row=row, column=0, stick=tk.W)
        self.accessory_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.accessory_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1

        self.ccd_label = tk.Label(self.settings_frame, text="Camera:")
        self.ccd_label.grid(row=row, column=0, stick=tk.W)
        self.ccd_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.ccd_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1

        self.ccd_gain_label = tk.Label(
            self.settings_frame, text="Gain, e-/ADU:")
        self.ccd_gain_label.grid(row=row, column=0, stick=tk.W)
        self.ccd_gain_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width)
        self.ccd_gain_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1

        self.date_obs_label = tk.Label(
            self.settings_frame, text="Date-Obs (JD):")
        self.date_obs_label.grid(row=row, column=0, stick=tk.W)
        self.date_obs_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.date_obs_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1

        self.object_name_label = tk.Label(
            self.settings_frame, text="Object Name:")
        self.object_name_label.grid(row=row, column=0, stick=tk.W)
        self.object_name_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.object_name_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1

        self.object_notes_label = tk.Label(self.settings_frame, text="Notes:")
        self.object_notes_label.grid(row=row, column=0, stick=tk.W)
        self.object_notes_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.object_notes_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1

        self.object_kref_label = tk.Label(self.settings_frame, text="Use Check Star (AAVSO Label):")
        self.object_kref_label.grid(row=row, column=0, stick=tk.W)
        self.object_kref_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.object_kref_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1

        self.object_sel_comp_label = tk.Label(self.settings_frame, text="Select Comp Stars (AAVSO Label):")
        self.object_sel_comp_label.grid(row=row, column=0, stick=tk.W)
        self.object_sel_comp_entry = tk.Entry(
            self.settings_frame, width=extended_settings_entry_width, background='pink')
        self.object_sel_comp_entry.grid(row=row, column=1, sticky=tk.EW)
        row = row + 1


        # Here we have full-width settings dropdowns in the right frame
        self.weighting_label = tk.Label(
            self.right_frame, text="Ensemble Fit Weighting:")
        self.weighting_label.grid(row=3, column=0, sticky=tk.W)
        self.weighting_stringvar = tk.StringVar()
        self.weighting_stringvar.set("None")
        self.weighting_dropdown = tk.OptionMenu(
            self.right_frame, self.weighting_stringvar, "None", "Raw Flux", "Instrumental Magnitude", "PSF Sigma")
        self.weighting_dropdown.grid(row=4, column=0, sticky=tk.EW)

        self.catalog_label = tk.Label(
            self.right_frame, text="Comparison Catalog:")
        self.catalog_label.grid(row=5, column=0, sticky=tk.W)
        self.catalog_stringvar = tk.StringVar()
        self.catalog_stringvar.set("AAVSO")
        self.catalog_dropdown = tk.OptionMenu(
            self.right_frame, self.catalog_stringvar, "AAVSO", "APASS DR9", "URAT1", "USNO-B1.0", "Gaia DR2", "VizieR Catalog")
        self.catalog_dropdown.grid(row=6, column=0, sticky=tk.EW)

        self.vizier_catalog_label = tk.Label(
            self.right_frame, text="AAVSO ChartID or VizieR Catalog Number:")
        self.vizier_catalog_label.grid(row=7, column=0, sticky=tk.W)
        self.vizier_catalog_entry = tk.Entry(self.right_frame)
        self.vizier_catalog_entry.grid(row=8, column=0, sticky=tk.EW)

        self.fitter_label = tk.Label(self.right_frame, text="PSF Fitter:")
        self.fitter_label.grid(row=9, column=0, sticky=tk.W)
        self.fitter_stringvar = tk.StringVar()
        self.fitter_stringvar.set("Levenberg-Marquardt")
        self.fitter_dropdown = tk.OptionMenu(self.right_frame, self.fitter_stringvar,
                                             "Levenberg-Marquardt", "Linear Least Square", "Sequential LS Programming", "Simplex LS")
        self.fitter_dropdown.grid(row=10, column=0, sticky=tk.EW)

        self.remove_vsx_var = tk.BooleanVar()
        self.remove_vsx_checkbox = tk.Checkbutton(
            self.right_frame, text="Ignore VSX Sources in Ensemble", variable=self.remove_vsx_var)
        self.remove_vsx_checkbox.grid(row=11, column=0, sticky=tk.EW)
        self.remove_vsx_var.set(True)

        self.nearby_vsx_var = tk.BooleanVar()
        self.nearby_vsx_checkbox = tk.Checkbutton(
            self.right_frame, text="Report VSX Sources Nearby", variable=self.nearby_vsx_var)
        self.nearby_vsx_checkbox.grid(row=12, column=0, sticky=tk.EW)
        self.nearby_vsx_var.set(True)

        self.batch_psf_var = tk.BooleanVar()
        self.batch_psf_checkbox = tk.Checkbutton(
            self.right_frame, text="PSF Photometry Batch Processing", variable=self.batch_psf_var)
        self.batch_psf_checkbox.grid(row=13, column=0, sticky=tk.EW)
        self.batch_psf_var.set(True)

        self.plate_solve_on_open = tk.BooleanVar()
        self.plate_solve_on_open_checkbox = tk.Checkbutton(
            self.right_frame, text="Plate solve when loading FITS file", variable=self.plate_solve_on_open)
        self.plate_solve_on_open_checkbox.grid(row=14, column=0, sticky=tk.EW)
        self.plate_solve_on_open.set(True)

        # Histogram stretch sliders
        self.stretch_label = tk.Label(
            self.right_frame, text="Histogram Stretch Low/High:")
        self.stretch_label.grid(row=15, column=0, sticky=tk.W)
        self.stretch_low = tk.Scale(
            self.right_frame, from_=0, to=100, orient=tk.HORIZONTAL, command=self.update_histogram_low)
        self.stretch_low.grid(row=16, column=0, sticky=tk.EW)
        self.stretch_high = tk.Scale(
            self.right_frame, from_=0, to=100, orient=tk.HORIZONTAL, command=self.update_histogram_high)
        self.stretch_high.set(5)
        self.stretch_high.grid(row=17, column=0, sticky=tk.EW)

        self.stretching_label = tk.Label(
            self.right_frame, text="Image Stretching:")
        self.stretching_label.grid(row=18, column=0, sticky=tk.W)
        self.stretching_stringvar = tk.StringVar()
        self.stretching_stringvar.set("None")
        self.stretching_dropdown = tk.OptionMenu(
            self.right_frame, self.stretching_stringvar, "None", "Square Root", "Log", "Asinh")
        self.stretching_dropdown.grid(row=18, column=0, sticky=tk.EW)
        self.stretching_stringvar.trace(
            "w", lambda name, index, mode, sv=self.stretching_stringvar: self.update_display())

        # Console below
        self.console = tk.Text(self.center, height=40,
                               bg='black', fg='white', width=200)
        self.console.grid(sticky=tk.N+tk.S+tk.E+tk.W, column=0, row=3)
        self.console_scrollbar = tk.Scrollbar(self.center)
        self.console_scrollbar.grid(
            sticky=tk.N + tk.S + tk.E + tk.W, column=1, row=3)

        self.console.config(yscrollcommand=self.console_scrollbar.set)
        self.console_scrollbar.config(command=self.console.yview)

        self.console_msg(self.program_full_name)
        self.console_msg("Ready")

        # self.initialize_debug()
        # self.plot_photometry()

        tk.mainloop()


myGUI = MyGUI()

