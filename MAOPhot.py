4# -*- coding: utf-8 -*-
"""
 #     #    #    ####### ######                      
 ##   ##   # #   #     # #     # #    #  ####  ##### 
 # # # #  #   #  #     # #     # #    # #    #   #   
 #  #  # #     # #     # ######  ###### #    #   #   
 #     # ####### #     # #       #    # #    #   #   
 #     # #     # #     # #       #    # #    #   #   
 #     # #     # ####### #       #    #  ####    #   


   #        #####        ###   
  ##       #     #      #   #  
 # #             #     #     # 
   #        #####      #     # 
   #   ### #       ### #     # 
   #   ### #       ###  #   #  
 ##### ### ####### ###   ###   
                                                                                                                                                                                          
Welcome to MAOPhot 1.2.0, a PSF Photometry tool using Astropy and Photutils.psf

    1.2.0 Revision

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
 
Note: This module was developed with AI assistance (approx. 10% of logic).
All code has been human-verified for safety and performance.

This program was derived from MetroPSF by Maxym Usatov.  
It has been redesigned for AAVSO reporting only and includes, but not limited 
to the following enhancements:

- Generation of Effective PSF model, and ability to create a ‘rejection list’
- option to use an Integrated Gaussian PRF (Pixel Response Function) as model 
- PSF Photometry using an iterative algorithm to perform point spread function 
photometry in crowded fields
- Photometry using an ensemble of comparison stars or a single comp star
- Generation of Multi Color Photometry (B,V), (V,R), (V, I), etc., and Single Image
 Photometry reports in AAVSO extended format 
- Use of telescope Transformation Coefficients (needed for Multi Color Photometry)
- Image display shows comp star AAVSO label number and name of any found VSX 
objects in image field
- Intermediate results are saved as .csv files 
- User can optionally enter a AAVSO Chart ID when retrieving comparison star
 data
- User can specify check star and list of comp stars to use

    More about Single Image Photometry

        Single Image Photometry does not utilize the Transformation
        coefficients. Simple differential photometry is used.
        If more than 1 comp star then CMAG=ENSEMBLE and ensemble photometry 
        is performed
        
        Var mag = Var IM - Comp IM + Comp (known) mag
        Check star mag = Check star IM - Comp IM + Comp (known) mag

        where IM is the instrumental magnitude 
            -2.5 * np.log10(self.results_tab_df["flux_fit"] / exptime)
        self.results_tab_df["flux_fit"] represents the fitted flux for the 
        star (in that row)

        

    Multi Color Photometry 
        The process uses the AAVSO Reommended iterative method
        that matches AAVSO TransformApplier (TA v2.7.1). 
        The AAVSO recommended transforms that MAOPhot performs require the following coefficients:
        for  B     V     R     I
        BVRI Tb_bv Tv_bv Tr_vi Ti_vi
        BVR  Tb_bv Tv_bv Tvr
        BVI  Tb_bv Tv_bv       Ti_vi
        VRI        Tv_vi Tr_vi Tvi
        BV   Tb_bv Tv_bv
        VR         Tv_vr Tr_vr
        VI         Tv_vi       Tvi
        

        For more information see function color_photometry (dispatched from the
        seven menu callbacks BV_/VR_/VI_/BVR_/BVI_/VRI_/BVRI_multi_color_photometry).
        To calculate then generate a report, the correct sets of results_tab_df in csv
        format must exist, one for each color, and must have been derived from the 
        B, V, R, and/or I images of the Var under analysis.
        
        E.g., when the BV_generate_aavso_report callback is invoked by the menu item
        'Generate AAVSO Report->Multi Color Photometry->(B,V)', it calls
        generate_aavso_report_color('BV'), which reads the previously-saved
        Master-Report CSV (produced by the (B,V) entry under 'Multi Color Photometry'
        -> color_photometry('BV')).
        
        From these files/Panda databases, the formulas are calculated, and 
		results are displayed. 
        
        Error Estimation :
        MAOPhot mimics VPhot when calculating error estimation. 
        From VPhot documentation:
            In an ensemble solution with more than two comp stars, 
            the magnitude is estimated as the average of the individual 
            comp stars estimate [of the check star], and the error is taken as 
            the standard deviation of this sample. 

       


    MAOPhot derived from original MetroPSF (C) Copyright 2021, Maxym Usatov 
    <maxim.usatov@bcsatellite.net> Refer to metropsf.pdf for license information.

    This research made use of Photutils, an Astropy package for
    detection and photometry of astronomical sources (Bradley et al. 2024).

"""
#Tell user it's coming
print("MAOPhot is loading...please wait for GUI")

#
# Constants
#
__version__ = "1.2.0"
__label_prefix__ = "comp " # prepended to comp stars label's; forces type to str
__empty_cell__ = "%" #this forces cell to be type string
__our_padding__ = 10

from ast import Assert
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling.fitting import TRFLSQFitter, SLSQPLSQFitter, SimplexLSQFitter
from astropy.nddata import Cutout2D
from astropy.nddata import NDData
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.stats import SigmaClip
from astropy.table import Table, QTable
from astropy.time import Time, TimeDelta
from astropy.visualization import SqrtStretch, LogStretch, AsinhStretch, simple_norm
from astropy.wcs import WCS
from astroquery.astrometry_net import AstrometryNet
from astroquery.vizier import Vizier
from matplotlib import cm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.mplot3d import Axes3D
from photutils.background import Background2D
from photutils.background import LocalBackground
from photutils.background import MADStdBackgroundRMS
from photutils.background import MedianBackground
from photutils.background import MMMBackground
from photutils.detection import DAOStarFinder
from photutils.detection import find_peaks
from photutils.psf import CircularGaussianPRF, MoffatPSF, SourceGrouper
from photutils.psf import EPSFBuilder,EPSFFitter
from photutils.psf import extract_stars
from photutils.psf import IterativePSFPhotometry, PSFPhotometry
from PIL import Image, ImageTk, ImageMath
from time import gmtime, strftime
from tkinter import filedialog as fd
from tkinter import simpledialog
from tkinter import ttk
from tkinter.messagebox import askokcancel, askyesno
from tqdm import tqdm
import csv
import datetime
import io
import logging
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd
import pickle
import requests
import sys
import tkinter as tk
import warnings
import zipfile
import tempfile

from astropy.modeling import models as astropy_models
from astropy.modeling import fitting as astropy_fitting

#############################################################################
#
#  My Utilities
#
#############################################################################

#
# is_number returns True if s (usually a str) can be converted to a float
#
def is_number(s):
    try:
        float(s)  # Try converting to float
        return True
    except ValueError:
        return False


#############################################################################
#
#  PSF Parameter Estimation Functions
#
#  Based on estimate_psf_params.py - estimates FWHM and Moffat beta from
#  stars in an image by fitting 2D Moffat profiles.
#
#############################################################################

def extract_star_cutout(data, x, y, size=25):
    """Extract a square cutout centered on (x, y)."""
    half = size // 2
    y_int, x_int = int(round(y)), int(round(x))
    
    y_min = max(0, y_int - half)
    y_max = min(data.shape[0], y_int + half + 1)
    x_min = max(0, x_int - half)
    x_max = min(data.shape[1], x_int + half + 1)
    
    cutout = data[y_min:y_max, x_min:x_max].copy()
    
    # Only return if we got a full cutout
    if cutout.shape[0] == size and cutout.shape[1] == size:
        return cutout, x_int - x_min, y_int - y_min
    return None, None, None


def fit_moffat_2d(cutout, x_center, y_center):
    """
    Fit a 2D Moffat function to a background-subtracted star cutout.
    
    Parameters
    ----------
    cutout : 2D array
        Background-subtracted star cutout
    x_center : float
        Initial x position of star center in cutout coordinates
    y_center : float
        Initial y position of star center in cutout coordinates
    
    Returns
    -------
    dict with fitted parameters or None if fit fails
    """
    size = cutout.shape[0]
    y, x = np.mgrid[0:size, 0:size]
    
    # Data is already background-subtracted
    data_sub = cutout
    
    # Initial estimates
    amplitude = np.max(data_sub)
    if amplitude <= 0:
        return None
    
    # Initial FWHM estimate from second moments
    total = np.sum(data_sub[data_sub > 0])
    if total <= 0:
        return None
    
    # Estimate gamma (related to FWHM) - start with ~2 pixels
    gamma_init = 2.0
    alpha_init = 2.5  # Typical Moffat beta
    
    # Create Moffat2D model
    moffat_init = astropy_models.Moffat2D(
        amplitude=amplitude,
        x_0=x_center,
        y_0=y_center,
        gamma=gamma_init,
        alpha=alpha_init
    )
    
    # Set bounds to keep parameters physical
    moffat_init.amplitude.bounds = (0, amplitude * 2)
    moffat_init.x_0.bounds = (x_center - 3, x_center + 3)
    moffat_init.y_0.bounds = (y_center - 3, y_center + 3)
    moffat_init.gamma.bounds = (0.5, 10.0)
    moffat_init.alpha.bounds = (1.0, 10.0)
    
    # Fit
    fitter = astropy_fitting.LevMarLSQFitter()
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            moffat_fit = fitter(moffat_init, x, y, data_sub, maxiter=500)
        except Exception:
            return None
    
    # Extract parameters
    gamma = moffat_fit.gamma.value
    alpha = moffat_fit.alpha.value  # This is beta in the Moffat literature
    amplitude_fit = moffat_fit.amplitude.value
    
    # Calculate FWHM from gamma and alpha (beta)
    # FWHM = 2 * gamma * sqrt(2^(1/alpha) - 1)
    fwhm = 2.0 * gamma * np.sqrt(2.0**(1.0/alpha) - 1.0)
    
    # Calculate fit quality (normalized residual)
    model_image = moffat_fit(x, y)
    residual = data_sub - model_image
    
    # Use pixels within FWHM for residual calculation
    mask = ((x - moffat_fit.x_0.value)**2 + (y - moffat_fit.y_0.value)**2) < (fwhm * 1.5)**2
    if np.sum(mask) > 0:
        fit_quality = np.sqrt(np.mean(residual[mask]**2)) / amplitude_fit
    else:
        fit_quality = np.sqrt(np.mean(residual**2)) / amplitude_fit
    
    # Sanity checks
    if fwhm < 0.5 or fwhm > 20:
        return None
    if alpha < 1.0 or alpha > 10:
        return None
    if fit_quality > 0.5:  # Very poor fit
        return None
    
    return {
        'fwhm': fwhm,
        'beta': alpha,
        'gamma': gamma,
        'amplitude': amplitude_fit,
        'x_fit': moffat_fit.x_0.value,
        'y_fit': moffat_fit.y_0.value,
        'fit_quality': fit_quality,
        'model': moffat_fit  # Store the fitted model for plotting
    }



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
settings_filename = ""
image_data = []
       


class MyGUI:

    zoom_level = 1
    linreg_error = 0
    photometry_results_plotted = False
    ePSF_samples_plotted = False
    results_tab_df = pd.DataFrame()
    image_bkg_value = 0
    fit_shape = 5
    error_raised = False
    histogram_slider_low = 0
    histogram_slider_high = 5
    last_clicked_x = 0
    last_clicked_y = 0
    ensemble_size = 0
    jd = 0
    image_file = ""
    inner_fits_name = None
    image_id = None
    settings_filename = ""
    photometry_circles = {}
    valid_parameter_list = {} # contains user's settings
    valid_config_list = {} # contains 'global' MAOPhot settings/config parameters
    ePSF_rejection_list = pd.DataFrame({'x':[],'y':[],"stale":[]})
    ePSF_pending_rejection_list = pd.DataFrame({'x':[],'y':[], "stale":[]})
    epsf_model = None
    stars_tbl = None
    isolated_stars_tbl = None
    fits_header_filter = ""
    is_inverted = False # Track if image is inverted
    current_oversampling = 2

    # Parameter declaration  and init 
    find_peaks_npeaks_entry = None
    fit_width_entry = None
    max_ensemble_magnitude_entry = None
    fwhm_entry = None
    oversampling_entry = None
    star_detection_threshold_factor_entry = None
    photometry_iterations_entry = None
    sharplo_entry = None
    matching_radius_entry = None
    aavso_obscode_entry = None
    telescope_entry = None
    tb_bv_entry = None
    tb_br_entry = None
    tb_bi_entry = None
    tv_bv_entry = None
    tv_vr_entry = None
    tr_vr_entry = None
    tr_ri_entry = None
    ti_ri_entry = None
    tv_vi_entry = None
    ti_vi_entry = None
    tr_vi_entry = None
    tbv_entry = None
    tbr_entry = None
    tbi_entry = None
    tvr_entry = None
    tri_entry = None
    tvi_entry = None
    tb_bv_err_entry = None
    tb_br_err_entry = None
    tb_bi_err_entry = None
    tv_bv_err_entry = None
    tv_vr_err_entry = None
    tr_vr_err_entry = None
    tr_ri_err_entry = None
    ti_ri_err_entry = None
    tv_vi_err_entry = None
    ti_vi_err_entry = None
    tr_vi_err_entry = None
    tbv_err_entry = None
    tbr_err_entry = None
    tbi_err_entry = None
    tvr_err_entry = None
    tri_err_entry = None
    tvi_err_entry = None
    linearity_limit_entry = None
    catalog_stringvar = None
    vizier_catalog_entry = None
    fitter_stringvar = None
    astrometrynet_entry = None
    astrometrynet_key_entry = None
    object_kref_entry = None
    object_sel_comp_entry = None
    object_min_snr_entry = None
    object_name_entry = None
    object_name_alpha_entry = None
    object_name_delta_entry = None
    object_notes_entry = None
    display_all_objects = None
    use_gaussian_prf_model = None
    generate_residual_image = None
    candidate_stars = None
    settings_filename_entry = None
    user_note_entry = None
    auto_behavior = None
    filter_entry = None
    max_qfit_entry = None
    moffat_beta_entry = None
    min_separation_bias_entry = None
    extinction_B_entry = None
    extinction_V_entry = None
    extinction_I_entry = None
    extinction_R_entry = None
    extinction_C_entry = None

    # The TopLoevel window containing the settings
    es_top = None

    # Some FIT files have DATE-END 
    date_end_jd = None

#######################################################################################
#
# console_msg
#
#######################################################################################

    def console_msg(self, MAOPhot_message, level=logging.INFO):
        # add a time stamp
        message = datetime.datetime.now().strftime("%d %b %Y %H:%M:%S")
        message += "      " + MAOPhot_message
        self.console.insert(tk.END, message+"\n")
        self.console.see(tk.END)
        self.window.update_idletasks()
        self.our_logger.log(level=level, msg=MAOPhot_message)

#######################################################################################
#
# check_image_parameters
#
# Validates and returns image analysis parameters from settings entries.
# Returns a dictionary of validated parameters or None if critical parameters invalid.
#
#######################################################################################

    def check_image_parameters(self):
        """
        Validates and returns image analysis parameters from settings entries.
        
        Returns
        -------
        dict or None
            Dictionary containing validated parameters:
            - fwhm: float, FWHM estimate in pixels
            - threshold_factor: int, star detection threshold factor
            - linearity_limit: int, saturation/linearity limit
            - n_peaks: int or None, maximum number of peaks
            Returns None if critical parameters are invalid.
        """
        params = {}
        
        # Test fwhm
        if self.fwhm_entry is not None and is_number(self.fwhm_entry.get().strip()):
            params['fwhm'] = float(self.fwhm_entry.get())
        else:
            self.console_msg("FWHM not numeric, using 3.0")
            params['fwhm'] = 3.0

        # Test oversampling
        if self.oversampling_entry is not None and is_number(self.oversampling_entry.get().strip()):
            params['oversampling'] = int(abs(float(self.oversampling_entry.get())))
        else:
            self.console_msg("oversampling not numeric, using 2, (most common)")
            params['oversampling'] = 2

        # Test star_detection_threshold_factor
        if self.star_detection_threshold_factor_entry is not None:
            thresh_val = self.star_detection_threshold_factor_entry.get().strip()
            if thresh_val.isnumeric():
                params['threshold_factor'] = int(thresh_val)
            else:
                self.console_msg("Detection threshold factor not numeric, using 10")
                params['threshold_factor'] = 10
        else:
            params['threshold_factor'] = 10

        # Test linearity_limit_entry 
        if self.linearity_limit_entry is not None:
            linearity_val = self.linearity_limit_entry.get().strip()
            if linearity_val and linearity_val.isnumeric():
                params['linearity_limit'] = int(linearity_val)
            else:
                self.console_msg("Linearity limit not valid, using 60000")
                params['linearity_limit'] = 60000
        else:
            params['linearity_limit'] = 60000

        # Test find_peaks_npeaks_entry (can be None for unlimited)
        if self.find_peaks_npeaks_entry is not None:
            npeaks_val = self.find_peaks_npeaks_entry.get().strip()
            if npeaks_val and npeaks_val.isnumeric():
                params['n_peaks'] = int(npeaks_val)
            else:
                params['n_peaks'] = None  # Unlimited
        else:
            params['n_peaks'] = None
            
        return params

#######################################################################################
#
# estimate_and_update_psf_params
#
# Estimates FWHM and Moffat beta from stars in the loaded image and updates
# the corresponding settings entries.
#
#######################################################################################

    def estimate_and_update_psf_params(self):
        """
        Estimates FWHM and Moffat beta from stars in the loaded image.
        Updates self.fwhm_entry and self.moffat_beta_entry with the results.
        Also generates a diagnostic plot.
        """
        global image_data
        
        try:
            self.console_msg("PSF Estimation: starting...")
            
            # Check that image is loaded
            if len(image_data) == 0:
                self.console_msg("PSF Estimation: no image loaded, skipping")
                return
            
            # Get validated parameters
            params = self.check_image_parameters()
            if params is None:
                self.console_msg("PSF Estimation: invalid parameters, skipping")
                return
            
            fwhm_estimate = params['fwhm']
            threshold_factor = params['threshold_factor']
            saturation_limit = params['linearity_limit']
            n_stars = params['n_peaks'] if params['n_peaks'] is not None else 30
            
            # Limit n_stars for PSF estimation (don't need hundreds)
            n_stars = min(n_stars, 100)
            if n_stars < 10:
                n_stars = 30
            
            cutout_size = 25  # Hardcoded as specified
            
            self.console_msg(f"PSF Estimation: using n_stars={n_stars}, fwhm_estimate={fwhm_estimate:.1f}, threshold_factor={threshold_factor}, saturation={saturation_limit}")
            
            # Estimate background
            mean, median, std = sigma_clipped_stats(image_data, sigma=2.0)
            threshold = threshold_factor * std
            
            self.console_msg(f"PSF Estimation: background median={median:.1f}, std={std:.1f}, threshold={threshold:.1f}")
            
            # Background-subtracted data
            data_sub = image_data - median
            
            # Detect stars using DAOStarFinder
            daofind = DAOStarFinder(
                fwhm=fwhm_estimate,
                threshold=threshold,
                sharplo=0.2,
                sharphi=1.0,
                roundlo=-1.0,
                roundhi=1.0,
                peakmax=float(saturation_limit) * 0.9,  # Stay below saturation
                exclude_border=True,
                min_separation=fwhm_estimate
            )
            
            sources = daofind(data_sub)
            
            if sources is None or len(sources) == 0:
                self.console_msg("PSF Estimation: no stars detected, skipping")
                return
            
            self.console_msg(f"PSF Estimation: detected {len(sources)} sources")
            
            # Compute minimum separation to avoid overlapping cutouts
            min_separation = cutout_size + fwhm_estimate
            
            # Sort by flux (descending) and select isolated stars
            sources.sort('flux', reverse=True)
            
            selected_sources = []
            for row in sources:
                x, y = row['xcentroid'], row['ycentroid']
                
                # Check isolation
                is_isolated = True
                for selected in selected_sources:
                    sx, sy = selected['xcentroid'], selected['ycentroid']
                    dist = np.sqrt((x - sx)**2 + (y - sy)**2)
                    if dist < min_separation:
                        is_isolated = False
                        break
                
                if is_isolated:
                    selected_sources.append(row)
                    if len(selected_sources) >= n_stars * 2:  # Get extra in case some fits fail
                        break
            
            self.console_msg(f"PSF Estimation: selected {len(selected_sources)} isolated stars for fitting")
            
            # Fit Moffat to each star
            results = []
            for row in selected_sources:
                x, y = row['xcentroid'], row['ycentroid']
                
                # Extract cutout from background-subtracted data
                cutout, x_local, y_local = extract_star_cutout(data_sub, x, y, size=cutout_size)
                if cutout is None:
                    continue
                
                fit_result = fit_moffat_2d(cutout, x_local, y_local)
                if fit_result is not None:
                    fit_result['x_image'] = x
                    fit_result['y_image'] = y
                    results.append(fit_result)
                
                if len(results) >= n_stars:
                    break
            
            if len(results) < 3:
                self.console_msg(f"PSF Estimation: only {len(results)} stars fit successfully, need at least 3")
                return
            
            self.console_msg(f"PSF Estimation: successfully fit {len(results)} stars")
            
            # Calculate statistics using sigma clipping
            fwhms = np.array([r['fwhm'] for r in results])
            betas = np.array([r['beta'] for r in results])
            qualities = np.array([r['fit_quality'] for r in results])

            # First pass: sigma clip to remove extreme outliers
            fwhm_clipped, fwhm_low, fwhm_high = sigma_clip(fwhms, sigma=2.5, masked=False, return_bounds=True)
            beta_clipped, beta_low, beta_high = sigma_clip(betas, sigma=2.5, masked=False, return_bounds=True)
            
            # Create mask for values that survived clipping
            fwhm_mask = (fwhms >= fwhm_low) & (fwhms <= fwhm_high)
            beta_mask = (betas >= beta_low) & (betas <= beta_high)
            combined_mask = fwhm_mask & beta_mask
            
            # Also exclude fits that hit parameter bounds (beta=1 is suspicious)
            bound_mask = betas > 1.05  # Exclude fits that landed near lower bound
            combined_mask = combined_mask & bound_mask
            
            # Apply mask
            fwhms_good = fwhms[combined_mask]
            betas_good = betas[combined_mask]
            qualities_good = qualities[combined_mask]
            
            if len(fwhms_good) < 5:
                # Fall back to simple sigma-clipped stats if too few good fits
                self.console_msg("Warning: Few high-quality fits (<5), using sigma-clipped stats")
                fwhm_mean, fwhm_median, fwhm_std = sigma_clipped_stats(fwhms, sigma=2.5)
                beta_mean, beta_median, beta_std = sigma_clipped_stats(betas, sigma=2.5)
            else:
                # Convert quality to weights (lower residual = higher weight)
                # Use inverse quality, normalized
                weights = 1.0 / (qualities_good + 0.001)  # Add small value to avoid div by zero
                weights = weights / np.sum(weights)  # Normalize to sum to 1
                
                # Weighted statistics
                fwhm_mean = np.average(fwhms_good, weights=weights)
                beta_mean = np.average(betas_good, weights=weights)
                
                # For weighted median, use weighted percentile approach
                def weighted_median(values, weights):
                    """Calculate weighted median."""
                    sorted_indices = np.argsort(values)
                    sorted_values = values[sorted_indices]
                    sorted_weights = weights[sorted_indices]
                    cumsum = np.cumsum(sorted_weights)
                    median_idx = np.searchsorted(cumsum, 0.5)
                    return sorted_values[min(median_idx, len(sorted_values)-1)]
                
                fwhm_median = weighted_median(fwhms_good, weights)
                beta_median = weighted_median(betas_good, weights)
                
                # Weighted standard deviation
                fwhm_std = np.sqrt(np.average((fwhms_good - fwhm_mean)**2, weights=weights))
                beta_std = np.sqrt(np.average((betas_good - beta_mean)**2, weights=weights))
            

            self.console_msg(f"PSF Estimation Results:")
            self.console_msg(f"  FWHM: median={fwhm_median:.3f}, mean={fwhm_mean:.3f}, std={fwhm_std:.3f} pixels")
            self.console_msg(f"  Beta: median={beta_median:.3f}, mean={beta_mean:.3f}, std={beta_std:.3f}")
            self.console_msg(f"  Median fit quality: {np.median(qualities):.4f}")
            
            # Update the settings entries
            if self.fwhm_entry is not None:
                self.set_entry_text(self.fwhm_entry, f"{fwhm_median:.2f}")
                self.console_msg(f"PSF Estimation: updated FWHM entry to {fwhm_median:.2f}")
            
            #Keep beta_median, resonable
            beta_median = max(beta_median, 1.5)

            if self.moffat_beta_entry is not None:
                self.set_entry_text(self.moffat_beta_entry, f"{beta_median:.2f}")
                self.console_msg(f"PSF Estimation: updated Moffat \u03B2 entry to {beta_median:.2f}")
            
            raw_fwhm_median = np.median(fwhms)
            raw_beta_median = np.median(betas)

            # Generate diagnostic plot (pass data_sub for example star profile)
            self.plot_psf_estimation_results(results, fwhm_median, beta_median,
                                            fwhm_std, beta_std, median_quality=None, data_sub=data_sub,
                                            raw_fwhm_median=raw_fwhm_median, raw_beta_median=raw_beta_median)

            self.console_msg("PSF Estimation: complete")
            self.console_msg("Ready")
            
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg(f"PSF Estimation Exception at line {exc_tb.tb_lineno}: {str(e)}", level=logging.ERROR)

#######################################################################################
#
# plot_psf_estimation_results
#
# Generates diagnostic plots for PSF estimation results
#
#######################################################################################

    def plot_psf_estimation_results(self, results, fwhm_median, beta_median, 
                                     fwhm_std, beta_std, median_quality, data_sub,
                                     raw_fwhm_median=None, raw_beta_median=None):
        """
        Generate diagnostic plots for PSF estimation results.
        
        Parameters
        ----------
        results : list of dict
            List of fit results for each star
        fwhm_median : float
            Median FWHM value (quality-weighted if weighting was applied)
        beta_median : float
            Median Moffat beta value (quality-weighted if weighting was applied)
        fwhm_std : float
            Standard deviation of FWHM
        beta_std : float
            Standard deviation of beta
        median_quality : float
            Median fit quality
        data_sub : 2D array
            Background-subtracted image data for extracting example cutout
        raw_fwhm_median : float, optional
            Unweighted median FWHM (for comparison, if weighting was applied)
        raw_beta_median : float, optional
            Unweighted median beta (for comparison, if weighting was applied)
        """
        try:
            import matplotlib.pyplot as plt
            
            # Create a new figure with a specific number to avoid interfering 
            # with MAOPhot's embedded canvases
            fig = plt.figure(num='PSF Estimation Results', figsize=(14, 9), clear=True)
            axes = fig.subplots(2, 3)
            fig.suptitle('PSF Parameter Estimation Results', fontsize=12)
            
            fwhms = np.array([r['fwhm'] for r in results])
            betas = np.array([r['beta'] for r in results])
            qualities = np.array([r['fit_quality'] for r in results])
            
            # Identify boundary hits (beta near lower bound of 1.0)
            boundary_mask = betas < 1.05
            n_boundary = np.sum(boundary_mask)
            
            # 1. FWHM histogram
            ax = axes[0, 0]
            ax.hist(fwhms, bins=15, edgecolor='black', alpha=0.7)
            ax.axvline(fwhm_median, color='red', linestyle='--', linewidth=2,
                       label=f"Median: {fwhm_median:.2f}")
            if raw_fwhm_median is not None and abs(raw_fwhm_median - fwhm_median) > 0.05:
                ax.axvline(raw_fwhm_median, color='orange', linestyle=':', linewidth=2,
                           label=f"Raw median: {raw_fwhm_median:.2f}")
            ax.set_xlabel('FWHM (pixels)')
            ax.set_ylabel('Count')
            ax.set_title('FWHM Distribution')
            ax.legend()
            
            # 2. Beta histogram - enhanced to show boundary hits
            ax = axes[0, 1]
            
            # Create histogram with custom coloring for boundary hits
            bin_edges = np.linspace(max(0.5, betas.min() - 0.2), min(8, betas.max() + 0.2), 16)
            n, bins, patches = ax.hist(betas, bins=bin_edges, edgecolor='black', alpha=0.7)
            
            # Color bars red if they contain boundary hits (beta < 1.05)
            for i, patch in enumerate(patches):
                bin_left = bins[i]
                bin_right = bins[i + 1]
                if bin_right <= 1.05:  # Bin entirely in boundary region
                    patch.set_facecolor('red')
                elif bin_left < 1.05:  # Bin partially in boundary region
                    patch.set_facecolor('salmon')
                else:
                    patch.set_facecolor('steelblue')
            
            ax.axvline(beta_median, color='green', linestyle='-', linewidth=2,
                       label=f"Weighted median: {beta_median:.2f}")
            if raw_beta_median is not None and abs(raw_beta_median - beta_median) > 0.1:
                ax.axvline(raw_beta_median, color='red', linestyle=':', linewidth=2,
                           label=f"Raw median: {raw_beta_median:.2f}")
            
            ax.set_xlabel('Moffat β')
            ax.set_ylabel('Count')
            ax.set_title('β Distribution')
            ax.legend(loc='upper right', fontsize=9)
            
            # Add warning text if boundary hits detected
            if n_boundary > 0:
                ax.text(0.05, 0.95, f'{n_boundary} fits hit bound (β≈1)',
                        transform=ax.transAxes, fontsize=9,
                        verticalalignment='top', color='red',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
            
            # 3. FWHM vs Beta scatter - color by fit quality
            ax = axes[0, 2]
            
            # Color points: red for boundary hits, otherwise by quality
            colors = np.where(boundary_mask, 'red', 'steelblue')
            scatter = ax.scatter(fwhms[~boundary_mask], betas[~boundary_mask], 
                                c=qualities[~boundary_mask], cmap='viridis_r',
                                alpha=0.7, edgecolors='black', linewidths=0.5,
                                label='Good fits')
            if n_boundary > 0:
                ax.scatter(fwhms[boundary_mask], betas[boundary_mask], 
                          c='red', marker='x', s=60, linewidths=2,
                          label=f'Boundary hits ({n_boundary})')
            
            ax.axhline(beta_median, color='green', linestyle='-', alpha=0.7, linewidth=1.5)
            ax.axvline(fwhm_median, color='green', linestyle='-', alpha=0.7, linewidth=1.5)
            ax.set_xlabel('FWHM (pixels)')
            ax.set_ylabel('Moffat β')
            ax.set_title('FWHM vs β (color = fit quality)')
            ax.legend(loc='upper right', fontsize=8)
            
            # Add colorbar for quality
            cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
            cbar.set_label('Fit quality', fontsize=8)
            
            # 4. Beta vs Fit quality scatter - key diagnostic
            ax = axes[1, 0]
            
            colors = np.where(boundary_mask, 'red', 'steelblue')
            ax.scatter(betas[~boundary_mask], qualities[~boundary_mask], 
                      c='steelblue', alpha=0.6, edgecolors='black', linewidths=0.5,
                      label='Used in statistics')
            if n_boundary > 0:
                ax.scatter(betas[boundary_mask], qualities[boundary_mask], 
                          c='red', marker='x', s=60, linewidths=2,
                          label=f'Excluded ({n_boundary})')
            
            ax.axvline(1.05, color='red', linestyle=':', alpha=0.5, 
                       label='Boundary threshold')
            ax.set_xlabel('Moffat β')
            ax.set_ylabel('Fit Quality (lower is better)')
            ax.set_title('β vs Fit Quality')
            ax.legend(loc='upper right', fontsize=8)
            
            # 5. Example star profile - choose from GOOD fits only
            ax = axes[1, 1]
            
            # Find the star with FWHM closest to median, excluding boundary hits
            good_results = [r for r in results if r['beta'] >= 1.05]
            if good_results:
                closest_idx = np.argmin([abs(r['fwhm'] - fwhm_median) for r in good_results])
                example = good_results[closest_idx]
            else:
                # Fall back to closest overall if no good fits
                closest_idx = np.argmin([abs(r['fwhm'] - fwhm_median) for r in results])
                example = results[closest_idx]
            
            # Extract cutout for this star from background-subtracted data
            cutout_size = 21
            cutout, x_local, y_local = extract_star_cutout(
                data_sub, example['x_image'], example['y_image'], size=cutout_size)
            
            if cutout is not None:
                size = cutout.shape[0]
                center = size // 2
                
                # Extract profile along the central row of the cutout
                profile_data = cutout[center, :]
                
                # X positions relative to cutout center
                x_plot = np.arange(size) - center
                ax.plot(x_plot, profile_data, 'ko', markersize=4, label='Data')
                
                # Use the stored model to evaluate along the central row
                # The model expects coordinates in cutout space
                x_fine = np.linspace(-center, center, 100)
                model = example['model']
                # Evaluate model at (x + x_fit, y_fit) to get profile through star center
                y_model = model(x_fine + example['x_fit'], 
                               np.full_like(x_fine, example['y_fit']))
                ax.plot(x_fine, y_model, 'r-', linewidth=2, label='Moffat fit')
                
                ax.set_xlabel('Offset from center (pixels)')
                ax.set_ylabel('Counts (background subtracted)')
                ax.set_title(f'Example Star Profile\nFWHM={example["fwhm"]:.2f}, β={example["beta"]:.2f}, qual={example["fit_quality"]:.4f}')
                ax.legend()
                ax.set_xlim(-8, 8)
            else:
                ax.text(0.5, 0.5, 'Could not extract example star', 
                        transform=ax.transAxes, ha='center', va='center')
                ax.set_title('Example Star Profile')
            
            # 6. Star positions on image - mark boundary hits differently
            ax = axes[1, 2]
            # Show image with stars marked
            norm = simple_norm(data_sub, 'sqrt', percent=99)
            ax.imshow(data_sub, origin='lower', cmap='gray', norm=norm)
            
            xs = np.array([r['x_image'] for r in results])
            ys = np.array([r['y_image'] for r in results])
            
            # Plot good fits in cyan, boundary hits in red
            ax.scatter(xs[~boundary_mask], ys[~boundary_mask], s=50, 
                      facecolors='none', edgecolors='cyan', linewidths=1.5,
                      label=f'Good fits ({np.sum(~boundary_mask)})')
            if n_boundary > 0:
                ax.scatter(xs[boundary_mask], ys[boundary_mask], s=50, 
                          facecolors='none', edgecolors='red', linewidths=2,
                          marker='s', label=f'Boundary hits ({n_boundary})')
            
            ax.set_title(f'Stars Used ({len(results)} total)')
            ax.set_xlabel('X (pixels)')
            ax.set_ylabel('Y (pixels)')
            ax.legend(loc='upper right', fontsize=8)
            
            plt.tight_layout()
            
            # Save the figure
            if self.image_file:
                plot_file = self.image_file + "_psf_estimation.png"
            else:
                plot_file = "psf_estimation.png"
            
            plt.savefig(plot_file, dpi=150, bbox_inches='tight')
            self.console_msg(f"PSF Estimation: diagnostic plot saved to {plot_file}")
            
            # Show only this specific figure
            fig.show()
            
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg(f"PSF Estimation Plot Exception at line {exc_tb.tb_lineno}: {str(e)}", level=logging.ERROR)


#######################################################################################
#
# display_image
#
#######################################################################################

    def display_image(self, verbose=False):
        if len(image_data) > 0:
            self.canvas.delete("all")
            global generated_image
            generate_FITS_thumbnail(self.histogram_slider_low,
                                    self.histogram_slider_high,
                                    self.zoom_level,
                                    self.stretching_stringvar.get())
            
            # Apply inversion if enabled
            if self.is_inverted:
                generated_image = ImageMath.eval("255 - a", a=generated_image)
                
            self.image = ImageTk.PhotoImage(generated_image)
            self.image_id = self.canvas.create_image(0, 0, anchor=tk.NW, image=self.image)
            self.canvas.config(scrollregion=self.canvas.bbox(tk.ALL))
            self.canvas.bind("<Button-1>", self.mouse_main_canvas_click)

            #
            # For dragging
            #
            self.canvas.tag_bind(self.image_id, "<Shift-Button-1>", self.on_drag_start)
            self.canvas.tag_bind(self.image_id, "<Shift-B1-Motion>", self.on_drag_move)

            #
            # For zooming & scrolling
            #
            self.canvas.bind("<MouseWheel>", self.on_canvas_mousewheel)
            self.canvas.bind("<Shift-MouseWheel>", self.on_canvas_shift_mousewheel)

            #
            # For centering the image
            #
            self.canvas.bind("<Button-3>", self.on_button_3_click)

            if self.ePSF_samples_plotted:
                self.display_ePSF_samples()
            self.plot_photometry(verbose=verbose)

#######################################################################################
#
# load_FITS
#
# 
#
#######################################################################################

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
                # check if this is a zip file
                # We need to get the inner_fits_name of the zip file
                # and use it later in case use want to execute a plain save.
                if image_file.endswith('.zip'):
                    with zipfile.ZipFile(image_file, 'r') as zf:
                        self.inner_fits_name = zf.namelist()[0]
                else:
                    self.inner_fits_name = None

                # Detect compressed FITS: primary HDU has no data, 
                # and there's a CompImageHDU in extension 1
                if image[0].data is None and len(image) > 1 and \
                    isinstance(image[1], fits.CompImageHDU):
                    self.console_msg("Detected CFITSIO compressed FITS image, using extension 1")
                    hdu_index = 1
                else:
                    hdu_index = 0

                header = image[hdu_index].header
                image_data = image[hdu_index].data
                self.image_file = image_file
                self.filename_label['text'] = "FITS: " + image_file
                self.canvas.delete("all")
                self.zoom_level = 1
                self.photometry_results_plotted = False
                self.ePSF_samples_plotted = False
                #Load previous work if it exists
                if os.path.isfile(self.image_file+".csv"):
                    self.results_tab_df = pd.read_csv(self.image_file + ".csv")
                else:
                    self.results_tab_df = pd.DataFrame()
                image_width = image_data.shape[1]
                image_height = image_data.shape[0]
                self.wcs_header = WCS(image[hdu_index].header)

                if not self.wcs_header.has_celestial:
                    self.console_msg(
                        "Note, Celestial coordinate information not in header.")

                FITS_minimum = np.min(image_data)
                FITS_maximum = np.max(image_data)
                self.console_msg("Width: " + str(image_width) +
                                 " Height: " + str(image_height))
                self.console_msg(
                    "FITS Minimum: " + str(FITS_minimum) + " Maximum: " +
                    str(FITS_maximum))
                if 'filter' in header:
                    self.fits_header_filter = str(header['filter'])
                    # User may want to use his own filter or filter not in Fits file
                    if self.auto_behavior.get():
                        self.filter = self.fits_header_filter
                        self.filter_entry.config(state='normal')
                        self.set_entry_text(self.filter_entry, self.filter)
                        self.filter_entry.config(state='disable')
                        self.console_msg("Filter: " + self.filter)
                    else:
                        self.filter_entry.config(state='normal')
                        self.filter = self.filter_entry.get().strip()
                        self.console_msg("Filter set by user: " + self.filter)
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
                        self.date_obs_entry.config(state='normal')
                        self.set_entry_text(self.date_obs_entry, str(self.jd))
                        self.date_obs_entry.config(state='readonly')
                        self.date_utc_entry.config(state='normal')
                        self.set_entry_text(self.date_utc_entry, str(Time(self.jd, format='jd', scale='utc').iso))
                        self.date_utc_entry.config(state='readonly')

                        
                    except Exception as e:
                        self.error_raised = True
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)

                if 'date-end' in header:
                    try:
                        date_end = Time(header['date-end'])
                        _date_end_jd = Time(date_end, format='jd')
                        self.date_end_jd = str(_date_end_jd)
                        self.console_msg("DATE-END: " + str(_date_end_jd.to_value('iso')) + " UTC; JD: " + self.date_end_jd)
                        
                    except Exception as e:
                        self.error_raised = True
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)
                else:
                    self.date_end_jd = None
                

                if 'jd' in header:
                    jd = header['jd']
                    self.console_msg(
                        "Julian date at the start of exposure (from JD): " + str(jd))
                    self.jd = jd
                    self.date_obs_entry.config(state='normal')
                    self.date_obs_entry.delete(0, tk.END)
                    self.date_obs_entry.insert(0, str(self.jd))
                    self.date_obs_entry.config(state='readonly')
                    self.date_utc_entry.config(state='normal')
                    self.date_utc_entry.delete(0, tk.END)
                    self.date_utc_entry.insert(0, str(Time(self.jd, format='jd', scale='utc').iso))
                    self.date_utc_entry.config(state='readonly')

                self.image_bkg_value = np.median(image_data)
                self.console_msg("Median background level, ADU: " + str(round(self.image_bkg_value, 2)))
                self.console_msg("Ready")

                
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) 
                +" "+str(e), level=logging.ERROR)

    ###############################################################
    #
    #
    #  open_FITS_file
    # 
    #
    ###############################################################
    def open_FITS_file(self):
        global header

        options = {}
        options['defaultextension'] = '.fit'
        options['title'] = 'Open FIT image file (compressed or not)...'

        fits_exts = '.fit .fits .fts'
        gz_exts = '.fit.gz .fits.gz .fts.gz'

        options['filetypes'] = [
            ('All FITS', f'{fits_exts} {gz_exts}'),
            ('FITS (uncompressed)', fits_exts),
            ('FITS (gzipped)', gz_exts),
            ('ZIP archives', '.zip'),
            ('All files', '*'),
        ]

        image_file = fd.askopenfilename(**options)
        
        if len(image_file) > 0:
            try:
                self.load_FITS(image_file)
                self.clear_ePSF()
                
                result = askyesno(title="Option to Estimate PSF Parameters", message="Do you want to estimate FWHM and Moffat \u03B2?")

                if result==True:
                    # Estimate PSF parameters (FWHM and Moffat beta) from the image
                    # and update the corresponding settings entries
                    self.estimate_and_update_psf_params()
                else:
                    self.console_msg("User not estimating and updating PSF parameters")
                    self.console_msg("Ready")
                
            except Exception as e:
                self.error_raised = True
                exc_type, exc_obj, exc_tb = sys.exc_info()
                self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
            

    ###############################################################
    #
    #
    #  save_FITS_file_as
    # 
    #  If a .gz or .zip file extension has been chosen, then the uncompressed
    #  file is compressed accordingly before saving.
    #
    ###############################################################
    def save_FITS_file_as(self):
        global image
        global image_data
        global header

        file_name = self.image_file

        options = {}
        options['title'] = 'Save FIT image file as...'

        # Build save dialog with compression options
        base, ext = os.path.splitext(file_name)  # .gz already stripped

        fits_exts = '.fit .fits .fts'

        options = {
            'defaultextension': ext,
            'filetypes': [
                (f'FITS ({ext})', ext),
                ('FITS (uncompressed)', fits_exts),
                (f'FITS gzipped ({ext}.gz)', f'{ext}.gz'),
                ('ZIP archive', '.zip'),
                ('All files', '*'),
            ],
            'initialfile': os.path.basename(file_name),
        }

        file_name = fd.asksaveasfilename(**options)

        if file_name != None and len(str(file_name)) > 0:
            if file_name.endswith('.zip'):
                # Save FITS first, then zip it
                fits_name = os.path.basename(base + ext)

                with tempfile.NamedTemporaryFile(suffix=ext, delete=False) as tmp:
                    temp_fits = tmp.name
                try:
                    self.console_msg("Saving FITS as " + str(file_name))
                    fits.writeto(temp_fits, image_data, header, overwrite=True)
                    with zipfile.ZipFile(file_name, 'w', zipfile.ZIP_DEFLATED) as zf:
                        zf.write(temp_fits, fits_name)  # archive name keeps the nice filename inside the zip
                except Exception as e:
                    self.error_raised = True
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  +" "+str(e), level=logging.ERROR)
                finally:
                    os.remove(temp_fits)

            else: #  .gz or .fits
                self.console_msg("Saving FITS as " + str(file_name))
                try:
                    fits.writeto(file_name, image_data, header, overwrite=True)
                except Exception as e:
                    self.error_raised = True
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  +" "+str(e), level=logging.ERROR)

            self.console_msg("Ready")


    ###############################################################
    #
    #
    #  save_FITS_file
    # 
    #  If file has a .gz or .zip file extension, then the uncompressed
    #  file is compressed accordingly before saving.
    #
    ###############################################################
    def save_FITS_file(self):
        global image
        global image_data
        global header
        file_name = self.image_file

        try:
            if file_name != None and len(str(file_name)) > 0:
                if file_name.endswith('.zip'):
                    fits_name = self.inner_fits_name or 'image.fits'
                    temp_fits = os.path.join(tempfile.gettempdir(), fits_name)
                    try:
                        self.console_msg("Saving FITS " + str(file_name))
                        fits.writeto(temp_fits, image_data, header, overwrite=True)
                        with zipfile.ZipFile(file_name, 'w', zipfile.ZIP_DEFLATED) as zf:
                            zf.write(temp_fits, fits_name)
                    finally:
                        if os.path.exists(temp_fits):
                            os.remove(temp_fits)
                else: #  .gz or .fits
                    self.console_msg("Saving FITS " + str(file_name))
                    try:
                        fits.writeto(file_name, image_data, header, overwrite=True)
                    except Exception as e:
                        self.error_raised = True
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  +" "+str(e), level=logging.ERROR)

                self.console_msg("Saved.")

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)   +" "+str(e), level=logging.ERROR)
            
##########################################################################################
#
# find_peaks
#
# Find peaks in loaded image above threshold (median + (std*10))
#
#
# Important table definitions:
#
# peaks_tbl, prelim_stars_tbl, isolated_stars_tbl, stars_tbl
# 
# peaks_tbl: initial set of stars from find_peaks
# prelim_stars_tbl: non_saturated_stars_tbl minus the ones near the edge
# isolated_stars_tbl: prelim_stars_tbl minus ones with close companions
# stars_tbl: isolated_stars_tbl minus ones rejected by user
#
#
# candidate_stars: Extracted stars from stars_tbl (cutouts from image). These
#                  are plotted in the selstars pages
#
#
##########################################################################################

    def find_peaks(self):
        global header
        try:
            self.console_msg("Find Peaks: starting find peaks...")

            #make sure an image is loaded
            if len(image_data) == 0:
                self.console_msg("Find Peaks: cannot proceed; an image must be loaded first; use File->Open...")
                return

            # test fit_width
            _shape = self.fit_width_entry.get().strip()
            if not _shape or not _shape.isnumeric():
                self.console_msg("Find Peaks: fitting Width not set (correctly); Set Fitting Width and Height in Setting Window")
                self.console_msg("Ready")
                return
            else:
                self.fit_shape = int(_shape)
                #if not odd (which IterativePSFPhotometry requires) make it so
                if self.fit_shape % 2 == 0:
                    self.fit_shape += 1


            # test fwhm
            if is_number(self.fwhm_entry.get().strip()):
                fwhm = float(self.fwhm_entry.get())
            else:
                self.console_msg("FWHM not numeric, using 3")
                fwhm = 3.0

            # test star_detection_threshold_factor
            if self.star_detection_threshold_factor_entry.get().isnumeric():
                star_detection_threshold_factor = int(self.star_detection_threshold_factor_entry.get())
            else:
                self.console_msg("DAOStarFinder threshold factor not numeric, using 10")
                star_detection_threshold_factor = 10

            # linearity_limit will get passed to DAOStarFinder
            # test linearity_limit_entry 
            linearity_limit = self.linearity_limit_entry.get().strip()
            if not linearity_limit or not linearity_limit.isnumeric():
                self.console_msg("Find Peaks: linearity limit is not valid....setting to 60000")
                linearity_limit = 60000

            # "Max num of Peaks" will get passed as brightest to DAOStarFinder
            # Test the user setting of "Max num of Peaks"
            user_npeaks = self.find_peaks_npeaks_entry.get().strip()
            if not user_npeaks or not user_npeaks.isnumeric():
                self.console_msg("Find Peaks: setting max num of peaks to 'unlimited'")
                user_npeaks = None
            else:
                self.console_msg("Find Peaks: limiting max num of Peaks to the user setting: " + str(int(user_npeaks)))
                user_npeaks = int(user_npeaks)


            """
            Determine the background using simple statistics
            ------------------------------------------------
            """
            # just for reference, lets looks at these stats first
            mean, median, std = sigma_clipped_stats(image_data, sigma=2.0)
            self.console_msg("Find Peaks: median sigma clipped level: " + str(round(median,2)))
            self.console_msg("Find Peaks: mean sigma clipped level: " + str(round(mean,2)))
            self.console_msg("Find Peaks: std sigma clipped level: " + str(round(std,2)))

            # now ready to find peaks
            star_find = DAOStarFinder(threshold = star_detection_threshold_factor*std,
                                       fwhm=fwhm,
                                        peakmax=float(linearity_limit),
                                        brightest=user_npeaks,
                                        exclude_border=True,
                                        min_separation=fwhm
                                       )

            peaks_tbl = QTable(star_find(data=image_data - median))

            if peaks_tbl == None or len(peaks_tbl) == 0:
                self.console_msg("Find Peaks: no peaks found!!!")
                self.console_msg("Ready")
                return
            else:
                peaks_tbl_len = len(peaks_tbl)
                self.console_msg("Find Peaks: found " + str(peaks_tbl_len) + " peaks.")

            #
            # Important table definitions
            #
            # peaks_tbl, prelim_stars_tbl, isolated_stars_tbl, stars_tbl
            # 
            # peaks_tbl: initial set of stars from find_peaks
            # prelim_stars_tbl: non_saturated_stars_tbl minus the ones near the edge
            # isolated_stars_tbl: prelim_stars_tbl minus ones with close companions
            # stars_tbl: isolated_stars_tbl minus ones rejected by user
            # candidate_stars: extracted stars from stars_tbl
            #

            #
            # mask out peaks near the edge
            #
        
            
            size = 2*self.fit_shape + 1 # Eg., 11
            hsize = (size - 1)/2 
            x = peaks_tbl['xcentroid']  
            y = peaks_tbl['ycentroid']  
            peak = peaks_tbl['peak']
            _image = Image.fromarray(image_data)
            width, height = _image.size
            mask = ((x > hsize) & (x < (width -1 - hsize)) &
                    (y > hsize) & (y < (height -1 - hsize)))  

            # prelim_stars_tbl are inbound stars (not to close to the edge)
            prelim_stars_tbl = Table()
            prelim_stars_tbl['x'] = x[mask]  
            prelim_stars_tbl['y'] = y[mask]  
            prelim_stars_tbl['peak'] = peak[mask]  
            prelim_stars_tbl['rejected'] = False #init

            prelim_stars_tbl_len = len(prelim_stars_tbl)
            self.console_msg("Find Peaks: found and removed " + str(peaks_tbl_len - prelim_stars_tbl_len) + " peaks on edge.")
            self.console_msg("Find Peaks: " + str(prelim_stars_tbl_len) + " peaks remain.")

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
            peak = prelim_stars_tbl['peak']
            reject_this = prelim_stars_tbl['rejected']

            mask = reject_this == False  # only keep ones we don't reject

            self.isolated_stars_tbl = Table()
            self.isolated_stars_tbl['x'] = x[mask]  
            self.isolated_stars_tbl['y'] = y[mask]  
            self.isolated_stars_tbl['peak'] = peak[mask]  
            self.isolated_stars_tbl['rejected'] = False #init

            isolated_stars_tbl_len = len(self.isolated_stars_tbl)
            self.console_msg("Find Peaks: found and removed " + str(prelim_stars_tbl_len - isolated_stars_tbl_len) + " close companions.")
            self.console_msg("Find Peaks: " + str(isolated_stars_tbl_len) + " peaks remain.")


            #re-count
            isolated_stars_tbl_len = len(self.isolated_stars_tbl)

            # now set 'rejected' to True for any stars that are proximate to a 
            # coordinate in ePSF_rejection_list
            # The 'x' and 'y' columns each do not necessariy contain unique values.
            # But the combination of multiple columns results in unique rows.
            self.isolated_stars_tbl.add_index(['x', 'y'])
            for isolated_index, isolated_row in enumerate(self.isolated_stars_tbl):
                psf_x = isolated_row['x']
                psf_y = isolated_row['y']
                for index, row in self.ePSF_rejection_list.iterrows():
                    reject_x = row['x']
                    reject_y = row['y']
                    if abs(reject_x - psf_x) <= hsize+3 and abs(reject_y - psf_y) <= hsize+3:
                        #user does not want this one
                        self.isolated_stars_tbl[isolated_index]['rejected'] = True
                        break
    
            x = self.isolated_stars_tbl['x']  
            y = self.isolated_stars_tbl['y']  
            reject_this = self.isolated_stars_tbl['rejected']

            mask = reject_this == False  # only keep ones we don't reject

            self.stars_tbl = Table()
            self.stars_tbl['x'] = x[mask]  
            self.stars_tbl['y'] = y[mask]  

            stars_tbl_len = len(self.stars_tbl)
            self.console_msg("Find Peaks: found and removed " + str(isolated_stars_tbl_len - stars_tbl_len) + " peaks rejected by user.")
            self.console_msg("Find Peaks: " + str(stars_tbl_len) + " peaks remain for EPSF Builder.")

            self.clear_selstars()

            # subtract background
            mean_val, median_val, std_val = sigma_clipped_stats(image_data, sigma=2.0)
            clean_image = image_data - median_val

            working_image = NDData(data=clean_image)

            #
            # 1 set of extracted stars is needed:  
            #  1) candidate_stars that only include non-rejected stars for the EPSFBuiler
            #

            self.candidate_stars = extract_stars(working_image, self.stars_tbl, size=size)  
            self.candidate_stars_index = 0
            
            for self.candidate_stars_index in range(min(len(self.candidate_stars),(self.nrows*self.ncols))):
                norm = simple_norm(self.candidate_stars[self.candidate_stars_index], 'log', percent=99.0)
                self.selstars_plot[self.candidate_stars_index].imshow(self.candidate_stars[self.candidate_stars_index],
                     norm=norm, origin='lower', cmap='viridis')

            #self.console_msg("candidate_stars index = " + str(self.candidate_stars_index), level=logging.DEBUG)

            plt.subplots_adjust(hspace=self.selstars_hspace, wspace=self.selstars_wspace)
            self.selstars_plot_canvas.draw()
            self.fig_selstars.canvas.mpl_connect('button_press_event', self.mouse_selstars_canvas_click)

            # display the rejected ones if any (with red circle) on main canvas
            self.ePSF_samples_plotted = True
            self.display_image()

            # update label with page (n of x) display E.g., Page: 1 of 6
            self.update_selstars_page_label(page_num=1)

            self.console_msg("Ready")

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  +" "+str(e), level=logging.ERROR)


##########################################################################################
#
# update_selstars_page_label
#
#
##########################################################################################

    def update_selstars_page_label(self, page_num):
        total_number_pages = math.ceil(len(self.candidate_stars)/(self.nrows*self.ncols))
        self.selstars_page_num_label['text'] = "Page: " + str(page_num) + " of " + str(total_number_pages)
        #Enable Forward and Back buttons accordingly
        # if on the first page, diable back
        if page_num == 1:
            self.back_selstars_button.config(state=tk.DISABLED)
        else:
            self.back_selstars_button.config(state=tk.NORMAL)

        # if on the last page disable Forward
        if page_num == total_number_pages:
            self.forward_selstars_button.config(state=tk.DISABLED)
        else:
            self.forward_selstars_button.config(state=tk.NORMAL)

        return


##########################################################################################
#
# plot_psf_model
#
# Plot model 
#
#
##########################################################################################

    def plot_psf_model(self, data, plotting_gaussian=False, plotting_loaded_psf=False, loaded_psf_file=None):

        model_length = len(data)
        x = np.arange(0, model_length, 1)
        y = np.arange(0, model_length, 1)
        x, y = np.meshgrid(x, y)

        #plot the psf model
        self.psf_plot.clear()
        self.psf_plot.plot_surface(x, y, data, cmap=cm.jet)
        self.psf_plot_canvas.draw()

        #plot ePSF samples
        self.clear_epsf_plot()
        if not plotting_gaussian and not plotting_loaded_psf:
            self.display_ePSF_samples()

        norm = simple_norm(data, 'log', percent=99.)
        im = self.ePSF_plot.imshow(data, norm=norm, origin='lower', cmap='viridis')
        if plotting_gaussian and self.use_gaussian_prf_model.get() != '1':
            self.ePSF_plot.set_title("Moffat")
        elif plotting_gaussian and self.use_gaussian_prf_model.get() != '0':
            self.ePSF_plot.set_title("Circular Gaussian")
        else:
            self.ePSF_plot.set_title("Effective PSF")

        if plotting_loaded_psf:
            self.ePSF_plotname_label.config(text="Effective PSF: " +loaded_psf_file)
        else:
            self.ePSF_plotname_label.config(text="PSF")


        self.fig_ePSF.colorbar(im, ax=self.ePSF_plot)

        self.ePSF_plot_canvas.draw()

        return

        


##########################################################################################
#
# create_ePSF
#
# Create and Effective PSF (Point Spread Function from peaks)
#
#
##########################################################################################

    def create_ePSF(self):
        global header
        self.console_msg("Starting Effective PSF building...")

        #make sure an image is loaded
        if len(image_data) == 0:
            self.console_msg("Cannot proceed; an image must be loaded first; use File->Open...")
            return

        if self.candidate_stars == None:
            self.console_msg("Cannot proceed; Run Find Peaks first; use File->Photometry->Find Peaks")
            return

        try:
            # test oversampling
            if is_number(self.oversampling_entry.get().strip()):
                self.current_oversampling = int(abs(float(self.oversampling_entry.get())))
            else:
                self.console_msg("oversampling not numeric, using 2, (most common)")
                self.current_oversampling = 2


            self.console_msg("Starting ePSF Builder...(check console progress bar)")
            epsf_builder = EPSFBuilder(oversampling=self.current_oversampling, maxiters=50, progress_bar=True) 

            self.epsf_model, fitted_stars = epsf_builder(stars=self.candidate_stars)  

            self.console_msg("self.epsf_model.data.shape="+str(self.epsf_model.data.shape))
            self.console_msg("ePSF Normalization Index = "+format(self.epsf_model.data.sum()/self.current_oversampling**2, '.2f'))

            self.plot_psf_model(self.epsf_model.data)
    
            self.ePSF_samples_plotted = True

            self.console_msg("Ready")

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  +" "+str(e), level=logging.ERROR)

##########################################################################################
#
# load_ePSF
#
# Load Effective PSF from file system
#
#
##########################################################################################

    def load_ePSF(self):
        try:
            options = {}
            options['defaultextension'] = '.epsf'
            options['filetypes'] = [('ePSF', 'epsf')]
            options['title'] = 'Load Effective PSF (ePSF)...'

            file_name = fd.askopenfilename(**options)

            if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
                self.console_msg("Loading ePSF from: " + str(file_name))
                with open(file_name, "rb") as file:
                    bundle = pickle.load(file)

                # check for old files
                if isinstance(bundle, dict):
                    self.epsf_model = bundle['epsf_model']
                    self.current_oversample = bundle['oversampling']
                else:
                    self.epsf_model = bundle
                    self.console_msg("WARNING: " + str(file_name) + " does not contain oversampling info")
        
                self.console_msg("Current oversampling = "+format(self.current_oversampling, 'd'))
                self.console_msg("self.epsf_model.data.shape="+str(self.epsf_model.data.shape))
                self.console_msg("ePSF Normalization Index = "+format(self.epsf_model.data.sum()/self.current_oversampling**2, '.4f'))

                self.plot_psf_model(self.epsf_model.data, plotting_loaded_psf=True,
                                     loaded_psf_file=os.path.basename(str(file_name)))
        
                self.ePSF_samples_plotted = True
                
            else:
                return
            
            self.console_msg("Ready")

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)

##########################################################################################
#
# save_as_ePSF
#
# Save As... Effective PSF to the file system(Point Spread Function from peaks)
#
#
##########################################################################################

    def save_as_ePSF(self):
        options = {}
        options['defaultextension'] = '.epsf'
        options['filetypes'] = [('ePSF', 'epsf')]
        options['title'] = 'Save Current Effective PSF (ePSF) as...'

        file_name = fd.asksaveasfile(**options)

        try:
            if len(str(file_name)) > 0:
                # bumdle current_oversampling along with the data
                bundle = {
                    'epsf_model': self.epsf_model,
                    'oversampling': self.current_oversampling
                }
                self.console_msg("Saving current ePSF as " + str(file_name.name))
                with open(str(file_name.name), 'wb') as file:
                    pickle.dump(bundle, file)
                self.ePSF_plotname_label.config(text="Effective PSF: "+os.path.basename(str(file_name.name)))
                self.console_msg("Saving ePSF as " + str(file_name.name))
                self.console_msg("Ready")
            else:   
                return
        
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)

##########################################################################################
#
# clear_ePSF
#
# Clear ePSF Plot and selstars
#
#
##########################################################################################

    def clear_ePSF(self):
        global header
        self.console_msg("Clearing ePSF model, Rejection List, plot...")
        #drop all the rows but keep the 'x' and 'y' column
        self.ePSF_pending_rejection_list.drop(self.ePSF_pending_rejection_list.index, inplace=True)
        self.ePSF_rejection_list.drop(self.ePSF_rejection_list.index, inplace=True)
        self.candidate_stars = None
        self.epsf_model = None #reset
        self.stars_tbl = None
        self.isolated_stars_tbl = None
        self.clear_psf_label()
        self.clear_epsf_plot()
        self.clear_selstars()
        self.ePSF_samples_plotted = False
        self.display_image()
        self.console_msg("Ready")
        return

##########################################################################################
#
# clear_selstars
#
# Clear selstars and reset label
#
#
##########################################################################################

    def clear_selstars(self):
        self.clear_selstars_plot()
        self.selstars_page_num_label['text'] = "Page:"
        return

    def load_ePSF_rejection_list(self):
        global header

        try:
            options = {}
            options['defaultextension'] = '.csv'
            options['filetypes'] = [('CSV', '.csv')]
            options['title'] = 'Choose the rejection list (*.csv)'

            file_name = fd.askopenfilename(**options)

            if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
                self.console_msg("Loading Rejection list from: " + str(file_name))
                self.ePSF_rejection_list = pd.read_csv(str(file_name))
                self.ePSF_rejection_list["stale"] = True #reset 
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

    def execute_iterative_psf_photometry(self):
        global header
        self.console_msg(
            "Starting iterative PSF photometry...")
        if len(image_data) == 0:
            self.console_msg("Cannot proceed; an image must be loaded first; use File->Open...")
            self.console_msg("Ready")
            return
        try:
            # test fwhm
            if is_number(self.fwhm_entry.get().strip()):
                fwhm = float(self.fwhm_entry.get())
            else:
                self.console_msg("FWHM not numeric, using 3")
                fwhm = 3.0


            # Confirm use of Circular Gaussian before DAOStarFinder
            # Why ask after finding stars if user then cancels...it wastes time
            if self.epsf_model != None:
                self.console_msg("Using derived Effective PSF Model")
                psf_model = self.epsf_model
            else:
                #Ask if using Gausian or Moffat OK
                if self.use_gaussian_prf_model.get() != '1':
                    # test moffat beta 
                    if is_number(self.moffat_beta_entry.get().strip()):
                        moffat_beta = abs(float(self.moffat_beta_entry.get()))
                    else:
                        moffat_beta = 1 # Lorentzian
                    self.console_msg("Moffat \u03B2 set to: " + str(moffat_beta))

                    result = askokcancel(title="Use Moffat Model?", message="Is it OK to use Moffat, \u03B2="+str(moffat_beta)+" model?")

                    if result==False:
                        self.console_msg("User cancelling iterative PSF photometry")
                        return

                    # Here user wants to use Moffat Model
                    alpha=fwhm/2/math.sqrt((2**(1/moffat_beta))-1)
                    self.console_msg("Using Moffat for model; fwhm = "+format(fwhm, '.1f')+"; alpha="+format(alpha, '.1f')+"; \u03B2="+str(moffat_beta))
                    psf_model = MoffatPSF(x_0=22.0, y_0=22.0, alpha=alpha, beta=moffat_beta)

                else:
                    result = askokcancel(title="Use Circular Gaussian Model?", message="Is it OK to use Circular Gaussian model?")

                    if result==False:
                        self.console_msg("User cancelling iterative PSF photometry")
                        return
                    
                    # Here user wants to use CircularPRF model

                    """
                        Create a PRF image model from a Circular Gaussian PRF.
                        In this case, we use the CircularGaussianPRF model 
                        directly as a PSF model
                    """
                    self.console_msg("Using Circular Gaussian PRF for model; fwhm = "+format(fwhm, '.1f'))
                    psf_model = CircularGaussianPRF(x_0=22.0, y_0=22.0, fwhm=fwhm)

                # Plot either CircularGaussian or Moffat
                yy, xx = np.mgrid[:45, :45]
                psf_data = psf_model(xx, yy)
                self.plot_psf_model(psf_data, plotting_gaussian=True)


            # test star_detection_threshold_factor
            if self.star_detection_threshold_factor_entry.get().isnumeric():
                star_detection_threshold_factor = int(self.star_detection_threshold_factor_entry.get())
            else:
                self.console_msg("DAOStarFinder threshold factor not numeric, using 10")
                star_detection_threshold_factor = 10

            # test iterations
            if self.photometry_iterations_entry.get().isnumeric():
                iterations = int(self.photometry_iterations_entry.get())
            else:
                self.console_msg("Photometry Iteration not positive int, using 3")
                iterations = 3

            # test sharplo
            if is_number(self.sharplo_entry.get().strip()):
                sharplo = float(self.sharplo_entry.get())
            else:
                self.console_msg("Lower Bound for Sharpness not numeric, using 0.2 (DAOStarFinder default)")
                sharplo = 0.2 # DAOStarFinder default

            # test min_separation_bias_entry
            if is_number(self.min_separation_bias_entry.get().strip()):
                min_separation_bias = float(self.min_separation_bias_entry.get().strip())
            else:
                self.console_msg("Min Separation: NaN; setting to 0")
                min_separation_bias = 0

            self.console_msg("Min Separation bias set to: " + format(min_separation_bias, '4.2f'))
            min_separation = fwhm/2 + (self.fit_shape-1)/2 + min_separation_bias
            grouper = SourceGrouper(min_separation)
            self.console_msg("Min Separation for Grouper set to: " + format(min_separation, '4.2f'))

            # test max qfit
            if is_number(self.max_qfit_entry.get().strip()):
                max_qfit = float(self.max_qfit_entry.get())
            else:
                max_qfit = 0.1
            self.console_msg("Max qfit set to: " + str(max_qfit))

            #
            # Determine the background using simple statistics
            #
            
            # just for reference, lets looks at these stats first
            mean, median_val, std = sigma_clipped_stats(image_data, sigma=2.0)
            self.console_msg("Median sigma clipped level: " + str(round(median_val,2)))
            self.console_msg("Mean sigma clipped level: " + str(round(mean,2)))
            self.console_msg("Std sigma clipped level: " + str(round(std,2)))

            # pass linearity_limit to DAOStarFinder
            # test linearity_limit_entry 
            linearity_limit = self.linearity_limit_entry.get().strip()
            if not linearity_limit or not linearity_limit.isnumeric():
                self.console_msg("Find Peaks: linearity limit is not valid....setting to 60000")
                linearity_limit = 60000

            # "Max num of Peaks" will get passed as brightest to DAOStarFinder
            # Test the user setting of "Max num of Peaks"
            user_npeaks = self.find_peaks_npeaks_entry.get().strip()
            if not user_npeaks or not user_npeaks.isnumeric():
                self.console_msg("Find Peaks: setting max num of peaks to 'unlimited'")
                user_npeaks = None
            else:
                self.console_msg("Find Peaks: limiting max num of Peaks to the user setting: " + str(int(user_npeaks)))
                user_npeaks = int(user_npeaks)

            star_find = DAOStarFinder(threshold = star_detection_threshold_factor*std,
                                      fwhm=fwhm,
                                      sharplo=sharplo,
                                      peakmax=float(linearity_limit),
                                      brightest=user_npeaks,
                                      exclude_border=True,  
                                      min_separation=fwhm  
                                      )

            # subtract background
            clean_image = image_data - median_val

            working_image = NDData(data=clean_image)

            # How many stars will it find?
            star_finder_result = star_find.find_stars(clean_image)
            if star_finder_result == None or len(star_finder_result) == 0:
                self.console_msg("DAOStarFinder did not find any stars, adjust lower bound of sharpness")
                self.console_msg("Ready")
                return

            self.console_msg("DAOStarFinder number of stars found : "  + str(len(star_finder_result)))

            bkgstat = MMMBackground()

            local_bkg = LocalBackground(inner_radius=fwhm*4, outer_radius=fwhm*8, bkg_estimator=bkgstat)

            if self.fitter_stringvar.get() == "Sequential LS Programming":
                self.console_msg(
                    "Setting fitter to Sequential Least Squares Programming")
                selected_fitter = SLSQPLSQFitter()

            elif self.fitter_stringvar.get() == "Simplex LS":
                self.console_msg(
                    "Setting fitter to Simplex and Least Squares Statistic")
                selected_fitter = SimplexLSQFitter()

            #default is TRF LS
            else:
                self.console_msg(
                    "Setting fitter to TRF and Least Squares Statistic")
                selected_fitter = TRFLSQFitter()

 
            photometry = IterativePSFPhotometry(
                                                psf_model = psf_model,
                                                fit_shape = self.fit_shape,
                                                finder = star_find,
                                                grouper = grouper,
                                                fitter = selected_fitter,
                                                #fitter_maxiters = 10,
                                                maxiters = iterations,
                                                localbkg_estimator = local_bkg,
                                                aperture_radius=1.5*fwhm,
                                                sub_shape=None, #defaults to model bounding box
                                                progress_bar=True
                                                )

            sys.setrecursionlimit(10000)
            self.console_msg("Starting Photometry...(check console progress bar)")
            result_tab = photometry(data=working_image)

            if 'message' in selected_fitter.fit_info:
                self.console_msg("Done. PSF fitter message(s): " + str(selected_fitter.fit_info['message']))
            else:
                self.console_msg("Done. PSF fitter; no message available")

            #get the residuals

            if self.generate_residual_image.get():
                """
                # Needs to be changed: psf_shape is now an optional keyword 
                # in the make_model_image and make_residual_image methods 
                # of PSFPhotometry and IterativePSFPhotometry. 
                # The value defaults to using the model bounding box to 
                # define the shape and is required only if the PSF model 
                # does not have a bounding box attribute.   
                """

                residual_image = photometry.make_residual_image(data=clean_image)

                #append current gm time to residual filename
                residual_file_name = self.image_file + "_residuals_" + strftime("%Y_%m_%d %H_%M_%S", gmtime()) + ".fit"
                fits.writeto(residual_file_name, residual_image, header, overwrite=True)
                self.console_msg("Residuals saved to: " + residual_file_name)

            # Remove any multidimensional columns
            goodnames = [name for name in result_tab.colnames if len(result_tab[name].shape) <=1]
            
            self.results_tab_df = result_tab[goodnames].to_pandas()

            """
            flags : bitwise flag values

                1 : one or more pixels in the fit_shape region were masked

                2 : the fit x and/or y position lies outside of the input data

                4 : the fit flux is less than or equal to zero

                8 : the fitter may not have converged. In this case, you can try increasing the maximum number of fit iterations using the fitter_maxiters keyword.

                16 : the fitter parameter covariance matrix was not returned

                32 : the fit x or y position is at the bounded value            
            
            """
            #
            # Remove any entries with flags != 0; 
            # 
            self.results_tab_df.drop(self.results_tab_df[self.results_tab_df['flags'] != 0].index, inplace=True)

            #
            # Remove poor fits qfit > max_qfit (defaults to .1)
            #
            self.results_tab_df.drop(self.results_tab_df[self.results_tab_df['qfit'] > max_qfit].index, inplace=True)

            self.results_tab_df["removed_from_ensemble"] = False

            _date_obs = self.date_obs_entry.get().strip()
            if not _date_obs or not is_number(_date_obs):
                self.console_msg("DATE-OBS not valid setting to 1 (4712 B.C.)")
                self.results_tab_df["date-obs"] = float(1)
            else:
                self.results_tab_df["date-obs"] = float(_date_obs)

            _date_end = self.date_end_jd
            if _date_end != None and is_number(_date_end):
                self.results_tab_df["date-end"] = float(_date_end)

            _airmass = self.airmass_entry.get().strip()
            if not _airmass or not is_number(_airmass):
                self.console_msg("AIRMASS not valid setting to na")
                self.results_tab_df["AMASS"] = "na"
            else:
                self.results_tab_df["AMASS"] = float(_airmass)

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

    def display_ePSF_samples(self):
        try:
            """
             Circle color (non inverted color/inverted color)
             white/black: stars that from find_peak and not rejected for ePSF Generation (isolated_stars_tbl)
             red/red: stars that rejected by user and in the isolated_stars_tbl (ePSF_rejection_list)
             yellow/yellow: stars in a loaded rejection list file that is not in the isolated_stars_tbl, they 
             have already been removed (ePSF_rejection_list)
             
            """
            ## make all the circles same as fit_shape; derive hsize (halfsize or radius)

            # test fit_width
            _shape = self.fit_width_entry.get().strip()
            if not _shape or not _shape.isnumeric():
                self.console_msg("Find Peaks: fitting Width not set (correctly); Set Fitting Width and Height in Setting Window")
                self.console_msg("Ready")
                return
            else:
                self.fit_shape = int(_shape)
                #if not odd (which IterativePSFPhotometry requires) make it so
                if self.fit_shape % 2 == 0:
                    self.fit_shape += 1

            size = 2*self.fit_shape + 1
            hsize = (size - 1)/2

            if len(self.isolated_stars_tbl) != 0:
                self.console_msg("Displaying ePSF samples; reject list size: " + str(len(self.ePSF_rejection_list)))

                #display the non-rejected stars as white, and rejected as red circles
                for psf_x, psf_y in self.isolated_stars_tbl.iterrows('x', 'y'):
                    color = ('white', "black")[self.is_inverted] # it is a white/black circle until a reject match is found
                    for index, row in self.ePSF_rejection_list.iterrows():
                        reject_x = row['x']
                        reject_y = row['y']

                        if abs(reject_x - psf_x) <= hsize and abs(reject_y - psf_y) <= hsize:
                            color = 'red'  #paint rejects red
                            self.ePSF_rejection_list.loc[index, "stale"] = False
                            break;
                    self.create_circle(x=psf_x * self.zoom_level, y=psf_y * self.zoom_level,
                                        r=hsize * self.zoom_level, canvas_name=self.canvas, outline=color)
                    
            #Always make a yellow circle for any stars that are "stale"
            for index, row in self.ePSF_rejection_list.iterrows():
                reject_x = row['x']
                reject_y = row['y']
                if row["stale"] == True:
                    self.create_circle(x=reject_x * self.zoom_level, y=reject_y * self.zoom_level,
                                        r=hsize * self.zoom_level, canvas_name=self.canvas, outline='yellow')

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) + " " + str(e), level=logging.ERROR)

    def plot_photometry(self, verbose=True):
        try:
            if verbose:
                self.console_msg("Plotting Photometry...")

            labels_in_photometry_table = "label" in self.results_tab_df
            vsx_ids_in_photometry_table = "vsx_id" in self.results_tab_df

            if os.path.isfile(self.image_file+".csv"):
                _shape = self.fit_width_entry.get().strip()
                if not _shape or not _shape.isnumeric():
                    self.console_msg("Cannot Plot Photometry with existing... ")
                    self.console_msg("...\"" + self.image_file + ".csv\"")
                    self.console_msg("...because fitting width is not recognized")
                    self.console_msg("Ready")
                    return
                else:
                    fit_shape = int(_shape)
                    #if not odd (which IterativePSFPhotometry requires) make it so
                    if fit_shape % 2 == 0:
                        fit_shape += 1

                self.fit_shape = fit_shape
                self.results_tab_df = pd.read_csv(self.image_file + ".csv")
                if "removed_from_ensemble" not in self.results_tab_df:
                    # This prefilling is required for backwards compatibility to read .phot
		            #(now called .csv) files from older versions.
                    self.results_tab_df["removed_from_ensemble"] = False

                self.photometry_results_plotted = True
                sel_comps = [] #init

                if labels_in_photometry_table and \
                   self.display_all_objects.get() != '1':
                    # if here we need to get the list of user selected comp stars
                    # ony the ones in this list will get displayed
                    #
                    # loop through all the selected comp stars and fill this sel_comps list
                    sel_comps = [] #init
                    sel_comps_to_use = self.object_sel_comp_entry.get()
                    sel_comps_to_use = [comp.strip() for comp in sel_comps_to_use.split(',')]
                    # Make unique, just in case
                    sel_comps_to_use = list(set(sel_comps_to_use))

                    #make array of int called sel_comps            
                    for comp in sel_comps_to_use:
                        sel_comps.append(comp.strip())

                    # add check star to the list so it can be displayed                       
                    check_star_to_use = self.object_kref_entry.get().strip()
                    if is_number(check_star_to_use):
                        # remove the check star from this list just in case
                        # this removes all of them 
                        sel_comps = [star for star in sel_comps if star != check_star_to_use]
                        #add it to the list of comps
                        sel_comps.append(check_star_to_use)

                for index, row in self.results_tab_df.iterrows():
                    outline = ("grey50", "black")[self.is_inverted]
                    if labels_in_photometry_table:
                        if str(row["label"]) != __empty_cell__: 
                               if self.display_all_objects.get() == '1' or \
                                  str(row["label"])[len(__label_prefix__):] in sel_comps: #ignore label prefix
                                # here if all comps are to be displayed or 
                                # if only user's comps are being displayed 
                                    outline = ("pink", "black")[self.is_inverted]
                                    self.create_circle(x=row["x_fit"] * self.zoom_level,
                                        y=row["y_fit"] * self.zoom_level,
                                        r=(self.fit_shape/2) * self.zoom_level,
                                        canvas_name=self.canvas, outline=outline)
                                    self.create_text(  x=row["x_fit"] * self.zoom_level,
                                        y=row["y_fit"] * self.zoom_level, 
                                        r=self.fit_shape * self.zoom_level,
                                        canvas_name=self.canvas,
                                        anchor=tk.CENTER,
                                        text=str(row["label"])[len(__label_prefix__):],
                                        fill=("pink", "black")[self.is_inverted])
                                    continue

                    if row["removed_from_ensemble"]:
                        assert False, "Found an entry 'removed from ensemble???!'"

                    if vsx_ids_in_photometry_table:
                        if str(row["vsx_id"]) != __empty_cell__:
                            if self.display_all_objects.get() == '1' or \
                               str(row["vsx_id"]) == self.object_name_entry.get().strip():
                                # here if all vsx objects to be displayed or 
                                # if only user's "Object Name" object is being displayed 
                                outline = ("yellow", "black")[self.is_inverted]
                                self.create_circle(x=row["x_fit"] * self.zoom_level,
                                    y=row["y_fit"] * self.zoom_level,
                                    r=(self.fit_shape/2) * self.zoom_level,
                                    canvas_name=self.canvas,
                                    outline=outline)
                                self.create_text(  x=row["x_fit"] * self.zoom_level,
                                    y=row["y_fit"] * self.zoom_level, 
                                    r=self.fit_shape * self.zoom_level,
                                    canvas_name=self.canvas,
                                    anchor=tk.CENTER,
                                    text=str(row["vsx_id"]).strip(),
                                    fill=('yellow', "black")[self.is_inverted])

            if verbose:
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
            return(True,  #indicate match
                   matched_objects.iloc[0]["x_fit"],
                   matched_objects.iloc[0]["y_fit"],
                   matched_objects.iloc[0]["flux_fit"],
                   matched_objects.iloc[0]["inst_mag"]
                   )
        else:
            return (False, # indicate no match
                     0, 0, 0, 0)
        
    ###############################################################
    #
    #
    #  mouse_selstars_canvas_click
    #
    #  callback when button press down in selstar caavas
    #  
    #  Place "Reject" title on Ax to indicate that this star to be rejected
    #  It can get added to the ePSF_rejection_list by the Submit button
    #
    #  "Reject" markers persist when going forward or backward
    #  by reading the ePSF_pending_rejection_list
    #
    #  If Ax already has a "Reject" then remove it
    #
    #
    ###############################################################
   
    def mouse_selstars_canvas_click(self,event):
        myax = event.inaxes
        if myax == None:
            return
        selected_index = myax.figure.axes.index(myax)
        self.console_msg("Ax number: " + str(selected_index))

        # Calculate index into self.candidate_stars that was just mouse clicked.
        # self.candidate_stars_index is pointing to the last displayed candidate in the displayed page.
        # Rewind the pointer to the first one in the page, then add selected_index to it.
        # Just subtracting (candidate_stars_index % (self.ncols*self.nrows)) from candidate_stars_index
        # will get the pointer to the first one in the page.
        #
        # First one in the page:
        candidate_stars_selected_index = self.candidate_stars_index - (self.candidate_stars_index % (self.ncols*self.nrows))
        # selected one with mouse click:
        candidate_stars_selected_index += selected_index
        norm = simple_norm(self.candidate_stars[candidate_stars_selected_index], 'log', percent=99.0)
        self.selstars_plot[selected_index].imshow(self.candidate_stars[candidate_stars_selected_index],
                     norm=norm, origin='lower', cmap='viridis')
        #Check if "Reject" is already there. If it is then remove it from Ax and removed from ePSF_pending_rejection_list
        (cand_x, cand_y) = self.candidate_stars[candidate_stars_selected_index].center
        if not ((self.ePSF_pending_rejection_list['x'] == cand_x) & (self.ePSF_pending_rejection_list['y'] == cand_y)).any():
            # No reject, 
            # Add it in 
            self.selstars_plot[selected_index].text(x=0,y=5, s="Reject")
            #update ePSF_pending_rejection_list
            self.ePSF_pending_rejection_list.loc[len(self.ePSF_pending_rejection_list.index)] = [cand_x, cand_y, True]
        else:
            # "Reject" already in; erase it
            self.selstars_plot[selected_index].clear()
            # remove from ePSF_pending_rejection_list (using mask)
            self.ePSF_pending_rejection_list = \
                self.ePSF_pending_rejection_list[~((self.ePSF_pending_rejection_list['x'] == cand_x) &
                                                     (self.ePSF_pending_rejection_list['y'] == cand_y))]

        self.selstars_plot[selected_index].imshow(self.candidate_stars[candidate_stars_selected_index],
                     norm=norm, origin='lower', cmap='viridis')
        plt.subplots_adjust(hspace=self.selstars_hspace, wspace=self.selstars_wspace)
        self.selstars_plot_canvas.draw()

        # Submit button is enabled only when there is something to submit
        if len(self.ePSF_pending_rejection_list) > 0:
            submit_button_state = tk.NORMAL
        else:
            submit_button_state = tk.DISABLED

        self.submit_rejects_selstars_button.config(state=submit_button_state)

        return

    ###############################################################
    #
    #
    #  on_auto_behavior_checkbox_checked
    #
    #  if checked, diable CCD Filter Text
    # 
    ###############################################################

    def on_auto_behavior_checkbox_checked(self):
        if self.auto_behavior.get():
            self.filter_entry.config(state='normal')
            self.set_entry_text(self.filter_entry, self.fits_header_filter)
            self.filter_entry.config(state='disabled')
        else:
            self.filter_entry.config(state='normal')

    ###############################################################
    #
    #
    #  on_canvas_mousewheel
    #
    #  Zoom in and out
    # 
    ###############################################################

    def on_canvas_mousewheel(self, event):
        if event.delta > 0:
            self.zoom_level *= 1.1
        else:
            self.zoom_level /= 1.1
        self.canvas.scale("all", 0, 0, self.zoom_level, self.zoom_level)
        self.display_image(verbose=False)

    ###############################################################
    #
    #  on_canvas_shift_mousewheel
    #
    #  scroll horizontally 
    # 
    ###############################################################

    def on_canvas_shift_mousewheel(self, event):
        self.canvas.xview_scroll(-1 * (event.delta // 120), "units")  # Scroll horizontally


    ###############################################################
    #
    #  on_drag_start
    #
    #  start dragging image 
    # 
    ###############################################################
    
    def on_drag_start(self, event):
        # Stores initial click position.
        self.canvas.scan_mark(event.x, event.y)
    
    ###############################################################
    #
    #  on_drag_move
    #
    #  drag image 
    # 
    ###############################################################
    
    def on_drag_move(self, event):
        # Moves the image with the mouse while Shift + Button-1 is held.
        self.canvas.scan_dragto(event.x, event.y, gain=1)
        self.plot_photometry(verbose=False)

    ###############################################################
    #
    #  on_button_3_click
    #
    #  center image
    # 
    ###############################################################
    
    def on_button_3_click(self, event):
        """Center the view accounting for the visible portion"""
        # For a canvas where scroll region is larger than visible area
        # The visible fraction is what fraction of total content is visible
        
        # Get visible fraction (returns tuple: (left, right) for x, (top, bottom) for y)
        x_view = self.canvas.xview()
        y_view = self.canvas.yview()
        
        # Calculate how much is visible
        x_visible_fraction = x_view[1] - x_view[0]
        y_visible_fraction = y_view[1] - y_view[0]
        
        # Center position is 0.5 minus half the visible fraction
        x_center = 0.5 - (x_visible_fraction / 2.0)
        y_center = 0.5 - (y_visible_fraction / 2.0)
        
        # Move to center
        self.canvas.xview_moveto(x_center)
        self.canvas.yview_moveto(y_center)
        self.plot_photometry(verbose=False)


    ###############################################################
    #
    #
    #  mouse_main_canvas_click
    #
    # 
    ###############################################################

    def mouse_main_canvas_click(self, event):
        global image_data

        # if care if shift key is not pressed; only use <Shift-Button-1> callback 
        if event.state & 0x01: 
            return
        
        x = int(self.canvas.canvasx(event.x) / self.zoom_level)
        y = int(self.canvas.canvasy(event.y) / self.zoom_level)
        self.last_clicked_x = x
        self.last_clicked_y = y
        ADU = image_data[y-1, x-1]
        sky = self.wcs_header.pixel_to_world(x, y)
        sky_coordinate_string = ""
        clicked_coordinate = None

        if hasattr(sky, 'ra'):
            clicked_coordinate = SkyCoord(ra=sky.ra, dec=sky.dec)
            sky_coordinate_string = "α δ: " + clicked_coordinate.to_string("hmsdms", precision=2)
            self.console_msg("Position X: "+str(x)+"\t Y: "+str(y) +
                            "\t ADU: "+str(ADU) + "\t\t\t" + sky_coordinate_string)
            alpha_delta = [] # used to keep alpha and delta, may end up in settings
            alpha_delta = clicked_coordinate.to_string("hmsdms", precision=2).split()

        if self.ePSF_samples_plotted:
            #add the selected coordinate into the ePSF_rejection_list
            #initally all assumed to be stale until ePSF builder is run
            self.ePSF_rejection_list.loc[len(self.ePSF_rejection_list.index)] = [x, y, True]
            #indicate the rejected ones
            self.display_image()

        elif self.photometry_results_plotted:
            vsx_ids_in_photometry_table = "vsx_id" in self.results_tab_df
            self.display_image()
            self.console_msg("")

            fit_matched, x_fit, y_fit, flux_fit, inst_mag = self.match_photometry_table(x, y)
            if fit_matched:
                sky = self.wcs_header.pixel_to_world(x_fit, y_fit)
                sky_coordinate_string = ""
                if hasattr(sky, 'ra'):
                    c = SkyCoord(ra=sky.ra, dec=sky.dec)
                    sky_coordinate_string = " α δ: " + c.to_string("hmsdms", precision=2)
                if x_fit != 0 and y_fit != 0:
                    psf_canvas_x = x_fit
                    psf_canvas_y = y_fit
                if str(x_fit)+str(y_fit) in self.photometry_circles:
                    self.canvas.delete(self.photometry_circles[str(x_fit)+str(y_fit)])

                self.canvas.create_line(x_fit*self.zoom_level, y_fit*self.zoom_level - 35*self.zoom_level, x_fit *
                                        self.zoom_level, y_fit*self.zoom_level - 10*self.zoom_level, fill="white")  # Draw "target" lines
                self.canvas.create_line(x_fit*self.zoom_level+35*self.zoom_level, y_fit*self.zoom_level,
                                        x_fit*self.zoom_level + 10*self.zoom_level, y_fit*self.zoom_level, fill="white")
                self.console_msg("Photometry fits, X: " + str(round(x_fit, 2)) + " Y: " + str(round(y_fit, 2)) + " Flux (ADU): " + str(
                    round(flux_fit, 2)) + " Instrumental magnitude: " + str(round(inst_mag, 3)) + " " + sky_coordinate_string)

                if "match_id" in self.results_tab_df:
                    matching_star_criterion = (self.results_tab_df["x_fit"] == x_fit) & (
                        self.results_tab_df["y_fit"] == y_fit)
                    if len(self.results_tab_df[matching_star_criterion]) > 0:
                        matching_star = self.results_tab_df[matching_star_criterion].iloc[0]
                        if vsx_ids_in_photometry_table and len(str(matching_star["vsx_id"])) > 1:
                            self.console_msg("Matching VSX Source: " + str(matching_star["vsx_id"]))
                                
                        elif type(matching_star["match_id"]) in (str, int):
                            self.console_msg(
                                "Matching catalog source ID: " + str(matching_star["match_id"]) + 
                                    "; label: " + str(matching_star["label"]) +
                                    " magnitude: " + str(matching_star["match_mag"]))
                        else:
                            #ask if user wants to name this object
                            result = simpledialog.askstring("Object Name", "Replace Object Name (and assign vsx status) with: ", 
                                                            initialvalue="my-user-obj")

                            if result:
                                user_name = result.strip()
                                """
                                Remove any pre-existing user_name(s) in results_tab_df first.
                                If not, user could populate df with multiple user_names.
                                """
                                self.results_tab_df = self.results_tab_df.loc[self.results_tab_df['vsx_id'] != user_name]
                                self.console_msg("Object Name is now: " + user_name)
                                self.set_entry_text(self.object_name_entry, user_name)
                                self.results_tab_df.loc[matching_star.name, "vsx_id"] = user_name

                                if clicked_coordinate != None:
                                    self.set_entry_text(self.object_name_alpha_entry, alpha_delta[0])
                                    self.set_entry_text(self.object_name_delta_entry, alpha_delta[1])
                                    self.results_tab_df.loc[matching_star.name, "RAJ2000"] = clicked_coordinate.ra.deg
                                    self.results_tab_df.loc[matching_star.name, "DEJ2000"] = clicked_coordinate.dec.deg

                                self.results_tab_df.to_csv(self.image_file + ".csv", index=False)
                                self.console_msg("Photometry saved to " + str(self.image_file + ".csv") + "; len = " + str(len(self.results_tab_df)))
                                self.display_image()
                                
            else:
                # These lines are "red" because object not in table
                self.canvas.create_line(x*self.zoom_level, y*self.zoom_level - 35*self.zoom_level,
                                        x*self.zoom_level, y*self.zoom_level - 10*self.zoom_level, fill="red")  # Draw "target" lines
                self.canvas.create_line(x*self.zoom_level+35*self.zoom_level, y*self.zoom_level,
                                        x*self.zoom_level + 10*self.zoom_level, y*self.zoom_level, fill="red")


            self.console_msg("Ready")

    ###############################################################
    #
    #
    #  clear_epsf_plot
    # 
    #
    ###############################################################
    def clear_epsf_plot(self):
        self.ePSF_plot.clear()
        self.fig_ePSF.clear()
        self.ePSF_plot_canvas.draw()
        self.ePSF_plotname_label.config(text="PSF")
        self.fig_ePSF, self.ePSF_plot = plt.subplots()
        self.ePSF_plot_canvas = FigureCanvasTkAgg(self.fig_ePSF, self.right_frame)
        self.ePSF_plot_canvas.draw()
        self.ePSF_canvas = self.ePSF_plot_canvas.get_tk_widget()
        self.ePSF_canvas.config(width=int(self.screen_width/8.5), height=int(self.screen_width/8.5))
        # Allocate small PSF canvas to a new grid inside the right_frame
        self.ePSF_canvas.grid(row=3, column=0)   #was row0

    ###############################################################
    #
    #
    #  clear_selstars_plot
    # 
    #
    ###############################################################
    def clear_selstars_plot(self):
        self.fig_selstars.clear()
        #plt.subplots_adjust(hspace=self.selstars_hspace, wspace=self.selstars_wspace)
        #self.selstars_plot_canvas.draw()
        self.fig_selstars, self.selstars_plot = plt.subplots(nrows=self.nrows,
                                                              ncols=self.ncols,
                                                              squeeze=False)
        self.selstars_plot = self.selstars_plot.ravel()
        self.selstars_plot_canvas = FigureCanvasTkAgg(self.fig_selstars, self.right_frame)
        plt.subplots_adjust(hspace=self.selstars_hspace, wspace=self.selstars_wspace)
        self.selstars_plot_canvas.draw()
        self.selstars_canvas = self.selstars_plot_canvas.get_tk_widget()
        self.selstars_canvas.config(width=int(self.screen_width/5), height=int(self.screen_width/5))
        # Allocate small PSF canvas to a new grid inside the right_frame
        self.selstars_canvas.grid(row=5, column=0)   #was row0


    def clear_psf_label(self):
        #clear plot label
        self.plotname_label['text'] = "Plot: "
        self.psf_plot.clear()
        self.psf_plot_canvas.draw()

    def zoom_in(self):
        self.zoom_level *= 1.1
        self.canvas.scale("all", 0, 0, self.zoom_level, self.zoom_level)
        self.console_msg("Zoom: "+str(self.zoom_level))
        self.display_image(verbose=False)

    def zoom_out(self):
        self.zoom_level /= 1.1
        self.canvas.scale("all", 0, 0, self.zoom_level, self.zoom_level)
        self.console_msg("Zoom: "+str(self.zoom_level))
        self.display_image(verbose=False)

    def zoom_100(self):
        self.canvas.scale("all", 0, 0, 1, 1)
        self.zoom_level = 1
        self.console_msg("Zoom: " + str(self.zoom_level))
        self.display_image()

    def invert_image(self):
        """Invert the displayed image and refresh the display"""
        try:
            global image_data
            if len(image_data) > 0:
                self.is_inverted = not self.is_inverted
                self.display_image()
                self.console_msg("Image display inverted" if self.is_inverted else "Image display restored")
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) + " " + str(e), level=logging.ERROR)

    ###############################################################
    #
    #
    #  solve_image
    # 
    #
    ###############################################################
    def solve_image(self):
        global generated_image
        global header

        if "flux_fit" not in self.results_tab_df:
            self.console_msg("Cannot plate solve before executing PSF Photometry; execute 'Photometry->Iterative PSF Photometry' first.")
            self.console_msg("Ready")
            return

        if self.wcs_header.has_celestial:
            #file is already plate solved
            result = askokcancel(title="Image already Plate Solved.", message="Is it OK to Plate Solve this Plate Solved image?")
            if result==False:
                self.console_msg("User cancelling Plate Solving")
                self.console_msg("Ready")
                return

        self.console_msg("Solving via Astrometry.Net...")
        try:
            # key must exist
            if len(self.astrometrynet_key_entry.get()) == 0:
                self.console_msg("Astrometry.net API Key required")
                self.console_msg("Ready")
                return

            # check if a photometry table exists.
            ast = AstrometryNet()
            ast.api_key = self.astrometrynet_key_entry.get()


            #ast.URL = "http://" + self.astrometrynet_entry.get()
            #ast.API_URL = "http://" + self.astrometrynet_entry.get() + "/api"

            sources_df = self.results_tab_df.sort_values("flux_fit", ascending=False)
            image_data = fits.getdata(self.image_file)
            image_width = image_data.shape[1]
            image_height = image_data.shape[0]


            self.wcs_header = ast.solve_from_source_list(sources_df['x_fit'], sources_df['y_fit'],
                                                         image_width, image_height,
                                                         solve_timeout=360, verbose=True)
            if self.wcs_header:
                self.console_msg("Astrometry.Net solution reference point RA: " + 
                                    str(self.wcs_header["CRVAL1"]) + " Dec: " + 
                                    str(self.wcs_header["CRVAL2"]))
                header = header + self.wcs_header
                self.wcs_header = WCS(header)
            else:
                self.console_msg("Astrometry.Net solution NOT FOUND")
            
            self.console_msg("Ready")

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
            

    #
    #   The seven menu-item callbacks below all dispatch to the single
    #   color_photometry() method, parameterized by input_color.
    #
    def BV_multi_color_photometry(self):
        self.color_photometry('BV')

    def VR_multi_color_photometry(self):
        self.color_photometry('VR')

    def VI_multi_color_photometry(self):
        self.color_photometry('VI')

    def BVR_multi_color_photometry(self):
        self.color_photometry('BVR')

    def BVI_multi_color_photometry(self):
        self.color_photometry('BVI')

    def VRI_multi_color_photometry(self):
        self.color_photometry('VRI')

    def BVRI_multi_color_photometry(self):
        self.color_photometry('BVRI')


    #######################################################################################
    #
    # color_photometry
    #
    # Unified multi-color photometry entry point.  Replaces the former
    # two_color_photometry / three_color_photometry / four_color_photometry
    # methods with a single algorithm, one code path, and one aux-report schema.
    #
    # Parameter input_color selects the filter set and the TA equation form for
    # each filter:
    #
    #     input_color   filters         filter equations (A = Type A, B = Type B)
    #     -----------   ---------       ---------------------------------------------
    #     'BV'          (B, V)          B:A(B-V) V:A(B-V)
    #     'VR'          (V, R)          V:A(V-R) R:A(V-R)
    #     'VI'          (V, I)          V:A(V-I) I:B(V->I)
    #     'BVR'         (B, V, R)       B:A(B-V) V:A(B-V) R:B(V->R)
    #     'BVI'         (B, V, I)       B:A(B-V) V:A(B-V) I:A(V-I)
    #     'VRI'         (V, R, I)       V:A(V-I) R:A(V-I) I:B(V->I)
    #     'BVRI'        (B, V, R, I)    B:A(B-V) V:A(B-V) R:A(V-I) I:A(V-I)
    #
    # Equation forms, matching the AAVSO TransformApplier (TA v2.71) convention:
    #
    #   Type A (standard):   f* = df + Tf_ci * D(CI) + f_comp
    #       where D(CI) = (f1* - f2*) - (f1c - f2c) depends on the unknowns,
    #       so the solution is iterated.  Tf_ci is a MAGNITUDE coefficient.
    #
    #   Type B (direct):     f2* = f1* - (f1c - f2c) - T_color * inst_color
    #       computed from the Type A result of its source filter f1*.
    #       T_color is a COLOR coefficient; inst_color = (f1s-f2s) - (f1c-f2c)
    #       is a constant (not iterated).
    #
    # The coupled system is solved by Gauss-Seidel iteration
    # (see _iterative_transform below).
    #
    # The TC coefficients consumed per input_color are:
    #
    #     coeff      B     V     R     I
    #     BV         Tb_bv Tv_bv
    #     VR               Tv_vr Tr_vr
    #     VI               Tv_vi        Tvi      (Tvi is a color coefficient)
    #     BVR        Tb_bv Tv_bv Tvr             (Tvr is a color coefficient)
    #     BVI        Tb_bv Tv_bv        Ti_vi
    #     VRI              Tv_vi Tr_vi Tvi       (Tvi is a color coefficient)
    #     BVRI       Tb_bv Tv_bv Tr_vi Ti_vi
    #
    #######################################################################################

    def color_photometry(self, input_color):

        #######################################################################
        # _iterative_transform  (helper for color_photometry)
        #
        # Gauss-Seidel iterative transformation for N filters, supporting both
        # Type A (standard) and Type B (color-coefficient) equations.
        #
        # Parameters:
        #   star_IM   - dict {filter: instrumental_mag} for the star
        #   comp_IM   - dict {filter: instrumental_mag} for the comp star
        #   comp_cat  - dict {filter: catalog_mag}      for the comp star
        #   filters   - tuple of filter names in iteration order
        #   eqns      - dict {filter: (type, ci_f1, ci_f2, coeff_value)}
        #
        #       Type A: ('A', ci_f1, ci_f2, mag_coeff)
        #           f* = delta_f + mag_coeff * D(CI) + f_comp
        #           where D(CI) = (est[ci_f1]-est[ci_f2]) - (comp_cat[ci_f1]-comp_cat[ci_f2])
        #
        #       Type B: ('B', ci_f1, ci_f2, color_coeff)
        #           f* = est[ci_f1] - (comp_cat[ci_f1]-comp_cat[ci_f2])
        #                - color_coeff * inst_color(ci_f1, ci_f2)
        #           ci_f1 is the "source" filter (solved by Type A first);
        #           ci_f2 is the Type B filter itself.
        #
        # Returns:
        #   (est, iterations)
        #######################################################################
        def _iterative_transform(star_IM, comp_IM, comp_cat, filters, eqns):
            """Gauss-Seidel iterative N-filter transformation (TA-matched)."""

            MAX_TRANSFORM_ITER = 20
            TRANSFORM_PRECISION = 1e-6

            # Instrumental differentials
            delta = {filt: star_IM[filt] - comp_IM[filt] for filt in filters}

            # Pre-compute catalog color indices and instrumental colors for each
            # (ci_f1, ci_f2) pair referenced by any filter's equation.
            comp_colors = {}  # (f1,f2) -> comp_cat[f1] - comp_cat[f2]
            inst_colors = {}  # (f1,f2) -> (star_IM[f1]-star_IM[f2]) - (comp_IM[f1]-comp_IM[f2])
            for filt in filters:
                _, ci_f1, ci_f2, _ = eqns[filt]
                key = (ci_f1, ci_f2)
                if key not in comp_colors:
                    comp_colors[key] = comp_cat[ci_f1] - comp_cat[ci_f2]
                    inst_colors[key] = (star_IM[ci_f1] - star_IM[ci_f2]) - (comp_IM[ci_f1] - comp_IM[ci_f2])

            # Initial guess: untransformed differential magnitudes
            est = {filt: delta[filt] + comp_cat[filt] for filt in filters}

            # Type A filters must be updated first each pass; Type B then reads
            # from the freshly-updated Type A source filter (ci_f1).
            type_a_filters = [f for f in filters if eqns[f][0] == 'A']
            type_b_filters = [f for f in filters if eqns[f][0] == 'B']

            iterations = 0
            for iterations in range(1, MAX_TRANSFORM_ITER + 1):
                prev = {f: est[f] for f in filters}

                # Type A: standard iterative equation
                for filt in type_a_filters:
                    _, ci_f1, ci_f2, mag_coeff = eqns[filt]
                    delta_color = (est[ci_f1] - est[ci_f2]) - comp_colors[(ci_f1, ci_f2)]
                    est[filt] = delta[filt] + mag_coeff * delta_color + comp_cat[filt]

                # Type B: direct from source filter (ci_f1), using color coeff
                for filt in type_b_filters:
                    _, ci_f1, ci_f2, color_coeff = eqns[filt]
                    est[filt] = est[ci_f1] - comp_colors[(ci_f1, ci_f2)] \
                                - color_coeff * inst_colors[(ci_f1, ci_f2)]

                # Convergence: every filter stable
                if all(abs(est[f] - prev[f]) < TRANSFORM_PRECISION for f in filters):
                    break

            return est, iterations

        # --------------------------------------------------------------------
        # Main entry: color_photometry
        # --------------------------------------------------------------------
        try:
            # ----------------------------------------------------------------
            # Configuration table: filters and equation entries per input_color.
            # Each filter entry: (type, ci_f1, ci_f2,
            #                     tk_coeff_entry, tk_coeff_err_entry,
            #                     coeff_text, coeff_err_text)
            # ----------------------------------------------------------------
            if input_color == 'BV':
                filters = ('B', 'V')
                eqn_entries = {
                    'B': ('A', 'B', 'V', self.tb_bv_entry, self.tb_bv_err_entry, "Tb_bv", "Tb_bvErr"),
                    'V': ('A', 'B', 'V', self.tv_bv_entry, self.tv_bv_err_entry, "Tv_bv", "Tv_bvErr"),
                }
            elif input_color == 'VR':
                filters = ('V', 'R')
                eqn_entries = {
                    'V': ('A', 'V', 'R', self.tv_vr_entry, self.tv_vr_err_entry, "Tv_vr", "Tv_vrErr"),
                    'R': ('A', 'V', 'R', self.tr_vr_entry, self.tr_vr_err_entry, "Tr_vr", "Tr_vrErr"),
                }
            elif input_color == 'VI':
                filters = ('V', 'I')
                eqn_entries = {
                    'V': ('A', 'V', 'I', self.tv_vi_entry, self.tv_vi_err_entry, "Tv_vi", "Tv_viErr"),
                    'I': ('B', 'V', 'I', self.tvi_entry,   self.tvi_err_entry,   "Tvi",   "TviErr"),
                }
            elif input_color == 'BVR':
                filters = ('B', 'V', 'R')
                eqn_entries = {
                    'B': ('A', 'B', 'V', self.tb_bv_entry, self.tb_bv_err_entry, "Tb_bv", "Tb_bvErr"),
                    'V': ('A', 'B', 'V', self.tv_bv_entry, self.tv_bv_err_entry, "Tv_bv", "Tv_bvErr"),
                    'R': ('B', 'V', 'R', self.tvr_entry,   self.tvr_err_entry,   "Tvr",   "TvrErr"),
                }
            elif input_color == 'BVI':
                filters = ('B', 'V', 'I')
                eqn_entries = {
                    'B': ('A', 'B', 'V', self.tb_bv_entry, self.tb_bv_err_entry, "Tb_bv", "Tb_bvErr"),
                    'V': ('A', 'B', 'V', self.tv_bv_entry, self.tv_bv_err_entry, "Tv_bv", "Tv_bvErr"),
                    'I': ('A', 'V', 'I', self.ti_vi_entry, self.ti_vi_err_entry, "Ti_vi", "Ti_viErr"),
                }
            elif input_color == 'VRI':
                filters = ('V', 'R', 'I')
                eqn_entries = {
                    'V': ('A', 'V', 'I', self.tv_vi_entry, self.tv_vi_err_entry, "Tv_vi", "Tv_viErr"),
                    'R': ('A', 'V', 'I', self.tr_vi_entry, self.tr_vi_err_entry, "Tr_vi", "Tr_viErr"),
                    'I': ('B', 'V', 'I', self.tvi_entry,   self.tvi_err_entry,   "Tvi",   "TviErr"),
                }
            elif input_color == 'BVRI':
                filters = ('B', 'V', 'R', 'I')
                eqn_entries = {
                    'B': ('A', 'B', 'V', self.tb_bv_entry, self.tb_bv_err_entry, "Tb_bv", "Tb_bvErr"),
                    'V': ('A', 'B', 'V', self.tv_bv_entry, self.tv_bv_err_entry, "Tv_bv", "Tv_bvErr"),
                    'R': ('A', 'V', 'I', self.tr_vi_entry, self.tr_vi_err_entry, "Tr_vi", "Tr_viErr"),
                    'I': ('A', 'V', 'I', self.ti_vi_entry, self.ti_vi_err_entry, "Ti_vi", "Ti_viErr"),
                }
            else:
                raise Exception("color_photometry: unknown input_color: " + str(input_color))

            n_filters = len(filters)
            method_label = {2: "Two", 3: "Three", 4: "Four"}.get(n_filters, str(n_filters))
            self.console_msg("Performing " + method_label + " Color Photometry (" +
                             input_color + ", Gauss-Seidel iterative method)...")

            # ----------------------------------------------------------------
            # Validate Object Name
            # ----------------------------------------------------------------
            variable_star = self.object_name_entry.get().strip()
            if variable_star is None or len(variable_star) == 0:
                self.console_msg(
                    "Color Photometry requires 'Object Name' to be filled; eg. 'V1117 Her'")
                self.console_msg("Ready")
                return

            # ----------------------------------------------------------------
            # Ask for one CSV file per filter
            # ----------------------------------------------------------------
            options = {}
            options['defaultextension'] = '.csv'
            options['filetypes'] = [('CSV', '.csv')]

            filter_data = {}  # filter -> DataFrame
            for filt in filters:
                options['title'] = 'Choose a file for filter ' + filt
                file_name = fd.askopenfilename(**options)
                if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
                    self.console_msg("Loading filter " + filt + " from " + str(file_name))
                    filter_data[filt] = pd.read_csv(str(file_name), dtype={'check_star': bool})
                else:
                    return
                if "label" not in filter_data[filt]:
                    self.console_msg("Cannot proceed; run 'Photometry->Get Comparison Stars' for " + filt + " first.")
                    return
            last_file_name = file_name  # remembered for the output directory

            # ----------------------------------------------------------------
            # Parse transformation coefficients
            # ----------------------------------------------------------------
            try:
                eqns = {}      # filter -> (type, ci_f1, ci_f2, coeff_value)
                eqn_info = {}  # filter -> (type, ci_f1, ci_f2, coeff_val, coeff_err,
                               #            coeff_text, coeff_err_text)  -- for aux report
                for filt in filters:
                    eqn_type, ci_f1, ci_f2, coeff_entry, coeff_err_entry, coeff_text, coeff_err_text = eqn_entries[filt]
                    coeff_val = float(coeff_entry.get().strip())
                    coeff_err = float(coeff_err_entry.get().strip())
                    eqns[filt] = (eqn_type, ci_f1, ci_f2, coeff_val)
                    eqn_info[filt] = (eqn_type, ci_f1, ci_f2, coeff_val, coeff_err, coeff_text, coeff_err_text)
            except:
                self.console_msg("Cannot proceed with Color Photometry: Missing or non-numeric transform coefficient(s)")
                return

            # ----------------------------------------------------------------
            # Check-star validation (must exist in all filter tables)
            # ----------------------------------------------------------------
            check_star_in_setting = self.object_kref_entry.get().strip()
            check_star_in_setting_with_prefix = __label_prefix__ + check_star_in_setting

            for filt in filters:
                if not filter_data[filt]["label"].isin([check_star_in_setting_with_prefix]).any():
                    self.console_msg("Settings check star: " + check_star_in_setting +
                                      " not found in " + filt + "; select another check star")
                    # Suggest comps common to all filter tables
                    common_labels = filter_data[filters[0]]["label"]
                    for other_filt in filters[1:]:
                        common_labels = common_labels[common_labels.isin(filter_data[other_filt]["label"])]
                    common_labels = common_labels[common_labels != "%"]
                    list_good_comps = ", ".join(
                        comp.split(maxsplit=1)[-1] for comp in common_labels.astype(str))
                    self.console_msg("Pick from following: " + list_good_comps)
                    return

            # ----------------------------------------------------------------
            # Variable-star validation (must exist in all filter tables)
            # ----------------------------------------------------------------
            for filt in filters:
                if not filter_data[filt]["vsx_id"].isin([variable_star]).any():
                    self.console_msg("Settings Object Name: " + variable_star +
                                      " not found in " + filt + "; select " +
                                      variable_star + " in " + filt +
                                      " image with mouse; make sure you click on the same star in each image")
                    return

            # ----------------------------------------------------------------
            # Reset and set the check_star flag in each filter table
            # ----------------------------------------------------------------
            check_star_rows = {}  # filter -> row Series
            for filt in filters:
                df = filter_data[filt]
                if df["check_star"].isin([True]).any():
                    df.loc[df[df["check_star"] == True].index, "check_star"] = False
                df.loc[df[df["label"] == check_star_in_setting_with_prefix].index, "check_star"] = True
                check_star_rows[filt] = df[df["label"] == check_star_in_setting_with_prefix].iloc[0]

            check_star_label = check_star_in_setting
            self.console_msg("Using check star " + check_star_label)

            # Extract per-filter check-star instrumental and catalog mags
            check_IM  = {filt: check_star_rows[filt]["inst_mag"]  for filt in filters}
            check_cat = {filt: check_star_rows[filt]["match_mag"] for filt in filters}

            # Find the variable star in each filter table
            var_star_rows = {}
            var_IM = {}
            for filt in filters:
                var_star_rows[filt] = filter_data[filt][filter_data[filt]["vsx_id"] == variable_star].iloc[0]
                var_IM[filt] = var_star_rows[filt]["inst_mag"]
            var_star_label = var_star_rows[filters[0]]["vsx_id"]

            # ----------------------------------------------------------------
            # Date-obs and airmass per filter.
            #   raw_jd[filt]   -> original shutter-open JD (aux "JD" column)
            #
            #   If "date-end" exists 
            #       date_obs[filt] -> (raw JD + date-end)/2    (aux "Date-Obs" column)
            #   else
            #       date_obs[filt] -> raw JD + EXPOSURE/2    (aux "Date-Obs" column)
            #       
            # ----------------------------------------------------------------
            date_obs = {}
            amass    = {}
            raw_jd   = {}
            for filt in filters:
                row = filter_data[filt][filter_data[filt]["check_star"] == True].iloc[0]
                raw_jd[filt] = row["date-obs"]
                if 'date-end' in row:
                    date_obs[filt] = (row["date-obs"] + row["date-end"])/2
                else:
                    half_exp = row["exposure"] / 2
                    _obs = Time(row["date-obs"], format='jd') + TimeDelta(half_exp, format='sec')
                    date_obs[filt] = _obs.jd

                amass[filt]    = row["AMASS"]

            # ----------------------------------------------------------------
            # Determine unique color indices referenced by the equation set.
            # These drive the diagnostic delta_* columns in the result tables.
            # ----------------------------------------------------------------
            ci_set = []
            ci_seen = set()
            for filt in filters:
                ci = (eqns[filt][1], eqns[filt][2])
                if ci not in ci_seen:
                    ci_set.append(ci)
                    ci_seen.add(ci)

            # ----------------------------------------------------------------
            # Build result-table column layout
            # ----------------------------------------------------------------
            result_columns = ["type", "name", "comp"]
            result_columns += ["IM" + f for f in filters]   # instrumental mag of comp
            result_columns += list(filters)                  # catalog mag of comp
            for ci_f1, ci_f2 in ci_set:
                result_columns.append("delta_" + ci_f1.lower() + "_minus_" + ci_f2.lower())  # inst
                result_columns.append("delta_" + ci_f1 + "_minus_" + ci_f2)                  # std
            result_columns += [f + "_star" for f in filters]  # transformed mags
            result_columns.append("iterations")

            check_columns = result_columns + ["outlier"]
            var_columns   = result_columns[:]

            result_check_star = pd.DataFrame(columns=check_columns)
            result_var_star   = pd.DataFrame(columns=var_columns)

            # ----------------------------------------------------------------
            # Loop through the selected comp stars
            # ----------------------------------------------------------------
            sel_comps_to_use = list(set(
                comp.strip() for comp in self.object_sel_comp_entry.get().split(',')))
            sel_comps = [(c, __label_prefix__ + c) for c in sel_comps_to_use]

            for comp_no_prefix, comp in sel_comps:

                # Skip the check star and any commented-out entries
                if comp_no_prefix == check_star_label:
                    continue
                if not comp_no_prefix or not comp_no_prefix[0].isdigit():
                    continue

                # Selected comp must be present in every filter table
                comp_star = {}
                skip = False
                for filt in filters:
                    if comp in filter_data[filt]["label"].values:
                        comp_star[filt] = filter_data[filt][filter_data[filt]["label"] == comp].iloc[0]
                    else:
                        self.console_msg("Comp star: " + comp_no_prefix + " not in " + filt + " table")
                        skip = True
                        break
                if skip:
                    continue

                comp_IM      = {f: comp_star[f]["inst_mag"]         for f in filters}
                comp_cat_mag = {f: float(comp_star[f]["match_mag"]) for f in filters}

                # ------------------------------------------------------------
                # Solve Gauss-Seidel for a star and build the result-row dict
                # ------------------------------------------------------------
                def _solve_and_build_row(star_IM, row_type, row_name):
                    est, n_iter = _iterative_transform(
                        star_IM, comp_IM, comp_cat_mag, filters, eqns)
                    row = {"type": row_type, "name": row_name, "comp": comp_no_prefix}
                    for f in filters:
                        row["IM" + f] = comp_IM[f]
                    for f in filters:
                        row[f] = comp_star[f]["match_mag"]
                    for ci_f1, ci_f2 in ci_set:
                        dinst = (star_IM[ci_f1] - star_IM[ci_f2]) - (comp_IM[ci_f1] - comp_IM[ci_f2])
                        dstd  = (est[ci_f1]     - est[ci_f2])     - (comp_cat_mag[ci_f1] - comp_cat_mag[ci_f2])
                        row["delta_" + ci_f1.lower() + "_minus_" + ci_f2.lower()] = dinst
                        row["delta_" + ci_f1 + "_minus_" + ci_f2] = dstd
                    for f in filters:
                        row[f + "_star"] = est[f]
                    row["iterations"] = n_iter
                    return row

                # Check star
                row = _solve_and_build_row(check_IM, "check", check_star_label)
                row["outlier"] = ''
                result_check_star.loc[len(result_check_star)] = row

                # Variable star
                row = _solve_and_build_row(var_IM, "var", var_star_label)
                result_var_star.loc[len(result_var_star)] = row

            # ----------------------------------------------------------------
            # Means and standard deviations (fallback to 1/SNR for a single comp)
            # ----------------------------------------------------------------
            means_check = {}
            stds_check  = {}
            means_var   = {}
            stds_var    = {}
            for filt in filters:
                col = filt + "_star"
                means_check[filt] = result_check_star[col].mean()
                means_var[filt]   = result_var_star[col].mean()
                if len(result_check_star[col]) > 1:
                    stds_check[filt] = result_check_star[col].std()
                    stds_var[filt]   = result_var_star[col].std()
                else:
                    stds_check[filt] = 1 / (check_star_rows[filt]['flux_fit'] / self.image_bkg_value)
                    stds_var[filt]   = 1 / (var_star_rows[filt]['flux_fit']   / self.image_bkg_value)

            # ----------------------------------------------------------------
            # IQR outlier detection (check-star table only)
            # ----------------------------------------------------------------
            iqr_limits = {}  # filt -> (lower, upper)   (kept for console output)
            for filt in filters:
                col = filt + "_star"
                q3, q1 = np.percentile(result_check_star[col], [75, 25])
                iqr = q3 - q1
                upper = q3 + (iqr * 1.5)
                lower = q1 - (iqr * 1.5)
                iqr_limits[filt] = (lower, upper)
                if (result_check_star[col] < lower).any():
                    result_check_star.loc[result_check_star[col] < lower, "outlier"] = "<--OUTLIER"
                if (result_check_star[col] > upper).any():
                    result_check_star.loc[result_check_star[col] > upper, "outlier"] = "<--OUTLIER"

            # ----------------------------------------------------------------
            # Console output
            # ----------------------------------------------------------------
            self.console_msg("\n")

            # Check-star summary
            check_msg = "Check Star Estimates using check star: " + str(check_star_label)
            for filt in filters:
                check_msg += " (" + filt + ": " + format(float(check_cat[filt]), ' >6.3f') + ")"
            check_msg += ";    (qfit--->"
            for filt in filters:
                check_msg += " " + filt + ": " + format(check_star_rows[filt]["qfit"], ' >5.3f') + ";"
            check_msg = check_msg.rstrip(';') + ')'
            check_msg += '\n' + result_check_star.sort_values(by="name").to_string() + '\n'
            ave_line = ""
            std_line = ""
            for filt in filters:
                ave_line += filt + "* Ave: " + format(means_check[filt], ' >6.3f') + "  "
                std_line += filt + "* Std: " + format(stds_check[filt], ' >6.3f') + "  "
            check_msg += ave_line.rstrip().rjust(137) + '\n'
            check_msg += std_line.rstrip().rjust(137)
            self.console_msg(check_msg)

            # IQR limits
            for filt in filters:
                lower, upper = iqr_limits[filt]
                self.console_msg(("Check Star IQR limit for " + filt + "*: " +
                                  format(lower, ' >6.3f') + ';' + format(upper, ' >6.3f')).rjust(123))
            self.console_msg('\n')

            # Variable-star summary
            var_msg = "Variable Star Estimates of Var: " + var_star_label
            var_msg += ";    (qfit--->"
            for filt in filters:
                var_msg += " " + filt + ": " + format(var_star_rows[filt]["qfit"], ' >5.3f') + ";"
            var_msg = var_msg.rstrip(';') + ')'
            var_msg += '\n' + result_var_star.sort_values(by="name").to_string() + '\n'
            ave_line = ""
            std_line = ""
            for filt in filters:
                ave_line += filt + "* Ave: " + format(means_var[filt], ' >6.3f') + "  "
                std_line += filt + "* Std: " + format(stds_var[filt], ' >6.3f') + "  "
            var_msg += ave_line.rstrip().rjust(137) + '\n'
            var_msg += std_line.rstrip().rjust(137)
            self.console_msg(var_msg)

            # ----------------------------------------------------------------
            # Aux report (generic schema for the AAVSO notes section)
            # ----------------------------------------------------------------
            aux_columns = ["color", "JD", "KMAGS", "KMAGINS", "KREFMAG",
                           "VMAGINS", "Date-Obs", "KNAME", "AMASS",
                           "coeff_val", "coeff_err", "coeff_text", "coeff_err_text"]
            result_aux_report = pd.DataFrame(columns=aux_columns)
            for filt in filters:
                eqn_type, ci_f1, ci_f2, coeff_val, coeff_err, coeff_text, coeff_err_text = eqn_info[filt]
                result_aux_report.loc[len(result_aux_report)] = {
                    "color":          filt,
                    "JD":             raw_jd[filt],         # orig
                    "KMAGS":          means_check[filt],
                    "KMAGINS":        check_IM[filt],
                    "KREFMAG":        check_cat[filt],
                    "VMAGINS":        var_IM[filt],
                    "Date-Obs":       date_obs[filt],       # orig + EXPOSURE/2
                    "KNAME":          check_star_label,
                    "AMASS":          amass[filt],
                    "coeff_val":      coeff_val,
                    "coeff_err":      coeff_err,
                    "coeff_text":     coeff_text,
                    "coeff_err_text": coeff_err_text,
                }

            # ----------------------------------------------------------------
            # Save master report
            # ----------------------------------------------------------------
            master_frames = [result_check_star, result_var_star, result_aux_report]
            master_report = pd.concat(master_frames, keys=['check', 'var', 'aux'])
            dir_path = os.path.dirname(os.path.realpath(last_file_name)) + "\\"
            report_name = self.object_name_entry.get() + "-" + input_color + "-Master-Report.csv"
            master_report.to_csv(dir_path + report_name, index=False)
            self.console_msg("Master Report saved to " + str(dir_path + report_name))
            self.console_msg("Ready")

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) + " " + str(e), level=logging.ERROR)


    #######################################################################################
    #
    # get_comparison_stars
    # 
    # This is callback after Photometry__>Get Comparison Stars is clicked
    # 
    #
    #  APASS DR9 uses designation with 'B" (e.g., 123B_2)  mag 12.3
    #  APASS DR10 uses designation with 'A" (e.g., 123A_2)
    #  APASS DR11 uses designation with 'C" (e.g., 123C_2)
    # 
    #
    ########################################################################################
    
    def get_comparison_stars(self):
        global image_width, image_height
        comp_stars_found = [] #init
        try:
            if not os.path.exists(self.image_file + ".csv"):
                # No photometry exists!
                self.console_msg("Image needs photometry first, run 'Photometry->Iterative PSF Photometry' first")
                self.console_msg("Ready")
                return

            self.filter = self.filter_entry.get()

            frame_center = self.wcs_header.pixel_to_world(
                int(image_width / 2), (image_height / 2))

            if not self.wcs_header.has_celestial:
                #file needs to be plate solved first
                self.console_msg("Image needs to be plate solved first! execute 'Photometry->Solve Image'")
                self.console_msg("Ready")
                return

            """
            Get Minimum SNR value from settings amd if not a number set to default 20.

            Later a found comp or check star gets its SNR evaluated and if it is less than min_snr
            that comp or check star is not used, and not displayed.

            """
            min_snr_to_use = self.object_min_snr_entry.get().strip()
            if is_number(min_snr_to_use):
                min_snr = float(min_snr_to_use)
            else:
                min_snr = 20
                self.set_entry_text(self.object_min_snr_entry, "20")
              

            frame_center_coordinates = SkyCoord(ra=frame_center.ra, dec=frame_center.dec)
            frame_edge = self.wcs_header.pixel_to_world(int(image_width), (image_height / 2))
            frame_edge_coordinates = SkyCoord(ra=frame_edge.ra, dec=frame_edge.dec)
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

            ## Init 
            using_aavso_catalog = False
            using_apass_dr9 = False
            using_apass_dr10 = False
    
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
                mag_error_column_name = "Mag Error"
            
            elif catalog_selection == "APASS DR10":
                using_apass_dr10 = True
                
                #APASS DR10 (not in Vizier but on a AAVSO server)
                self.console_msg("Getting APASS DR10 Comparison Stars...")
                comparison_stars = \
                    self.APASS_DR10_get_comparison_stars(frame_center, filter=str(self.filter_entry.get()),
                    #divide fov_horizontal by 60 because in arcminutes add 10%; api is limited to 2 degrees
                    field_of_view=fov_horizontal/60*1.1,
                    maglimit=self.max_ensemble_magnitude_entry.get())
                
                if type(comparison_stars) == int:
                    self.console_msg("Get Comparison Stars failed (APASS DR10)", level=logging.ERROR)
                    self.console_msg("Ready")
                    return

                ra_column_name = "RA"
                dec_column_name = "Dec"
                mag_column_name = "Mag"
                mag_error_column_name = "Mag Error"
            
            else:

                # Default value
                # Changes for some catalogs
                mag_column_name = self.filter + "mag"
                mag_error_column_name = self.filter + "mag error"  #<--- needs to be verified/fixed

                #Use VizierR catalog
                if catalog_selection == "APASS DR9":
                    catalog = "II/336"
                    ra_column_name = "RAJ2000"
                    dec_column_name = "DEJ2000"
                    catalog_columns = ["RAJ2000", "DEJ2000", "Vmag", "Bmag", "r'mag", "i'mag"]
                    DR9_Mag_name = {
                        "V" : "Vmag",
                        "B" : "Bmag",
                        "R" : "r_mag", # Sloan
                        "I" : "i_mag", # Sloan
                        }
                    mag_column_name = DR9_Mag_name[self.filter]
                    #mag_error_column_name = DR9_Mag_Error_name[self.filter]   #<--- needs to be fixed
                    using_apass_dr9 = True
    
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
                    sourceId_column_name = "DR2Name"
                    ra_column_name = "RA_ICRS"
                    dec_column_name = "DE_ICRS"
                    catalog_columns = ["DR2Name", "RA_ICRS", "DE_ICRS", "Plx", "Gmag", "RPmag", "BPmag", "Lum"]
    
                else:
                    self.console_msg("UNEXPECTED ERROR UNKOWN CATALOG SELECTION!!!!!", level=logging.ERROR)
                    return

                self.console_msg(
                    "Inquiring VizieR Catalog " + catalog_selection + ", center α δ " +
                        frame_center.to_string("hmsdms", precision=2) +
                        ", radius " + str(frame_radius))
    
                # [0] implies I/345/gaia2; [1] implies I/345/varres
                comparison_stars = Vizier(catalog=catalog,
                                          columns=catalog_columns,
                                          row_limit=-1).query_region(frame_center, radius=frame_radius)[0] 


                # print(comparison_stars)
                if mag_column_name not in comparison_stars.colnames:
                    self.console_msg(
                        "Catalog " + self.catalog_stringvar.get() + " does not list " + self.filter + " magnitudes.")
                    return

            if len(comparison_stars) == 0:
                self.console_msg(
                    "NO Comparison stars found in the field; make sure filter, chart Id (optional), and comp names are correct.")
                self.console_msg("Ready")
                return
            else:                
                self.console_msg(
                    "Found " + str(len(comparison_stars)) + " comp stars in the field.")


            self.console_msg("Matching image to catalog...")
            matching_radius = float(
                self.matching_radius_entry.get()) / 3600.0  # arcsec to degrees
            
            if using_aavso_catalog:
                catalog_comparison = SkyCoord(
                    comparison_stars[ra_column_name],
                    comparison_stars[dec_column_name],
                    unit=(u.hourangle, u.deg))
                
            elif using_apass_dr10:
                catalog_comparison = SkyCoord(
                    comparison_stars[ra_column_name],
                    comparison_stars[dec_column_name],
                    unit=(u.deg, u.deg))
            else:
                catalog_comparison = SkyCoord(
                    comparison_stars[ra_column_name],
                    comparison_stars[dec_column_name])
                #here not using aavso cat so we need check star info
                check_star_to_use = self.object_kref_entry.get().strip()
                
            for index, row in self.results_tab_df.iterrows():
                photometry_star_coordinates = SkyCoord(
                    ra=row["ra_fit"] * u.deg, dec=row["dec_fit"] * u.deg)
                match_index, d2d_match, d3d_match = photometry_star_coordinates.match_to_catalog_sky(
                    catalog_comparison)

                #Catalog specific hacks
                if using_aavso_catalog:
                    match_id = comparison_stars.iloc[match_index]["AUID"]

                if using_apass_dr10:
                    match_id = comparison_stars.iloc[match_index]["Label"] # DR10 doesn't have AUID

                if using_aavso_catalog or using_apass_dr10:
                    match_label = comparison_stars.iloc[match_index]["Label"]
                    match_ra = catalog_comparison[match_index].ra.degree
                    match_dec= catalog_comparison[match_index].dec.degree
                    match_mag = comparison_stars.iloc[match_index][mag_column_name]
                    match_mag_error = comparison_stars.iloc[match_index][mag_error_column_name]
                    match_is_check = comparison_stars.iloc[match_index]["Check Star"]

                elif using_apass_dr9:
                    match_id = "RA" + format(comparison_stars[match_index][ra_column_name], '.2f') +\
                                     "+" + \
                                    "DE" + format(comparison_stars[match_index][dec_column_name], '.2f')
                    # Create a label similar to APASS' 
                    match_label = format(comparison_stars[match_index][mag_column_name] * 10, '3.0f') + "B"

                    match_ra = comparison_stars[match_index][ra_column_name]
                    match_dec = comparison_stars[match_index][dec_column_name]
                    match_mag = comparison_stars[match_index][mag_column_name]
                    match_mag_error = comparison_stars[match_index][mag_error_column_name]
                    match_is_check = (match_label == check_star_to_use)

                else:
                    match_id = comparison_stars[match_index][0]
                    match_label = comparison_stars[match_index][sourceId_column_name]
                    match_ra = comparison_stars[match_index][ra_column_name]
                    match_dec = comparison_stars[match_index][dec_column_name]
                    match_mag = comparison_stars[match_index][mag_column_name]
                    match_mag_error = comparison_stars[match_index][mag_error_column_name]
                    
                match_coordinates = SkyCoord(ra=match_ra * u.deg, dec=match_dec * u.deg)
                separation = photometry_star_coordinates.separation(match_coordinates)

                # Look for minimu separation
                if separation < (matching_radius * u.deg):
                    # Check for minimum SNR
                    comp_snr = self.results_tab_df.loc[index, "flux_fit"] / self.image_bkg_value
                    if comp_snr < min_snr:
                        self.console_msg("Comp star: " + str(match_label) + 
                                        " SNR too low " + 
                                        " (" + format(comp_snr, '.1f') + 
                                        " < " + str(min_snr) + "); skipping")
                        #Here if low snr
                        self.results_tab_df.loc[index, "match_id"] = ""
                        self.results_tab_df.loc[index, "label"] = __empty_cell__ # prevent nan
                        self.results_tab_df.loc[index, "match_ra"] = ""
                        self.results_tab_df.loc[index, "match_dec"] = ""
                        self.results_tab_df.loc[index, "match_mag"] = ""
                        self.results_tab_df.loc[index, "check_star"] = False # all values in this column must be of boolean 
                    else:
                        #For a given comp star and matching_radius, there could be 
                        #more than one match. Choose the entry with the lowest qfit (Quality of Fit Score)
                        if "label" in self.results_tab_df:
                            already_gotten = self.results_tab_df.loc[self.results_tab_df["label"] == (__label_prefix__ + str(match_label))]    
                            if not already_gotten.empty:
                                # Here if we already got this comp star 
                                # Choose the entry with the lowest qfit (Quality of Fit Score)
                                if already_gotten['qfit'].iloc[0] < self.results_tab_df.iloc[index]['qfit']:
                                    # The the new match is NOT a better fit. ignore it
                                    continue
                                else:
                                    # The new match IS a better fit, move the label to it
                                    # Erase the old
                                    self.results_tab_df.loc[already_gotten.index, "match_id"] = ""
                                    self.results_tab_df.loc[already_gotten.index, "label"] = __empty_cell__ 
                                    self.results_tab_df.loc[already_gotten.index, "match_ra"] = ""
                                    self.results_tab_df.loc[already_gotten.index, "match_dec"] = ""
                                    self.results_tab_df.loc[already_gotten.index, "match_mag"] = ""
                                    self.results_tab_df.loc[already_gotten.index, "match_mag_error"] = ""
                                    self.results_tab_df.loc[already_gotten.index, "check_star"] = False 

                        #Found a match within matching_radius
                        self.results_tab_df.loc[index, "match_id"] = \
                            str(self.catalog_stringvar.get()) + \
                                " " + str(match_id)
                        self.results_tab_df.loc[index, "label"] = __label_prefix__ + str(match_label)
                        self.results_tab_df.loc[index, "match_ra"] = match_ra
                        self.results_tab_df.loc[index, "match_dec"] = match_dec
                        self.results_tab_df.loc[index, "match_mag"] = match_mag
                        self.results_tab_df.loc[index, "match_mag_error"] = match_mag_error

                        if using_aavso_catalog or using_apass_dr9 or using_apass_dr10:
                            self.results_tab_df.loc[index, "check_star"] = match_is_check
                            #record comp stars used for console if AAVSO comp stars
                            comp_stars_found.append((str(match_label), match_is_check, self.results_tab_df.loc[index, "qfit"]))
                        else:
                            self.results_tab_df.loc[index, "check_star"] = False
                else:
                    #Here if separation >= matching_radius
                    self.results_tab_df.loc[index, "match_id"] = ""
                    self.results_tab_df.loc[index, "label"] = __empty_cell__ # prevent nan
                    self.results_tab_df.loc[index, "match_ra"] = ""
                    self.results_tab_df.loc[index, "match_dec"] = ""
                    self.results_tab_df.loc[index, "match_mag"] = ""
                    self.results_tab_df.loc[index, "check_star"] = False # all values in this column must be of boolean 
            
            self.console_msg("Inquiring VizieR (B/vsx/vsx) for VSX variables in the field...")


            # B/vsx : AAVSO International Variable Star Index VSX 
            # B/vsx/vsx : Variable Star indeX,
            vsx_result = Vizier(catalog="B/vsx/vsx", row_limit=-1).query_region(frame_center,
                                                                                 radius=frame_radius)


            # Prep object_name
            object_name = self.object_name_entry.get().strip()
            object_name_exist =  object_name != None and len(object_name) > 0

            #Look for any and all VSX stars
            #first init this flag 
            found_user_object = False
            if len(vsx_result) > 0:

                vsx_stars = vsx_result[0]
                self.console_msg(
                    "Found " + str(len(vsx_stars)) + " VSX sources in the field. Matching...")
                catalog_vsx = SkyCoord(
                    vsx_stars["RAJ2000"], vsx_stars["DEJ2000"])
                for index, row in self.results_tab_df.iterrows():
                    # First check if this entry is a comp star. Comp stars have higher 
                    # priority over any vsx names ( a comp could also be a GAIA entry ).
                    # So ignore any vsx matching to a comp star entry
                    if self.results_tab_df.loc[index, "label"] != __empty_cell__:
                        continue
                    photometry_star_coordinates = SkyCoord(
                        ra=row["ra_fit"] * u.deg, dec=row["dec_fit"] * u.deg)
                    match_index, d2d_match, d3d_match = photometry_star_coordinates.match_to_catalog_sky(
                        catalog_vsx)
                  
                    match_id = vsx_stars[match_index]["Name"]
                    match_ra = vsx_stars[match_index]["RAJ2000"]
                    match_dec = vsx_stars[match_index]["DEJ2000"]
                    match_coordinates = SkyCoord(ra=match_ra * u.deg, dec=match_dec * u.deg)
                    alpha_delta = [] # used to keep alpha and delta, may end up in settings
                    alpha_delta = match_coordinates.to_string("hmsdms", precision=2).split()

                    separation = photometry_star_coordinates.separation(match_coordinates)
                    
                    if separation < (matching_radius * u.deg):
                        if "vsx_id" in self.results_tab_df:
                            already_gotten = self.results_tab_df.loc[self.results_tab_df["vsx_id"] == match_id]    
                            if not already_gotten.empty:
                                # Here if we already got this vsx
                                # Choose the entry with the lowest qfit (Quality of Fit Score)
                                if already_gotten['qfit'].iloc[0] < self.results_tab_df.iloc[index]['qfit']:
                                    # The the new match is NOT a better fit. ignore it
                                    continue
                                else:
                                    # The new match IS a better fit, move the label to it
                                    # Erase the old
                                    self.results_tab_df.loc[already_gotten.index, "match_id"] = ""
                                    self.results_tab_df.loc[already_gotten.index, "label"] = __empty_cell__ 
                                    self.results_tab_df.loc[already_gotten.index, "match_ra"] = ""
                                    self.results_tab_df.loc[already_gotten.index, "match_dec"] = ""
                                    self.results_tab_df.loc[already_gotten.index, "match_mag"] = ""
                                    self.results_tab_df.loc[already_gotten.index, "check_star"] = False 

                        # Found a match within matching_radius
                        self.results_tab_df.loc[index, "vsx_id"] = str(match_id)
                        self.results_tab_df.loc[index, "RAJ2000"] = str(match_ra)
                        self.results_tab_df.loc[index, "DEJ2000"] = str(match_dec)
                        self.results_tab_df.loc[index, "separation"] = str(separation)
                        self.console_msg("Match VSX source:" +\
                                          " (qfit:" + format(self.results_tab_df.loc[index, "qfit"], '0.4f') +") " +\
                                          str(match_id))
                        # See if user has specified an Object Name and if so, then 
                        # if there is a match and (separation < matching_radius), then 
                        # insert alpha and delta into Settings.
                        if object_name_exist and object_name == str(match_id):
                            self.set_entry_text(self.object_name_alpha_entry, alpha_delta[0])
                            self.set_entry_text(self.object_name_delta_entry, alpha_delta[1])
                            found_user_object_qfit = self.results_tab_df.loc[index, "qfit"]
                            found_user_object = True

                    else:
                        self.results_tab_df.loc[index, "vsx_id"] = __empty_cell__ # prevent nan
            else:
                self.console_msg("Found no VSX sources in the field.")

            # If necessary search using Object Alpha and Delta
            object_alpha = self.object_name_alpha_entry.get().strip()
            object_alpha_exist =  object_alpha != None and len(object_alpha) > 0
            object_delta = self.object_name_delta_entry.get().strip()
            object_delta_exist =  object_delta != None and len(object_delta) > 0
            if object_name_exist and object_alpha_exist and object_delta_exist and not found_user_object:
                user_object = SkyCoord(object_alpha, object_delta, frame='icrs')
                for index, row in self.results_tab_df.iterrows():
                    photometry_star_coordinates = SkyCoord(ra=row["ra_fit"] * u.deg, dec=row["dec_fit"] * u.deg, frame='icrs')
                    separation = photometry_star_coordinates.separation(user_object)
                    if separation < (matching_radius * u.deg):
                        # We found it 
                        self.results_tab_df.loc[index, "vsx_id"] = object_name
                        self.results_tab_df.loc[index, "RAJ2000"] = user_object.ra.deg
                        self.results_tab_df.loc[index, "DEJ2000"] = user_object.dec.deg
                        found_user_object_qfit = self.results_tab_df.loc[index, "qfit"]
                        found_user_object = True
                        self.console_msg("Found Object using α:" +  object_alpha + "; δ:" + object_delta)
                        break

            if object_name_exist and not found_user_object:
                self.console_msg("DID NOT FIND Object Name:" + object_name + "!!!!")
            else:
                # Reprint at the bottom of list for convenience
                self.console_msg("Match VSX source:" +\
                                " (qfit:" + format(found_user_object_qfit, '0.4f') +") " + object_name)

            if "label" in self.results_tab_df:
                # Now clean up results_tab_df pf rows that were never accessed. 
                # Remove NaN values
                self.results_tab_df = self.results_tab_df.dropna(subset=['label'])  
                # Remove rows with unknown label
                self.results_tab_df.drop(self.results_tab_df[(self.results_tab_df['label'] != "%") &
                                                        (~self.results_tab_df['label'].str.contains(r'^'+__label_prefix__, regex=True))].index, inplace=True)


            self.results_tab_df.to_csv(self.image_file + ".csv", index=False)
            self.console_msg("Photometry table saved to " + str(self.image_file + ".csv"))

            if using_aavso_catalog:
                comp_list = ''
                #output comp_stars_found
                found_check = False #init
                for comp in comp_stars_found:
                    (label, ischeck, qfit) = comp
                    found_check |= ischeck == True
                    comp_list += str(label) + " (" + format(qfit, '0.4f') + "); " 
                self.console_msg("AAVSO comp stars (qfit) found: " + comp_list)
                    
                check_star = self.object_kref_entry.get().strip()
                if not found_check and check_star != '':
                    self.console_msg("WARNING!! The requested check_star: " 
                                     + check_star + " was NOT FOUND in image! Choose another.", level=logging.WARNING)

            # No need to show ePSF data
            self.ePSF_samples_plotted = False
            
            self.display_image()
            self.console_msg("Ready")
        
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)

    def set_entry_text(self, entry, text):
        entry.delete(0, tk.END)
        entry.insert(0, text)

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
        self.display_image()

    def update_histogram_high(self, value):
        self.histogram_slider_high = int(value)
        self.display_image()


    ############################################################################################
    #
    # save_settings_as
    #
    # function callback for Save As.. button in Settings Window
    #
    ############################################################################################

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
                self.settings_filename = str(file_name.name)
                
                '''
                Use the valid_parameter list which contains the official list
                of user parameters
                '''

                # But first check if auto_behavior is set (CCD Filter from FITS), if set
                # then dont save any Fits data, erase filter_entry, put back later
                if self.auto_behavior.get():
                    filter_entry = self.filter_entry.get()
                    self.filter_entry.config(state='normal')
                    self.set_entry_text(self.filter_entry, "")
                    self.filter_entry.config(state='disable')

                settings = {}
                mao_parameters = list(self.valid_parameter_list)
                for param in mao_parameters:
                    settings.update({param : self.valid_parameter_list[param].get()})

                with open(str(file_name.name), 'w') as f:
                    w = csv.DictWriter(f, settings.keys())
                    w.writeheader()
                    w.writerow(settings)

                #special case for settings filename
                self.set_entry_text(self.settings_filename_entry, self.settings_filename)
                self.settings_filename_entry.xview_scroll(len(self.settings_filename), tk.UNITS)
                self.console_msg("Saved.")

                # Put back that filter_entry if auto_behavior on (the From FITS checkbox in settings)
                if self.auto_behavior.get():
                    self.filter_entry.config(state='normal')
                    self.set_entry_text(self.filter_entry, filter_entry)
                    self.filter_entry.config(state='disable')
            else:
                return
        
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
        finally:
            # bring the settings window back to the front                
            self.es_top.lift()


    ############################################################################################
    #
    # load_config
    #
    # load configuration file (./.config)
    #
    ############################################################################################
    def load_config(self):
        try:
            # Check if .config file has been created
            # if there then load existing config parameters
            if os.path.exists(self.config_file):
                configs = {}
                with open(str(self.config_file)) as f:
                    r = csv.DictReader(f)
                    for row in r:
                        # dict from OrderedDict, required for Python < 3.8 as DictReader behavior changed
                        row = dict(row)
                        # append configs dictionary with the read row
                        configs.update(row)
                    for key in configs:
                        if key not in self.valid_config_list:
                            continue
                        if hasattr(self, key):  # Check if the attribute exists
                            setattr(self, key, configs[key])
                    pass
        
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            #Logger not up yet do print to console
            print("Exception at line no: " + str(exc_tb.tb_lineno) )
            os._exit(1)

    ############################################################################################
    #
    # load_settings
    #    
    # command for "Load..." button in Settings window
    #
    ############################################################################################
    
    def load_settings(self):
        try:
            options = {}
            options['defaultextension'] = '.txt'
            options['filetypes'] = [('TXT', '.txt')]
            #options['initialfile'] = ''
            options['title'] = 'Load MAOPhot settings...'

            file_name = fd.askopenfilename(**options)
            if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
                self.open_settings(file_name)

                self.console_msg("Loaded settings from " + str(file_name))
                self.console_msg("Ready")
            else:
                return

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
        finally:
            # bring the settings window back to the front                
            self.es_top.lift()


    ############################################################################################
    #
    #
    #      open_settings  
    #
    #
    ############################################################################################

    def open_settings(self, file_name):
        try:            
            if len(str(file_name)) > 0:
                self.settings_filename = str(file_name)
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
                    #special case for settings filename
                    self.set_entry_text(self.settings_filename_entry, self.settings_filename)
                    self.settings_filename_entry.xview_scroll(len(self.settings_filename), tk.UNITS)

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            pass

    ############################################################################################
    #
    # show_settings; From File->Edit Settings.... 
    #
    ############################################################################################

    def show_settings(self):
        self.es_top.deiconify()

    ############################################################################################
    #
    # launch_settings  
    #
    # command for "File->Edit Settings..." menu item which creates popup window
    #
    ############################################################################################

    def launch_settings(self):
        try:
            height_factor_ = .7

            es_ = tk.Toplevel(self.window, padx=15, pady=15, takefocus=True)

            es_.title("Settings")
            self.es_top = es_

            self.es_top.protocol("WM_DELETE_WINDOW", self.es_top.withdraw)

            tk.Grid.columnconfigure(self.es_top, 0, weight=1)

            settings_entry_width = 6
            settings_entry_pad = 0
            extended_settings_entry_width = 30
            extended_settings_entry_pad = 0


            #
            #
            #
            #          Settings Left [Side] Frame
            #
            #
            settings_left_frame = tk.Frame(self.es_top, padx=__our_padding__, pady=__our_padding__)
            settings_left_frame.grid(row=0, column=0, sticky=tk.NSEW)

            row = 0

            """
                    ePSF and PSF Photometry Parameters
            """

            separator_ = ttk.Separator(settings_left_frame, orient='horizontal')
            separator_.grid(row=row, columnspan=5, pady=5, sticky=tk.EW)
            row += 1

            _label_ = tk.Label(settings_left_frame, text="ePSF and PSF Photometry Parameters")
            _label_.grid(row=row, columnspan=5, sticky=tk.EW)
            row += 1

            separator_ = ttk.Separator(settings_left_frame, orient='horizontal')
            separator_.grid(row=row, columnspan=5, pady=5, sticky=tk.EW)
            row += 1

            find_peaks_npeaks_label = tk.Label(settings_left_frame, text="Max Number of Peaks:")
            find_peaks_npeaks_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.find_peaks_npeaks_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.find_peaks_npeaks_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            qfit_label = tk.Label(settings_left_frame, text="Max qfit:")
            qfit_label.grid(row=row, column=3, sticky=tk.E)
            self.max_qfit_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.max_qfit_entry.grid(row=row, column=4, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            fit_width_label = tk.Label(settings_left_frame, text="Fitting Width/Height, px:")
            fit_width_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.fit_width_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.fit_width_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            min_separation_bias_label = tk.Label(settings_left_frame, text="Min Separation Bias,(added to FWHM/2+Fitting Width/2):")
            min_separation_bias_label.grid(row=row, column=3, sticky=tk.E)
            self.min_separation_bias_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.min_separation_bias_entry.grid(row=row, column=4, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            max_ensemble_magnitude_label = tk.Label(settings_left_frame, text="Maximum Ensemble Magnitude:")
            max_ensemble_magnitude_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.max_ensemble_magnitude_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.max_ensemble_magnitude_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            separator_ = ttk.Separator(settings_left_frame, orient='horizontal')
            separator_.grid(row=row, column=3, columnspan=2, pady=5, sticky=tk.EW)
            row += 1

            fwhm_label = tk.Label(settings_left_frame, text="FWHM, px:")
            fwhm_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.fwhm_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.fwhm_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            oversampling_label = tk.Label(settings_left_frame, text="oversampling:")
            oversampling_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.oversampling_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.oversampling_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            default_model_label = tk.Label(settings_left_frame, text="Default PSF Model:")
            default_model_label.grid(row=row, column=3, columnspan=2, sticky=tk.W)
            row += 1

            star_detection_threshold_factor_label = tk.Label(settings_left_frame, text="DAOStarFinder Threshold Factor (x std):")
            star_detection_threshold_factor_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.star_detection_threshold_factor_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.star_detection_threshold_factor_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            use_gaussian_prf_model = ttk.Radiobutton(settings_left_frame, text="Circular Gaussian PRF",
                                                         variable=self.use_gaussian_prf_model, value=1)
            use_gaussian_prf_model.grid(row=row, column=3, columnspan=2, sticky=tk.W)
            row += 1

            photometry_iterations_label = tk.Label(settings_left_frame, text="Photometry Iterations:")
            photometry_iterations_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.photometry_iterations_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.photometry_iterations_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            use_moffat_model = ttk.Radiobutton(settings_left_frame, text="Moffat PSF",
                                                  variable=self.use_gaussian_prf_model,  value=0)
            use_moffat_model.grid(row=row, column=3, sticky=tk.W)

            beta_label = tk.Label(settings_left_frame, text="\u03B2:") # β: 
            beta_label.grid(row=row, column=3, sticky=tk.E)
            self.moffat_beta_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.moffat_beta_entry.grid(row=row, column=4, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            sharplo_label = tk.Label(settings_left_frame, text="Lower Bound for Sharpness:")
            sharplo_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.sharplo_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.sharplo_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            separator_ = ttk.Separator(settings_left_frame, orient='horizontal')
            separator_.grid(row=row, column=3, columnspan=2, pady=5, sticky=tk.EW)
            row += 1

            matching_radius_label = tk.Label(settings_left_frame, text="Matching Radius, arcsec:")
            matching_radius_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.matching_radius_entry = tk.Entry(settings_left_frame, width=settings_entry_width)
            self.matching_radius_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            
            checkbox_ = tk.Checkbutton(settings_left_frame, text="Generate Residual Image", variable=self.generate_residual_image)
            checkbox_.grid(row=row, column=3, columnspan=2, sticky=tk.W)
            row += 1

            fitter_label = tk.Label(settings_left_frame, text="PSF Fitter:")
            fitter_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            fitter_dropdown = tk.OptionMenu(settings_left_frame, self.fitter_stringvar,
                                                "TRF LS", "Sequential LS Programming", "Simplex LS")
            fitter_dropdown.grid(row=row, column=2, sticky=tk.EW)
            row += 1
            
            separator_telescope = ttk.Separator(settings_left_frame, orient='horizontal')
            separator_telescope.grid(row=row, columnspan=5, pady=5, sticky=tk.EW)
            row += 1

            #
            #        Telescope Name and Parameters (only used for convenience)
            #

            telescope_label = tk.Label(settings_left_frame, text="Telescope:")
            telescope_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.telescope_entry = tk.Entry(settings_left_frame, width=extended_settings_entry_width)
            self.telescope_entry.grid(row=row, column=2, sticky=tk.W)
            row += 1
            
            #
            # Filter Band Coefficients
            #
            # 
            # Tb_bv,
            # Tb_br,
            # Tb_bi, 
            # Tv_bv,
            # Tv_vr, 
            # Tr_vr,
            # Tr_ri
            # Ti_ri, 
            # Tv_vi,
            # Ti_vi, 
            # Tr_ri,
            # Tbv,
            # Tbr,
            # Tbi,
            # Tvr,
            # Tri,
            # Tvi,
            # 
            #   ... same order as seen in VPhot Telescope Profiles
            #

            #
            # Extinction Coefficients
            #
            # B, V, R, I, C
            #
            #

            _label = tk.Label(settings_left_frame, text="Transformation Coefficients:")
            _label.grid(row=row, column=0, columnspan=3, sticky=tk.W)
            _label = tk.Label(settings_left_frame, text="Extinction Coefficients:")
            _label.grid(row=row, column=2, sticky=tk.E)

            row += 1

            """/* Tb_bv */"""
            tb_bv_label = tk.Label(settings_left_frame, text="Tb_bv:")
            tb_bv_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tb_bv_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tb_bv_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tb_bv_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tb_bv_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            # B extinction
            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=2, sticky=tk.E)
            _extinction_label = tk.Label(f_helper, text="B:")
            _extinction_label.grid(row=0, column=0, sticky=tk.E)
            self.extinction_B_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.extinction_B_entry.grid(row=0, column=1, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tb_br */"""
            tb_br_label = tk.Label(settings_left_frame, text="Tb_br:")
            tb_br_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tb_br_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tb_br_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tb_br_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tb_br_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            # V extinction
            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=2, sticky=tk.E)
            _extinction_label = tk.Label(f_helper, text="V:")
            _extinction_label.grid(row=0, column=0, sticky=tk.E)
            self.extinction_V_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.extinction_V_entry.grid(row=0, column=1, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tb_bi */"""
            tb_bi_label = tk.Label(settings_left_frame, text="Tb_bi:")
            tb_bi_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tb_bi_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tb_bi_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tb_bi_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tb_bi_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            # R extinction
            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=2, sticky=tk.E)
            _extinction_label = tk.Label(f_helper, text="R:")
            _extinction_label.grid(row=0, column=0, sticky=tk.E)
            self.extinction_R_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.extinction_R_entry.grid(row=0, column=1, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tv_bv */"""
            tv_bv_label = tk.Label(settings_left_frame, text="Tv_bv:")
            tv_bv_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tv_bv_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tv_bv_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tv_bv_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tv_bv_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            # I extinction
            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=2, sticky=tk.E)
            _extinction_label = tk.Label(f_helper, text="I:")
            _extinction_label.grid(row=0, column=0, sticky=tk.E)
            self.extinction_I_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.extinction_I_entry.grid(row=0, column=1, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tv_vr */"""
            tv_vr_label = tk.Label(settings_left_frame, text="Tv_vr:")
            tv_vr_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tv_vr_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tv_vr_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tv_vr_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tv_vr_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            # C extinction
            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=2, sticky=tk.E)
            _extinction_label = tk.Label(f_helper, text="C:")
            _extinction_label.grid(row=0, column=0, sticky=tk.E)
            self.extinction_C_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.extinction_C_entry.grid(row=0, column=1, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tr_vr */"""
            tr_vr_label = tk.Label(settings_left_frame, text="Tr_vr:")
            tr_vr_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tr_vr_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tr_vr_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tr_vr_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tr_vr_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tr_ri */"""
            tr_ri_label = tk.Label(settings_left_frame, text="Tr_ri:")
            tr_ri_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tr_ri_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tr_ri_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tr_ri_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tr_ri_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Ti_ri */"""
            ti_ri_label = tk.Label(settings_left_frame, text="Ti_ri:")
            ti_ri_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.ti_ri_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.ti_ri_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.ti_ri_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.ti_ri_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tv_vi */"""
            tv_vi_label = tk.Label(settings_left_frame, text="Tv_vi:")
            tv_vi_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tv_vi_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tv_vi_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tv_vi_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tv_vi_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Ti_vi */"""
            ti_vi_label = tk.Label(settings_left_frame, text="Ti_vi:")
            ti_vi_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.ti_vi_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.ti_vi_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.ti_vi_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.ti_vi_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tr_vi */"""
            tr_vi_label = tk.Label(settings_left_frame, text="Tr_vi:")
            tr_vi_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tr_vi_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tr_vi_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tr_vi_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tr_vi_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tbv */"""
            tbv_label = tk.Label(settings_left_frame, text="Tbv:")
            tbv_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tbv_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tbv_entry.grid(row=0, column=0, ipadx=settings_entry_pad, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tbv_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tbv_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tbr */"""
            tbr_label = tk.Label(settings_left_frame, text="Tbr:")
            tbr_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tbr_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tbr_entry.grid(row=0, column=0, ipadx=settings_entry_pad, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tbr_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tbr_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tbi */"""
            tbi_label = tk.Label(settings_left_frame, text="Tbi:")
            tbi_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tbi_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tbi_entry.grid(row=0, column=0, ipadx=settings_entry_pad, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tbi_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tbi_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tvr */"""
            tvr_label = tk.Label(settings_left_frame, text="Tvr:")
            tvr_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tvr_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tvr_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tvr_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tvr_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            """/* Tri */"""
            tri_label = tk.Label(settings_left_frame, text="Tri:")
            tri_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tri_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tri_entry.grid(row=0, column=0, ipadx=settings_entry_pad, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tri_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tri_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1


            """/* Tvi */"""
            tvi_label = tk.Label(settings_left_frame, text="Tvi:")
            tvi_label.grid(row=row, column=0, sticky=tk.E)

            f_helper = tk.Frame(settings_left_frame)
            f_helper.grid(row=row, column=1, sticky=tk.W)
            self.tvi_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tvi_entry.grid(row=0, column=0, sticky=tk.W)
            _err_label = tk.Label(f_helper, text="+/-")
            _err_label.grid(row=0, column=1, sticky=tk.W)
            self.tvi_err_entry = tk.Entry(f_helper, width=settings_entry_width, background='pink')
            self.tvi_err_entry.grid(row=0, column=2, ipadx=settings_entry_pad, sticky=tk.W)

            linearity_limit_label = tk.Label(settings_left_frame, text="Linearity Limit (ADU):")
            linearity_limit_label.grid(row=row, column=2, sticky=tk.E)
            self.linearity_limit_entry = tk.Entry(settings_left_frame, width=settings_entry_width, background='pink')
            self.linearity_limit_entry.grid(row=row, column=3, sticky=tk.W)
            row += 1

            #
            #
            #
            #          Settings Right [Side] Frame
            #
            #

            settings_right_frame = tk.Frame(self.es_top, padx=__our_padding__, pady=__our_padding__)
            settings_right_frame.grid(row=0, column=1, sticky=tk.NSEW)

            row = 0


            separator_ = ttk.Separator(settings_right_frame, orient='horizontal')
            separator_.grid(row=row, columnspan=3, pady=5, sticky=tk.EW)
            row += 1

            """
                    AAVSO Report Settings
            """
            _label_ = tk.Label(settings_right_frame, text="AAVSO Report Settings")
            _label_.grid(row=row, columnspan=3, sticky=tk.EW)
            row += 1

            separator_ = ttk.Separator(settings_right_frame, orient='horizontal')
            separator_.grid(row=row, columnspan=3, pady=5, sticky=tk.EW)
            row += 1

            aavso_obscode_label = tk.Label(
                settings_right_frame, text="Observer Code:")
            aavso_obscode_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.aavso_obscode_entry = tk.Entry(settings_right_frame, width=extended_settings_entry_width, background='pink')
            self.aavso_obscode_entry.grid(row=row, column=2, ipadx=extended_settings_entry_pad)
            row += 1

            exposure_label = tk.Label(
                settings_right_frame, text="Exposure Time:")
            exposure_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.exposure_entry = tk.Entry(
                settings_right_frame, width=settings_entry_width, background='pink')
            self.exposure_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            airmass_label = tk.Label(settings_right_frame, text="Airmass:")
            airmass_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.airmass_entry = tk.Entry(
                settings_right_frame, width=extended_settings_entry_width)
            self.airmass_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            date_obs_label = tk.Label(
                settings_right_frame, text="Date-Obs (JD):")
            date_obs_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.date_obs_entry = tk.Entry(
                settings_right_frame, width=extended_settings_entry_width)
            self.date_obs_entry.config(state='readonly')
            self.date_obs_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            date_utc_label = tk.Label(
                settings_right_frame, text="Date-Obs (UTC):")
            date_utc_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.date_utc_entry = tk.Entry(
                settings_right_frame, width=extended_settings_entry_width)
            self.date_utc_entry.config(state='readonly')
            self.date_utc_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            object_notes_label = tk.Label(settings_right_frame, text="Notes:")
            object_notes_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.object_notes_entry = tk.Entry(
                settings_right_frame, width=extended_settings_entry_width)
            self.object_notes_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            catalog_label = tk.Label(settings_right_frame, text="Comparison Catalog:")
            catalog_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            catalog_dropdown = tk.OptionMenu(
                settings_right_frame, self.catalog_stringvar, "AAVSO", "Gaia DR2", "APASS DR9", "APASS DR10")
            #, "URAT1", "USNO-B1.0", "VizieR Catalog") <--not supporting now
            catalog_dropdown.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            vizier_catalog_label = tk.Label(settings_right_frame, text="AAVSO ChartID:")
            vizier_catalog_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.vizier_catalog_entry = tk.Entry(settings_right_frame,width=extended_settings_entry_width)
            self.vizier_catalog_entry.grid(row=row, column=2, sticky=tk.E)
            row += 1

            separator_ = ttk.Separator(settings_right_frame, orient='horizontal')
            separator_.grid(row=row, columnspan=3, pady=5, sticky=tk.EW)
            row += 1

            from_fits_header_check = tk.Checkbutton(settings_right_frame, text="From FITS",
                                                     variable=self.auto_behavior,
                                                     command=self.on_auto_behavior_checkbox_checked)
            from_fits_header_check.grid(row=row, column=0, sticky=tk.W)
            
            filter_label = tk.Label(settings_right_frame, text="CCD Filter:")
            filter_label.grid(row=row, column=1, sticky=tk.E)
            self.filter_entry = tk.Entry(settings_right_frame,
                                          width=extended_settings_entry_width, background='pink')
            self.filter_entry.grid(row=row, column=2, ipadx=settings_entry_pad, sticky=tk.W)
            row += 1

            separator_ = ttk.Separator(settings_right_frame, orient='horizontal')
            separator_.grid(row=row, columnspan=3, pady=5, sticky=tk.EW)
            row += 1

            """
                    Miscellaneous Settings
            """
            _label_ = tk.Label(
                settings_right_frame, text="Miscellaneous Settings")
            _label_.grid(row=row, columnspan=3, sticky=tk.EW)
            row += 1

            separator_ = ttk.Separator(settings_right_frame, orient='horizontal')
            separator_.grid(row=row, columnspan=3, pady=5, sticky=tk.EW)
            row += 1

            object_name_label = tk.Label(
                settings_right_frame, text="Object Name:")
            object_name_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.object_name_entry = tk.Entry(
                settings_right_frame, width=extended_settings_entry_width, background='pink')
            self.object_name_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            object_name_alpha_label = tk.Label(
                settings_right_frame, text="α:")
            object_name_alpha_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.object_name_alpha_entry = tk.Entry(
                settings_right_frame, width=extended_settings_entry_width)
            self.object_name_alpha_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            object_name_delta_label = tk.Label(
                settings_right_frame, text="δ:")
            object_name_delta_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.object_name_delta_entry = tk.Entry(
                settings_right_frame, width=extended_settings_entry_width)
            self.object_name_delta_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            object_kref_label = tk.Label(settings_right_frame, text="Use Check Star (AAVSO Label):")
            object_kref_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.object_kref_entry = tk.Entry(
                settings_right_frame, width=extended_settings_entry_width, background='pink')
            self.object_kref_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1

            object_sel_comp_label = tk.Label(settings_right_frame, text="Select Comp Stars (AAVSO Label):")
            object_sel_comp_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.object_sel_comp_entry = tk.Entry(
                settings_right_frame, width=extended_settings_entry_width, background='pink')
            self.object_sel_comp_entry.grid(row=row, column=2, sticky=tk.EW)
            row += 1
    
            object_min_snr_label = tk.Label(settings_right_frame, text="Minimum SNR:")
            object_min_snr_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.object_min_snr_entry = tk.Entry(
                settings_right_frame, width=settings_entry_width)
            self.object_min_snr_entry.grid(row=row, column=2, sticky=tk.W)
            row += 1
    
            display_users_objects_only = ttk.Radiobutton(settings_right_frame, text="Display selected objects only",
                                                         variable=self.display_all_objects, value=0)
            display_users_objects_only.grid(row=row, column=2, columnspan=2, sticky=tk.W)
            row += 1

            display_all_objects_rb = ttk.Radiobutton(settings_right_frame, text="Display all objects",
                                                  variable=self.display_all_objects,  value=1)
            display_all_objects_rb.grid(row=row, column=2, columnspan=2, sticky=tk.W )
            row += 1

            astrometrynet_label = tk.Label(
                settings_right_frame, text="Astrometry.net Server:")
            astrometrynet_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.astrometrynet_entry = tk.Entry(
                settings_right_frame, width=extended_settings_entry_width)
            self.astrometrynet_entry.grid(row=row, column=2, ipadx=extended_settings_entry_pad)
            row += 1

            astrometrynet_key_label = tk.Label(
                settings_right_frame, text="Astrometry.net API Key:")
            astrometrynet_key_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.astrometrynet_key_entry = tk.Entry(
                settings_right_frame, width=extended_settings_entry_width)
            self.astrometrynet_key_entry.grid(row=row, column=2, ipadx=extended_settings_entry_pad)
            self.astrometrynet_key_entry.config(show="*")
            row += 1

            separator_ = ttk.Separator(settings_right_frame, orient='horizontal')
            separator_.grid(row=row, columnspan=3, pady=5, sticky=tk.EW)
            row += 1

            settingsfile_key_label = tk.Label(settings_right_frame, text="Settings Filename:")
            settingsfile_key_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.settings_filename_entry = tk.Entry(settings_right_frame, width=extended_settings_entry_width)
            self.settings_filename_entry.grid(row=row, column=2, ipadx=extended_settings_entry_pad)
            self.set_entry_text(self.settings_filename_entry, self.settings_filename)
            #self.settings_filename_entry.xview_scroll(len(self.settings_filename), tk.UNITS)
            row += 1

            #
            #        User Note (only used for convenience)
            #

            user_note_label = tk.Label(settings_right_frame, text="User Note:")
            user_note_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
            self.user_note_entry = tk.Entry(settings_right_frame, width=extended_settings_entry_width)
            self.user_note_entry.grid(row=row, column=2, ipadx=extended_settings_entry_pad)
            row += 1



            # Separator and Buttons across the bottom
            row=1
            separator_ = ttk.Separator(self.es_top, orient='horizontal')
            separator_.grid(row=row, columnspan=2, pady=5, sticky=tk.EW)
            row += 1

            #
            # Buttons 
            #
            f_helper = tk.Frame(self.es_top)
            f_helper.grid(row=row, columnspan=2, sticky=tk.EW)

            f_helper.grid_columnconfigure(0, weight=1)  # Left-align button in column 0
            f_helper.grid_columnconfigure(1, weight=2)  # Center button in column 1
            f_helper.grid_columnconfigure(2, weight=1)  # Right-align button in column 2

            load_settings_button = tk.Button(f_helper, text="Load...", command=self.load_settings)
            load_settings_button.grid(row=0, column=0, padx=20, sticky=tk.W)

            save_settings_button = tk.Button(f_helper, text="Save As...", command=self.save_settings_as)
            save_settings_button.grid(row=0, column=1, padx=20)

            close_settings_button = tk.Button(f_helper, text="  OK/Hide  ", command=self.es_top.withdraw)
            close_settings_button.grid(row=0, column=2, padx=20, sticky=tk.E)
            row += 1

            # Update layout to calculate dimensions
            self.es_top.update_idletasks()

            # Automatically adjust window size to fit contents
            self.es_top.geometry(f"{self.es_top.winfo_reqwidth()}x{self.es_top.winfo_reqheight()}")

            self.es_top.resizable(False, False)

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)

##############################################################################
#
# forward_selstars_list
#
# Callback for "Forward" button. This invokes display of next page of 
# candidate stars.
# 
# 
##############################################################################

    def forward_selstars_list(self):
        if self.candidate_stars_index < len(self.candidate_stars) - 1:
            self.clear_selstars_plot()
            # look forward
            self.candidate_stars_index += 1
            i = 0
            selstars_plot_index = 0
            resolve_index = False 
            for self.candidate_stars_index in range(self.candidate_stars_index, len(self.candidate_stars)):
                if i == (self.ncols*self.nrows):
                    resolve_index = True
                    break
                i += 1
                norm = simple_norm(self.candidate_stars[self.candidate_stars_index], 'log', percent=99.0)
                self.selstars_plot[selstars_plot_index].imshow(self.candidate_stars[self.candidate_stars_index],
                        norm=norm, origin='lower', cmap='viridis')
                # check if this has already been rejected
                (cand_x, cand_y) = self.candidate_stars[self.candidate_stars_index].center
                if ((self.ePSF_pending_rejection_list['x'] == cand_x) & (self.ePSF_pending_rejection_list['y'] == cand_y)).any():
                    self.selstars_plot[selstars_plot_index].text(x=0,y=5, s="Reject")
                selstars_plot_index += 1
            plt.subplots_adjust(hspace=self.selstars_hspace, wspace=self.selstars_wspace)
            self.selstars_plot_canvas.draw()
            # need to resolve index so that when we forward again we dont skip current index
            # (un)fortunately index never reaches the stop value in the last frame so 
            if resolve_index:
                self.candidate_stars_index -= 1

            # update label with page (n of x) display E.g., Page: 1 of 6
            new_page_number = math.ceil(self.candidate_stars_index/(self.nrows*self.ncols))
            self.update_selstars_page_label(page_num=new_page_number)
        else:
            # nothing past this
            self.console_msg("There are no more selected stars to show.");
        
        self.fig_selstars.canvas.mpl_connect('button_press_event', self.mouse_selstars_canvas_click)
        self.console_msg("candidate_stars index = " + str(self.candidate_stars_index), level=logging.DEBUG)
        return

##############################################################################
#
# back_selstars_list
#
# Callback for "Back" button. This invokes display of previous page of 
# candidate stars.
# 
# 
##############################################################################

    def back_selstars_list(self):
        if self.candidate_stars_index >= (self.ncols*self.nrows): 
            # not displaying the first set
            self.clear_selstars_plot()
            # 
            # We want the index to point to beginning of the last frame (frame being
            # a set of (self.ncols*self.nrows) subplots),
            # so we have to account for it now pointing to somewhere in the middle.
            # This would be true if (len(self.candidate_stars) % self.ncols*self.nrows) != 0
            # so if index+ 1 is not a multiple of self.ncols*self.nrows, then we are at the last frame.
            # If at the last frame, then index must be subtracted accordingly
            self.candidate_stars_index += 1
            candidate_stars_remainder = self.candidate_stars_index % (self.ncols*self.nrows)
            if candidate_stars_remainder != 0:
                self.candidate_stars_index -= (self.ncols*self.nrows + candidate_stars_remainder)
            else:
                self.candidate_stars_index -= 2*self.ncols*self.nrows  #backup

            i = 0
            selstars_plot_index = 0
            for self.candidate_stars_index in range(self.candidate_stars_index, len(self.candidate_stars)):
                if i == (self.ncols*self.nrows):
                    break
                i += 1
                norm = simple_norm(self.candidate_stars[self.candidate_stars_index], 'log', percent=99.0)
                self.selstars_plot[selstars_plot_index].imshow(self.candidate_stars[self.candidate_stars_index],
                        norm=norm, origin='lower', cmap='viridis')
                # check if this has already been rejected
                (cand_x, cand_y) = self.candidate_stars[self.candidate_stars_index].center
                if ((self.ePSF_pending_rejection_list['x'] == cand_x) & (self.ePSF_pending_rejection_list['y'] == cand_y)).any():
                    self.selstars_plot[selstars_plot_index].text(x=0,y=5, s="Reject")
                selstars_plot_index += 1
            plt.subplots_adjust(hspace=self.selstars_hspace, wspace=self.selstars_wspace)
            self.selstars_plot_canvas.draw()
            # need to resolve index so that when we forward again we dont skip current index
            self.candidate_stars_index -= 1
            # update label with page (n of x) display E.g., Page: 1 of 6
            new_page_number = math.ceil(self.candidate_stars_index/(self.nrows*self.ncols))
            self.update_selstars_page_label(page_num=new_page_number)
        else:
            # nothing past this
            self.console_msg("There are no more selected stars to show.")
        
        self.fig_selstars.canvas.mpl_connect('button_press_event', self.mouse_selstars_canvas_click)
        self.console_msg("candidate_stars index = " + str(self.candidate_stars_index), level=logging.DEBUG)
        return

############################################################
#
# submit_rejects_selstars_list
#
############################################################
    def submit_rejects_selstars_list(self):
        #update ePSF_rejection_list with ePSF_pending_rejection_list
        # then update main canvas 
        self.ePSF_rejection_list = pd.concat([self.ePSF_rejection_list, self.ePSF_pending_rejection_list], ignore_index=True)
        #reset pending
        self.ePSF_pending_rejection_list.drop(self.ePSF_pending_rejection_list.index, inplace=True)
        # display the rejected ones (red circle) on main canvas
        self.ePSF_samples_plotted = True
        self.display_image()
        # display an updated selstars area
        self.find_peaks()
        # Now nothing to submit so...
        self.submit_rejects_selstars_button.config(state=tk.DISABLED)
        return


############################################################
#
# clear_rejects_selstars_list
#
# Clears the rejected stars in the selstars list and 
# re-displays them
#
############################################################

    def clear_rejects_selstars_list(self):
        #reset pending
        self.ePSF_pending_rejection_list.drop(self.ePSF_pending_rejection_list.index, inplace=True)
        # display the rejected ones (red circle) on main canvas
        self.ePSF_samples_plotted = True
        self.display_image()
        # display an updated selstars area
        self.find_peaks()
        # Now nothing to submit so...
        self.submit_rejects_selstars_button.config(state=tk.DISABLED)
        return


############################################################
#
#  aavso_get_comparison_stars
#
############################################################

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
            sel_comps_to_use = self.object_sel_comp_entry.get()
            sel_comps_to_use = [comp.strip() for comp in sel_comps_to_use.split(',')]
            # Make unique, just in case
            sel_comps_to_use = list(set(sel_comps_to_use))
                
            for comp in sel_comps_to_use:
                # ignore comp numbers that are not numbers like '#120' or '-120'
                # this way you can 'remove' comp stars from list by just prepending 
                # them with a '-' or some non-numeric char
                if comp[0].isdigit() and int(comp.strip()) > 0: 
                    sel_comps.append(int(comp.strip()))
                
            check_star_to_use = self.object_kref_entry.get().strip()
            if is_number(check_star_to_use):
                # remove the check star from this list just in case
                # this removes all of them 
                sel_comps = [star for star in sel_comps if star != check_star_to_use]
                #add it to the list of comps
                check_star = int(check_star_to_use)
                sel_comps.append(check_star)
            else:
                check_star = -1 #this value never matches a label

            chart_id = r.json()['chartid']
            self.console_msg(
                'Downloaded AAVSO Comparison Star Chart ID ' + str(chart_id))
            
            result = pd.DataFrame(columns=["AUID", "RA", "Dec", "Mag", "Mag Error", "Label", "Chart ID", "Check Star"])

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
                        mag_error = band['error']

                if mag == None:
                    self.console_msg("label: " + str(label) + " has no mag for " + filter_band + "..skipping")
                    continue #skip this one
                        
                result.loc[len(result)] = {"AUID": auid,
                    "RA": ra,
                    "Dec": dec,
                    "Mag": mag,
                    "Mag Error": mag_error,
                    "Label": label,
                    "Chart ID": chart_id,
                    "Check Star": is_check_star}

            return result
        
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)

############################################################
#
#  APASS_DR10_get_comparison_stars
#
############################################################

    def APASS_DR10_get_comparison_stars(self, frame_center, filter='V', field_of_view=2, maglimit=20):
        try:
           
            ra = frame_center.ra
            dec = frame_center.dec

            if not is_number(maglimit):
                maglimit = 20
            else:
                maglimit = float(maglimit)

            r = requests.get('https://www.aavso.org/cgi-bin/apass_dr10_download.pl?ra='+str(ra)+
                             '&dec='+str(dec)+
                             '&radius='+str(field_of_view)+
                             '&maglimit='+str(maglimit)+
                             '&outtype=1')

            if r.status_code != 200:
                self.console_msg("Fetch of APASS DR10 Comparison Stars Failed! status: "+ str(r.status_code))
                return -1
            else:
                self.console_msg("Downloaded APASS DR10 Comparison Stars")

               
            # Back to the result
            result_df = pd.read_csv(io.StringIO(r.text))

            #
            #  KEEP UNTIL cgi-bin/apass_dr10_download.pl IS FIXED!!!!!!
            #
            #remove any with Johnson V > maglimit (until cgi-bin/apass_dr10_download.pl is fixed!)
            self.console_msg("Removing entries where 'Johnson_V (V)' > maglimit REMOVE after cgi-bin/apass_dr10_download.pl is fixed")
            result_df.drop(result_df[result_df['Johnson_V (V)'] > maglimit].index, inplace=True)

            """
            If we are using R or I We have to convert the Sloan Filters to Johnson equivalents
            
            From code George sent me:
            sr: Sloan_r (SR)
            si: Sloan_i (SI)
            sr_e: SRerr
            si_e: SIerr


            Formula

            If sr > 0 And si > 0 Then
                r = v - 1.09 * (sr - si) - 0.22
                r_e = Math.Sqrt(v_e ^ 2 + 1.09 ^ 2 * (sr_e ^ 2 + si_e ^ 2))
                i = r - (sr - si) - 0.21
                i_e = Math.Sqrt(r_e ^ 2 + (sr_e ^ 2 + si_e ^ 2))

            Existing columns:
                       
            radeg,raerr("),decdeg,decerr("),Johnson_V (V),Verr,Vnobs,Johnson_B (B),
            Berr,Bnobs,Sloan_u (SU),SUerr,SUnobs,Sloan_g (SG),SGerr,SGnobs,Sloan_r (SR),
            SRerr,SRnobs,Sloan_i (SI),SIerr,SInobs,Sloan_z (SZ),SZerr,SZnobs,PanSTARRS_Y (Y),Yerr,Ynobs

            """

            """
            filter_scheme = {
                'V' : 'Johnson_V (V)',
                'B' : 'Johnson_B (B)',
                'R' : 'Sloan_r (SR)',
                'I' : 'Sloan_i (SI)'
            }

            # update renaming_scheme with proper filter to use
            filter_in_apass = filter_scheme[filter]

            # add this filter_in_apass to list of columns we want to keep
            # this becomes the Mag we are interrested in
            renaming_list[filter_in_apass] = 'Mag'

            if 'Johnson_V (V)' not in renaming_list: 
                # This column always has to be in the table because is derives values under the Label
                renaming_list['Johnson_V (V)'] = 'V'

            
            """

            # rename columns 
            renaming_list = {
                'radeg' : 'RA',
                'decdeg' : 'Dec',
                'Johnson_V (V)' : 'V',
                'Johnson_B (B)' : 'B',
                'Sloan_r (SR)' : 'sr',
                'Sloan_i (SI)' : 'si',
                'Verr' : 'V_e',
                'SRerr' : 'sr_e',
                'SIerr' : 'si_e'
            }

            # build a list of columns to drop, the ones not in the renaming_list
            drop_list = []
            for col in result_df.columns:
                if col not in renaming_list:
                    drop_list.append(col)
            # drop them
            result_df.drop(columns=drop_list, inplace=True)

            # rename using the renaming_list
            result_df = result_df.rename(columns=renaming_list)

            if filter == 'I' or filter == 'R':
                #we need to convert 
                """
                If sr > 0 And si > 0 Then
                    R = V - 1.09 * (sr - si) - 0.22
                    R_e = Math.Sqrt(V_e ^ 2 + 1.09 ^ 2 * (sr_e ^ 2 + si_e ^ 2))
                    I = R - (sr - si) - 0.21
                    I_e = Math.Sqrt(r_e ^ 2 + (sr_e ^ 2 + si_e ^ 2))

                """ 

                #Define some helper functions
                def calc_R(v, sr, si):
                    return v - 1.09 * (sr - si) -.22

                def calc_R_e(v_e, sr_e, si_e):
                    return np.sqrt(v_e**2 + 1.09**2 * (sr_e**2 + si_e**2))

                def calc_I(r, sr, si):
                    return r - (sr - si) -.21

                def calc_I_e(r_e, sr_e, si_e):
                    return np.sqrt(r_e**2 + (sr_e**2 + si_e**2))


                # Remove NaN's
                result_df.dropna(subset=['V', 'sr', 'si'], inplace=True)

                #remove any sr<0 or si<0
                result_df.drop(result_df[result_df['sr'] < 0].index, inplace=True)
                result_df.drop(result_df[result_df['si'] < 0].index, inplace=True)


                # Calculate Johnson equivalents
                result_df['R'] = result_df.apply(lambda row: calc_R(row['V'], row['sr'], row['si']), axis=1)
                result_df["R_e"] = result_df.apply(lambda row: calc_R_e(row['V_e'], row['sr_e'], row['si_e']), axis=1)
                result_df['I'] = result_df.apply(lambda row: calc_I(row['R'], row['sr'], row['si']), axis=1)
                result_df["I_e"] = result_df.apply(lambda row: calc_I_e(row['R_e'], row['sr_e'], row['si_e']), axis=1)

            else:
                # Here if V or B so Drop NaN's for these irst get rid of any NaN in 'V' column
                result_df = result_df.dropna(subset=['V', 'B'])


            # Add the Label column E.g., Vmag = 13.2, then Label is 132A

            result_df["Label"] = result_df['V'].apply(lambda x: str(int(round(x * 10.0, 0)))+'A')

            # Rename the filter of interest (filter) to just Mag
            result_df = result_df.rename(columns={filter : 'Mag'})

            # Now append any duplicate labels with _N, where N is the nth duplicate

            # Identify duplicates and create a counter
            result_df['dup_count'] = result_df.groupby('Label').cumcount()

            # Append counter only to duplicates
            result_df['Label'] = result_df.apply(lambda x: f"{x['Label']}_{x['dup_count']}" if x['dup_count'] > 0 else x['Label'], axis=1)

            # Drop the helper column
            result_df.drop(columns=['dup_count'], inplace=True)


            #Add the "Check Star" column where it is set to True if it is the user's check star
            check_star_to_use = self.object_kref_entry.get().strip()

            result_df['Check Star'] = result_df['Label'] == check_star_to_use

            return result_df
        
        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)

    #####################################################################################
    #
    #   generate_aavso_report_1image
    # 
    #   Generates a AAVSO report in extended format for a single filter.
    #   This is done for a single comp or an ensemble.
    #   In the ensemble case, CMAGINS and CREFMAG are calculated for use by 
    #   VPhot's Transform Applier (TA). 
    #
    #   For single comp:
    #       CDEC: dec of comp star in degrees
    #       CRA:  ra of comp star in degrees
    #       CMAGINS: the comp star instrumental magnitude
    #       CREFMAG: the reference magnitudes of comp  (from chart)
    #       CREFERR: the error of the comp reference magnitude (from chart)
    #
    #       KDEC: dec of check star in degrees
    #       KRA:  ra of check star in degrees
    #       KMAGINS: the check star instrumental magnitude
    #       KREFMAG: the reference magnitudes of check star (from chart)
    #       KREFERR: the error of the check reference magnitude (from chart)
    #       KMAGSTD: the calculated check star standard mag 
    #
    #       CREFMAG - CMAGINS = KMAGSTD - KMAGINS   
    #
    #
    #   For ensemble:
    #       CMAGINS: the average of the ensemble comp star instrumental magnitudes
    #       CREFMAG: the average of the ensemble comp star reference magnitudes
    #       CRA: average comp's RAs  (Unknown definition)
    #       CDEC: average comp's DECs
    #
    #       KDEC: dec of check star in degrees
    #       KRA:  ra of check star in degrees
    #       KMAGINS: the check star instrumental magnitude
    #       KREFMAG: the reference magnitudes of check star (from chart)
    #       KREFERR: the error of the check reference magnitude (from chart)
    #       KMAGSTD: the calculated check star standard mag 
    #
    #       CREFMAG - CMAGINS = KMAGSTD - KMAGINS   
    #
    #       ENSTYPE = 1 (Unknown definition)
    #
    #
    #   For single comp and ensemble:
    #       VMAGINS: the target star instrumental magnitude
    #       
    #
    #####################################################################################

    def generate_aavso_report_1image(self):
        global image_width, image_height

        """
        Typical Single Image Reports shown below.

        Requirements to run:
         - Fits file loaded, plate solved and comp and vsx star shown and entered in settings
         If there are more than 1 comp stars in the list 'Select Comp Stars', then CNAME=ENSEMBLE
         and the average Target Estimate is used. 
        """

        self.console_msg("Beginning Generate AAVSO Single Image Photometry Report...")

        var_star_name = self.object_name_entry.get().strip()
        if len(var_star_name) == 0:
            self.console_msg(
                "'Object Name' must be specified; eg. 'W Her'")
            self.console_msg("Ready")
            return

        check_star_name = self.object_kref_entry.get().strip()
        if len(check_star_name) == 0:
            self.console_msg(
                "'Use Check Star (AAVSO Label)' must be specified; eg. '144'")
            return

        comp_star_list = self.object_sel_comp_entry.get()
        comp_star_list = [comp.strip() for comp in comp_star_list.split(',')]
        # Make unique, just in case
        comp_star_list = list(set(comp_star_list))
        # remove the check star from this list just in case
        comp_star_list = [star for star in comp_star_list if star != check_star_name]
        # remove any commented out comps
        comp_star_list = [star for star in comp_star_list if star[0].isdigit()]
        #check if anything left
        if len(comp_star_list) == 0:
            self.console_msg(
                "'Select Comp Stars (AAVSO Label)' must be specified; eg. '144'")
            return

        #if there is a list of comp stars, then perform ENSEMBLE photometry
        if len(comp_star_list) == 1:
            # single comp star
            comp_star_name = comp_star_list[0]
            ensemble = False
        else:
            ensemble = True

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
            self.console_msg("Cannot proceed; run 'Photometry->Get Comparison Stars' first .")
            self.console_msg("Hint: Check if proper photometry file was selected.")
            return

        try:
            #Check if the Var to report on has been measured
            if not var_star_name in self.results_tab_df_color["vsx_id"].values:
                self.console_msg("Var not found in table; "+ str(var_star_name) + " not found!", level=logging.ERROR)
                return

            if not ensemble:
                #Check if comp star has been measured
                if not (__label_prefix__ + comp_star_name) in self.results_tab_df_color["label"].values:
                    self.console_msg("Comp star not found in table; "+ comp_star_name + " not found!", level=logging.ERROR)
                    return
            else:
                #Check if any comp star has been measured
                comp_star_list_with_prefix = [__label_prefix__ + x for x in comp_star_list]
                if not self.results_tab_df_color["label"].isin(comp_star_list_with_prefix).any():
                    self.console_msg("No comp stars found in table!", level=logging.ERROR)
                    return
            

            with open(report_filename, mode='w') as f:
    
                decimal_places = 4
                decimal_places_for_ra_dec = 5
                
    
                f.write("#TYPE=Extended\n")
                f.write("#OBSCODE="+self.aavso_obscode_entry.get()+"\n")
                f.write("#SOFTWARE=Self-developed; " + self.program_full_name + "\n") 
                f.write("#DELIM=,\n")
                f.write("#DATE=JD\n")
                f.write("#OBSTYPE=CCD\n")
                f.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")

                #var_star_name, check_star_name were determined above
                var_star = self.results_tab_df_color[self.results_tab_df_color["vsx_id"] == var_star_name].iloc[0]
                var_IM = var_star["inst_mag"]

                check_star = self.results_tab_df_color[self.results_tab_df_color["label"] == (__label_prefix__ + check_star_name)].iloc[0]
                check_IM = check_star["inst_mag"]
                check_star_mag_ref = check_star["match_mag"]
                check_star_mag_ref_error = check_star["match_mag_error"]

                check_ra = float(check_star["match_ra"])
                check_dec = float(check_star["match_dec"])

                if not ensemble:
                    """
                    Example Single Image Report (Single comp example)
                    
                        #TYPE=EXTENDED
                        #OBSCODE=Zzzz
                        #SOFTWARE=Self-developed; MAOPhot 1.2.0 using photutils.psf
                        #DELIM=,
                        #DATE=JD
                        #OBSTYPE=CCD
                        #NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
                        S Lyn,2461127.587720,14.869,0.015,B,NO,STD,101,-10.992,124,-8.475,1.1730,na,X42057BFL,Mittelman ATMoB Observatory|
                        CDEC=57.669976|CMAGINS=-10.992|CRA=101.262311|CREFERR=0.438|CREFMAG=10.807|KDEC=57.921511|KMAGINS=-8.475|
                        KMAGSTD=13.324|KRA=101.285689|KREFERR=0.015|KREFMAG=13.089|VMAGINS=-6.930
                    
                    """

                    #comp_star_name was determined above
                    comp_star = self.results_tab_df_color[self.results_tab_df_color["label"] == (__label_prefix__ + comp_star_name)].iloc[0]
                    comp_star_mag = float(comp_star["match_mag"])
                    comp_star_mag_error = float(comp_star["match_mag_error"])
                    comp_star_ra = float(comp_star["match_ra"])
                    comp_star_dec = float(comp_star["match_dec"])
                    comp_IM = comp_star["inst_mag"]

                    #differential photometry
                    var_star_mag = var_IM - comp_IM + comp_star_mag
                    check_star_mag = check_IM - comp_IM + comp_star_mag

                    #error estimates
                    # MERR
                    snr_var_star = var_star['flux_fit']/self.image_bkg_value
                    snr_check_star = check_star['flux_fit']/self.image_bkg_value

                    var_star_err = 1/snr_var_star
                    check_star_err = 1/snr_check_star

                    starid = var_star_name
                    date = format(float(self.date_obs_entry.get()), '.6f') 
                    mag = str(round(var_star_mag, decimal_places))
                    merr = str(round(var_star_err, decimal_places))
                    filt = self.filter_entry.get().strip()
                    trans = "NO"
                    mtype = "STD"
                    cname = comp_star_name
                    cmag = str(round(comp_IM, decimal_places))
                    kname = check_star_name
                    kmag = str(round(check_IM, decimal_places)) #not same as KMAG in notes
                    amass = self.airmass_entry.get().strip() if len(self.airmass_entry.get()) > 0 else "na"
                    group = "na"
                    chart = self.vizier_catalog_entry.get().strip()
                    notes = self.object_notes_entry.get().strip()
                    notes += "|CMAGINS=" + cmag + \
                            "|CRA=" + str(round(comp_star_ra, decimal_places_for_ra_dec)) +\
                            "|CDEC=" + str(round(comp_star_dec, decimal_places_for_ra_dec)) +\
                            "|CREFERR=" + str(round(comp_star_mag_error, decimal_places)) +\
                            "|CREFMAG=" + str(round(comp_star_mag, decimal_places)) +\
                            "|KMAGINS=" + str(round(check_IM, decimal_places)) +\
                            "|KMAGSTD=" + str(round(check_star_mag, decimal_places)) +\
                            "|KREFERR=" + str(round(check_star_mag_ref_error, decimal_places)) +\
                            "|KREFMAG=" + str(round(float(check_star_mag_ref), decimal_places)) +\
                            "|KRA=" + str(round(check_ra, decimal_places_for_ra_dec)) +\
                            "|KDEC=" + str(round(check_dec, decimal_places_for_ra_dec)) +\
                            "|VMAGINS=" + str(round(var_IM, decimal_places))

                    # Print the results to console (this is the Not Ensemble case)
                    self.console_msg("Check Star Estimates using check star: " + check_star_name +
                                     " ("+filt+": " + str(round(float(check_star_mag_ref), decimal_places)) + ")\n" +
                                     (filt+"* Mag: " + format(check_star_mag, ' >6.3f')).rjust(79) +
                                     '\n' +
                                     (filt+"* Err: " + format(check_star_err, ' >6.3f')).rjust(79))

                    self.console_msg(var_star_name+" Variable Star Estimates:\n" +
                                    (filt+"* Mag: " + format(var_star_mag, ' >6.3f')).rjust(79)  +
                                    '\n' +
                                    (filt+"* Err: " + format(var_star_err, ' >6.3f')).rjust(79))
                else:
                    # ENSEMBLE case

                    """
                    Example Single Image Report (Ensemble example)
                                
                        #TYPE=Extended
                        #OBSCODE=Zzzz
                        #SOFTWARE=Self-developed; MAOPhot 1.2.0 using photutils.psf
                        #DELIM=,
                        #DATE=JD
                        #OBSTYPE=CCD
                        #NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
                        S Lyn,2461127.58764,14.722,0.114,B,NO,STD,ENSEMBLE,na,134,14.063,1.1733,na,X42057BFO,Mittelman ATMoB Observatory|
                        KDEC=57.86925|KMAGINS=-7.493|KMAGSTD=14.063|KRA=101.05321|KREFERR=0.017|KREFMAG=14.047|VMAGINS=-6.834 

                    """
                    comp_data = self.results_tab_df_color[self.results_tab_df_color["label"].isin(comp_star_list_with_prefix)]

                    check_star_mag = [check_IM - comp_data["inst_mag"] + comp_data["match_mag"] for _, row in comp_data.iterrows()]
                    ave_check_star_mag = np.mean(check_star_mag)
                    stdev_check_star_mag = np.std(check_star_mag)

                    var_star_mag = [var_IM - comp_data["inst_mag"] + comp_data["match_mag"] for _, row in comp_data.iterrows()]
                    ave_var_star_mag = np.mean(var_star_mag)
                    ave_comp_match_mag = np.mean(comp_data["match_mag"])
                    ave_comp_inst_mag = np.mean(comp_data["inst_mag"])
                    ave_comp_ra = np.mean(comp_data["match_ra"])
                    ave_comp_dec = np.mean(comp_data["match_dec"])

                    starid = var_star_name
                    date = format(float(self.date_obs_entry.get()), '.6f') 
                    mag = str(round(ave_var_star_mag, decimal_places))
                    merr = str(round(stdev_check_star_mag, decimal_places))
                    filt = self.filter_entry.get().strip()
                    trans = "NO"
                    mtype = "STD"
                    cname = "ENSEMBLE"
                    cmag = "na"
                    kname = check_star_name
                    kmag = str(round(ave_check_star_mag, decimal_places)) # same as KMAGSTD
                    amass = self.airmass_entry.get().strip() if len(self.airmass_entry.get()) > 0 else "na"
                    group = "na"
                    chart = self.vizier_catalog_entry.get().strip()
                    notes = self.object_notes_entry.get().strip()
                    # Unknown def (ave DEC?) "|CDEC=" + str(round(ave_comp_dec, decimal_places))
                    # Unknown def (ave RA?)  "|CRA=" + str(round(ave_comp_ra, decimal_places)) 
                    # Unknown def (ENSTYPE)  "|ENSTYPE=1"
                    # Not included in VPhot but required for TA:
                    #  "|CMAGINS=" + str(round(ave_comp_inst_mag, decimal_places))
                    #  "|CREFMAG=" + str(round(ave_comp_match_mag, decimal_places))
                    notes += \
                            "|CMAGINS=" + str(round(ave_comp_inst_mag, decimal_places)) +\
                            "|CREFMAG=" + str(round(ave_comp_match_mag, decimal_places)) +\
                            "|KDEC=" + str(round(check_dec, decimal_places_for_ra_dec)) +\
                            "|KMAGINS=" + str(round(check_IM, decimal_places)) +\
                            "|KMAGSTD=" + str(round(ave_check_star_mag, decimal_places)) +\
                            "|KRA=" + str(round(check_ra, decimal_places_for_ra_dec)) +\
                            "|KREFERR=" + str(round(float(check_star_mag_ref_error), decimal_places)) +\
                            "|KREFMAG=" + str(round(float(check_star_mag_ref), decimal_places)) +\
                            "|VMAGINS=" + str(round(var_IM, decimal_places))
                
                    # make a table for "pretty" (for printing) Comparison Stars data
                    # Rename some of the columns of interest:
                    pretty_comp_data = pd.DataFrame()
                    pretty_comp_data["Star"] = comp_data["label"].str.split().str[1] # We just want the 3-digit number
                    pretty_comp_data["IM"] = comp_data["inst_mag"]
                    pretty_comp_data["X"] = comp_data["x_fit"]
                    pretty_comp_data["Y"] = comp_data["y_fit"]
                    pretty_comp_data[filt+"-mag"] = comp_data["match_mag"]
                    pretty_comp_data["TargetEstimate"] = round(var_IM - comp_data["inst_mag"] + comp_data["match_mag"], decimal_places)

                    #determine the OUTLIERS, IQR, Interquartile Range 
                    q3, q1 = np.percentile(pretty_comp_data["TargetEstimate"], [75 ,25])
                    iqr = q3 - q1
                    upper_limit = q3 + (iqr * 1.5)
                    lower_limit = q1 - (iqr * 1.5)
                    pretty_comp_data["outlier"] = np.where((pretty_comp_data["TargetEstimate"] < lower_limit) | (pretty_comp_data["TargetEstimate"] > upper_limit), "<--OUTLIER", "")

                    # Print the Comparison Stars Data (this is the Ensemble case)
                    self.console_msg("Variable Star Estimates of Var: " + var_star_name + "\n" +
                                      pretty_comp_data.sort_values(by="Star").to_string())

                    # Print the results to console (this is the Ensemble case)
                    self.console_msg("Check Star Estimates using check star: " + check_star_name +
                                     " ("+filt+": " + str(round(float(check_star_mag_ref), decimal_places)) + ")\n" +
                                     (filt+"* Mag: " + format(ave_check_star_mag, ' >6.3f')).rjust(79) +
                                     '\n' +
                                     (filt+"* Err: " + format(stdev_check_star_mag, ' >6.3f')).rjust(79))

                    self.console_msg(var_star_name+" Variable Star Estimates:\n" +
                                    (filt+"* Mag: " + format(ave_var_star_mag, ' >6.3f')).rjust(79)  +
                                    '\n' +
                                    (filt+"* Err: " + format(stdev_check_star_mag, ' >6.3f')).rjust(79))

                # Add " " after notes, because TA clobbers last char
                f.write(starid+","+date+","+mag+","+merr+","+filt+","+trans+","+mtype+"," +
                        cname+","+cmag+","+kname+","+kmag+","+amass+","+group+","+chart+","+notes+" \n")

            self.console_msg("AAVSO Photometry report saved to " + str(report_filename))
            self.console_msg("Ready")

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno)  + " " + str(e), level=logging.ERROR)

        return
        
    #####################################################################################
    #
    #   AAVSO-report menu callbacks.  All seven dispatch to the single
    #   generate_aavso_report_color(input_color) method below.
    #
    #####################################################################################
    def BV_generate_aavso_report(self):
        self.generate_aavso_report_color('BV')

    def VR_generate_aavso_report(self):
        self.generate_aavso_report_color('VR')

    def VI_generate_aavso_report(self):
        self.generate_aavso_report_color('VI')

    def BVR_generate_aavso_report(self):
        self.generate_aavso_report_color('BVR')

    def BVI_generate_aavso_report(self):
        self.generate_aavso_report_color('BVI')

    def VRI_generate_aavso_report(self):
        self.generate_aavso_report_color('VRI')

    def BVRI_generate_aavso_report(self):
        self.generate_aavso_report_color('BVRI')


    #####################################################################################
    #
    #   generate_aavso_report_color
    #
    #   Unified AAVSO extended-format report generator.  Replaces the former
    #   generate_aavso_report_2color and generate_aavso_report_3color methods,
    #   and adds handling for input_color == 'BVRI'.
    #
    #   Reads the Master-Report CSV produced by color_photometry(), extracts
    #   per-filter means/stds and aux rows, prints a console summary, and writes
    #   one AAVSO extended-format data row per filter to a .txt file in the
    #   aavso_reports/ subdirectory.
    #
    #   The aux rows now carry the generic schema
    #       color, JD, KMAGS, KMAGINS, KREFMAG, VMAGINS, Date-Obs, KNAME, AMASS,
    #       coeff_val, coeff_err, coeff_text, coeff_err_text
    #   which is what color_photometry writes for ALL input_color values.
    #
    #   Supported input_color: 'BV', 'VR', 'VI', 'BVR', 'BVI', 'VRI', 'BVRI'.
    #
    #   Typical AAVSO extended-format row (one per filter):
    #
    #     #TYPE=EXTENDED
    #     #OBSCODE=FPIA
    #     #SOFTWARE=Self-developed; MAOPhot ...
    #     #DELIM=,
    #     #DATE=JD
    #     #OBSTYPE=CCD
    #     #NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES
    #     Z Tau,2461117.520799,9.883,0.002,V,YES,STD,ENSEMBLE,na,118,11.678,1.2010,0101,,
    #       TEST Mittelman ATMoB Observatory|KMAGINS=-9.842|KMAGSTD=11.745|KREFMAG=11.810
    #       |Tv_vi=0.0250|Tv_viErr=0.0040|VMAGINS=-11.726
    #
    #####################################################################################

    def generate_aavso_report_color(self, input_color):
        global image_width, image_height

        # ----------------------------------------------------------------
        # Map input_color -> ordered filter tuple
        # ----------------------------------------------------------------
        filters_map = {
            'BV':   ('B', 'V'),
            'VR':   ('V', 'R'),
            'VI':   ('V', 'I'),
            'BVR':  ('B', 'V', 'R'),
            'BVI':  ('B', 'V', 'I'),
            'VRI':  ('V', 'R', 'I'),
            'BVRI': ('B', 'V', 'R', 'I'),
        }
        if input_color not in filters_map:
            raise Exception("generate_aavso_report_color: unknown input_color: " + str(input_color))
        filters = filters_map[input_color]
        n_filters = len(filters)
        method_label = {2: "Two", 3: "Three", 4: "Four"}[n_filters]

        self.console_msg("Beginning Generate AAVSO " + method_label +
                         " Color Ensemble Report (" + input_color + ")...")

        var_star_name = self.object_name_entry.get().strip()
        if len(var_star_name) == 0:
            self.console_msg("'Object Name' must be specified; eg. 'V1117 Her'")
            self.console_msg("Ready")
            return

        # ----------------------------------------------------------------
        # Ask user for <ObjectName>-<input_color>-Master-Report.csv
        # ----------------------------------------------------------------
        options = {}
        options['defaultextension'] = '.csv'
        options['filetypes'] = [('CSV', '.csv')]
        options['title'] = ('Choose the ' + var_star_name + '-' +
                            input_color + '-Master-Report.csv')

        file_name = fd.askopenfilename(**options)
        if len(str(file_name)) > 0 and os.path.isfile(str(file_name)):
            self.console_msg("Loading Master-Report from " + str(file_name))
            master_report = pd.read_csv(str(file_name))
        else:
            return

        # ----------------------------------------------------------------
        # Create aavso_reports subdir
        # ----------------------------------------------------------------
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
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +
                             " " + str(e), level=logging.ERROR)

        _basename = os.path.basename(file_name)
        report_filename = os.path.join(
            report_dir, "AAVSO " + os.path.splitext(_basename)[0] +
                        " " + str(self.object_name_entry.get()) + ".txt")

        decimal_places = 3   # AAVSO standard

        # ----------------------------------------------------------------
        # Extract per-filter means and stds from check/var rows.
        # NOTE: MERR for each filter uses the CHECK-star scatter (stds_check),
        # matching the original two/three-color reports.
        # ----------------------------------------------------------------
        result_check_star = master_report[master_report["type"] == "check"]
        result_var_star   = master_report[master_report["type"] == "var"]

        means_check = {f: result_check_star[f + "_star"].mean() for f in filters}
        means_var   = {f: result_var_star  [f + "_star"].mean() for f in filters}
        stds_check  = {f: result_check_star[f + "_star"].std()  for f in filters}
        stds_var    = {f: result_var_star  [f + "_star"].std()  for f in filters}

        # ----------------------------------------------------------------
        # Per-filter aux rows (keyed by "color" column, which holds 'B','V',...)
        # ----------------------------------------------------------------
        aux = {f: master_report[master_report["color"] == f] for f in filters}

        # KNAME may be a float (AAVSO chart comps, e.g. 150.0 -> "150")
        # or a string (APASS comps, e.g. "150A_1"). Apply once here.
        first_kname_raw = aux[filters[0]]["KNAME"].iloc[0]
        if isinstance(first_kname_raw, str):
            check_kname_str = first_kname_raw
        elif isinstance(first_kname_raw, (np.float64, int)):
            check_kname_str = str(int(first_kname_raw))
        else:
            check_kname_str = '???'
            self.console_msg("Unrecognised KNAME type in Master-Report")

        # ----------------------------------------------------------------
        # Console summary (check star, then variable star)
        # ----------------------------------------------------------------
        # Check-star line: name + one (FILT: KREFMAG) per filter
        check_hdr = "Check Star Estimates using check star: " + check_kname_str
        for f in filters:
            check_hdr += " (" + f + ": " + \
                str(round(float(aux[f]['KREFMAG']), decimal_places)) + ")"
        ave_line = ""
        std_line = ""
        for f in filters:
            ave_line += f + "* Ave: " + format(means_check[f], ' >6.3f') + "  "
            std_line += f + "* Std: " + format(stds_check[f],  ' >6.3f') + "  "
        self.console_msg(check_hdr + '\n' +
                         ave_line.rstrip().rjust(72) + '\n' +
                         std_line.rstrip().rjust(72))

        var_hdr = var_star_name + " Variable Star Estimates"
        ave_line = ""
        std_line = ""
        for f in filters:
            ave_line += f + "* Ave: " + format(means_var[f], ' >6.3f') + "  "
            std_line += f + "* Std: " + format(stds_var[f],  ' >6.3f') + "  "
        self.console_msg(var_hdr + '\n' +
                         ave_line.rstrip().rjust(72) + '\n' +
                         std_line.rstrip().rjust(72))

        # ----------------------------------------------------------------
        # Write AAVSO extended-format .txt
        # ----------------------------------------------------------------
        try:
            with open(report_filename, mode='w') as f_out:
                f_out.write("#TYPE=Extended\n")
                f_out.write("#OBSCODE=" + self.aavso_obscode_entry.get() + "\n")
                f_out.write("#SOFTWARE=Self-developed; " + self.program_full_name + "\n")
                f_out.write("#DELIM=,\n")
                f_out.write("#DATE=JD\n")
                f_out.write("#OBSTYPE=CCD\n")
                f_out.write("#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,"
                            "KNAME,KMAG,AMASS,GROUP,CHART,NOTES\n")

                starid = str(self.object_name_entry.get())
                for filt in filters:
                    aux_row = aux[filt]

                    date  = format(float(aux_row['Date-Obs'].iloc[0]), '.6f')
                    mag   = str(round(means_var[filt],  decimal_places))
                    merr  = str(round(stds_check[filt], decimal_places))
                    trans = "YES"
                    mtype = "STD"
                    cname = "ENSEMBLE"
                    cmag  = "na"
                    kname = check_kname_str
                    kmag  = str(round(means_check[filt], decimal_places))

                    amass_val = aux_row["AMASS"].iloc[0]
                    amass = (format(float(amass_val), '.3f')
                             if type(amass_val) == np.float64 else "na")

                    group = "na"
                    chart = self.vizier_catalog_entry.get().strip()
                    notes = self.object_notes_entry.get().strip()
                    notes += (
                        "|KMAGINS=" + str(round(float(aux_row["KMAGINS"].iloc[0]), decimal_places)) +
                        "|KMAGSTD=" + str(round(means_check[filt], decimal_places)) +
                        "|KREFMAG=" + str(round(float(aux_row["KREFMAG"].iloc[0]), decimal_places)) +
                        "|" + str(aux_row["coeff_text"].iloc[0]) + "=" +
                              str(round(float(aux_row["coeff_val"].iloc[0]), decimal_places)) +
                        "|" + str(aux_row["coeff_err_text"].iloc[0]) + "=" +
                              str(round(float(str(aux_row["coeff_err"].iloc[0])), decimal_places)) +
                        "|VMAGINS=" + str(round(float(aux_row["VMAGINS"].iloc[0]), decimal_places))
                    )

                    # Trailing space after notes: TA clobbers the last char.
                    f_out.write(",".join([
                        starid, date, mag, merr, filt, trans, mtype,
                        cname, cmag, kname, kmag, amass, group, chart, notes
                    ]) + " \n")

            self.console_msg("AAVSO Photometry report saved to " + str(report_filename))
            self.console_msg("Ready")

        except Exception as e:
            self.error_raised = True
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +
                             " " + str(e), level=logging.ERROR)
            

    ##########################################################################
    # 
    # 
    #  exit_app
    # 
    # 
    # Ask for confirmation before saving conig items exiting; 
    #
    ##########################################################################        
    
    def exit_app(self):
        if tk.messagebox.askokcancel("Quit", "Do you really want to exit?"):
            # save config items
            try:
                config = {}
                '''
                Use the valid_config list which contains the official list
                of config parameters
                '''
                mao_config = list(self.valid_config_list)
                for param in mao_config:
                    config.update({param : getattr(self, param)})

                with open(str(self.config_file), 'w') as f:
                    w = csv.DictWriter(f, config.keys())
                    w.writeheader()
                    w.writerow(config)
                os._exit(0)
            
            except Exception as e:
                self.error_raised = True
                exc_type, exc_obj, exc_tb = sys.exc_info()
                self.console_msg("Exception at line no: " + str(exc_tb.tb_lineno) +" "+str(e), level=logging.ERROR)
                os._exit(1)

    ##########################################################################
    # 
    # 
    #  exit_fullscreen
    # 
    # 
    # Add a way to exit fullscreen (E.g., press 'Esc')
    #
    ##########################################################################        
    
    def exit_fullscreen(self, event):
        self.window.wm_attributes('-fullscreen', False)

    ##########################################################################
    # 
    # 
    #  __init__
    # 
    # 
    ##########################################################################        

    def __init__(self):
 
        ############################################################################
        #
        #  valid_config_list list of config parameters loaded from .config file
        # 
        ############################################################################
        self.valid_config_list = {
            'settings_filename': self.settings_filename
        }

        
        #Wie heißen Sie?
        self.program_name = "MAOPhot"
        self.program_version = __version__
        self.program_name_note = "using photutils.psf"
        self.program_full_name = self.program_name + " " + self.program_version + " " + self.program_name_note
        self.config_file = os.getcwd() + "\\.config"

        # Check if there is a ./log dir
        self.logging_dir = os.getcwd() + "\\logs\\"

        # Check if directory exists
        if not os.path.exists(self.logging_dir):
            # Create the directory
            os.makedirs(self.logging_dir)

        #set the logger up
        self.our_logger = logging.getLogger(self.program_name + self.program_version + ".log")
        self.our_fh = logging.FileHandler(self.logging_dir + self.program_name + self.program_version + ".log", encoding='utf-8')
        self.our_logger.setLevel(logging.INFO)

        # create formatter and add it to the handlers
        self.our_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        self.our_fh.setFormatter(self.our_formatter)
        self.our_logger.addHandler(self.our_fh)
        
        self.load_config()

        #
        #
        #
        #               Set up GUI
        #
        #

        self.window = tk.Tk()
        
        # Maximize
        self.window.state('zoomed')

        self.screen_width = self.window.winfo_screenwidth()
        self.screen_height = self.window.winfo_screenheight()


        # Matplotlib settings
        matplotlib.rc('xtick', labelsize=7)
        matplotlib.rc('ytick', labelsize=7)

        self.window.bind('<Escape>', self.exit_fullscreen)
        
        # Bind the "X" button to the custom close function
        self.window.protocol("WM_DELETE_WINDOW", self.exit_app)

        self.window.title(self.program_full_name)

        #
        #
        #
        #               Menu
        #
        #

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
        self.viewmenu.add_separator()
        self.viewmenu.add_command(label="Invert", command=self.invert_image)
        self.viewmenu.add_separator()
        self.viewmenu.add_command(label="Refresh", command=self.display_image)
        self.menubar.add_cascade(label="View", menu=self.viewmenu)

        self.photometrymenu = tk.Menu(self.menubar, tearoff=0)
        self.photometrymenu.add_command(label="Find Peaks", command=self.find_peaks)
        self.ePSF_sub_menu = tk.Menu(self.photometrymenu, tearoff=False)
        self.ePSF_sub_menu.add_command(label="Create", command=self.create_ePSF)
        self.ePSF_sub_menu.add_command(label="Load...", command=self.load_ePSF)
        self.ePSF_sub_menu.add_command(label="Save As...", command=self.save_as_ePSF)
        self.ePSF_sub_menu.add_command(label="Clear", command=self.clear_ePSF)
        self.photometrymenu.add_cascade(label="Effective PSF", menu=self.ePSF_sub_menu)
        self.photometrymenu.add_separator()
        self.photometrymenu.add_command(label="Load Rejection List...", command=self.load_ePSF_rejection_list)
        self.photometrymenu.add_command(label="Save Rejection List As...", command=self.save_as_ePSF_rejection_list)
        self.photometrymenu.add_separator()

        self.photometrymenu.add_command(
            label="Iterative PSF Photometry", command=self.execute_iterative_psf_photometry)
        self.photometrymenu.add_separator()

        self.photometrymenu.add_command(
            label="Solve Image", command=self.solve_image)
        self.photometrymenu.add_separator()

        self.photometrymenu.add_command(
            label="Get Comparison Stars", command=self.get_comparison_stars)
        self.menubar.add_cascade(label="Photometry", menu=self.photometrymenu)
        
        self.multi_color_photo_menu = tk.Menu(self.menubar, tearoff=0)
        #self.multi_color_photo_menu.add_command(
        self.multi_color_photo_menu.add_command(label = "(B,V)", command=self.BV_multi_color_photometry)
        self.multi_color_photo_menu.add_command(label = "(V,R)", command=self.VR_multi_color_photometry)
        self.multi_color_photo_menu.add_command(label = "(V,I)", command=self.VI_multi_color_photometry)
        self.multi_color_photo_menu.add_separator()
        self.multi_color_photo_menu.add_command(label = "(B,V,R)", command=self.BVR_multi_color_photometry)
        self.multi_color_photo_menu.add_command(label = "(B,V,I)", command=self.BVI_multi_color_photometry)
        self.multi_color_photo_menu.add_command(label = "(V,R,I)", command=self.VRI_multi_color_photometry)
        self.multi_color_photo_menu.add_separator()
        self.multi_color_photo_menu.add_command(label="(B,V,R,I)", command=self.BVRI_multi_color_photometry)
        self.menubar.add_cascade(label="Multi Color Photometry", menu=self.multi_color_photo_menu)

        self.reportmenu = tk.Menu(self.menubar, tearoff=0)
        self.reportmenu.add_command(label="Single Image Photometry", command=self.generate_aavso_report_1image)
        self.multi_color_sub_menu = tk.Menu(self.reportmenu, tearoff=False)
        self.multi_color_sub_menu.add_command(label = "(B,V)",     command=self.BV_generate_aavso_report)
        self.multi_color_sub_menu.add_command(label = "(V,R)",     command=self.VR_generate_aavso_report)
        self.multi_color_sub_menu.add_command(label = "(V,I)",     command=self.VI_generate_aavso_report)
        self.multi_color_sub_menu.add_separator()
        self.multi_color_sub_menu.add_command(label = "(B,V,R)",   command=self.BVR_generate_aavso_report)
        self.multi_color_sub_menu.add_command(label = "(B,V,I)",   command=self.BVI_generate_aavso_report)
        self.multi_color_sub_menu.add_command(label = "(V,R,I)",   command=self.VRI_generate_aavso_report)
        self.multi_color_sub_menu.add_separator()
        self.multi_color_sub_menu.add_command(label = "(B,V,R,I)", command=self.BVRI_generate_aavso_report)
        self.reportmenu.add_cascade(label="Multi Color Photometry", menu=self.multi_color_sub_menu)

        self.menubar.add_cascade(label="Generate AAVSO Report", menu=self.reportmenu)

        self.window.config(menu=self.menubar)

        #
        # Layout left, center, and right frames
        #


        #
        #
        #
        #               Left [Side] Frame
        #
        #

        # We will lay image stretching sliders into the left_frame
        self.left_frame = tk.Frame(self.window, padx=__our_padding__, pady=__our_padding__)  # Left half of the window
        self.left_frame.grid(row=0, column=0, sticky=tk.NSEW)

        row = 0

        self.stretching_label = tk.Label(
            self.left_frame, text="Image Stretching:")
        self.stretching_label.grid(row=row, column=0, sticky=tk.NSEW)

        row += 1
        self.stretching_stringvar = tk.StringVar()
        self.stretching_stringvar.set("Asinh")
        self.stretching_dropdown = tk.OptionMenu(
            self.left_frame, self.stretching_stringvar, "None", "Square Root", "Log", "Asinh")
        self.stretching_dropdown.grid(row=row, column=0, sticky=tk.NW)


        row += 1
        self.stretching_stringvar.trace_add("write", self.display_image())
        # Histogram stretch sliders
        self.stretch_label = tk.Label(
            self.left_frame, text="Histogram Stretch Low/High:")
        self.stretch_label.grid(row=row, column=0, sticky=tk.NW)
        row += 1
        self.stretch_low = tk.Scale(
            self.left_frame, from_=0, to=100, orient=tk.HORIZONTAL, command=self.update_histogram_low)
        self.stretch_low.grid(row=row, column=0, columnspan=2, sticky=tk.NSEW)
        row += 1
        self.stretch_high = tk.Scale(
            self.left_frame, from_=0, to=100, orient=tk.HORIZONTAL, command=self.update_histogram_high)
        self.stretch_high.set(5)
        self.stretch_high.grid(row=row, column=0, columnspan=2, sticky=tk.NSEW)


        #
        #
        #
        #               Center [middle] Frame
        #
        #

        self.center_frame = tk.Frame(self.window, padx=__our_padding__, pady=__our_padding__)  # Center of the window
        self.center_frame.grid(row=0, column=1, sticky=tk.NSEW)

        row = 0
        self.filename_label = tk.Label(self.center_frame, text="FITS:" + image_file)
        self.filename_label.grid(row=row, column=0)  # Place label

        row += 1
        self.canvas = tk.Canvas(self.center_frame, bg='black')  # Main canvas

        # Place main canvas, sticky to occupy entire cell
        self.canvas.grid(row=row, column=0, sticky=tk.NSEW)

        # Expand main canvas column to fit whole  cell
        tk.Grid.columnconfigure(self.center_frame, 0, weight=1)

        # Give the canvas the most weight, it will do all the stretching
        tk.Grid.rowconfigure(self.center_frame, 1, weight=1)

        self.canvas_scrollbar_V = tk.Scrollbar(
            self.center_frame, orient=tk.VERTICAL)  # Main canvas scrollbars
        self.canvas_scrollbar_V.grid(sticky=tk.NSEW, row=row, column=1)

        row += 1
        self.canvas_scrollbar_H = tk.Scrollbar(
            self.center_frame, orient=tk.HORIZONTAL)
        self.canvas_scrollbar_H.grid(row=row, column=0)
        self.canvas_scrollbar_H.grid(sticky=tk.NSEW, row=row, column=0)
        self.canvas_scrollbar_H.config(command=self.canvas.xview)
        self.canvas_scrollbar_V.config(command=self.canvas.yview)
        self.canvas.config(xscrollcommand=self.canvas_scrollbar_H.set)
        self.canvas.config(yscrollcommand=self.canvas_scrollbar_V.set)

        row += 1
        # Console below
        self.console = tk.Text(self.center_frame, #height=40,
                               bg='black', fg='white', wrap='none')
        self.console.grid(sticky=tk.NSEW, row=row, column=0)


        self.console_scrollbar_V = tk.Scrollbar(
            self.center_frame, orient=tk.VERTICAL)  # Main canvas scrollbars
        self.console_scrollbar_V.grid(sticky=tk.NSEW, row=row, column=1)

        row += 1
        self.console_scrollbar_H = tk.Scrollbar(
            self.center_frame, orient=tk.HORIZONTAL)
        self.console_scrollbar_H.grid(row=row, column=0)
        self.console_scrollbar_H.grid(sticky=tk.NSEW, row=row, column=0)
        self.console_scrollbar_H.config(command=self.console.xview)
        self.console_scrollbar_V.config(command=self.console.yview)
        self.console.config(xscrollcommand=self.console_scrollbar_H.set)
        self.console.config(yscrollcommand=self.console_scrollbar_V.set)

        self.console_msg(self.program_full_name)
        self.console_msg("Ready")

        #
        #
        #
        #               Right [Side] Frame
        #
        #

        # Place right_frame 
        self.right_frame = tk.Frame(self.window, padx=__our_padding__, pady=__our_padding__)  # Right half of the window
        self.right_frame.grid(row=0, column=2, sticky=tk.NSEW)

        # Place label
        row = 0
        self.plotname_label = tk.Label(self.right_frame, text="Plot:")
        self.plotname_label.grid(row=row, column=0)  #row=0

        row += 1
        self.fig_psf = Figure()
        self.psf_plot = self.fig_psf.add_subplot(111, projection='3d')
        # PSF 3D plot canvas - Matplotlib wrapper for Tk
        self.psf_plot_canvas = FigureCanvasTkAgg(self.fig_psf, self.right_frame)
        self.psf_plot_canvas.draw()
        self.psf_canvas = self.psf_plot_canvas.get_tk_widget()
        self.psf_canvas.config(width=int(self.screen_width/8.5), height=int(self.screen_width/8.5))
        # Allocate small PSF canvas to a new grid inside the right_frame
        self.psf_canvas.grid(row=row, column=0) #row=1
        
        #
        #make another canvas for 2D plot of effectivePSF
        #
        # Place label
        row += 1
        self.ePSF_plotname_label = tk.Label(self.right_frame, text="PSF:")
        self.ePSF_plotname_label.grid(row=row, column=0)

        row += 1
        self.fig_ePSF, self.ePSF_plot = plt.subplots()
        self.ePSF_plot_canvas = FigureCanvasTkAgg(self.fig_ePSF, self.right_frame)
        self.ePSF_plot_canvas.draw()
        self.ePSF_canvas = self.ePSF_plot_canvas.get_tk_widget()
        self.ePSF_canvas.config(width=int(self.screen_width/8.5), height=int(self.screen_width/8.5))
        # Allocate small PSF canvas to a new grid inside the right_frame
        self.ePSF_canvas.grid(row=row, column=0) #row=3

        #
        #make another canvas for selected stars
        #
        # Place label
        row += 1
        self.selstars_title_label = tk.Label(self.right_frame, text="Selected Stars")
        self.selstars_title_label.grid(row=row, column=0) #row=4

        self.nrows = 5
        self.ncols = 5
        self.selstars_hspace = .55
        self.selstars_wspace = .55
        self.fig_selstars, self.selstars_plot = plt.subplots(nrows=self.nrows,
                                                              ncols=self.ncols,
                                                              squeeze=False)
        self.fig_selstars.canvas.mpl_connect('button_press_event', self.mouse_selstars_canvas_click)
        
        self.selstars_plot = self.selstars_plot.ravel()
        #

        row += 1
        self.selstars_plot_canvas = FigureCanvasTkAgg(self.fig_selstars, self.right_frame)
        plt.subplots_adjust(hspace=self.selstars_hspace, wspace=self.selstars_wspace)
        self.selstars_plot_canvas.draw()
        self.selstars_canvas = self.selstars_plot_canvas.get_tk_widget()

        self.selstars_canvas.config(width=int(self.screen_width/5), height=int(self.screen_width/5))
        self.selstars_canvas.config
        
        # Allocate small PSF canvas to a new grid inside the right_frame
        self.selstars_canvas.grid(row=row, column=0)#row=5

        #
        #make another canvas for selected stars
        #
        # Place label
        row += 1
        self.selstars_page_num_label = tk.Label(self.right_frame, text="Page:")
        self.selstars_page_num_label.grid(row=row, column=0)  #row=6

        row += 1
        separator_reject_buttons = ttk.Separator(self.right_frame, orient='horizontal')
        separator_reject_buttons.grid(row=row, pady=5, sticky=tk.EW) #row=7

        row += 1 #only for right_subframe
        self.right_subframe = tk.Frame(self.right_frame)
        self.right_subframe_sub0 = tk.Frame(self.right_subframe)
        self.right_subframe_sub1 = tk.Frame(self.right_subframe)
        self.right_subframe_sub2 = tk.Frame(self.right_subframe)

        self.right_subframe.columnconfigure(0, minsize=175)
        self.right_subframe.columnconfigure(2, minsize=175)

        self.submit_rejects_selstars_button = tk.Button(self.right_subframe_sub0, text="Submit", command=self.submit_rejects_selstars_list)
        self.submit_rejects_selstars_button.config(state=tk.DISABLED)
        self.submit_rejects_selstars_button.grid()

        self.back_selstars_button_label = tk.Label(self.right_subframe_sub1, text="<-----:")
        self.back_selstars_button_label.grid(row=0, column=0, sticky=tk.E)  # Place label

        self.back_selstars_button = tk.Button(self.right_subframe_sub1, text="Back", command=self.back_selstars_list)
        self.back_selstars_button.config(state=tk.DISABLED)
        self.back_selstars_button.grid(row=0, column=1, ipadx=12, padx=0, sticky=tk.E)

        self.forward_selstars_button = tk.Button(self.right_subframe_sub1, text="Forward", command=self.forward_selstars_list)
        self.forward_selstars_button.config(state=tk.DISABLED)
        self.forward_selstars_button.grid(row=0, column=2, padx=0, sticky=tk.W)

        self.forward_selstars_button_label = tk.Label(self.right_subframe_sub1, text=":----->")
        self.forward_selstars_button_label.grid(row=0, column=3, sticky=tk.W)  # Place label

        clear_rejects_selstars_button = tk.Button(self.right_subframe_sub2, text="Clear", command=self.clear_rejects_selstars_list)
        clear_rejects_selstars_button.grid()

        self.right_subframe_sub0.grid(row=0, column=0, sticky=tk.W)
        self.right_subframe_sub1.grid(row=0, column=1)
        self.right_subframe_sub2.grid(row=0, column=2, sticky=tk.E)

        self.right_subframe.grid(row=row, column=0, sticky=tk.N)#row=8

        # when shrinking, don't let rows with "Page" label, seperator, and buttons 
        # get cut off
        tk.Grid.rowconfigure(self.right_frame, 0, weight=1) 
        tk.Grid.rowconfigure(self.right_frame, 1, weight=1)
        tk.Grid.rowconfigure(self.right_frame, 2, weight=1)
        tk.Grid.rowconfigure(self.right_frame, 3, weight=1)
        tk.Grid.rowconfigure(self.right_frame, 4, weight=1)
        tk.Grid.rowconfigure(self.right_frame, 5, weight=1)

        #
        # end of right (side) frame layout
        #
        #

        # Assign weights so that only column 1 shrinks
        tk.Grid.columnconfigure(self.window, 1, weight=1)

        # Row 1 (console) expands the most
        tk.Grid.rowconfigure(self.window, 0, weight=1)

        # Update layout to calculate dimensions
        self.window.update_idletasks()

        #
        # laumch_settings; pops up settings window and initializes all settings
        #
        #
        # But some parameter settings need initial values first
        #
        self.catalog_stringvar = tk.StringVar()
        self.catalog_stringvar.set("AAVSO")
        self.fitter_stringvar = tk.StringVar()
        self.fitter_stringvar.set("TRF LS")
        self.display_all_objects = tk.StringVar(None, 0) #init to display user objects only
        self.use_gaussian_prf_model = tk.StringVar(None, 1) #init to gaussian prf  
        self.generate_residual_image = tk.BooleanVar()
        self.generate_residual_image.set(False) #init to not generate residual
        self.auto_behavior = tk.BooleanVar()
        self.auto_behavior.set(False)

        

        self.launch_settings()

        # Now that the widjets are defined we can create this dict for saving data:
        
        ############################################################################
        #
        #  valid_parameter_list facilitates loading from and saving to a settings file
        #  buttons in popup window: Settings
        # 
        ############################################################################
        self.valid_parameter_list = {
            'find_peaks_npeaks_entry': self.find_peaks_npeaks_entry,
            'fit_width_entry': self.fit_width_entry,
            'max_ensemble_magnitude_entry': self.max_ensemble_magnitude_entry,
            'fwhm_entry': self.fwhm_entry,
            'oversampling_entry': self.oversampling_entry,
            'star_detection_threshold_factor_entry': self.star_detection_threshold_factor_entry,
            'photometry_iterations_entry': self.photometry_iterations_entry,
            'sharplo_entry': self.sharplo_entry,
            'matching_radius_entry': self.matching_radius_entry,
            'aavso_obscode_entry': self.aavso_obscode_entry,
            'telescope_entry': self.telescope_entry,
            'user_note_entry': self.user_note_entry,
            'tb_bv_entry': self.tb_bv_entry,
            'tb_br_entry': self.tb_br_entry,
            'tb_bi_entry': self.tb_bi_entry,
            'tv_bv_entry': self.tv_bv_entry,
            'tv_vr_entry': self.tv_vr_entry,
            'tr_vr_entry': self.tr_vr_entry,
            'tr_ri_entry': self.tr_ri_entry,
            'ti_ri_entry': self.ti_ri_entry,
            'tv_vi_entry': self.tv_vi_entry,
            'ti_vi_entry': self.ti_vi_entry,
            'tr_vi_entry': self.tr_vi_entry,
            'tbv_entry': self.tbv_entry,
            'tbr_entry': self.tbr_entry,
            'tbi_entry': self.tbi_entry,
            'tvr_entry': self.tvr_entry,
            'tri_entry': self.tri_entry,
            'tvi_entry': self.tvi_entry,
            'tb_bv_err_entry': self.tb_bv_err_entry,
            'tb_br_err_entry': self.tb_br_err_entry,
            'tb_bi_err_entry': self.tb_bi_err_entry,
            'tv_bv_err_entry': self.tv_bv_err_entry,
            'tv_vr_err_entry': self.tv_vr_err_entry,
            'tr_vr_err_entry': self.tr_vr_err_entry,
            'tr_ri_err_entry': self.tr_ri_err_entry,
            'ti_ri_err_entry': self.ti_ri_err_entry,
            'tv_vi_err_entry': self.tv_vi_err_entry,
            'ti_vi_err_entry': self.ti_vi_err_entry,
            'tr_vi_err_entry': self.tr_vi_err_entry,
            'tbv_err_entry': self.tbv_err_entry,
            'tbr_err_entry': self.tbr_err_entry,
            'tbi_err_entry': self.tbi_err_entry,
            'tvr_err_entry': self.tvr_err_entry,
            'tri_err_entry': self.tri_err_entry,
            'tvi_err_entry': self.tvi_err_entry,
            'linearity_limit_entry': self.linearity_limit_entry,
            'catalog_stringvar': self.catalog_stringvar,
            'vizier_catalog_entry': self.vizier_catalog_entry,
            'fitter_stringvar': self.fitter_stringvar,
            'astrometrynet_entry': self.astrometrynet_entry,
            'astrometrynet_key_entry': self.astrometrynet_key_entry,
            'object_kref_entry': self.object_kref_entry,
            'object_sel_comp_entry': self.object_sel_comp_entry,
            'object_min_snr_entry': self.object_min_snr_entry,
            'object_name_entry': self.object_name_entry,
            'object_name_alpha_entry': self.object_name_alpha_entry,
            'object_name_delta_entry': self.object_name_delta_entry,
            'object_notes_entry': self.object_notes_entry,
            'display_all_objects': self.display_all_objects,
            'use_gaussian_prf_model': self.use_gaussian_prf_model,
            'generate_residual_image': self.generate_residual_image,
            'auto_behavior': self.auto_behavior,
            'filter_entry': self.filter_entry,
            'max_qfit_entry': self.max_qfit_entry,
            'moffat_beta_entry': self.moffat_beta_entry,
            'min_separation_bias_entry': self.min_separation_bias_entry,
            'extinction_B_entry': self.extinction_B_entry,
            'extinction_V_entry': self.extinction_V_entry,
            'extinction_I_entry': self.extinction_I_entry,
            'extinction_R_entry': self.extinction_R_entry,
            'extinction_C_entry': self.extinction_C_entry
            }

        # if .config had a valid settings_filename, then load that one in
        if os.path.exists(self.settings_filename):
            self.open_settings(self.settings_filename)
            self.console_msg("Loaded settings from " + str(self.settings_filename))
            
        self.console_msg("Ready")

        tk.mainloop()

myGUI = MyGUI()

