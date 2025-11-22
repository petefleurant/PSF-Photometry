# PSF-Photometry
Welcome to MAOPhot 1.1.8, a PSF Photometry tool using Astropy 7.0.1 and Photutils 2.3.0

## Version 1.1.8 Changes

1. Button-2 click centers image
2. Added View → Invert; this inverts image
3. Added Effective PSF → Load... and Effective PSF → Save As...
4. ENSEMBLE calculated for Single Image Photometry
5. ENSEMBLE outliers displayed for Single Image Photometry
6. Button-1 + Shift drags image
7. Mouse wheel zooms image in and out
8. DAOStarFinder replaces find_peaks, replaces IRAFStarFinder
9. SourceGrouper is used in IterativePSFPhotometry
10. In APASS DR10, remove any entries with Johnson (V) > maglimit until cgi-bin/apass_dr10_download.pl is fixed
11. In "Iterative PSF Photometry" remove rows with flags != 0, qfit > 0.07; APASS DR10 only keep lowest qfit match; Replace IRAFStarFinder with DAOStarFinder
12. Added Moffat with beta parameter as option for PSF model
13. Added option to generate residual image

**Note:** Version 1.1.7 has been merged into 1.1.8

## Version 1.1.6 Changes

1. Ensemble is supported in Single Image Photometry
2. Select comp stars input list is stripped of whitespace and uniquefied
3. Added Transformation Coefficients Error figures to AAVSO report
4. Added Effective PSF → Load... and Effective PSF → Save As...

## Version 1.1.5 Changes

1. Image can be dragged by Shift + Left-mouse-button
2. Image can be zoomed in and zoomed out with mouse wheel

## Version 1.1.4 Changes

1. Support for APASS DR10 (along with GAIA DR2, APASS DR9)
2. Added Setting: Max qfit ('qfit' is quality of PSF fit, lower number is better fit)
3. Added Setting: Min Separation Factor (this × FWHM = the minimum distance in pixels such that any two sources separated by less than this distance will be placed in the same group)
4. Added Setting: "From Fits" checkbox for CCD Filter; when unchecked user can manually override FILTER value in FITS header
5. "Find Peaks" and "Iterative PSF Photometry" now use DAOStarFinder
6. In APASS DR10, remove any entries with Johnson (V) > maglimit until cgi-bin/apass_dr10_download.pl is fixed
7. Use mouse wheel and shift-mouse wheel to scroll

## Overview

PSF photometry models the star's light distribution as a mathematical function called the Point Spread Function (PSF), which describes how the star's light spreads across the detector. The PSF is then fitted to the star to determine its total flux, accounting for overlaps with nearby stars.

This program was derived from "MetroPSF" by Maxym Usatov. It has been renamed and extended to produce AAVSO reports exclusively and to facilitate generating an effective PSF for PSF photometry.

MAOPhot calculates stellar magnitudes from FITS formatted digital photographs using PSF photometry. It produces an extended AAVSO (American Association of Variable Star Observers) report (https://www.aavso.org/aavso-extended-file-format) which can be submitted to AAVSO using their online tool WebObs (https://www.aavso.org/webobs).

MAOPhot uses PSF (Point Spread Function) Photometry exclusively.

MAOPhot is written in Python using Astropy (a common core package for astronomy). MAOPhot also uses Photutils. See "PSF Photometry" documentation which describes many of the classes and methods used in MAOPhot: https://photutils.readthedocs.io/en/stable/user_guide/psf.html

## Key Features

MAOPhot has been redesigned for AAVSO reporting only and includes, but is not limited to, the following enhancements:

- Uses Astropy 7.0.1 and Photutils 2.3.0

- Generation of an Effective PSF model (ePSF model), and the ability to create a 'rejection list' of stars that the user can select that will not be part of the ePSF model generated

- Option to use a Gaussian PRF (Pixel Response Function) or Moffat as model

- Uses Iterative PSF Photometry, an iterative version of PSF Photometry where new sources are detected in the residual image after the fit sources are subtracted

- PSF Photometry using an ensemble of comparison stars or a single comp star

- Generation of Two-Color Photometry (B, V), (V, R) or (V, I), and Single Image Photometry reports in AAVSO extended format

- Use of telescope Transformation Coefficients (needed for Two Color Photometry)

- User can specify check star and list of comp stars

- Manually select a star for measurement

- Intermediate results are saved as .csv files

- Optionally enter an AAVSO Chart ID when retrieving comparison star data

## Documentation

See HelpFile.pdf for more information.