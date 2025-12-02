# PSF-Photometry

Welcome to MAOPhot 1.1.8, a PSF Photometry tool using Astropy 7.1.1 and Photutils 2.3.0

## Version 1.1.8 Changes

1)	Intermediate results display qfit metrics
2)	Button-2 click centers image
3)	Button-1 + Shift drags image
4)	Mouse wheel zooms image in and out
5)	Added View-->Invert; this inverts image
6)	Added Effective PSF-->Load... and Effective PSF-->Save As...
7)	ENSEMBLE calculated for Single Image Photometry
8)	ENSEMBLE outliers displayed for Single Image Photometry
9)	Replace IRAFStarFinder with DAOStarFinder
10)	SourceGrouper is used in IterativePSFPhotometry
11)	In APASS DR10, remove any entries with Johnson (V) > maglimit until cgi-bin/apass_dr10_download.pl is fixed
12)	Added Moffat with beta parameter as option for PSF model
13)	Added option to generate residual image
14)	MAOPhot 1.1.7 has been merged into MAOPhot 1.1.8


**Note:** Version 1.1.7 has been merged into 1.1.8


## Overview

MAOPhot calculates stellar magnitudes from FITS formatted digital photographs using PSF photometry. It produces an extended AAVSO (American Association of Variable Star Observers) report (https://www.aavso.org/aavso-extended-file-format) which can be submitted to AAVSO using their online tool WebObs (https://www.aavso.org/webobs).
MAOPhot uses PSF (point spread function) photometry exclusively.

MAOPhot is written in Python using Astropy (a common core package for astronomy). MAOPhot also uses Photutils. See "PSF Photometry" which describes many of the classes and methods used in MAOPhot: https://photutils.readthedocs.io/en/stable/user_guide/psf.html

## Key Features

MAOPhot has been redesigned for AAVSO reporting only and includes, but is not limited to, the following enhancements:
•	Uses Astropy 7.1.1 and Photutils 2.3.0
•	Generation of an Effective PSF model (EPSF model) following the prescription of Anderson and King (2000; PASP 112, 1360), with the ability to create a 'rejection list' of stars that the user can select that will not be part of the EPSF model generated
•	Option to use a Circular Gaussian PRF (Pixel Response Function) as a model
•	Option to use a Moffat PSF model
•	Uses Iterative PSF Photometry - an iterative version of PSF Photometry (IterativePSFPhotometry class from Photutils) where new sources are detected in the residual image after the fit sources are subtracted. This is particularly useful for crowded fields where faint sources are very close to bright sources and may not be detected in the first pass
•	PSF Photometry using an ensemble of comparison stars or a single comp star
•	Generation of Two-Color Photometry (B, V), (V, R) or (V, I), and Single Image Photometry reports in AAVSO extended format
•	Use of telescope Transformation Coefficients (needed for Two Color Photometry)
•	User can specify check star and list of comp stars
•	Manually select a star for measurement
•	Intermediate results are saved as .csv files
•	Optionally enter an AAVSO Chart ID when retrieving comparison star data
## Documentation
See HelpFile.pdf for more information.
