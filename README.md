# PSF-Photometry

MAOPhot 1.2.0
Welcome to MAOPhot 1.2.0, a PSF photometry tool built with Astropy 7.1.1 and Photutils 2.3.0.
Version 1.2.0 Changes
1)	New feature: Multi-color photometry
Adds support for 3-color and 4-color transformed reports. Analysis uses the AAVSO-recommended iterative solution, following the approach used by AAVSO TransformApplier (v2.7.1).
MAOPhot applies the AAVSO-recommended transformation coefficients.
MAOPhot requires the following coefficients:
coeff	B		V		R		I
-------------------------------------
BVRI	Tb_bv	Tv_bv	Tr_vi	Ti_vi
BVR		Tb_bv	Tv_bv	Tvr	
BVI		Tb_bv	Tv_bv			Ti_vi
VRI				Tv_vi	Tr_vi	Tvi
BV		Tb_bv	Tv_bv		
VR				Tv_vr	Tr_vr	
VI				Tv_vi			Tvi

2)	Upgraded packages:
a.	Astropy 7.1.1   7.2
b.	NumPy 2.3.4   2.4
c.	Pandas 2.3.3   3.0.1
d.	Pillow 11.0.0   12.1.1
e.	SciPy 1.16.2   1.17.1
f.	tqdm 4.67.1   4.67.3

3)	Settings changes:
a.	Added Oversampling, used when building the ePSF model (default: 2).
b.	Added additional transformation coefficients (e.g., Tb_br, Tb_bi, Tr_ri).
c.	Added User Note.
d.	Added Date-Obs (UTC).
e.	Date-Obs (JD) and Date-Obs (UTC) are read-only and are obtained from the FITS file.
4)	Saved the current oversampling value in the exported ePSF file.
5)	Corrected the single-comparison-star error calculation (uses 1/SNR).
6)	Added Minimum SNR to Settings. If any comparison star or check star has an SNR below the minimum, that star is excluded.
7)	Corrected a case where the check star was not labeled.
8)	Reports now record the correct transformation coefficients in the Notes section of the AEFF (AAVSO Extended File Format) report.
9)	AAVSO report names are now derived from the Master Report name rather than the loaded FITS file.


## Overview

MAOPhot calculates stellar magnitudes from FITS-formatted images using PSF photometry. It produces an AAVSO Extended File Format (AEFF) report that can be submitted to the AAVSO via WebObs.
MAOPhot uses PSF (point spread function) photometry exclusively.
MAOPhot is written in Python using Astropy (a common core package for astronomy) and Photutils. For background on the Photutils classes and methods used by MAOPhot, see the Photutils “PSF Photometry” user guide.
Key Features
MAOPhot is designed for AAVSO reporting and includes the following features:
•	Built with Astropy and Photutils
•	Generates an effective PSF (ePSF) model following Anderson and King (2000; PASP 112, 1360), with optional creation of a rejection list of stars to exclude from the ePSF build
•	Option to use a Circular Gaussian PRF (Pixel Response Function) as a model
•	Option to use a Moffat PSF model
•	Supports iterative PSF photometry (Photutils IterativePSFPhotometry), where new sources are detected in the residual image after fitted sources are subtracted—useful for crowded fields
•	PSF Photometry using an ensemble of comparison stars or a single comp star
•	Generates multi-color photometry reports (BV, VR, VI, BVR, BVI, VRI, and BVRI) using transformation coefficients (TCs) in AAVSO extended format; single-image photometry reports are also supported
•	Use of telescope Transformation Coefficients (needed for Color Photometry)
•	User can specify check star and list of comp stars
•	Manually select a star for measurement
•	Intermediate results are saved as .csv files
•	Optionally enter an AAVSO Chart ID when retrieving comparison star data
•	When image file is loaded, optionally elect to find weighted median Moffat β value and weighted median FWHM value

See HelpFile.pdf for more information.
