# PSF-Photometry
Welcome to MAOPhot 1.1.2, a PSF Photometry tool using Astropy 6.1.6 and Photutils 2.0.2

PSF photometry models the star's light distribution as a mathematical function called the Point Spread Function (PSF), which describes how the star's light spreads across the detector. The PSF is then fitted to the star to determine its total flux, accounting for overlaps with nearby stars.

This program was derived from “MetroPSF” by Maxym Usatov.  It has been renamed and extended to produce AAVSO reports exclusively and to facilitate generating an effective PSF for PSF photometry. 

MAOPhot calculates stellar magnitudes from FIT formatted digital photographs using PSF photometry. It produces an extended AAVSO (American Association of Variable Star Observers)
report (https://www.aavso.org/aavso-extended-file-format) which can be submitted to AAVSO using their online tool WebObs (https://www.aavso.org/webobs).

MAOPhot uses the PSF (point spread function) Photometry exclusively. 

MAOPhot is written in Python using Astropy (a common core package for astronomy). MAOPhot also uses Photutils. See "PSF Photometry" which describes many of the classes and methods used in MAOPhot. (https://photutils.readthedocs.io/en/stable/user_guide/psf.html)

MAOPhot has been redesigned for AAVSO reporting only and includes, but is not limited to the following enhancements:

o	uses Astropy 6.1.6 and Photutils 2.0.2 

o	generation of an Effective PSF model (EPSF model), and the ability to create a ‘rejection list’ of stars that the user can select that will not be part of the EPSF model generated 

o	option to use a Gaussian PRF (Pixel Response Function) as model

o	Uses Iterative PSF Photometry, an iterative version of PSF Photometry where new sources are detected in the residual image after the fit sources are subtracted

o	PSF Photometry using an ensemble of comparison stars or a single comp star

o	generation of Two-Color Photometry (B, V), (V, R) or (V, I), and Single Image Photometry reports in AAVSO extended format 

o	use of telescope Transformation Coefficients (needed for Two Color Photometry)

o	user can specify check star and list of comp stars 

o	manually select a star for measurement

o	intermediate results are saved as .csv files

o	optionally enter an AAVSO Chart ID when retrieving comparison star data

See HelpFile.pdf for more information.
