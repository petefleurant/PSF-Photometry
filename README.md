# PSF-Photometry

Welcome to MAOPhot 0.1, a PSF Photometry tool using Astropy and Photutils.psf

MAOPhot calculates stellar magnitudes from Fit formatted digital photographs. 
It produces an extended AAVSO (American Association of Variable Star Observers)
report which can be submitted to the AAVSO using the online tool WebObs 
(http://www.aavso.org/webobs).
There are many photometry measuring programs available such as VPhot 
(http://www.aavso.org/vphot) and AstroImageJ (University of Louisville). 
VPhot uses the aperture photometry method.
MAOPhot uses the PSF Photometry method exclusively. PSF (point spread function)
modeling is well suited for measuring stellar magnitudes in crowded fields, or
the magnitude of a star that has a close companion, e.g., Z Tau. 
(See https://www.aavso.org/lpv-double-trouble-campaign-0)

MAOPhot is written in Python. It uses many Python 'astropy' 
(https://www.astropy.org/) libraries. The astropy package contains key 
functionality and common tools for performing astronomy and astrophysics with
Python. Included in the package is Photutils.psf. See "PSF Photometry" 
(https://photutils.readthedocs.io/en/stable/psf.html) which describes many of 
the classes and methods used in MAOPhot.

This program was derived from MetroPSF by Maxym Usatov.
It has been redesigned for AAVSO reporting only and includes, but not limited
 to the following enhancements:
- Generation of Effective PSF model, and ability to create a ‘rejection list.’
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

See HelpFile.pdf for more information.
