# PSF-Photometry
Welcome to MAOPhot 1.0.0, a PSF Photometry tool using Astropy 6.0.0 and 
photutils.psf 1.10.0

This program was derived from MetroPSF by Maxym Usatov. 

MAOPhot calculates stellar magnitudes from Fit (*.fit, *.fits, *.fts) formatted
digital photographs using PSF photometry. It produces an extended AAVSO 
(American Association of Variable Star Observers)
report which can be submitted to the AAVSO using their online tool WebObs
(http://www.aavso.org/webobs).

MAOPhot uses the PSF (point spread function) Photometry method exclusively. 

PSF modeling is well suited for measuring stellar magnitudes in crowded fields,
or the magnitude of a star that has a close companion, e.g., Z Tau. 
(See https://www.aavso.org/lpv-double-trouble-campaign-0)

MAOPhot is written in Python using Astropy [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
(a common core package for astronomy). MAOPhot also uses photutils.psf.
See "PSF Photometry" (https://photutils.readthedocs.io/en/stable/psf.html)
which describes many of the classes and methods used in MAOPhot.

MAOPhot has been redesigned for AAVSO reporting only and includes, but not
limited to the following enhancements:

o	uses Astropy 6.0.0 and photutils.psf 1.10.0 (all packages required have
    version specified in requirements.txt file)

o	generation of Effective PSF model (EPSF model), and ability to create a
    ‘rejection list’ of stars that the user can select that will not be part of
    the EPSF model 

o	option to use an Integrated Gaussian PRF (Pixel Response Function) as model

o	PSF Photometry using an iterative algorithm to perform point spread
    function photometry

o	PSF Photometry using a non-iterative algorithm to perform point spread
    function

o	PSF Photometry using an ensemble of comparison stars or a single comp star

o	generation of Two-Color Photometry (B, V), (V, R) or (V, I), and Single
    Image Photometry reports in AAVSO extended format 

o	use of telescope Transformation Coefficients (needed for Two Color
    Photometry)

o	user can specify check star and list of comp stars to use

o	a Radio Button option to display all AAVSO comp stars as AAVSO label
    numbers and any found VSX objects in image field or only VSX and AAVSO comp
    stars specified in settings 

o	intermediate results are saved as .csv files

o	user can optionally enter a AAVSO Chart ID when retrieving comparison star
    data

o	if intermediate results include found comp stars or VSX objects associated
    with an image being loaded, then these objects are displayed

See HelpFile.pdf for more information.
