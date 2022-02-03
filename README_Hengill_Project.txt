============================ README Hengill Project===========================
A. Duran 
Feb 2021

Here you will find the main descriptions and codes for the Hengill project.
For the explicit codes of SIMULPS14 and FDtimes ask directly to Tobias or Edi Kissling. 

Have a look to the files "Magnitude_Comput_Iceland_AD.pdf" and "TOMO_GRAL_PROCEDURE.pdf" for a general overview of the problems faced and the procedures. 


Description of the FOLDERS:

COSEISMIQ_DOC
Documents about the coseismiq project (the project from where this work forms part). 

ISOR_catalog
Seismic catalogue from Icelandic Geosurvey (ISOR)

Python_Obspy
Various codes (they are more like exercises) to download the data from the server, plot waveforms, apply filters, compute the local magnitude of the events, among others.

Tomography
Codes to perform the tomography and synthetic tests. But, for the SIMULPS14 code you better request it to Tobias or Edi as mentioned.  

Min1D_Model
Tests performed with VELEST software to obtain the minimum 1D velocity model (Vp and Vs) and station corrections for the Hengill region. 

Automatic_catalogue_analysis
Study of the seismicity and b-values for the automatic catalogue on Hengill. The catalogue used for location of events the min 1D models (from simultaneous P+S inversion) obtained in the previous folder. Here you find various catalogues according to different quality criteria set by Tobias.  Check the folder QC_AUT_Data/Data_Tobias/COSEISMIQ_Automatic_Catalog_for_Alejandros_manuscript

The criteria are:

1) High-Quality HQ subset:
gap < 180deg, nobs >= 10, RMS<=0.09 s (3xminimum 1D RMS), distance to closest station <= 7.5 km, event-score >= -1.0
 
2) Medium quality MQ subset:
gap < 250deg, nobs >= 10, RMS<=0.20 s, distance to closest station <= 10 km, event-score >= -50 