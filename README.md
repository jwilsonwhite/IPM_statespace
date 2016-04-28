# IPM_statespace
State-space integral projection models

This is the code described in White et al., Ecological Applications (2016)
Fitting state-space integral projection models to time series data.
All code is written for Matlab.

There are two runme files that are the primary wrappers for all the code.

runme_IPM_mockdata.m calls all of the necessary files for running the simulated blue rockfish data fits reported in the manuscript.
The simulated data is actually generated by rockfish_mockdata.m

runme_IPM.m calls all of the necessary files for fitting the IPM to the Pt Lobos rockfish datasets.

The original PISCO data (for the selected species and sites) is in directory PISCO_data_2015
The fits to those data are in runs_for_publication_pre2007_July2015

The simulated data used in the manuscript is in mockdata_May2015
The fits to those simulated data are in mockdata_fits_June2015

The data from Rick Starr's hook & line survey data used to generate the observation ogive is in Starr_hookline_survey_data

Post-processing and figure-generating code includes:
gridsize_figs.m
mockdata_plots.m
postproc_pisco_rockfish_fit.m


