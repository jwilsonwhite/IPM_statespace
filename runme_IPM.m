function runme_IPM(Site,Species,Type,Fprior)
% Run Central Coast IPM model

% This master file reads in the data and parameters, and calls the
% model-fitting routine

% As currently set up this code is designed to be used for a single species
% at a time, and to be used for fitting model to pre-2007 data or to
% post-2007 data (not shown in White et al. PLoS ONE paper)
% There is some prelim code for forecasting post-2007 dynamics, but this is
% not fully implemented.

% Site: can be an MPA region (Pt_Lobos, White_Rock, or Big_Creek) or, for
% post-2007 runs, an individual site

% Species: 
% SMYS (Sebastes mystinus, blue rockfish) 
% SCAR (S. carnatus, gopher rockfish)
% SATR (S. atrovirens, kelp rockfish [not used in White et al. PLoS One])
% SMEL (S. melanops, black rockfish - this species no longer used bc central coast is the southern edge of distribution, few large fish are observed)

% Type: pre2007, post2007 % describes the type of fitting to be done.
% For forward runs (not fitting) starting in 2007, use pre2007

% For post2007, have to input the prior on F (presumably either close to
% zero (for MPA) or derived from pre2007 estimate)

% A metadata structure file (Meta) will be associated with particular runs and
% keep track of the various filenames, species, sites, etc. used in
% particular runs.  This facilitates passing information between m-files.

% Chain: An integer indicating which MCMC chain is being run in this
% particular call.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate parameter file and read it in.
Meta = IPM_parameters(Species,Site,Type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Read in the data
% PISCO data is in the directory PISCO_data_2015

read_PISCOdata(Meta.savename)

% reads in the data, saves it in a structure file.  Also counts up number
% of transects in each survey.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Fit the data to obtain historical values of F
pauseflag = false; % change this if you want to pause the fitting routine for debugging

Chain = NaN; % if you are calling rockfish_fit_pisco directly from this file, 
             % it will run with the default number of chains specified in
             % the IPM_parameters file. If you want to manually parallelize
             % the chains, call rockfish_fit_pisco separately on the
             % command line and specify which Chain you are running.
rockfish_fit_pisco(Meta.savename,Chain,pauseflag)

% calls: IPM_histo.m (to convert count data into histogram)
%        do_IPM.m (to run the IPM)
%        makeSimpVec (to make the Simpson vector for integration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



