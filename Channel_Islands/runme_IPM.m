function runme_IPM(Site,Species,Type,F_initial)
% Run state-space IPM model for fitting to PISCO datasets
% Methods originally documented in White et al. (2016) Ecological
% Applications 26:2675-2692

% This master file reads in the data and parameters, and calls the
% model-fitting routine

% As currently set up this code is designed to be used for a single species
% at a time, and to be used for fitting model to pre-MPA data or to
% post-MPA data (the latter not shown in White et al. 2016)

% This version of the code was forked from the version used by White et al.
% 2016 in 2019, and modified to operate on PISCO subtidal fish data from
% the Santa Barbara Channel Islands

%----------------------------------------
% Input arguments:

% Site: The name of the PISCO study site

% Species: 
% SMYS (Sebastes mystinus, blue rockfish) 
% SCAR (S. carnatus, gopher rockfish)
% SATR (S. atrovirens, kelp rockfish [not used in White et al. PLoS One])
% SMEL (S. melanops, black rockfish - this species no longer used bc central coast is the southern edge of distribution, few large fish are observed)
% PCLA (Paralabrax clathratus, kelp bass)
% SPUL (Semicossyphus pulcher, California sheephead)

% Type: pre2007, post2007 % describes the type of fitting to be done.
% For forward runs (not fitting) starting in 2007, use pre2007 (these
% options operate on the Central Coast data in White et al. 2016)
% For Channel Islands data, use 'post2003' 

% F: For post2007 or post2003, there is an option to input the prior on F.
% If not used the code will draw the prior from the parameters file.

% A metadata structure file (Meta) will be associated with particular runs and
% keep track of the various filenames, species, sites, etc. used in
% particular runs.  This facilitates passing information between m-files.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate parameter file and read it in.
if ~exist('F_initial','var'); F_initial=NaN; end
Meta = IPM_parameters_ChannelIslands(Species,Site,Type,F_initial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Read in the data
% PISCO data is in the directory PISCO_data_2015

read_PISCOdata_ChannelIslands(Meta.savename)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: post-processing, plotting results
postproc_pisco_rockfish_fit(Meta.savename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



