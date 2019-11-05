function[]=runme_IPM_SMYS_AM(varargin)
% Run Central Coast IPM model

% This master file reads in the data and parameters, and calls the
% model-fitting routine

%% As currently set up this code is designed to be used for a single species
% at a time, and to be used for fitting model to pre-2007 data

% Site: can be an MPA region (Andrew Molera)

% Species:
% SMYS (Sebastes mystinus, blue rockfish)

% Type: pre2007 % describes the type of fitting to be done.

% A metadata structure file (Meta) will be associated with particular runs and
% keep track of the various filenames, species, sites, etc. used in
% particular runs.  This facilitates passing information between m-files.

% Chain: An integer indicating which MCMC chain is being run in this
% particular call.

if(nargin>0)
Species = varargin{1};
else
Species = 'SMYS';
Site = 'Andrew_Molera';
Type = 'pre2007';
end
disp(Species)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate parameter file and read it in.
Meta = IPM_parameters_fit(Species,Site,Type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Read in the data
% PISCO data is in the directory PISCO_data_2015

read_PISCOdata_cluster(Meta)

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
rockfish_fit_pisco_cluster(Meta,Chain,0)

% calls: IPM_histo.m (to convert count data into histogram)
%        do_IPM.m (to run the IPM)
%        makeSimpVec (to make the Simpson vector for integration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



