function runme_IPM_Tomales(Modelstart,Modelstop,Cluster)

% This file describes all the code needed to run the IPM model for Tomales
% oyster data. This code was modified from that in White et al. (2016)
% Ecological Applications, so there are some legacy features that were
% necessary for that implementation that are not used here.

% The following structure file (Meta) will be associated with particular runs and
% keep track of the various filenames, species, sites, etc. used in
% particular runs.  This way you only need to tell each set of code which
% metadata file to use.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Model types:
Model_str = define_model_types;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Assorted Metadata:

Meta.T = 1:10; % full # years for model
Meta.Tdatafull = 5:10; % full timespan of data
Meta.Tpre =   1:4; % time prior to data that model will run
Meta.Tdata =  [5 6 7 10]; % years with data       
    
% Sites: 
Meta.Sites = {'E1','E2','E3','E4','W2','W3','W4'}; % cell array.  Insert NaN if you want all the sites
Meta.Nsites = 7;
Meta.Distances = [ 6.25, 8.25, 12.25, 15.5, 8, 12, 16]; % Distance of each site from the mouth of the bay

% Parameters that are fixed for the estimation step:
Meta.fixparm = []; % currently none?
              
% IPM parameters
Meta.IPM.meshsize = 100;
Meta.IPM.meshmin = 0;
Meta.IPM.meshmax = 100;
Meta.IPM.x = linspace(Meta.IPM.meshmin,Meta.IPM.meshmax,Meta.IPM.meshsize);
Meta.IPM.meshdiff = diff(Meta.IPM.x(1:2));
Meta.IPM.edges = Meta.IPM.x - Meta.IPM.meshdiff/2;

% MCMC parameters
Meta.MCMC.M = 5e4; % length of chain
Meta.MCMC.chains = 3; % number of chains
Meta.Q = 100; % number of particles (based on plotting std(L) vs. Q - see code in fit_tomales_IPM)

%Meta.fit_savename = 'SMYS_mock_fit_19July2013.mat'; % string
Meta.savename = strcat('Tomales_20Jan2018_metadata.mat');
Meta.fit_name = '20Jan2018_fit.mat';
%Meta.savename = 'SMYS_Mock2_19July2013_metadata.mat';

% Load in data & scale correctly:
[N,Meta.Sample_size,Meta.SiteOrder] = sort_oysterdata_IPM(Meta.IPM);
% N is the dataset for fitting. It has dimensions length x site x year. 

% Expand the Sample_size matrix to include years with no data
SS = nan(Meta.Nsites,length(Meta.T));
SS(:,Meta.Tdata) = Meta.Sample_size;
Meta.Sample_size = SS;

save(Meta.savename,'Meta') % save the metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over each model to find fit

% Make this a parfor loop for actual modeling
%PoolObj = parpool(Cluster,3);
for k = Modelstart:Modelstop

    disp(strcat('Model ',num2str(k)))
    
    pauseflag = false;
fit_tomales_IPM(N,Meta.savename,Model_str(k),pauseflag)

% calls: do_IPM.m (to run the IPM)
%        makeSimpVec (to make the Simpson vector for integration)

% Save results of each model in a separate file

end
%delete(PoolObj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing and model selection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



