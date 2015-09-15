function runme_IPM_mockdata(scenario)
% Run central coast rockfish IPM model on simulated blue rockfish datasets for validation

% Forked from runme_IPM...formats the mockdata to interact with the
% rockfish_fit_pisco code

% Load in the metadata used to create the mockdata
load('SMYS_Pt_Lobos_pre2007_13Dec2013_metadata')

data_savename = strcat('mockdata_May2015/SMYS_mockdata_',scenario,'_28May2015'); % where the mock data are
savename = strcat('SMYS_mockdata_',scenario,'_28May2015'); % where the fit will be stored


for i = 1:10 % loop over reps

% Read in data
load(strcat(data_savename,'_R',num2str(i),'.mat'))

% Now make some adjustments reflecting the mock-ness
Meta.Tdata = Years; % the years of actual data collection
Meta.MPAstatus = false;
Meta.MPAnew = false; % just for completeness...this won't actually be used
Meta.Sites = Site_Names;

Meta.Fprior = [1e-5 10; 1e-5 10]; % very very flat prior that allows small values of F
Meta.Rprior = NaN; % for now these will be NaNs (for pre-2007 runs)

Meta.ogive = [Inf 0]; % remove observation ogive for this dataset

% IPM parameters
Meta.IPM.meshsize = 100;
Meta.IPM.meshmin = 0;
Meta.IPM.meshmax = Meta.fixparm(1)*2;

% MCMC parameters
Meta.MCMC.M = 5e3; % length of chain
Meta.MCMC.chains = 4; % number of chains
Meta.Q = 100; % number of particles
    
% Savenames:
Meta.PISCO_filename = 'null';
Meta.data_savename = strcat(data_savename,'_R',num2str(i),'.mat');
Meta.fit_savename =  strcat(savename,'_R',num2str(i),'_fit.mat');
Meta.forward_savename = 'null';
Meta.savename = strcat(savename,'_R',num2str(i),'_meta.mat'); 

save(Meta.savename,'Meta')

rockfish_fit_pisco(Meta.savename,0)

end % end loop over reps