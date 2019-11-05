function runme_extra_SMYS_ss(fitnum, data_savename, savename, samplesize)
% Read in data
load('SMYS_Pt_Lobos_pre2007_31Aug2015_metadata');

load(strcat(data_savename,'_R',num2str(fitnum),'.mat'))
disp(strcat(data_savename,'_R',num2str(fitnum),'.mat'))
disp('it loaded mockdata!')
samplesize
% Now make some adjustments reflecting the mock-ness
Meta.Tdata = Years; % the years of actual data collection

%load file that has size data in it from PISCO
% Adjust samplesize:
for y = 1:length(Years)
Ctmp = D_str.SMYS.NULL(y).data.Count; % actual counts (using a large value of Rmean)
Cnew = Resamp(Ctmp,samplesize); % resample to get the desired total sample size
D_str.SMYS.NULL(y).data.Count = Cnew; %  re-insert this
end % end loop over years

Meta.MPAstatus = false;
Meta.MPAnew = false; % just for completeness...this won't actually be used
Meta.Sites = Site_Names;

Meta.Fprior = [1e-5 10; 1e-5 10]; % very very flat prior that allows small values of F
Meta.Rprior = NaN; % for now these will be NaNs (for pre-2007 runs)

Meta.ogive = [Inf 0]; % remove observation ogive for this dataset

% IPM parameters
Meta.IPM.meshsize  = 100;
Meta.IPM.meshmin = 0;
Meta.IPM.meshmax = Meta.fixparm(1)*2;

% MCMC parameters
Meta.MCMC.M = 5e3; % length of chain
Meta.MCMC.chains = 4; % number of chains
Meta.Q = 100; % number of particles

% Savenames:
Meta.PISCO_filename = 'null';
Meta.data_savename = strcat(data_savename,'_R',num2str(fitnum),'.mat');
Meta.fit_savename =  strcat(savename,'_R',num2str(fitnum),'_fit.mat');
Meta.forward_savename = 'null';
Meta.savename = strcat(savename,'_R',num2str(fitnum),'_meta.mat');

save(Meta.savename,'Meta')
save(Meta.data_savename,'D_str','RF','Site_Names','Species_Names','Years')
disp(Meta)
rockfish_fit_pisco_mock_ss(Meta,0,D_str,RF,Site_Names,Species_Names,Years)
end
