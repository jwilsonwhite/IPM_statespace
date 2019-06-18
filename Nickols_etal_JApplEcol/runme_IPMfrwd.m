function runme_IPMfrwd(Site,Species,Type)
%Run Central Coast IPM model to look at responses to MPA implementation

%This master file reads in the data and parameters (from Ecol Apps paper)
%and runs forward in time.

%Code designed to be used for single species and site at a time.

%Site is MPA region for recruitment data as well as single site for looking
%at predictions at sites that changed from non-MPA to MPA site. 

% Species: 
% SMYS (Sebastes mystinus, blue rockfish) 

% A metadata structure file (Meta) will be associated with particular runs. 
% This keeps track of the various filenames, species, sites, etc. used in
% particular runs.  This facilitates passing information between m-files.

% Call metadata file 
savename = strcat(Species,'_',Site,'_',Type,'_31Aug2015_metadata.mat');

load(savename)
Meta.Ffit_savename = strcat(Species,'_',Site,'_',Type,'_fit_31Aug2015_postproc.mat');
%The above file comes from White et al. 2016 fits
Meta.SSDfit_savename = strcat(Species,'_',Site,'_07Sep2016_SSDfit.mat');
%Can change this name to whatever you want, it is the file that will be
%saved
Meta.frwdfit_savename = strcat(Species,'_',Site,'_07Sep2016_frwdfit.mat');
%Can change this name to whatever you want, it is the file that will be
%saved
Meta.savename_new = strcat(Species,'_',Site,'_',Type,'_07Sep2016_metadata.mat');
%Can change this name to whatever you want, it is the file that will be
%saved

switch Site
    case 'White_Rock'
        Meta.MPAnew=logical([0 1 0 0]); %need site White_Rock to be only new
    otherwise
        Meta.MPAnew=Meta.MPAnew;
end

save(Meta.savename_new, 'Meta') %save updated metadata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Read in the data
% PISCO data is in the directory PISCO_data_2015

read_PISCOdata(Meta.savename)

% reads in the data, saves it in a structure file.  Also counts up number
% of transects in each survey.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Find SSD 
findSSD(Meta.savename_new)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Forward Projections from SSD 
% Run IPM forward from SSD with fitted F and zero fishing using recruitment
% estimates from model fits (constant and random)
rockfish_frwdproject_SSD(Meta.savename_new)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Forward Projections from MPA monitoring data
% Run with fitted F and zero fishing for random recruitment 
% estimated from modeling fits and starting abundance equal to site that 
% switches from fished to reserve during year of implementation

rockfish_frwdproject(Meta.savename_new)


