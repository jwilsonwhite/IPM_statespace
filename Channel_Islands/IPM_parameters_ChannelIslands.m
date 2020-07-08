function Meta = IPM_parameters_ChannelIslands(Species,Site,Type,F_initial)

%**************************************************************************
% Create all parameters used in the IPM statespace model.
% Forked from 'IPM_parameters.m' in July 2019 for use in the M. Yamane
% Channel Islands project.
% 
% 'Species' is a string giving 4-letter species code
% 'Type' is pre2003, post2003 
% 'Site' is a string giving the PISCO sampling site to be modeled
%**************************************************************************

Meta.Species = {Species};
Meta.Sites = {Site};

switch Type
    
    case 'pre2003'
    % fit sites in a given region with common parameters (except F)
        
    Meta.Tdatafull = 1999:2017; % full timespan of data
    Meta.Tpre =  1990:1998; % time prior to data that model will run
    Meta.Tdata =  1999:2003; % years with data before MLPA MPAs    
    Meta.Rprior = NaN; % no prior on R, estimate it from the data    
 

    case 'post2003'
    % fit to individual sites post-2003
    
    Meta.Tdatafull = 1999:2017; % full timespan of data
    Meta.Tpre =   1990:1998; % time prior to data that model will run
    % (NaN b/c start with initial distribution in 2003)
    
    Meta.Tdata =  1999:2017; % years with data to fit
    % (include 2003 as year 1, data not fit in that year though) 
    
    Meta.Rprior = NaN; % no prior on R, estimate it from the data  

    Meta.MPAstatus = false(size(Meta.Sites)); % assumes no MPA, fits an F
    Meta.MPAnew = false(size(Meta.Sites)); % this will not get used
    
    Meta.Tpre2 = 1999:2002; % years with data but before MPA
    Meta.Tpost = 2003:2017; % years with data  + MPA
      
end % end switch Type

% Parameters that are fixed for the estimation step:
%              [Linf   k   x0   M   Lfish   Lmat   Lvar];

% Note that some von Bert params are given in terms of t0 not x0. To
% convert, x0 = Linf*(1-exp(K*t0))
% Lvar is value of CV from stock assessment.

switch Species
    
    case 'SMYS'
        % Blue rockfish (SMYS):
        Meta.fixparm = [36 0.105 3.308 0.13 21 29 0.087];
        
        % Prior of Fs. Based on stock assessment values.  
        % Format: 1st row is pre-1999; 2nd row is post-1999
        % 1st column is mean, 2nd column is SD (use 1 as a good value on
        % lognormal scale?)
        switch Type
            case 'pre2003'
                Meta.Fprior = [0.25 0.23; 0.03 1]; % mean & SD of prior
            case 'post2003'
                Meta.Fprior = [0.25 0.23; 0.03 1];
        end
        
    case 'SATR'
        % Kelp rockfish (SATR): (t0 = -0.7, so x0 = 5.62)
        Meta.fixparm = [37.8 0.23 5.62 0.2 22.7 28 0.23];
        
        switch Type
            case 'pre2003'
                Meta.Fprior = [0.24 0.44; 0.07 1]; % mean & SD of prior
                % (F1: derived from S. carnatus values)
            case 'post2003'
                Meta.Fprior = [0.24 0.44; 0.07 1]; % mean & SD of prior
        end
        
    case 'PCLA'
        % Kelp Bass (PCLA):
        Meta.fixparm = [69.8 0.06 13.22 0.178 35.56 22.3 0.24];
        
        switch Type
            case 'pre2003'
                Meta.Fprior = [0.371 1; 0.101 1]; % mean & SD of prior
            case 'post2003'
                Meta.Fprior = [0.371 1; 0.101 1]; % mean & SD of prior 
        end
        
    case 'SPUL'
        % California Sheephead (SPUL):
        Meta.fixparm = [66.27 0.158 11.963 0.2 30.48 25.24 0.11]; % M adjusted from 0.1 to 0.2 on 1 April 2020
        
        switch Type
            case 'pre2003'
                Meta.Fprior = [0.19 .113; 0.07 0.23]; % mean & SD of prior
            case 'post2003'
                Meta.Fprior = [0.19 .113; 0.07 0.23]; % mean & SD of prior
        end

    case 'EJAC'
        % Black Surfperch (EJAC):
        Meta.fixparm = [25.745 0.742 7.458 0.183 21 23.3 0.127];
        
        switch Type
            case 'pre2003'
                Meta.Fprior = [0.01 1; 0.01 1]; % mean & SD of prior
            case 'post2003'
                Meta.Fprior = [0.01 1; 0.01 1]; % mean & SD of prior
        end
        
    case 'ELAT'
        % Striped Surfperch (ELAT):
        Meta.fixparm = [39.7 0.24 4.83 0.42 21 22.8 0.5]; % from fishbase.se
        
        switch Type
            case 'pre2003'
                Meta.Fprior = [0.01 .5; 0.01 0.5]; % mean & SD of prior
            case 'post2003'
                Meta.Fprior = [0.01 .5; 0.01 0.5]; % mean & SD of prior
        end
end % end switch species

% Recruit statistics - estimated empirically from the data
switch Species
    case 'SMYS'
        Meta.recruits.meansize = 7.75;% mean of 'recruits size 
        Meta.recruits.sdsize = 1.15; % sd of 'recruits' size
        Meta.recruits.Rsize = 10; % max size of recruits in data (YOY)
        Meta.ogive = [28 7]; % normal cdf parameters specifying
        
        % ogive giving probability of actually being 
        % observed in the kelp forest.  For blue
        % rockfish this is from onto_ogive.m based
        % on Rick Starr's data.                                
        % *** For species without an ogive, set these parameters to be: 
        % Meta.ogive = [NaN, NaN];
        
    case 'SATR'
        Meta.recruits.meansize = 4.23;% mean of 'recruits size
        Meta.recruits.sdsize = 1.04; % sd of 'recruits' size
        Meta.recruits.Rsize = 8; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
        
    case 'PCLA'
        Meta.recruits.meansize = 5.34;% mean of 'recruits size
        Meta.recruits.sdsize = 1.83; % sd of 'recruits' size
        Meta.recruits.Rsize = 9; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
        
    case 'SPUL'
        Meta.recruits.meansize = 14.82;% mean of 'recruits size
        Meta.recruits.sdsize = 3.09; % sd of 'recruits' size
        Meta.recruits.Rsize = 21; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters

    case 'EJAC'
        Meta.recruits.meansize = 8.85;% mean of 'recruits size
        Meta.recruits.sdsize = 1.08; % sd of 'recruits' size
        Meta.recruits.Rsize = 11; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
        
    case 'ELAT'
        Meta.recruits.meansize = 8.75;% mean of 'recruits size
        Meta.recruits.sdsize = 1.13; % sd of 'recruits' size
        Meta.recruits.Rsize = 11; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
end % end switch species

% IPM parameters
Meta.IPM.meshsize = 100;
Meta.IPM.meshmin = 0;
Meta.IPM.meshmax = Meta.fixparm(1)*2;

% MCMC parameters (start small to test, then expand)
Meta.MCMC.M = 1e4; % length of chain
Meta.MCMC.chains = 3; % number of chains
Meta.Q = 200; % number of particles

Meta.Fs = {'0','fit'}; % cell array.
% What F values to use in forward projections? (NOT USED)


Meta.F_initial=F_initial; % Fix the value of F to be used in the 'pre-data' period if it is not possible to estimate from data
% This is used in the 2003-onwards Channel Islands analysis to specify
% pre-2003 F for sites where data collection started well after 2003.


          
% give a more detailed PISCO_filename...
Meta.PISCO_filename = strcat('data/PISCO_fish_imported_',...
    Species,'_',Site,'.mat');
Meta.data_savename = strcat('data/',Species,'_',Site,'_',Type,...
    '_data_sorted.mat');
Meta.fit_savename = strcat(Species,'_',Site,'_',Type,'_fit_March2020.mat');
Meta.forward_savename = strcat(Species,'_',Site,'_forward_March2020.mat');
Meta.savename = strcat(Species,'_',Site,'_',Type,'_March2020_metadata.mat');

save(Meta.savename,'Meta') % save the metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


