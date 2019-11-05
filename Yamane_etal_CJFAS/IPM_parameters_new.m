%function [] = IPM_parameters_new(varargin)

%if(nargin>0)
%scenario = varargin{1};
%else
%scenario = 'PCLA','Pt_Lobos','pre2007';
%end

% Create all parameters used in the IPM statespace model

% 'Species' is a string giving 4-letter species code
% 'Type' is pre2007, post2007MPA (MPA sites) or post2007fished (fished
% sites)
% 'Site' is a string giving the MPA region to be modeled

Species = 'SATR';
Site = 'Pt_Lobos';
Type = 'pre2007';

%Meta = str[];
Meta.Species = {Species};

switch Type
    
    case 'pre2007' % strategy: fit sites in a given region with common parameters (except F)
        
    Meta.Tdatafull = 1999:2014; % full timespan of data
    Meta.Tpre =  1990:1998; % time prior to data that model will run
    Meta.Tdata =  1999:2007; % years with data before MLPA MPAs    
    Meta.Rprior = NaN; % no prior on R, estimate it from the data    
    
switch Site
    case 'Pt_Lobos'

% Sites: Pt Lobos:
Meta.Sites = {'BLUEFISH','MALPASO','MONASTERY',...
              'PALO_COLORADO','SOBERANES','WESTON'}; % cell array.  Insert NaN if you want all the sites
Meta.MPAstatus = logical([1 0 0 0 0 1 0]); % true if it has always been an MPA
Meta.MPAnew = logical([0 0 1 0 0 0 0]); % true if it became an MPA in 2007
   
    case 'Big_Creek'

% Sites: Big Creek:
Meta.Sites = {'BIG_CREEK','DOLAN','ESALEN',...
              'LOPEZ'}; % cell array.  Insert NaN if you want all the sites
Meta.MPAstatus = logical([1 0 0 0]); % true if it has always been an MPA
Meta.MPAnew = logical([0 1 0 0]); % true if it became an MPA in 2007

    case 'White_Rock'
% Sites: White Rock:
switch Species
    case 'SATR' % only observed at a few of the sites
Meta.Sites = {'WHITE_ROCK','CAYUCOS'}; % cell array.  Insert NaN if you want all the sites
Meta.MPAstatus = logical([0 0]); % true if it has always been an MPA
Meta.MPAnew = logical([1 0]); % true if it became an MPA in 2007        
    otherwise
Meta.Sites = {'ESTERO'...
              'WHITE_ROCK','CAYUCOS','HARMONY'}; % cell array.  Insert NaN if you want all the sites
Meta.MPAstatus = logical([0 0 0 0 0 0]); % true if it has always been an MPA
Meta.MPAnew = logical([1 1 1 1 0 0]); % true if it became an MPA in 2007
end % end switch Species

end % end switch Site

    case {'post2007MPA', 'post2007fished'}  % fit to individual sites post-2007
    
    Meta.Tdatafull = 1999:2014; % full timespan of data
    Meta.Tpre =   NaN; % time prior to data that model will run (NaN b/c start with initial distribution
    Meta.Tdata =  2007:2014; % years with data to fit   (include 2007 as year 1, data not fit in that year though) 
        
    switch Type
        case 'post2007MPA';
            switch Site
                case 'Pt_Lobos'
                Meta.Sites = {'BLUEFISH','MONASTERY','WESTON' };
                case 'Big_Creek'
                Meta.Sites = {'BIG_CREEK'};
                case 'White_Rock'
                    switch Species
                        case 'SATR' % SATR not observd at Estero
                            Meta.Sites = {'WHITE_ROCK'};
                        otherwise
                            Meta.Sites = {'WHITE_ROCK','ESTERO'}; 
                    end
            end % end switch site
        case 'post2007fished'
            switch Site
                case 'Pt_Lobos'
                Meta.Sites = {'MALPASO','PALO_COLORADO','SOBERANES'};
                case 'Big_Creek'
                Meta.Sites = {'ESALEN','LOPEZ'};
                case 'White_Rock'
                    switch Species
                        case 'SATR' % SATR not observed at Harmony
                            Meta.Sites = {'CAYUCOS'};
                        otherwise
                        Meta.Sites = {'CAYUCOS','HARMONY'};
                    end
            end % end switch site
    end % end switch Type
    
    Meta.MPAstatus = false(size(Meta.Sites)); % assume no MPA & just fit an F
    Meta.MPAnew = false(size(Meta.Sites)); % this will not actually get used
  
    % Get original fit from 2007, based on an existing pre-2007 fit
    switch Species
        case 'SATR'
    fname = strcat(Species,'_',Site,'_pre2007_fit_8July2015_postproc.mat');
        case 'SCAR'
     fname = strcat(Species,'_',Site,'_pre2007_fit_11June2015_postproc.mat'); 
        case 'SMYS'
            switch Site
                case 'Pt_Lobos'
                  fname = strcat(Species,'_',Site,'_pre2007_fit_1July2015_postproc.mat');  
                otherwise
                  fname = strcat(Species,'_',Site,'_pre2007_fit_8July2015_postproc.mat');
            end
    end
    load(fname);
    
    Meta.N_init = Post(1).predicted; % note this is a structure with fields for each site, then a size x year matrix 
    
    Rs = Post(1).posterior(:,4:12);
    Meta.Rprior = [mean(Rs(:)), std(mean(Rs))];
    Meta.Rfact= Post.Rfact;
    
    if ~exist('Fprior','var') % if Fprior is not specified, use the one from the fit
        switch Type
            case 'post2007fished'
    Fprior = [Post.F_mean, Post.F_std];
            case 'post2007MPA'
    Fprior = [0 1];
        end
    end

end % end switch Type



% Parameters that are fixed for the estimation step:
%              [Linf      k   x0  M  Lfish Lmat Lvar];

% Note that some von Bert params are given in terms of t0 not x0.  To
% convert, x0 = Linf*(1-exp(K*t0))
% Lvar is value of CV from stock assessment.

% For notes on parameter sources, see
% Life_History_Paramter_for_RF_models.xlsx (Created by K. Nickols)

switch Species
    
    case 'SMYS'
% For blue rockfish tzero: -0.95 for males, -1.34 for females (mean =
% -1.145

% Blue rockfish (SMYS):
Meta.fixparm = [38.150142 0.172 6.2533 0.14 21.0295 27.086 0.1];

% Prior of Fs. Based on stock assessment values (see rockfish_priors.m).  

% Format: 1st row is pre-1999; 2nd row is post-1999
%         1st column is mean, 2nd column is SD (use 1 as a good value on
%         lognormal scale?)
switch Type
    case {'pre2007', 'pre2012', 'pre2003'}
        Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior
    case {'post2007MPA','post2007fished'}
        Meta.Fprior = [Fprior; Fprior];
end

    case 'SMEL'
% Black rockfish (SMEL): (t0 = 0.75)
Meta.fixparm = [45.1109 0.17 -12.668 0.14 29.0404 40.2293 0.07];  %k previously 0.33

% Prior from stock assessment
switch Type
    case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior
    case {'post2007MPA','post2007fished'}
        Meta.Fprior = [Fprior; Fprior];
end

    case 'SATR'
% Kelp rockfish (SATR): (t0 = -0.7, so x0 = 5.6369)
Meta.fixparm = [33.15 0.2307 0.1668 0.2 20 18 0.1]; 
switch Type
    case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior (use same values as S. carnatus)
    case {'post2007MPA','post2007fished'}
        Meta.Fprior = [Fprior; Fprior];
end

    case 'SCAR'
% Gopher rockfish (SCAR): (t0 = -0.5, so x0 = 3.6375)
Meta.fixparm = [34.1 0.2256 3.6375 0.2 25.5 17 0.1];
% F prior based on stock assessment - very high F in mid-to-late 90s, much
% lower in 2000s
switch Type
    case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior
    case {'post2007MPA','post2007fished'}
        Meta.Fprior = [Fprior; Fprior];
end
                                                               
case 'SCAU'
% Copper rockfish (SCAU): (t0 = -1, so x0 = 7.15)
Meta.fixparm = [56.5 0.1354 7.15 0.09 29.72 32 0.1];
% F prior based on stock assessment - very high F in mid-to-late 90s, much
% lower in 2000s
switch Type
    case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior
    case {'post2007MPA','post2007fished'}
        Meta.Fprior = [Fprior; Fprior];
end
                                                               
case 'SFLA'
% Yellowtail rockfish (SFLA): (t0 = -1, so x0 = 7.15)
Meta.fixparm = [49.9 0.1783 9.67 0.11 28.5 37 0.1];
% F prior based on stock assessment - very high F in mid-to-late 90s, much
% lower in 2000s
switch Type
    case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior
    case {'post2007MPA','post2007fished'}
        Meta.Fprior = [Fprior; Fprior];
end

case 'SMIN'
% Vermilion rockfish (SMIN): (t0 = -0.178, so x0 = 1.56)
Meta.fixparm = [53.9 0.1643 1.56 0.1 22 36 0.07]; %Lmat was 36
% F prior based on stock assessment - very high F in mid-to-late 90s, much
% lower in 2000s
switch Type
    case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior
    case {'post2007MPA','post2007fished'}
        Meta.Fprior = [Fprior; Fprior];
end

case 'SAUR'
% Brown rockfish (SAUR): (t0 = -0.55, so x0 = 4.33)
Meta.fixparm = [51.4 0.16 4.33 0.14 27.5 27.5 0.1]; % Lmat previously 27.5
% F prior based on stock assessment - very high F in mid-to-late 90s, much
% lower in 2000s
switch Type
    case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior
    case {'post2007MPA','post2007fished'}
        Meta.Fprior = [Fprior; Fprior];
end
                                                               
                                                               
    case 'STRFRAAD' 
% red urchin(STRFRAAD): (t0 = 0 so x0 = 0)
Meta.fixparm = [11.8 0.22 0 0.08 8.9 6 0.3];
switch Type
case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior (use same values as S. carnatus)
case {'post2007MPA','post2007fished'}
Meta.Fprior = [Fprior; Fprior];
end

%[Linf  k   x0  M  Lfish Lmat Lvar]
case 'OELO'
% Lingcod(OELO): (t0 = -1.56, so x0 = 22.9)
Meta.fixparm = [96.74 0.1734327 22.9 0.25 60.1 49.3 0.09];
switch Type
case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior (use same values as S. carnatus)
case {'post2007MPA','post2007fished'}
Meta.Fprior = [Fprior; Fprior];
end
                                                               
case 'PCLA'
% Kelp Bass(PCLA): (t0 = -3.5, so x0 = 13.2)
Meta.fixparm = [69.8 0.06 13.2 0.178 30.48 22.3 0.1]; % k previously 0.06
switch Type
case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior
case {'post2007MPA','post2007fished'}
Meta.Fprior = [Fprior; Fprior];
end
                                                               
case 'SCHR'
% Black-and-Yellow (SCHR): (t0 = -0.3837722, so x0 = 2.10565492)
Meta.fixparm = [24.9 0.23 2.1 0.20 24 17.5 0.1]; %
switch Type
case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior (use same values as S. carnatus)
case {'post2007MPA','post2007fished'}
Meta.Fprior = [Fprior; Fprior];
end
                                                               
    case 'SNEB'
% China (SNEB): (t0 = -0.26, so x0 = 1.95)
Meta.fixparm = [33.62 0.23 1.95 0.056 30.48 27 0.1]; %
switch Type
case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior (use same values as S. carnatus)
case {'post2007MPA','post2007fished'}
Meta.Fprior = [Fprior; Fprior];
end
                                                               
case 'HDEC'
% Kelp Greenling (HDEC):
Meta.fixparm = [41.14805 0.2447207 15.36 0.301 30.48 30 0.1]; %
switch Type
case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior (use same values as S. carnatus)
case {'post2007MPA','post2007fished'}
Meta.Fprior = [Fprior; Fprior];
end
                                                               
case 'SPAU'
% Bocaccio (SPAU):
Meta.fixparm = [70 0.22 -0.69 0.15 38 35.5 0.1]; %
switch Type
case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % mean & SD of prior (only F's from 00s available)
case {'post2007MPA','post2007fished'}
Meta.Fprior = [Fprior; Fprior];
end

case 'SMAR'
% Cabezon (SMAR):
Meta.fixparm = [49.9 0.28 14.53 0.275 38.1 34 0.1]; %
switch Type
case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % same as bocaccio
case {'post2007MPA','post2007fished'}
Meta.Fprior = [Fprior; Fprior];
end

case 'SSER'
% Olive (SSER):
Meta.fixparm = [47.7 0.21 11.8 0.14 28.5 34 0.1]; %
switch Type
case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % same as yellowtail
case {'post2007MPA','post2007fished'}
Meta.Fprior = [Fprior; Fprior];
end

case 'SGUT'
% CA Scorpionfish (SGUT):
Meta.fixparm = [40.29 0.126 11.6 0.25 25.4 17.5 0.1]; %
switch Type
case {'pre2007', 'pre2012', 'pre2003'}
Meta.Fprior = [1E-5 10; 1E-5 10]; % same as cabezon
case {'post2007MPA','post2007fished'}
Meta.Fprior = [Fprior; Fprior];
end
                                                               
end

% Recruit statistics - estimated empirically from the data
switch Species
    case 'SMYS'
Meta.recruits.meansize = 7.75;% mean of 'recruits size (estimated empirically from the data)
Meta.recruits.sdsize = 1.15; % sd of 'recruits' size
Meta.recruits.Rsize = 10; % max size of recruits in data (YOY)
Meta.ogive = [28 7]; % normal cdf parameters specifying
                                % ogive giving probability of actually being 
                                % observed in the kelp forest.  For blue
                                % rockfish this is from onto_ogive.m based
                                % on Rick Starr's data.                                
% *** For species without an ogive, set these parameters to be: 
%Meta.ogive = [NaN, NaN]; 
    
    case 'SMEL'
        
    Meta.recruits.meansize = 7.43;% mean of 'recruits size (estimated empirically from the data)
    Meta.recruits.sdsize = 1.32; % sd of 'recruits' size
    Meta.recruits.Rsize = 10; % max size of recruits in data (YOY)
    Meta.ogive = [NaN,NaN]; % normcdf parameters
    
    case 'SATR'
        
    Meta.recruits.meansize = 4.23;% mean of 'recruits size (estimated empirically from the data)
    Meta.recruits.sdsize = 1.04; % sd of 'recruits' size
    Meta.recruits.Rsize = 8; % max size of recruits in data (YOY)
    Meta.ogive = [NaN,NaN]; % normcdf parameters
        
    case 'SCAR'
        
        Meta.recruits.meansize = 5.08;% mean of 'recruits size (estimated empirically from the data)
        Meta.recruits.sdsize = 0.70; % sd of 'recruits' size
        Meta.recruits.Rsize = 8; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
    
    case 'SCAU'
        
        Meta.recruits.meansize = 8.13;% mean of 'recruits size (estimated empirically from the data)
        Meta.recruits.sdsize = 1.26; % sd of 'recruits' size
        Meta.recruits.Rsize = 13; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
                                                               
   case 'SFLA'
        
        Meta.recruits.meansize = 4.55;% mean of 'recruits size (not estimated empirically)
        Meta.recruits.sdsize = 1.26; % sd of 'recruits' size (this is a guess!)
        Meta.recruits.Rsize = 10; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
                                                               
    case 'SMIN'
        
        Meta.recruits.meansize = 3;% mean of 'recruits size (not estimated empirically)
        Meta.recruits.sdsize = 1.21; % sd of 'recruits' size (this is a guess!)
        Meta.recruits.Rsize = 9; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
    
    case 'SAUR'
        
        Meta.recruits.meansize = 3;% mean of 'recruits size (not estimated empirically)
        Meta.recruits.sdsize = 1.26; % sd of 'recruits' size (this is a guess!)
        Meta.recruits.Rsize = 10; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters

    case 'OELO'
        
        Meta.recruits.meansize = 14.48;% mean of 'recruits size (averaged estimates from stock assessment)
        Meta.recruits.sdsize = 4.16; % sd of 'recruits' size, average from 5 studies' estimates
        Meta.recruits.Rsize = 20; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
                                                               
    case 'PCLA'
        
        Meta.recruits.meansize = 2; %7.85;% mean of 'recruits size (averaged estimates from stock assessment)
        Meta.recruits.sdsize = 0.2; %6.4; % sd of 'recruits' size, average from 5 studies' estimates
        Meta.recruits.Rsize = 3; %25; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
                                                                    
   case 'SCHR'
        Meta.recruits.meansize = 5.08;% mean of 'recruits size (used SCAR recruit sizes)
        Meta.recruits.sdsize = 0.7; % sd of 'recruits' size
        Meta.recruits.Rsize = 8; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
                                                       
   case 'HDEC'
        Meta.recruits.meansize = 3;% mean of 'recruits size (used SCAR recruit sizes)
        Meta.recruits.sdsize = 0.7; % sd of 'recruits' size
        Meta.recruits.Rsize = 19; % max size of recruits in data (YOY)%previously 7
        Meta.ogive = [NaN,NaN]; % normcdf parameters
                    
   case 'SPAU'
        Meta.recruits.meansize = 2.5;% mean of 'recruits size (educated guesses)
        Meta.recruits.sdsize = 5; % sd of 'recruits' size (except this was a total guess)
        Meta.recruits.Rsize = 20; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
        
   case 'SMAR'
        Meta.recruits.meansize = 14; % tried 4 before,mean of 'recruits size (educated guesses)
        Meta.recruits.sdsize = 4.5; % 1 before, sd of 'recruits' size (except this was a total guess)
        Meta.recruits.Rsize = 27; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
        
   case 'SNEB'
        Meta.recruits.meansize = 4;% mean of 'recruits size (guesses)
        Meta.recruits.sdsize = 1; % sd of 'recruits' size 
        Meta.recruits.Rsize = 6; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
  
    case 'SSER'
        Meta.recruits.meansize = 3.4;% mean of 'recruits size 
        Meta.recruits.sdsize = 1; % sd of 'recruits' size (guess)
        Meta.recruits.Rsize = 17; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
        
    case 'SGUT'
        Meta.recruits.meansize = 1.1;% mean of 'recruits size 
        Meta.recruits.sdsize = 1; % sd of 'recruits' size (guess)
        Meta.recruits.Rsize = 12; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
                                                        
                                                        
case 'STRFRAAD'
        Meta.recruits.meansize = 1;% mean of 'recruits size
        Meta.recruits.sdsize = 0.5; % sd of 'recruits' size (guess)
        Meta.recruits.Rsize = 2; % max size of recruits in data (YOY)
        Meta.ogive = [NaN,NaN]; % normcdf parameters
                                                        
    % previously:
   %Meta.recruits.meansize = 0.5;% mean of 'recruits size
   %Meta.recruits.sdsize = 0.2; % sd of 'recruits' size (guess)
   %Meta.recruits.Rsize = 1; % max size of recruits in data (YOY)
   %Meta.ogive = [NaN,NaN]; % normcdf parameters
                                                        
end % end switch species

% IPM parameters
Meta.IPM.meshsize = 100;
Meta.IPM.meshmin = 0;
Meta.IPM.meshmax = Meta.fixparm(1)*2;

% MCMC parameters
Meta.MCMC.M = 1e4; % length of chain
Meta.MCMC.chains = 4; % number of chains
Meta.Q = 200; % number of particles

Meta.Fs = {'0','fit'}; % cell array.  What F values to use in forward projections?
          
% give a more detailed PISCO_filename...
Meta.PISCO_filename = strcat('data/PISCO_fish_imported_',Species,'_',Site,'.mat');
Meta.data_savename = strcat('data/',Species,'_',Site,'_',Type,'_data_sorted.mat'); % string
Meta.fit_savename = strcat(Species,'_',Site,'_',Type,'_fit_31Aug2015.mat'); % string
Meta.forward_savename = strcat(Species,'_',Site,'_forward_31Aug2015.mat'); % string
Meta.savename = strcat(Species,'_',Site,'_',Type,'_31Aug2015_metadata.mat');

save(Meta.savename,'Meta') % save the metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


