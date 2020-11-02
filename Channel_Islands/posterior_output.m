function posterior_output


% List of sites, in desired order:
Sites = {'SMI_TYLER_BIGHT','SMI_CROOK_POINT','SMI_HARRIS_PT_RESERVE','SMI_CUYLER',...
         'SRI_SOUTH_POINT','SRI_JOHNSONS_LEE_SOUTH','SCI_FORNEY','SCI_PAINTED_CAVE',...
         'SCI_HAZARDS','SCI_GULL_ISLE','SCI_VALLEY','SCI_YELLOWBANKS','SCI_PELICAN',...
         'SCI_CAVERN_POINT','SCI_SCORPION','ANACAPA_WEST_ISLE','ANACAPA_MIDDLE_ISLE',...
         'ANACAPA_EAST_ISLE','ANACAPA_LIGHTHOUSE_REEF'};
isMPA = [0 0 1 0  1 0 0 1 0 1 0 0 0 1 1 1 1 1 0];
       
Species = {'SATR','SMYS','SPUL','PCLA'};

Priors = [0.07, 0.03, 0.08, 0.1]; % priors on F
Priors_sd = [1.0, 1.0, 0.23, 1.0]; % sds on priors

% Which species to plot for each site?
doSpecies = [1 1 1 0 1 1 1 0 1 1 0 0 1 0 0 0 0 0 0; % SATR
             1 1 1 1 1 1 0 1 1 1 0 0 1 0 0 0 0 0 0; % SMYS
             0 0 0 0 1 1 0 1 0 1 1 0 1 1 0 1 0 1 1; % SPUL
             0 0 0 0 0 0 0 1 1 0 0 1 1 1 1 1 0 0 1];% PCLA
         
         
    for s = length(Sites)
    
    for p = 1:length(Species)
        
        if doSpecies(p,s)
        fname = strcat(Species{p},'_',Sites{s},'_post2003_July2019_metadata.mat');
        
        postproc_pisco_rockfish_fit(fname);
        end
    end
    end
    
 
         
         