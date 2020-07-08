function plot_posterior_distributions

% Plots the posterior distributions for all of the sites and years.

% List of sites, in desired order:
Sites = {'SMI_TYLER_BIGHT','SMI_CROOK_POINT','SMI_HARRIS_PT_RESERVE','SMI_CUYLER',...
         'SRI_SOUTH_POINT','SRI_JOHNSONS_LEE_SOUTH','SCI_FORNEY','SCI_PAINTED_CAVE',...
         'SCI_HAZARDS','SCI_GULL_ISLE','SCI_VALLEY','SCI_YELLOWBANKS','SCI_PELICAN',...
         'SCI_CAVERN_POINT','SCI_SCORPION','ANACAPA_WEST_ISLE','ANACAPA_MIDDLE_ISLE',...
         'ANACAPA_EAST_ISLE','ANACAPA_LIGHTHOUSE_REEF'};
       
Species = {'SATR','SMYS','SPUL','PCLA'};

% Which species to plot for each site?
doSpecies = [1 1 1 0 1 1 1 0 1 1 0 0 1 0 0 0 0 0 0; % SATR
             1 1 1 1 1 1 0 1 1 1 0 0 1 0 0 0 0 0 0; % SMYS
             0 0 0 0 1 1 0 1 0 1 1 0 1 1 0 1 0 1 1; % SPUL
             0 0 0 0 0 0 0 1 1 0 0 1 1 1 1 1 0 0 1];% PCLA
         
for s = 1:length(Sites)
    for p = 1:length(Species)
        
        if doSpecies(p,s)
        % Assemble meta-savename
        Meta_savename = strcat(Species{p},'_',Sites{s},'_post2003_July2019_metadata.mat');
        load(Meta_savename,'Meta')
        load(Meta.data_savename,'D_str','Species_Names','Years','Site_Names')
        postproc_pisco_rockfish_fit(Meta_savename)
        end
    end
end

         