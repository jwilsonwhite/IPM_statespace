function [N,Sample_size,Site_order] = sort_oysterdata_IPM(Meta)
% This reads in Kimbro's Actual_Size_Data.xls and creates appropriate
% histograms for IPM analysis, accounting for differences in sampling
% effort across sites/times
% Uses input argument Meta to define edges of histogram, based on IPM grid
% size

Input = importdata('Data/Actual_Size_Data.xls'); % read in data

% Put the different columns into different variables.

% Need to trim out "NoFace" observations.  In early years, Kimbro did not
% record data for rockfaces that were missing on a given rock.  In 2009,
% "NoFace" was entered for those cases.  The code is set up to such that if
% there is a row entry for a rockface, that face was sampled.  So we need
% to trim out those rows.

NoFace = strcmp(Input.textdata(2:end,7),'NoFace');
NoFace_cell = [true; NoFace]; % add an extra row at the top for the cell arrays in Input

% Note that numeric data are in vectors, text data go into cell arrays
Year = Input.data(~NoFace,1); 
Site = strtrim(Input.textdata(~NoFace_cell,2));% strtrim will remove leading and trailing whitespaces
Transect = Input.data(~NoFace,3);
Rock = strtrim(Input.textdata(~NoFace_cell,4));
RockFace = Input.textdata(~NoFace_cell,6);
Length = Input.data(~NoFace,7);

% Get the number of unique types of each level of sampling
Uni_Year = unique(Year);
Uni_Site = unique(Site);
Uni_Trans = unique(Transect);
Uni_Rock = unique(Rock);

% Now loop over years, sites, transects, and rocks
for y = 1:length(Uni_Year)
    Y = Uni_Year(y);
    
    OKyear = Year == Y;
    
    for s = 1:length(Uni_Site)
        S = Uni_Site{s};
        
        OKsite = strcmp(Site,S) & OKyear;
        for t = 1:length(Uni_Trans)
            T = Uni_Trans(t);
            
            OKtrans = (Transect == T) & OKsite;
            
            for r = 1:length(Uni_Rock)
                R = Uni_Rock{r};
                
                % These are the indices of the dataset that are in this
                % combination of year, site, transect, and rock
                OKrock = strcmp(Rock,R) & OKtrans;
                % Now figure out how many faces there were on that rock
                NumFaces(y,s,t,r) = length(unique(RockFace(OKrock)));
                
                % How many individuals were counted?
                TotalOyster_rock(y,s,t,r) = sum(~isnan(Length(OKrock))); % number of actual oysters counted on this rock.
                
                %keyboard
            end % end loop over rock
        end % end loop over transect
        
        % Now get histogram for this site:
        edges = Meta.edges;
        Ntmp = histc(Length(OKsite),[0 Meta.edges(2:end-1) Inf]); % histogram
        
        % Correct for sample effort:
        NumFace_tmp = NumFaces(y,s,:,:);
        Sample_size(s,y) = sum(NumFace_tmp(:)); % total number of rock faces sampled at each site-year
        
        N(:,s,y) = Ntmp;
        
        
    end % end loop over site
end % end loop over year

        % Total sample size per site/year:
        
% Site order: currently alphabetical
Site_order = Uni_Site;
     


        
        
        
    


