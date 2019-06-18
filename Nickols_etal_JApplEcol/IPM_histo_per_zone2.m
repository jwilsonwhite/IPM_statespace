function N = IPM_histo_per_zone2(D_str,Years,Sites,edges)

% Create histogram of size distributions at each site in D_str
% Input argument D_str is a list of each fish observed at each
% site/year/month

% March 2016: forked from IPM_histo to make calculations on a per-zone basis

% Output is array of length x site x year x side
N = nan(length(edges),length(Sites),length(Years),4); % 4 zones at most sites??
for y = 1:length(Years)
    for s = 1:length(Sites)
        
        Ntmp = D_str.(Sites{s})(y).data.Zone;
        
        if isnan(Ntmp(1).Count)
        N(:,s,y,:) = NaN; % if missing data
        % elseif nansum(Ntmp(1).Count) == 0
        %  N(:,s,y) = 0; % if data were collected, but no fish
        else % else there were some data
            
            Nzones = length(Ntmp); % number of zones sampled
            for zn = 1:Nzones
            
            Ntmp2 = [];
        for n = 1:length(Ntmp(zn).Count) % expand counts of multiple fish
        Ntmp2 = [Ntmp2; repmat(Ntmp(zn).TL(n),[Ntmp(zn).Count(n),1])]; 
        end
        
        % actually do the histogram:
        Ntmp_histo = histc(Ntmp2,[0 edges(2:end-1) Inf]);
        %Ntmp_tot=sum(Ntmp_histo); %added by KJN
        %Ntmp_histo = Ntmp_histo/Ntmp_tot; %added by KJN
        
        %N(:,s,y) = Ntmp_histo*nansum(Ntmp.Count);   %added by KJN
        N(:,s,y,zn) = Ntmp_histo;
            end % end loop over sides
        end % end 
        
    end 
end

%numbers, site, year, zone