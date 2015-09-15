function N = IPM_histo(D_str,Years,Sites,edges)

% Create histogram of size distributions at each site in D_str
% Input argument D_str is a list of each fish observed at each
% site/year/month

% Output is (for now) array of length x site x year
N = nan(length(edges),length(Sites),length(Years));
for y = 1:length(Years)
    for s = 1:length(Sites)
        
        Ntmp = D_str.(Sites{s})(y).data;
        
        if isnan(Ntmp.Count)
        N(:,s,y) = NaN; % if missing data
        elseif nansum(Ntmp.Count) == 0
            N(:,s,y) = 0; % if data were collected, but no fish
        else % else there were some data
            Ntmp2 = [];
        for n = 1:length(Ntmp.Count) % expand counts of multiple fish
        Ntmp2 = [Ntmp2; repmat(Ntmp.TL(n),[Ntmp.Count(n),1])]; 
        end
        
        % actually do the histogram:
        Ntmp_histo = histc(Ntmp2,[0 edges(2:end-1) Inf]);
       
        N(:,s,y) = Ntmp_histo;
        end % end 
        
    end 
end
