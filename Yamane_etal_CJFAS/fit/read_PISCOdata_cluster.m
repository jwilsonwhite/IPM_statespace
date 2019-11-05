function read_PISCOdata_cluster(Meta)

% Ths function reads in PISCO data for state-space IPM

% Read in PISCO dataset, convert to array with size structure by year &
% site
if ~exist(Meta.PISCO_filename,'file') % if the file doesn't exist
Dir = 'PISCO_data_2015/';
File = 'PISCO_FISH.csv'; % this version lacks headers and only has data relevant to the PLoS One paper 

fid = fopen(strcat([Dir,File]));
D = textscan(fid, strcat(['%s %s ',...% campus,  method
                          '%f %f %f ',... % year month day
                           '%s %s %s %s %f',... % site side zone level transect
                           '%s %f %f %f %f',... % classcode, count fish_tl min_tl max_tl
                           '%s %f %f %f %s %s %f',... % observer, depth, vis, temp, surge, windwave, pctcnpy
                           '%f %f %f %f']),... % macpyr1_4 5_15 >15 nerlue stipes, notes
    'Delimiter',',','MultipleDelimsAsOne',false,'HeaderLines',0);
fclose(fid);

% Column names
ColNames = {'campus','method','year','month','day','site','side','zone','level','transect',...
            'classcode','count','fish_tl','min_tl','max_tl','observer','depth','vis','temp','surge',...
            'windwave','pctcnpy','macpyr1_4','macpyr5_15','>15','nerlue','stipes','notes'};

save(Meta.PISCO_filename,'D','ColNames')
else
load(Meta.PISCO_filename,'D','ColNames')
end        
        
Species_Names = Meta.Species; % from the metadata

if isnan(Meta.Sites{1}) % just take all of the sites
    Site_Names = unique(D{strcmp(ColNames,'site')}); 
else
    Site_Names = Meta.Sites;
end

Years = unique(D{strcmp(ColNames,'year')});
Yearcol = D{strcmp(ColNames,'year')};
Sitecol = D{strcmp(ColNames,'site')};
Classcol = D{11};
Countcol = D{12};
Fish_TLcol = D{13};
Monthcol = D{strcmp(ColNames,'month')};
Sidecol = D{strcmp(ColNames,'side')};
Zonecol = D{strcmp(ColNames,'zone')};
Transcol = D{strcmp(ColNames,'transect')};

D_str = struct([]);
for s = 1:length(Species_Names)
    for g = 1:length(Site_Names)
    for y = 1:length(Years)
        
    % Find unique months in each year
    OKgy = strcmp(Sitecol,Site_Names{g}) & Yearcol == Years(y); % this site & year was sampled
    OKsp = strcmp(Classcol,Species_Names{s}); % rows with this species
    Months = unique(Monthcol(OKgy & OKsp));
    
    Sides = unique(Sidecol(OKgy));
    Zones = unique(Zonecol(OKgy));
    Transects = unique(Transcol(OKgy));
    
    d_str = struct([]);
    if sum(OKgy) == 0 % no sampling in this site/year
    
        d_str(1).Count = NaN;
        d_str(1).TL = NaN;
        d_str(1).Month = NaN;
        d_str(1).numtrans = 0;
        
    else % else if there was sampling
   
    
    % If there were no fish
    if isempty(Months) % this will only happen if no fish observed
        Months2 = unique(Monthcol(OKgy));
        d_str(1).Count = 0;
        d_str(1).TL = NaN;
        d_str(1).Month = Months2;
        d_str(1).numtrans = 1; % doesn't matter...
    else
    
        transects = zeros(length(Months),1);
    for m = 1:length(Months)
        OKsgym = OKgy & OKsp & Monthcol == Months(m);

        % count up number of transects
        Side = unique(Sidecol(OKsgym)); % total number of sides in this month
        for sd = 1:length(Side)
            OKsgyms = OKsgym & strcmp(Sidecol,Side{sd});
            zone = unique(Zonecol(OKsgyms));
            for zn = 1:length(zone)
                OKsgymsz = OKsgyms & strcmp(Zonecol, zone{zn});
                transects(m) = transects(m) + length(unique(Transcol(OKsgymsz)));
            end
        end
        
        if m == 1
        d_str(1).Count = Countcol(OKsgym);
        d_str(1).TL = Fish_TLcol(OKsgym);
        d_str(1).Month = Months(m);
        else
          %  keyboard
            extraCount = Countcol(OKsgym);
            extraTL = Fish_TLcol(OKsgym);
            extraMonth = Months(m);
            d_str(1).Count = vertcat(d_str(1).Count,extraCount);
            d_str(1).TL = vertcat(d_str(1).TL,extraTL);
            d_str(1).Month = vertcat(d_str(1).Month,extraMonth);
        end
    
    
    end % end loop over months
    
    
    d_str.numtrans = sum(transects);
    
    end % end if isempty(Months)
    
    
    end % if sum(OKgy)
    
    D_str(1).(Species_Names{s}).(Site_Names{g})(y).data = d_str;
    D_str(1).(Species_Names{s}).(Site_Names{g})(y).year = Years(y);
    
    end % end loop over years
    end % end loop over sites
end % end loop over spp

    save(Meta.data_savename,'D_str','Species_Names','Site_Names','Years')
    
    

