function read_PISCOdata_zone2(meta_savename)

load(meta_savename); % load in metadata needed to get species, sites, etc.

% Read in PISCO dataset, convert to array with size structure by year &
% site
if ~exist(Meta.PISCO_filename,'file') % if the file doesn't exist
Dir = '../PISCO_data_2015/';
File = 'PISCO_FISH_v2.csv'; % V2 lacks headers, the original has headers 

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
        d_str(1).Zone(1).Count = NaN;
        
    else % else if there was sampling
   
    
    % If there were no fish
    if isempty(Months) % this will only happen if no fish observed
        Months2 = unique(Monthcol(OKgy));
        d_str(1).Count = 0;
        d_str(1).TL = NaN;
        d_str(1).Month = Months2;
        d_str(1).numtrans = 1; % doesn't matter...
        %d_str(1).Zone(1:4).Count = 0;
        d_str(1).Zone(:).Count = 0;
    else
    
        transects = zeros(length(Months),1);
        transects_perzone = zeros(length(Months),1); % 'zone' on the columns 
        Zone_str = struct([]); % new structure to hold per-zone data
        
       z_cnt = 0; % counter needed bc we are looping over sides but ignoring side when tallying up zones

    for m = 1:length(Months)
        OKsgym = OKgy & OKsp & Monthcol == Months(m);
        
        % count up number of transects & do tally at 'zone' scale
        Side = unique(Sidecol(OKsgym)); % total number of sides in this month
        for sd = 1:length(Side)
            OKsgyms = OKsgym & strcmp(Sidecol,Side{sd});
            Zone = unique(Zonecol(OKsgyms));
            
         %   transects_perside_tmp = 0; % restart this count
            for zn = 1:length(Zone)
                OKsgymsz = OKsgyms & strcmp(Zonecol, Zone{zn});
                transects_perzone(m,z_cnt + zn) = length(unique(Transcol(OKsgymsz)));
                transects(m) = transects(m) + length(unique(Transcol(OKsgymsz)));
              %  transects_perside_tmp = transects_perside_tmp + length(unique(Transcol(OKsgymsz)));
            
          %  transects_perside(m,sd) = transects_perside_tmp;
            
          %  if m == 1
                if isempty(Countcol(OKsgymsz));
                    Zone_str(zn+z_cnt).Count = 0;
                    Zone_str(zn+z_cnt).TL = NaN;
                else
                Zone_str(zn+z_cnt).Count = Countcol(OKsgymsz);
                Zone_str(zn+z_cnt).TL = Fish_TLcol(OKsgymsz);
                end
            Zone_str(zn+z_cnt).Month = Months(m);
            Zone_str(zn+z_cnt).Side_name = Sides{sd};
            Zone_str(zn+z_cnt).Zone_name = Zones{zn};
          %  else
          %  extraCount = Countcol(OKsgymsz);
          %  extraTL = Fish_TLcol(OKsgymsz);
          %  extraMonth = Months(m);   
          %  Zone_str(zn+z_cnt).Count = vertcat(Zone_str(sd).Count,extraCount); 
          %  Zone_str(zn+z_cnt).TL = vertcat(d_str(1).TL,extraTL);
          %  Zone_str(zn+z_cnt).Month = vertcat(d_str(1).Month,extraMonth);
          % [don't just add on the count like we do when aggregating by
          % site, because there is no guarantee that the same combination
          % of zones was surveyed in the second month. Just add on as
          % additional zones.
          %  Zone_str(zn+z_cnt).Count = Countcol(OKsgymsz);
          %  Zone_str(zn+z_cnt).TL = Fish_TLcol(OKsgymsz);
          %  Zone_str(zn+z_cnt).Month = Months(m);
          %  end
            
            end % end loop over zones
            
            z_cnt = z_cnt + zn; % so that in the next loop we just add on to the end instead of replacing
            
        end % end loop over sides
        
    
        
        if m == 1
        d_str(1).Count = Countcol(OKsgym);
        d_str(1).TL = Fish_TLcol(OKsgym);
        d_str(1).Month = Months(m);
        else
            
            extraCount = Countcol(OKsgym);
            extraTL = Fish_TLcol(OKsgym);
            extraMonth = Months(m);
            d_str(1).Count = vertcat(d_str(1).Count,extraCount);
            d_str(1).TL = vertcat(d_str(1).TL,extraTL);
            d_str(1).Month = vertcat(d_str(1).Month,extraMonth);
        end
    
    
    end % end loop over months
    
    
    d_str.numtrans = sum(transects);
    d_str.transects_perzone = transects_perzone; % this is a matrix with months on the rows and zones on the columns
    d_str.total_transects_perzone = sum(transects_perzone); % if there are multiple months surveyed per year, this gives the total # transects surveyed in each of the zones.
    d_str.Zone = Zone_str;
    
    end % end if isempty(Months)
    
    
    end % if sum(OKgy)
    
    D_str(1).(Species_Names{s}).(Site_Names{g})(y).data = d_str;
    D_str(1).(Species_Names{s}).(Site_Names{g})(y).year = Years(y);
    
    end % end loop over years
    end % end loop over sites
end % end loop over spp

    save(Meta.data_savename,'D_str','Species_Names','Site_Names','Years')
    
    

