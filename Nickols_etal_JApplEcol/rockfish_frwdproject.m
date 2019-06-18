function rockfish_frwdproject(Meta_savename)

%use with runme_IPMfrwd.m

%recruitment is randomly selected from a distribution of recruitment that 
%results from model fits to data, forward running simulations based on an 
%initial size distribution from PISCO data

%STEP 1: Get all parameters in
load(Meta_savename) %load in metadata
load(Meta.Ffit_savename)
F=Post.F_mean;

%timescale to run model
T = 31;

% IPM parameters:
meshsize = Meta.IPM.meshsize;
meshmin = Meta.IPM.meshmin;
meshmax = Meta.IPM.meshmax;

fixparm=Meta.fixparm;

x = linspace(meshmin,meshmax,meshsize);
meshdiff = diff(x(1:2));
edges = x - meshdiff/2;

%load file that has size data in it from PISCO
load(Meta.data_savename,'D_str','Species_Names','Site_Names','Years')

% Which sites and species:
Site_Names = Meta.Sites;
Species_Names = Meta.Species;

%get data into a histogram for site running projection based on
N = IPM_histo(D_str.(Species_Names{1}),Years,Site_Names,edges);

%data then need to be corrected for the number of transects carried out per
%year. 

%divide total fish numbers by number of transects
for i = 1:length(Years)
    for j = 1:length(Site_Names)
        NT = D_str.(Species_Names{1}).(Site_Names{j})(i).data.numtrans;
        N(:,j,i) = N(:,j,i)./NT;
    end
end

%start with starting distribution from 2007 for site that is a new MPA
N_init=N(:,Meta.MPAnew,9);

%specify max size of recruits in data (YOY)
Rsize = Meta.recruits.Rsize; 

%add part to use recruit from model fits instead of data
%Assemble fit data for each site into one matrix with site as second
%dimension and year as third dimension (organize like data)

sitecount = length(Meta.Sites);
yearuse = 10:18;
N_predfits = NaN(meshsize,sitecount,length(yearuse));

for i = 1:sitecount
    for j=1:9
        N_predfits(:,i,j)=Post.predicted.(Meta.Sites{i}).Npred(:,j+9);
    end
end

%simulate recruitment from distribution based on data

%get sum of recruits for each year per site. 
R1 = squeeze(nansum(N_predfits(x<Rsize,:,:),1));
%find mean and sd for each year across sites
R1m = nanmean(R1(:,:),1);
R1s = nanstd(R1(:,:));
Ryr = nanmean(R1m);

%make vector of recruits to add in each year, this will be multiplied by
%the magnitude calculated above. Specify mean of pdf for size according to
%specific species, specify sd too.

R = normpdf(x,Meta.recruits.meansize,Meta.recruits.sdsize)'; %Average size 
                                                             %of recruit, 
                                                             %with sd = 1
                         
%scale R so integrates appropriately
dy=x(2)-x(1);
Sy=makeSimpVec(dy,meshsize);
Symat= repmat(Sy(:)',[length(Sy),1]); % make Simpsons matrix for later 
                                      % array operations of projection runs
                                     
% Determine size classes to be included (account for size selectivity of sampling)
% see White et al. 2016 Ecological Applications
Ogive_b = Meta.ogive;
if isnan(Ogive_b(1)) % if there is no ogive
    Ogive = ones(size(x));
else
Ogive = 1-normcdf(x,Ogive_b(1),Ogive_b(2)); % probability of observation in the kelp forest
end
OKlen = Ogive; % probability of being observed

RR=1000; % # of runs

%set up matrix to hold model runs 
N_R = nan(length(x),T,RR); %new reserve 
N_F = N_R;  %continued fishing median of F distribution from fits
sigma_p=1e-1; %add process error parameter
for r=1:RR
    N_R(:,1,r)=N_init./OKlen';
    N_F(:,1,r)=N_R(:,1,r);

    for t = 2:T %loop over years
        [kmatR] = kernmatSimp(x, 0, fixparm, 1);
        [kmatF] = kernmatSimp(x, F, fixparm, 1);
       
        kmatR = Symat.*kmatR;
        kmatF = Symat.*kmatF;
        
        Radd=(R.*poissrnd(Ryr))+normrnd(0,sigma_p,size(N_R(:,1,1)));
        
        %RESERVE
        N_R(:,t,r) = kmatR*N_R(:,t-1,r)+Radd;
        N_R(:,t,r) = max(0,N_R(:,t,r));

        %FISHING
        N_F(:,t,r) = kmatF*N_F(:,t-1,r)+Radd;
        N_F(:,t,r) = max(0,N_F(:,t,r));

    end
end

%save output

save(Meta.frwdfit_savename)

