function findSSD(Meta_savename)

%find the stable size distribution for Central California Coast blue 
%rockfish, uses recruitment magnitude from IPM fits from White et al. 2016

load(Meta_savename)
load(Meta.Ffit_savename)
F=Post.F_mean;

fixparm = Meta.fixparm;
meshsize = Meta.IPM.meshsize;
meshmin = Meta.IPM.meshmin;
meshmax = Meta.IPM.meshmax;
x = linspace(meshmin,meshmax,meshsize);
meshdiff = diff(x(1:2));
edges = x - meshdiff/2;

sitecount = length(Meta.Sites);
yearuse = 10:18;
N_predfits = NaN(meshsize,sitecount,length(yearuse));

for i = 1:sitecount
    for j=1:9
        N_predfits(:,i,j)=Post.predicted.(Meta.Sites{i}).Npred(:,j+9);
    end
end

%specify max size of recruits in data (YOY)
Rsize = Meta.recruits.Rsize; 

%simulate recruitment from distribution based on data

%get sum of recruits for each year per site. 
R1 = squeeze(nansum(N_predfits(x<Rsize,:,:),1));
%find mean and sd for each year across sites
R1m = nanmean(R1(:,:),1);
R1s = nanstd(R1(:,:));
Ryr = nanmean(R1m);

%Generate 'unfished' kernel
%F=0
R = normpdf(x,Meta.recruits.meansize,Meta.recruits.sdsize)'; 

dy=x(2)-x(1);
Sy=makeSimpVec(dy,meshsize);
Symat= repmat(Sy(:)',[length(Sy),1]); % make Simpsons matrix for later 
                                      % array operations of projection runs                                     

Ogive_b = Meta.ogive;

if isnan(Ogive_b(1)) % if there is no ogive
    Ogive = ones(size(x));
else
Ogive = 1-normcdf(x,Ogive_b(1),Ogive_b(2)); % probability of observation in the kelp forest
end
OKlen = Ogive; % probability of being observed

% Create the kernel:
kmat = kernmatSimp(x,0,fixparm,1);
kmatF = kernmatSimp(x,F,fixparm,1);

% Simpson's integration:
kmat = Symat.*kmat;
kmatF = Symat.*kmatF;

% Do a quick iteration to get the initial size distribution.
N0 = ones(size(kmat,1),100);
for t = 2:100
N0(:,t) = kmat*N0(:,t-1) + R.*Ryr;
end
N_init = N0(:,end); 

%Fish for 30 years constant recruitment 
T=30; %timescale to run model
N1 = ones(size(kmatF,1),T);
N1R = ones(size(kmat,1),T);
N1(:,1)=N_init;
N1R(:,1)=N_init;
for t=2:T
    N1(:,t) = kmatF*N1(:,t-1) +R.*Ryr;
    N1R(:,t) = kmat*N1R(:,t-1) +R.*Ryr;
end
N2_init = N1(:,end);
N2R_init = N1R(:,end);

%Stop fishing for 30 years constant recruitment 
N2_rc = ones(size(kmat,1),T);
N2_rc(:,1)=N2_init;
N2R = ones(size(kmat,1),T);
N2F = ones(size(kmatF,1),T);
N2R(:,1)=N2R_init;
N2F(:,1)=N2_init;
for t=2:T
    N2_rc(:,t)=kmat*N2_rc(:,t-1)+R.*Ryr;
    N2F(:,t)=kmatF*N2F(:,t-1) +R.*Ryr;
    N2R(:,t) = kmat*N2R(:,t-1) +R.*Ryr;
end

%Add recruitment variability

RR=1000;
sigma_p=1e-1; %add process error parameter, only used for random runs

%Fish for 30 years with random recruitment and stop fishing with
%random recruitment 30 years
TT=60; %total time
N1_rr = ones(size(kmatF,1),TT,RR);
N1F_rr = ones(size(kmatF,1),TT,RR);

for r = 1:RR
    
    N1F_rr(:,1,r)=N_init;
        for t=2:30
            N1F_rr(:,t,r) = kmatF*N1F_rr(:,t-1,r) +R.*poissrnd(Ryr)+normrnd(0,sigma_p,size(N1F_rr(:,1,1)));
            N1F_rr(:,t,r) = max(0,N1F_rr(:,t,r));
            N1_rr(:,t,r) = N1F_rr(:,t,r);
        end
end
            
for r=1:RR            
        
        for t=31:60
            Radd=(R.*poissrnd(Ryr))+normrnd(0,sigma_p,size(N1F_rr(:,1,1)));
            N1F_rr(:,t,r) = kmatF*N1F_rr(:,t-1,r)+Radd;
            N1F_rr(:,t,r)=max(0,N1F_rr(:,t,r));
            
            N1_rr(:,t,r) = kmat*N1_rr(:,t-1,r)+Radd;
            N1_rr(:,t,r)=max(0,N1_rr(:,t,r));
            
        end
end



%save output

save(Meta.SSDfit_savename)