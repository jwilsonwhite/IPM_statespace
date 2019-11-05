function rockfish_fit_pisco_cluster(Meta,Chain,pauseflag)

%updated code as generalized as possible for running the IPM from PISCO
%rockfish data forward to make projections from PISCO data.
% KJ Nickols and JW White 2015

% Load in metadata

%STEP 1: Get all parameters in
%fixparm = [Linf k M Lfish Lmat Lvar  v_grow]; 2.5 is variance around Lmat
%and Lfish
fixparm = Meta.fixparm;

%timescale to run model
Tdata = Meta.Tdata; % timespan of data

Tpre = Meta.Tpre; % time prior to data that model will run
if ~isnan(Tpre) % for pre-2007 runs
T = 1:length([Tpre(:); Tdata(:)]);
else
T = 1:length(Tdata(:));
end % end if ~isnan(Tpre)

% format integration mesh
meshsize = Meta.IPM.meshsize;
meshmin = Meta.IPM.meshmin;
meshmax = Meta.IPM.meshmax;
x = linspace(meshmin,meshmax,meshsize);
meshdiff = diff(x(1:2));
edges = x - meshdiff/2;

% get weighting vector for Simpson's integration
dy=x(2)-x(1);
Sy=makeSimpVec(dy,meshsize);

%load file that has size data in it from PISCO
load(Meta.data_savename,'D_str','Species_Names','Years')

Site_Names = Meta.Sites;

%get data into a histogram
N = IPM_histo(D_str.(Species_Names{1}),Years,Site_Names,edges);
% this will only have dimension (lengths, sites, years)

%data then need to be corrected for the number of transects carried out per
%year. 

% calculate total number of transects
%divide total fish numbers by number of transects
for i = 1:length(Years) 
    for j = 1:length(Site_Names)
        NT(j,i) = D_str.(Species_Names{1}).(Site_Names{j})(i).data.numtrans;
    end
end
if isnan(Tpre) % if it is a post-2007 run
    NT = NT(:,end-length(Tdata):end); 
    % only take the later years. include 2007 as t = 1; also include 2006
    % here because the code assumes there is 1 pre-data year
end
%LY: NT will be #rows = #sites, #columns = #years

% strip out the years we don't need
OKyears = false(length(Years),1);
for i = 1:length(Years)
    if any(Meta.Tdata == Years(i))
        OKyears(i) = true;
    end
end
N = N(:,:,OKyears);

% specify ogive for observations - don't see fish that are too large
% b/c of ontogenetic migration
Ogive_b = Meta.ogive;
if isnan(Ogive_b(1)) % if there is no ogive
    Ogive = ones(size(x));
else
Ogive = 1-normcdf(x,Ogive_b(1),Ogive_b(2)); % probability of observation in the kelp forest
end

%specify max size of recruits in data (YOY)
Rsize = Meta.recruits.Rsize;

% Obtain recruitment statistics

%get sum of recruits for each year per site. 
if isnan(Meta.Rprior) % if this is a pre-2007 run
R1 = squeeze(nansum(N(x<Rsize,:,:),1));

%find mean for each year across sites
if size(R1,1)>1 && size(R1,2) >1 % if more than one site
R1m = nanmean(R1(:,:),1);
else
R1m = R1;    % just one site, so no need to take mean
end
%log transform time series
Rfact = min(1,min(R1m(R1m>0))); % addition factor to avoid log(0)
R1lnm = log(R1m+Rfact);
%find mean and std of the time series
mu1 = mean(R1lnm);
sig1 = std(R1lnm);
else % if it's a post-2007 run, then we have priors already
mu1 = Meta.Rprior(1);
sig1 = Meta.Rprior(2);
Rfact = Meta.Rfact;
end % end if isnan Meta.Rprior

%make size density function of recruits to add in each year, this will be multiplied by
%the magnitude calculated above. Specify mean of pdf for size according to
%specific species, specify sd too.
R = normpdf(x,Meta.recruits.meansize,Meta.recruits.sdsize)'; 
                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do MCMC
M = Meta.MCMC.M;

if ~exist('Chain','var')
chains = Meta.MCMC.chains;
elseif isnan(Chain)
chains = Meta.MCMC.chains;
else
chains = 1;
end

mc_str = struct([]); % structure to hold results

for c = 1:chains

% Initial parameter vector
% F1, F2, Recruitment1, Recruitment2, sigma_p
prior_vec = [Meta.Fprior(1,1),Meta.Fprior(2,1), mu1, repmat(mu1,1,length(Tdata)), 1e-1];
parm_vec = prior_vec;
parmSD = [1,1, sig1, repmat(sig1,1,length(Tdata)), NaN]; % these are just used for the candidate generation
CV = repmat(1./([0.5 1 5:3:40])',[1,length(parm_vec)]); % CVs used for delayed rejection
            stepmax = size(CV,1);
            
            parm_gam1 = 2; % inverse gamma hyperparameter for process error. Low confidence.
            
            k = 1; % parameter counter
            kk = 1; % delayed rejection counter
            rej = 0; % count number of rejections
            rej2= 0 ; % count number of stepped rejections
            acc = 0; % count number of acceptances
            
            L = [];
            
            % calculate likelihood separately for each site
            Ltmp = nan(1,length(Site_Names));
            
            for s = 1:length(Site_Names)
            Site_Names{s};
            % switch off fishing for pre-existing MPAs
            parm_vec_tmp = parm_vec;
            parm_vec_tmp(1:2) = parm_vec_tmp(1:2).*(1-Meta.MPAstatus(s));
            
            Ltmp(s) = do_IPM(meshsize,x,dy,Sy,parm_vec_tmp,fixparm,Rfact,R,squeeze(N(:,s,:)),T,Tpre,Ogive,Meta,NT(s,:),Site_Names{s},pauseflag);
            
            end % end loop over s
            
            % This is a little subroutine for estimating the value of Q required:
            doQruns = false;
            if doQruns
            
            Q = [5 50 100:100:2000];
            LQ = nan(10,length(Q));
            for q = 1:length(Q)
            Meta.Q = Q(q);
            for i = 1:100
            tic
            LQ(i,q) = do_IPM(meshsize,x,dy,Sy,parm_vec_tmp,fixparm,Rfact,R,squeeze(N(:,s,:)),T,Tpre,Ogive,Meta,NT(s,:),Site_Names{s});
            Time(q)= toc
            end
            end
            
            keyboard
            
            end % end if doQruns
            
            % Add prior calculation:
            PriorL(1) = log(max(realmin,normpdf(log(parm_vec_tmp(1)),log(Meta.Fprior(1,1)),Meta.Fprior(1,2)))); % F1
            PriorL(2) = log(max(realmin,normpdf(log(parm_vec_tmp(2)),log(Meta.Fprior(2,1)),Meta.Fprior(2,2)))); % F2
            PriorL(3) = sum(log(max(realmin,normpdf(parm_vec_tmp(3:(3+length(Tdata))),mu1,sig1)))); % R
            PriorL(4) = sum(log(max(realmin,gampdf(1./parm_vec_tmp(end-1:end),parm_gam1,1./prior_vec(end-1:end)./(parm_gam1-1))))); % process error
            
            L(1) = sum(Ltmp)+sum(PriorL);
            
            % MCMC iterations
            for m = 1:M
            
            if mod(m,1000)==0
            disp(strcat(['m = ',num2str(m)])) % counter
            end
            
            next_step = false;
            step = 1;
            
            while ~next_step
            
            % simulate candidate values
            cand_vec = parm_vec(kk,:);
            if ~isnan(parmSD(k))
            cand_vec(k) = normrnd(cand_vec(k),parmSD(k)*CV(step,k));
            else
            cand_vec(k) = 1./gamrnd(parm_gam1/CV(step,k),1/cand_vec(k)/(parm_gam1/CV(step,k)));
            end
            
            cand_vec(1:end-2) = abs(cand_vec(1:end-2)); % can't have a negative fishing rate or recruitment!
            
            % Loop over each site and calculate likelihood
            Ltmp = nan(1,length(Site_Names));
            for s = 1:length(Site_Names)
            
            % switch off fishing for pre-existing MPAs
            cand_vec_tmp = cand_vec;
            cand_vec_tmp(1:2) = cand_vec_tmp(1:2).*(1-Meta.MPAstatus(s));
            
            Ltmp(s) = do_IPM(meshsize,x,dy,Sy,cand_vec_tmp,fixparm,Rfact,R,squeeze(N(:,s,:)),T,Tpre,Ogive,Meta,NT(s,:),Site_Names{s});
            
            end % end loop over s
            
            % Add prior calculation:
            PriorL(1) = log(max(realmin,normpdf(log(parm_vec_tmp(1)),log(Meta.Fprior(1,1)),Meta.Fprior(1,2)))); % F1
            PriorL(2) = log(max(realmin,normpdf(log(parm_vec_tmp(2)),log(Meta.Fprior(2,1)),Meta.Fprior(2,2)))); % F2
            PriorL(3) = sum(log(max(realmin,normpdf(parm_vec_tmp(3:(3+length(Tdata))),mu1,sig1)))); % R
            PriorL(4) = sum(log(max(realmin,gampdf(1./parm_vec_tmp(end-1:end),parm_gam1,1./prior_vec(end-1:end)./(parm_gam1-1))))); % Process error
            
            L_cand = sum(Ltmp)+sum(PriorL); % Candidate likelihood
            
            % Metropolis-Hastings step
            if ~isnan(L_cand) && ~isinf(L_cand)
            MH_prob = min(1,exp(L_cand-L(kk)));
            else
            MH_prob = 0;
            end
            K = MH_prob > rand;
            
            % Delayed rejection step
            if K % accept proposal
            
            parm_vec(kk+1,:) = cand_vec;
            L(kk+1) = L_cand;
            next_step = true;
            step = 1;
            acc = acc + 1;
            
            elseif ~K && step < stepmax % move to alternate proposal
            
            step = step + 1;
            rej2 = rej2 + 1;
            
            else % if actually reject
            
            parm_vec(kk+1,:) = parm_vec(kk,:);
            L(kk+1) = L(kk);
            next_step = true;
            rej = rej + 1;
            
            end % end if K
            end % end if delayed rejection while loop
            
            % advance counters
            kk = kk+1;
            k = k + 1;
            if k > length(cand_vec); k = 1; end
            
            end % end loop over M
            
mc_str(c).L = L;
mc_str(c).parm_vec = parm_vec;
mc_str(c).acc = acc;
mc_str(c).rej = rej;
mc_str(c).rej2 = rej2;
mc_str(c).NT = NT;
            
end % end loop over chains
            
save(Meta.fit_savename,'mc_str')
            
% If you want to pause the code to look at the fit
if ~exist('pauseflag','var')
pauseflag = false;
end
if pauseflag
for s = 1:length(Site_Names)
[L1(s), Npred,Nact] = do_IPM(meshsize,x,dy,Sy,mean(parm_vec(1:end,:)),fixparm,Rfact,R,squeeze(N(:,s,:)),T,Tpre,Ogive,Meta,NT(s,:),Site_Names{s});
            
end
keyboard
end




