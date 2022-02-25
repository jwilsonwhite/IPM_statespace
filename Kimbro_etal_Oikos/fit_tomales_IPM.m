function fit_tomales_IPM(N,Meta_savename,Model_str,pauseflag)
% Code for running the oyster IPM for Tomales Bay data
% JWWhite 01.25.15

% Load in metadata
load(Meta_savename)

% Savename for this run:
Fit_name_tmp = Meta.fit_name;
Fit_name = strcat(['Tomales_IPM_fits/Tomales_Model_',num2str(Model_str.model_num),'_',Fit_name_tmp]);

%Get all parameters in

% Some things needed for integration:
x = Meta.IPM.x;
dy=x(2)-x(1);
Sy=makeSimpVec(dy,Meta.IPM.meshsize);

%fixparm [none right now??]
fixparm = Meta.fixparm;

%timescale to run model
Tdata = Meta.Tdata; % timespan of data
Tpre = Meta.Tpre; % time prior to data that model will run

% Determine model type to set up parameter vectors:
Prior = get_priors(Model_str,Meta);
Par(1).R = Prior.R_prior;
Par(1).aM = Prior.aM_prior;
Par(1).jM = Prior.jM_prior;
Par(1).pe = Prior.pe_prior;
Par(1).oe = Prior.oe_prior;

ParSD(1) = Par; % this will hold the SDs of each parameter for the candidate generating function. Will be populated below...
ParSD(1).R = ones(size(Par.R));
ParSD(1).aM = ones(size(Par.aM));
ParSD(1).jM = ones(size(Par.jM));
ParSD(1).pe = Prior.pe_hyper;
ParSD(1).oe = Prior.oe_hyper;

Npar = length(Par.R) + length(Par.aM) + length(Par.jM) + 1; % total # params
CV = 1./([0.1 0.5 1 5:3:20])'; % CVs for delayed rejection

% Get information on recruitment for constructing IPM:
D = importdata('Data/SizevsAge.csv'); 
D = D(D(:,2)==0,:); % only age-0 information

%specify max size of recruits in data (YOY)
Rsize = max(D(:,3)+1.96*D(:,4)); % upper limit of distribution, mean + 1.96 SD
switch Model_str.growth
    case 'Var'
        %get sum of recruits for each year per site. 
        R1 = squeeze(nansum(N(x<Rsize,:,:),1));
        R1 = R1./Meta.Sample_size(:,Meta.Tdata); % scale by sample size
        R1m = nanmean(R1(:)); % mean over all sites & years
        %log transform time series
        Rfact = min(1,min(R1m(R1m>0))); % addition factor to avoid log(0)
        R1lnm = log(R1+Rfact);
        %find mean and std of the time series (over time)
        mu1 = mean(R1lnm,2);
        sig1 = std(R1lnm,[],2);
        %
        % In this case, mu1 and sig1 are vectors (1 entry per site)
        
        Meta.recruits.meansize = D(:,3);
        Meta.recruits.sdsize = D(:,4);
        
    case 'Con'
        %get sum of recruits for each year per site. 
        R1 = squeeze(nansum(N(x<Rsize,:,:),1));
        R1 = R1./Meta.Sample_size(:,Meta.Tdata); % scale by sample size
        R1m = nanmean(R1(:)); % mean over all sites & years
        %log transform time series
        Rfact = min(1,min(R1m(R1m>0))); % addition factor to avoid log(0)
        R1lnm = log(mean(R1,1)+Rfact); % mean across all sites
        %find mean and std of the time series (over time)
        mu1 = mean(R1lnm,2);
        sig1 = std(R1lnm,[],2);
        % Here mu1 and sig1 are scalars
        
        Meta.recruits.meansize = repmat(mean(D(:,3)),[Meta.Nsites,1]);
        Meta.recruits.sdsize = repmat(mean(D(:,4)),[Meta.Nsites,1]);
        
end % end switch Var

%make vector of recruits to add in each year (R), this will be multiplied by
%the magnitude calculated above. Specify mean of pdf for size according to
%specific species, specify sd too.
% also create vector isJuv that indicates whether an individual is likely a
% juvenile or not
for k = 1:length(Meta.recruits.meansize)
Rvec(1:length(x),k) = normpdf(x,Meta.recruits.meansize(k),Meta.recruits.sdsize(k))';
isJuv(1:length(x),k) = 1-normcdf(x,Meta.recruits.meansize(k)+2*Meta.recruits.sdsize(k),1)';
end
                   
% Get growth parameters:
    load('Tomales_vB_fits.mat','X')
switch Model_str.growth
    case 'Var'
        Par(1).Growth = X;
    case 'Con'
        Par(1).Growth = repmat(mean(X),[7,1]); % just repeat the mean over several sites
end % end switch growth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do MCMC
M = Meta.MCMC.M;
chains = Meta.MCMC.chains;
stepmax = length(CV);

mc_str = struct([]); % structure to hold results

for c = 1:chains
    
    
    % Book-keeping for DR process:
    k = 1; % parameter counter
    kk = 1; % delayed rejection counter
    rej = 0; % count number of rejections
    rej2= 0 ; % count number of stepped rejections
    acc = 0; % count number of acceptances
    
    % reset parameter vector at the beginning of each chain
    Par = Par(1);
    
    % Initialize likelihood 
    L = [];
    Ltmp = nan(Meta.Nsites,1);
    for s = 1:Meta.Nsites
    Ltmp(s) = do_Tomales_IPM(dy,Sy,Par(1),Rvec(:,s),isJuv(:,s),squeeze(N(:,s,:)),Meta,Model_str,s);
    end
    
    % if we do need to do the runs for checking Q (will need to fix this if needed) 
    doQruns = false;
    if doQruns
        
        Q = [5 10 30 50 100:100:500];
        LQ = nan(10,length(Q));
        for q = 1:length(Q)
            Meta.Q = Q(q);
        for i = 1:100
            tic
            s = 4;
            LQ(i,q) = do_Tomales_IPM(dy,Sy,Par(1),Rvec(:,s),isJuv(:,s),squeeze(N(:,s,:)),Meta,Model_str,s);
            Time(q)= toc
        end
        end
        
        keyboard
        
    end % end if doQruns
    
    % Add prior calculation:
    PriorL(1) = sum(log(max(realmin,normpdf(Par.R,Prior.R_prior,Prior.R_hyper))));
    PriorL(2) = sum(log(max(realmin,normpdf(Par.aM,Prior.aM_prior,Prior.aM_hyper))));
    PriorL(3) = sum(log(max(realmin,normpdf(Par.jM,Prior.jM_prior,Prior.jM_hyper))));
    PriorL(4) = sum(log(max(realmin,gampdf(1./Par.pe,Prior.pe_hyper,1/Prior.pe_prior./(Prior.pe_hyper-1)))));
    %PriorL(5) = sum(log(max(realmin,gampdf(1./Par.oe,Prior.oe_hyper,1./Prior.oe_prior./(Prior.oe_hyper-1)))));
    
    
    L(1) = sum(Ltmp)+sum(PriorL);
    
   % keyboard
    
    for m = 1:M
        
        if mod(m,1e4)==0
            disp(strcat(['m = ',num2str(m)]))
            %keyboard
        end
        
        next_step = false;
        step = 1;
        
        while ~next_step
            
        % simulate candidate values 
        
        Par_cand = Par(kk); % start with existing parameter state
        Par_cand = cand_get_fxn(Par_cand,ParSD,CV,step,k);
        
        % Run the model:
        for s = 1:Meta.Nsites
        Ltmp(s) = do_Tomales_IPM(dy,Sy,Par_cand,Rvec(:,s),isJuv(:,s),squeeze(N(:,s,:)),Meta,Model_str,s);
        end
        
          % Add prior calculation:
          PriorL(1) = sum(log(max(realmin,normpdf(Par_cand.R,Prior.R_prior,Prior.R_hyper))));
          PriorL(2) = sum(log(max(realmin,normpdf(Par_cand.aM,Prior.aM_prior,Prior.aM_hyper))));
          PriorL(3) = sum(log(max(realmin,normpdf(Par_cand.jM,Prior.jM_prior,Prior.jM_hyper))));
          PriorL(4) = sum(log(max(realmin,gampdf(1./Par_cand.pe,Prior.pe_hyper,1./Prior.pe_prior./(Prior.pe_hyper-1)))));
          %PriorL(5) = sum(log(max(realmin,gampdf(1./Par_cand.oe,Prior.oe_hyper,1./Prior.oe_prior./(Prior.oe_hyper-1)))));
          
        L_cand = sum(Ltmp)+sum(PriorL);
        
        % M-H
        if ~isnan(L_cand) && ~isinf(L_cand) && isreal(L_cand)
            MH_prob = min(1,exp(L_cand-L(kk)));
        else
            MH_prob = 0;
        end
        
        K = MH_prob > rand;
        
        if kk>10 && all(L(end-5:end)==L(end)) % get out of trapping states caused by a rogue high L
            Par_cand = Par(kk); % start with existing parameter state
        
            % Run the model:
            for s = 1:Meta.Nsites
            Ltmp(s) = do_Tomales_IPM(dy,Sy,Par_cand,Rvec(:,s),isJuv(:,s),squeeze(N(:,s,:)),Meta,Model_str,s);
            end
        
          % Add prior calculation:
          PriorL(1) = sum(log(max(realmin,normpdf(Par_cand.R,Prior.R_prior,Prior.R_hyper))));
          PriorL(2) = sum(log(max(realmin,normpdf(Par_cand.aM,Prior.aM_prior,Prior.aM_hyper))));
          PriorL(3) = sum(log(max(realmin,normpdf(Par_cand.jM,Prior.jM_prior,Prior.jM_hyper))));
          PriorL(4) = sum(log(max(realmin,gampdf(1./Par_cand.pe,Prior.pe_hyper,1./Prior.pe_prior./(Prior.pe_hyper-1)))));
        %  PriorL(5) = sum(log(max(realmin,gampdf(1./Par_cand.oe,Prior.oe_hyper,1./Prior.oe_prior./(Prior.oe_hyper-1)))));
          
            L_cand = sum(Ltmp)+sum(PriorL);
            K = 1;
        end
        
        % DR
        if K
            
        Par(kk+1) = Par_cand;
        L(kk+1) = L_cand;
        next_step = true;
        step = 1;
        acc = acc + 1;
        
        elseif ~K && step < stepmax % move to alternate proposal
        
        step = step + 1;
        rej2 = rej2 + 1;
        
        else % if actually reject
    
        Par(kk+1) = Par(kk);
        L(kk+1) = L(kk);
        next_step = true;
        rej = rej + 1;
        
        end % end if K
        end % end if delayed rejection while loop
            
    % advance counters
    kk = kk+1;
    k = k + 1;
    if k > Npar; k = 1; end
     
    end % end loop over M
    
mc_str(c).L = L;
mc_str(c).Par = Par;
mc_str(c).acc = acc;
mc_str(c).rej = rej;
mc_str(c).rej2 = rej2;

% A more user-friendly way to store the final parameter values:
R_vec = nan(M,length(Par(1).R));
aM_vec = nan(M,length(Par(1).aM));
jM_vec = nan(M,length(Par(1).jM));
%oe_vec = nan(M,1);
pe_vec = nan(M,1);
for i = 1:M
    R_vec(i,:) = Par(i).R;
    aM_vec(i,:) = Par(i).aM;
    jM_vec(i,:) = Par(i).jM;
   % oe_vec(i) = Par(i).oe;
    pe_vec(i) = Par(i).pe;
end

mc_str(c).R = R_vec;
mc_str(c).aM = aM_vec;
mc_str(c).jM = jM_vec;
%mc_str(c).oe = oe_vec;
mc_str(c).pe = pe_vec;
    
mc_str(c).R_vec = Rvec;
mc_str(c).Rfact = Rfact;
mc_str(c).isJuv = isJuv;
mc_str(c).N = N;

mc_str(c).Model_str = Model_str;

end % end loop over chains



save(Fit_name,'mc_str')

% If you want to pause the code to look at the fit
if ~exist('pauseflag','var')
    pauseflag = false;
end
if pauseflag

[L1, Npred,Nact] = do_Tomales_IPM(dy,Sy,Par_cand,Rvec(:,s),isJuv(:,s),squeeze(N(:,s,:)),Meta,Model_str,s);

keyboard
end


%-------------------------------------------------------------------------
% Candidate-generating function:
function Par_out = cand_get_fxn(Par,ParSD,CV,step,k)

% Determine which parameter to vary
L = cumsum([length(Par.R),length(Par.aM),length(Par.jM),1]);
kL = k - L;

if kL(1) <= 0 && all(kL(2:4)<0) % Par.R
    Par.R(k) = normrnd(Par.R(k),ParSD.R(k)*CV(step));
    
elseif kL(2) <= 0 && all(kL(3:4) < 0) && kL(1) > 0 % 
    Par.aM(kL(1)) = normrnd(Par.aM(kL(1)),ParSD.aM(kL(1))*CV(step));
    
elseif kL(3) <= 0 && all(kL(4:4) < 0) && all(kL(1:2)>0) % 
    Par.jM(kL(2)) = normrnd(Par.jM(kL(2)),ParSD.jM(kL(2))*CV(step));
    
elseif kL(4) == 0
    Par.pe = sqrt( 1/gamrnd(1/CV(step),1/Par.pe.^2/(max(2,1/CV(step))-1)));
    
%elseif kL(5) == 0
%    Par.oe = sqrt( 1/gamrnd(1/CV(step),1/Par.oe.^2/(max(2,1/CV(step))-1)));
    
end % end if

Par_out = Par;


