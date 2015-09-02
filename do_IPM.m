%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run the IPM
function [L, Npred_out, Nact_out] = do_IPM(meshsize,x,dy, Sy,cand_vec,fixparm,Rfact,Rvec,data,T,Tpre,Ogive,Meta,NT,SiteName,pauseflag)

if ~exist('pauseflag','var')
    pauseflag = false;
end

Symat= repmat(Sy(:)',[length(Sy),1]);

% extract the parameters
sigma_o = cand_vec(end-1);
sigma_p = cand_vec(end);
F1 = cand_vec(1); % F in pre-data period
F2 = cand_vec(2); % F in data period
F = repmat(F2,length(T),1); % expand into vector
if ~isnan(Tpre)
F(1:length(Tpre)) = F1; % if Tpre = NaN, then this is a post 2007 run and we ignore this step
end 
R1 = exp(cand_vec(3))-Rfact; % recruitment in pre-data period
R2 = exp(cand_vec(4:end-2))-Rfact; % Rfact is log(0) correction factor
R = max(0,[repmat(R1,1,length(Tpre)), R2]); % constant recruitment in pre-data period, then R2. Be sure to restrict to > 0


Q = Meta.Q; % number of particles

% Initialize variables
N = nan(meshsize,length(T)); % population size
P = nan(meshsize,meshsize,length(T)); % variance


% If data are missing (i.e. no observations) in a given year, we cannot
% calculate a likelihood.  So just skip over those.
% do this by inserting a bunch of NaNs for pre-data years
if ~isnan(Meta.Tpre) % pre2007 run
Nact = [nan(size(data,1),length(Tpre)), data]; 
else % post2007 run
Nact = data;
end


% Do projections & calculate likelihoods
% Following deValpine & Hastins for likelihood calculation

% Determine size classes to be included in likelihood (account for size
% selectivity of sampling)
OKlen = Ogive; % probability of being observed
% if there is also a < 1 probability of being observed as a recruit, that
% needs to be included too.

% 1) Calculate P(y1) = P(n1)*P(y1|n1) (summed over all n1 values)
% (assume flat prior for P(n1)).  

% Create the kernel:
kmat = kernmatSimp(x,F(1),fixparm,1);

% Simpson's integration:
kmat = Symat.*kmat;

% Do a quick iteration to get the initial size distribution.
% This is necessary if we are starting in pre-data years (pre2007 runs).
% If starting in 2007, need to load in the starting distribution from
% the pre2007 fit.

if ~isnan(Tpre(1))
N0 = ones(size(kmat,1),100);
for t = 2:100
N0(:,t) = kmat*N0(:,t-1) + mean(R)*Rvec;
end
N(:,1) = N0(:,end); 

else

N(:,1) = Meta.N_init.(SiteName).Npred(:,end);
    
end


% Particle filter: (following Knape & deValpine 2012)
% Generate Q particles (independent simulations of N):
%Q = 500;
Nf = repmat(N,[1,1,Q]);
Nf(:,1,:) = max(0,Nf(:,1,:) + normrnd(0,sigma_p,size(Nf(:,1,:))));

% Use midpoint rule here to convert distribution to predicted number of
% fish, then assume normal distribution (could make So proportional to mean
% if desired)
%ftmp(1,:) = sum(log(max(realmin,normpdf(dy*Nf(:,1,:)-repmat(Nact(:,1),[1,1,Q]),0,sigma_o))));

ftmp(1,:) = ones(1,Q); % for year 1, there will never be any data, so use flat distribution

% resample for accuracy
Wgt = cumsum(squeeze(ftmp(1,:)/sum(ftmp(1,:))));
Rnd = rand(Q,1);
Wgt = repmat(Wgt(:),[1,Q]);
Rnd = repmat(Rnd(:)',[Q,1]);
Ind = Q - sum(Rnd < Wgt) + 1;

Nf(:,1,:) = Nf(:,1,Ind); % replace with resampled values

Ind2 = randperm(Q,1);
N(:,1) = Nf(:,1,Ind2); % pick one randomly to be *the* distribution

for t = 2:length(T)
    
    % Create the kernel:
    kmat = kernmatSimp(x,F(t),fixparm,1);

    % Simpson's integration:
    kmat = Symat.*kmat;
    
   % Advance the model (note - uses mmat for multidimensional matrix
   % multiplication)
 %  Nf(:,t,:) = max(0,mmat(kmat,Nf(:,t-1,:)) + repmat(Rvec*R(t),[1,1,Q]) + normrnd(0,sigma_p,size(Nf(:,t-1,:))));

%Nf(:,t,:) = max(0,mmat(kmat,Nf(:,t-1,:)) + repmat(Rvec*R(t),[1,1,Q]) + normrnd(0,sigma_p,size(Nf(:,t-1,:))));
    Nrand = normrnd(0,sigma_p,size(Nf(:,t-1,:)));
    RR = Rvec*R(t);
    for q = 1:Q
      
      Nf(:,t,q) = kmat*Nf(:,t-1,q) + RR + Nrand(:,:,q);  
      
    end

 
   Nf(:,t,:) = max(0,Nf(:,t,:));
   
    % Weights:
    if ~isnan(Nact(1,t)) % if there are data for that year
        % compare the integrals; so they have to be using the same
        % integration mesh
   % ftmp(t,:) = sum(log(max(realmin,normpdf((dy*Nf(:,t,:)).*repmat(OKlen(:),[1,1,Q])-repmat(dy*Nact(:,t),[1,1,Q]),0,sigma_o))));
    ftmp(t,:) = sum(log(max(realmin,poisspdf(repmat(Nact(:,t),[1,1,Q]),(NT(t-length(Tpre))*dy*Nf(:,t,:)).*repmat(OKlen(:),[1,1,Q])))));
    else
    ftmp(t,:) = 1; % if no data, just insert a bunch of 1s
    end
        
    if any(isinf(ftmp(t,:)))
        keyboard
    end
    
    % Resample
    Wgt = cumsum(squeeze(ftmp(t,:)/sum(ftmp(t,:))));
    Rnd = rand(Q,1);
    Wgt = repmat(Wgt(:),[1,Q]);
    Rnd = repmat(Rnd(:)',[Q,1]);
    Ind = Q - sum(Rnd < Wgt) + 1;

    
    Nf(:,t,:) = Nf(:,t,Ind); % replace with resampled values
    Ind2 = randperm(Q,1);
    N(:,t) = Nf(:,t,Ind2);
end


if pauseflag
    keyboard
end

%L = -( sum(Lo) + sum(Lp) );

L = nansum(mean(ftmp,1));

if exist('Nact','var')
    Npred_out = N;
    Nact_out = Nact;
else
    Npred_out = N;
    Nact_out = NaN;
end
%also output the last month's fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
