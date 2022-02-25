%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run the IPM
function [L, Npred_out, Nact_out] = do_Tomales_IPM(dy,Sy,Par,Rvec,isJuv,data,Meta,Model_str,Si,pauseflag)

if ~exist('pauseflag','var')
    pauseflag = false;
end

if ~exist('T','var')
    T = 1; % no. of years to run the model. Not currently used
end

Symat= repmat(Sy(:)',[length(Sy),1]);

% Initialize variables
N = nan(Meta.IPM.meshsize,length(Meta.T)); % population size
P = nan(Meta.IPM.meshsize,Meta.IPM.meshsize,length(T)); % variance

% If data are missing (i.e. no observations) in a given year, we cannot
% calculate a likelihood.  So just skip over those.
% do this by inserting a bunch of NaNs for non-data years
Nact = N;

Nact(:,Meta.Tdata) = data; % insert data for the correct years

% Expand recruitment information, to be contained in 1 x years vector R:
switch Model_str.rect_time
    case 'Con'  
        switch Model_str.rect_space
            case 'Con'
                Rt = exp(Par.R);%+Rfact;
            case 'Grad'
                Rt = exp(Par.R(1) + Meta.Distances(Si).*Par.R(2));%+Rfact;
            case 'Unim'
                Rt = exp(Par.R(1) + Par.R(2)*(Meta.Distances(Si)-Par.R(3)).^2);%+Rfact;
        end % end rect_space switch
                R = repmat(Rt,[1,length(Meta.T)]); % repeat over years
               
                
    case 'Var' % different value for each year
        switch Model_str.rect_space
            case 'Con'
                R = exp(Par.R);%+Rfact;
            case 'Grad'
                R = exp(Par.R(1:end-1) + Meta.Distances(Si).*Par.R(end));%+Rfact;
            case 'Unim'
                R = exp(Par.R(1:end-2) + Par.R(end-1)*(Meta.Distances(Si)-Par.R(end)).^2);%+Rfact;
        end % end rect_space switch
    
end % end rect_time switch

% Do projections & calculate likelihoods
% Following deValpine & Hastins for likelihood calculation

% 1) Calculate P(y1) = P(n1)*P(y1|n1) (summed over all n1 values)
  
% Create the kernel:
kmat = kernmatSimp_Tomales(Meta.IPM.x,Par,isJuv,Model_str,Meta,1,Si);

% Simpson's integration:
kmat = Symat.*kmat;

% Do a quick iteration to get the initial size distribution.
N0 = ones(size(kmat,1),50);

for t = 2:50
N0(:,t) = kmat*N0(:,t-1) + mean(R)*Rvec;
end
N(:,1) = N0(:,end); 



% Particle filter: (following Knape & deValpine 2012)
% Generate Q particles (independent simulations of N):
%Q = 500;
Nf = repmat(N,[1,1,Meta.Q]);
Nf(:,1,:) = max(0,Nf(:,1,:) + normrnd(0,Par.pe,size(Nf(:,1,:))));

% Weighting function
ftmp(1,:) = ones(1,Meta.Q); % for year 1, there will never be any data, so use flat distribution

% resample for accuracy
Wgt = cumsum(squeeze(ftmp(1,:)/sum(ftmp(1,:))));
Rnd = rand(Meta.Q,1);
Wgt = repmat(Wgt(:),[1,Meta.Q]);
Rnd = repmat(Rnd(:)',[Meta.Q,1]);
Ind = Meta.Q - sum(Rnd < Wgt) + 1;

Nf(:,1,:) = Nf(:,1,Ind); % replace with resampled values

Ind2 = randperm(Meta.Q,1);
N(:,1) = Nf(:,1,Ind2); % pick one randomly to be *the* distribution

for t = 2:length(Meta.T)
    
    % if time is varying for relevant parameters, recreate the kernel:
    if strcmp(Model_str.aM_time,'Var') || strcmp(Model_str.jM_time,'Var')
    kmat = kernmatSimp_Tomales(Meta.IPM.x,Par,isJuv,Model_str,Meta,t-1,Si);
    kmat = Symat.*kmat;
    % otherwise the kernel does not change
    end
    
   % Advance the model 
    Nrand = normrnd(0,Par.pe,size(Nf(:,t-1,:)));
    RR = Rvec*R(t);
    for q = 1:Meta.Q
      
      Nf(:,t,q) = kmat*Nf(:,t-1,q) + RR + Nrand(:,:,q);  
      
    end

 
   Nf(:,t,:) = max(0,Nf(:,t,:));
   
    % Weights:
    if ~isnan(Nact(1,t)) % if there are data for that year
    ftmp(t,:) = sum(log(max(realmin,poisspdf(repmat(Nact(:,t),[1,1,Meta.Q]),Meta.Sample_size(Si,t).*dy*Nf(:,t,:)))));
    %ftmp(t,:) = sum(log(max(realmin,normpdf(dy*Nf(:,t,:)-repmat(Nact(:,t),[1,1,Meta.Q]),0,Par.oe))))
    else
    ftmp(t,:) = 1; % if no data, just insert 1 as a dummy
    end
        

    % Resample
    Wgt = cumsum(squeeze(ftmp(t,:)/nansum(ftmp(t,:))));
    Rnd = rand(Meta.Q,1);
    Wgt = repmat(Wgt(:),[1,Meta.Q]);
    Rnd = repmat(Rnd(:)',[Meta.Q,1]);
    Ind = Meta.Q - sum(Rnd < Wgt) + 1;
  
    Nf(:,t,:) = Nf(:,t,Ind); % replace with resampled values
    Ind2 = randperm(Meta.Q,1);
    N(:,t) = Nf(:,t,Ind2);
end


if pauseflag
    keyboard
end


L = nansum(nanmean(ftmp,2));

if ~isreal(L)
    keyboard
end


if exist('Nact','var')
    Npred_out = N;
    Nact_out = Nact;
else
    Npred_out = N;
    Nact_out = NaN;
end
%also output the last month's fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
