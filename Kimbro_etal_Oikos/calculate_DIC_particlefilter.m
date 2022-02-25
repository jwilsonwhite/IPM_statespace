function [DIC, Pd, L_corrected] = calculate_DIC_particlefilter(L,R,aM,jM,pe,Rvec,isJuv,data,Meta,Model_str)
                

% Calculate Deviance Information Criterion for Tomales oyster model
% This version uses the particle filter IPM model (2015)

% Calculates both original DIC and DIC3 variant
% See Celeux et al. (2006) Bayesian Analysis 1:651-674

% Note: if raw likelihoods are quite small, then DIC3 will return Inf
% (because log(0) is -Inf).  To avoid this, consider adding some constant
% value to the likelihoods (not log-likelihoods!).  But in general the
% original should work fine for many applications

% Note: be sure that B and L have had burnin removed

% DevianceMeans is deviance calculated from posterior mean parameter values
    Par(1).R = mean(R);
    Par(1).aM = mean(aM);
    Par(1).jM = mean(jM);
 %   Par(1).oe = mean(oe);
    Par(1).pe = mean(pe);
    
    % First calculate muhat, the mean likelihood across the chains.
    % In the original IPM code, likelihood is combined with prior prior to
    % storing (bc we stored the value used in the Metropolis rule). So now
    % subtract that effect:
    Prior = get_priors(Model_str,Meta);
    n = size(R,1);
    PriorL(:,1) = (sum(log(max(realmin,normpdf(R,repmat(Prior.R_prior,[n,1]),repmat(Prior.R_hyper,[n,1])))),2));
    PriorL(:,2) = (sum(log(max(realmin,normpdf(aM,repmat(Prior.aM_prior,[n,1]),repmat(Prior.aM_hyper,[n,1])))),2));
    PriorL(:,3) = (sum(log(max(realmin,normpdf(jM,repmat(Prior.jM_prior,[n,1]),repmat(Prior.jM_hyper,[n,1])))),2));
    PriorL(:,4) = (sum(log(max(realmin,gampdf(1./pe,Prior.pe_hyper,1./Prior.pe_prior./(Prior.pe_hyper-1)))),2));
  %  PriorL(:,5) = (sum(log(max(realmin,gampdf(1./oe.^2,Prior.oe_hyper,1./Prior.oe_prior.^2./(Prior.oe_hyper-1)))),2));
    
    PriorLs = sum(PriorL,2);
  %  L_corrected = (L(:)-PriorLs); 
  L_corrected = L(:);  % !!!! nevermind, you are supposed to use the posterior after all
    
% muhat is the mean deviance during the Markov Chain
% deviance = -2.*log likelihood
% This is based on Gelman et al. (2004) Bayesian Data Analysis. This form
% seems to work better than the original muhat - D(mean(Theta))
%muhat = -2.*nanmean(L_corrected);
Pd = 0.5*var(L_corrected);


%if isnan(DIC)
%    keyboard
%end

do_mu3hat = false;
if do_mu3hat
% mu3hat is slightly different: 
% take log of mean likelihood instead of mean log likelihood
% This will have numerical problems if the likelihood is large
%mu3hat = -2.*log(nanmean(exp(L_corrected)));
    
    % Get growth parameters:
    load('Tomales_vB_fits.mat','X')
    switch Model_str.growth
    case 'Var'
        Par(1).Growth = X;
    case 'Con'
        Par(1).Growth = repmat(mean(X),[7,1]); % just repeat the mean over several sites
    end % end switch growth
    
    % Some things needed for integration:
    x = Meta.IPM.x;
    dy=x(2)-x(1);
    Sy=makeSimpVec(dy,Meta.IPM.meshsize);

    % calculate likelihood
    L_tmp = nan(1,7);
    for s = 1:7
    L_tmp(s) = do_Tomales_IPM(dy,Sy,Par,Rvec(:,s),isJuv(:,s),squeeze(data(:,s,:)),Meta,Model_str,s); 
    end
    L_mean = sum((L_tmp)); % log-likelihood

end
%DevianceMeans = -2*L_mean;
  muhat = -2*mean(L_corrected);  

  
DIC = muhat+Pd;

%Pd = muhat-DevianceMeans;
%Pd3 = muhat - mu3hat;


%DIC = muhat + Pd;
%DIC3 = muhat + Pd3;
%end % end if do_mu3hat
