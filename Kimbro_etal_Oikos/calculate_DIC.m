function [DIC, DIC3, Pd, Pd3] = calculate_DIC(B,L,j_ind, a_ind,r_ind,pj,y)

% Calculate Deviance Information Criterion for Tomales oyster model


% Calculates both original DIC and DIC3 variant
% See Celeux et al. (2006) Bayesian Analysis 1:651-674

% Note: if raw likelihoods are quite small, then DIC3 will return Inf
% (because log(0) is -Inf).  To avoid this, consider adding some constant
% value to the likelihoods (not log-likelihoods!).  But in general the
% original should work fine for many applications

% Note: be sure that B and L have had burnin removed

% muhat is the mean deviance during the Markov Chain
% deviance = -2.*log likelihood
muhat = -2.*nanmean(L);

% mu3hat is slightly different: 
% take log of mean likelihood instead of mean log likelihood
mu3hat = -2.*log(mean(exp(L)));

% DevianceMeans is deviance calculated from posterior mean parameter values
B_mean = mean(B);

% calculate likelihood
    m_j = B_mean(j_ind);
    m_a = B_mean(a_ind);
    R = B_mean(r_ind);
    R = reshape(R,pj,y); % make R into a p-by-y matrix
    J = reshape(m_j,pj,y); % p-by-y matrix
    A = reshape(m_a,pj,y); % p-by-y matrix

    L_mean = popmodel_v2(J,A,R); % this has to be a log likelihood!
    
DevianceMeans = -2*L_mean;
    
Pd = muhat-DevianceMeans;
Pd3 = muhat - mu3hat;


DIC = muhat + Pd;
DIC3 = muhat + Pd3;

