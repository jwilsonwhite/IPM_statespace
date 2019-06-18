function [fp, tp] = ROC(X1,X2)

% create ROC curve for two distributions
% Input: distribution of values X1 & X2. 
% Hypothesis: X2 > X1

% Output: 
% fp: Probability of false positive
% tp: Probability of true positive

% To obtain ROC curve:
% plot(fp,tp)

Z = 1e3;
z = linspace(min(X1),max(X2),Z); % possible cutoffs
z1 = repmat(z,[length(X1),1]);
z2 = repmat(z,[length(X2),1]);
X1a = repmat(X1(:),[1,Z]);
X2a = repmat(X2(:),[1,Z]);

fp = mean(X1a>z1); % false positive
tp = mean(X2a>z2); % true positive



