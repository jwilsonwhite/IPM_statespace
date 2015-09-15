function kxy = mkkern(x,y,F,fixparm,T)

% This function creates the IPM kernel.
% Adapted from Easterling et al. (2001) for use with the central coast
% rockfish model


%Define which sizes can be fished:
isjuv = 1 - normcdf(x,fixparm(5),fixparm(7)); 

%define pr(reproductive) using a maturity ogive
%not needed now because open population
%ismat = 1 - normcdf(x,fixparm(5),diff(x(1,1:2))/2);

%SURVIVAL PART OF KERNEL
M = fixparm(4); % natural mortality rate
m = ones(size(x)).*M + (1-isjuv).*F; %this is a matrix size x, mortality for each size 
                                              
p1 = exp(-m*T); % convert mortality rate to survivorship, iterate over time steps

%GROWTH PART OF KERNEL
Linf = fixparm(1);
k = fixparm(2);
pmean1=Linf - (Linf - x).*exp(-k); % (do not add in x0 for the one-step growth)
%add variability around von Bertalanffy growth
psig1 = pmean1*fixparm(7); 

%evaluate growth part of kernel
p2 = normpdf(y, pmean1, psig1);

p1 = max(0,p1);    %to make sure no negatives
p2 = max(0,p2);    %to make sure no negatives

kxy = p1.*p2;
%add fecundity to kernel if a closed population

