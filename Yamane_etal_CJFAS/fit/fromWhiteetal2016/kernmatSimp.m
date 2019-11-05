function [kmat] = kernmatSimp(x,F,fixparm,T)

% Set up the integration mesh kernel using Simpsons
% For MLPA monitoring model

%adapted from Easterling et al. (2001). Evenly spaced grid, now add weights so
%use Simpson's rule using Marissa's code.

y = x;
%this creates a vector (y) that is equal to x

[x,y] = meshgrid(x,y); % Matlab built in function
%x is an array (original x by original x) with each row equal to original
%vector x
%y is an array (original y by original y) with each column equal to
%original vector y
%X corresponds to size at time t
%Y corresponds to size at time t+1

% Make the kernel from the grid
kmat = mkkern(x,y,F,fixparm,T);

kmat = max(0,kmat); 


