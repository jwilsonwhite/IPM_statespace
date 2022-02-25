function mcmc_iterator

% A script to iterate loops of the oyster model mcmc code

% Open matlabpools for parallel processing
if matlabpool('size') > 0
    matlabpool close force local
end
matlabpool open 3


% Model types to be passed to the mcmc code
% Mortality type
Mn = {'Con','Con','Con','Con','Region','Region','Region','Region',...
       'RegionYear','RegionYear','RegionYear','RegionYear'};
% Recruitment type   
Rn = {'Con','RegionYear','Region','Year','Con','RegionYear','Region','Year','Con','RegionYear','Region','Year',};
% Model number (see Excel spreadsheet with results summary)
Num = [1,3,5,6,13,15,17,18,25,27,29,30];

% The loop
parfor m =  1:length(Mn)
   run_mcmc_kalman(Mn{m},Rn{m},Num(m),'12Jan2013'); 
end


matlabpool close
