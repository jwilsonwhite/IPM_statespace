function [] = runme_IPM_mockdata_PCLA_F15(varargin)

if(nargin>0)
scenario = varargin{1};
else
scenario = 'F0.15';
end
disp(scenario)

data_savename = strcat('mockdata_March2018/PCLA_mockdata_',scenario,'_8Aug2018'); % where the mock data are

savename = strcat('PCLA_mockdata_',scenario,'_8Aug2018'); % where the fit will be stored

% Dole out mock data to workers
poolobj = gcp;
addAttachedFiles(poolobj,{'PCLA_Pt_Lobos_pre2007_31Aug2015_metadata.mat'})
for j = 1:10
addAttachedFiles(poolobj,{strcat(data_savename,'_R',num2str(j),'.mat')})
end

parfor i = 1:10 % loop over reps
fitnum = i;

runme_extra_PCLA(fitnum, data_savename, savename)

%disp('it made it thru rockfish_fit')
end
