function [] = runme_IPM_mockdata_samplesize_SMYS_F02ss100_12bJuly(varargin) %(scenario,samplesize)
% Run CC IPM model on mock datasets for validation

% Forked from runme_IPM...formats the mockdata to interact with the
% rockfish_fit_pisco_mock_ss code

% March 2017 - forked from runme_IPM_mockdata to include the option of
% specifying particular sample sizes

% Load in the metadata used to create the mockdata
if(nargin>0)
scenario = varargin{1};
else
scenario = 'F0.2ss100';
end

data_savename = strcat('mockdata_March2018/SMYS_mockdata_',scenario,'_12July2018'); % where the mock data are
savename = strcat('SMYS_mockdata_',scenario,'_12July2018'); % where the fit will be stored

% Dole out the mockdata to the workers
poolobj = gcp;
addAttachedFiles(poolobj,{'SMYS_Pt_Lobos_pre2007_31Aug2015_metadata.mat'})
for j = 51:100
addAttachedFiles(poolobj,{strcat(data_savename,'_R',num2str(j),'.mat')})
end

parfor i = 1:50 % loop over reps
fitnum = i + 50;

runme_extra_SMYS_ss(fitnum, data_savename, savename, 100)

%disp('it made it thru rockfish_fit')

end % end loop over reps
