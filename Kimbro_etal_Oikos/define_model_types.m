function Ms = define_model_types
% Define which models will be run for Tomales simulations

%%% Read in from IPM_model_types.csv
M = importdata('IPM_model_types.csv');

Ms = struct([]);

for i = 1:length(M.data)
    
    Ms(i).rect_time = M.textdata{i+1,1};
    Ms(i).rect_space = M.textdata{i+1,2};
    Ms(i).aM_time = M.textdata{i+1,3};
    Ms(i).aM_space = M.textdata{i+1,4};
    Ms(i).jM_time = M.textdata{i+1,5};
    Ms(i).jM_space = M.textdata{i+1,6};
    Ms(i).growth = M.textdata{i+1,7};
    Ms(i).model_num = M.data(i);
    
end % loop over i






