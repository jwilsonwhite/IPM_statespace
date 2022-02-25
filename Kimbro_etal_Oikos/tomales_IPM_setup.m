function tomales_IPM_setup

% Create a bunch of useful variables for the IPM and other model functions
% Save them in tomales_IPM_setup.mat

% Read in some data
load Size_mat % growth trajectories
%%% NEED TO BUILD THIS INTO IPM KERNELS

load Sizeclass_actual.mat Sizeclass_actual Density_actual Density_SD_actual  Density_SD_log_actual
% this is the observed size distribution & total density of each population (survey data)

% Obsyears: vector of which of list of simulated years actually correspond to observed data
Obsyears = [5 6 7 10]; % (assume that the first X years are going to be ignored)
Nyears = Obsyears(end);  %the total number of years to be simulated (the final years have data)
Age_vec = 0:10; % a vector of 0, 1, 2, ... A, 0 included to account for YOY age class
%Age_vec also equals the number of annual recruitment events

% To do the bookkeeping for separate juvenile & adult mortality rates, it
% will be useful to have two additional age vectors.
% This assumes that juveniles turn into adults between age 0 & age 1.
aj = 1; % the age after which the juvenile-to-adult transition occurs

Age_vec_j = Age_vec; % 'juvenile' age vector
Age_vec_j(Age_vec_j > aj) = aj;  % so all the values after age 1 become 1

% Logical vector of mature/not mature
mat_vec = Age_vec > aj;

Age_vec_a = Age_vec; % 'adult' age vector
Age_vec_a(Age_vec_a <= aj) = 0; % so all values less than age 4 go to zero
Age_vec_a(Age_vec_a > aj) = Age_vec_a(Age_vec_a > aj) - aj; % so this vector now counts age past the juveline stage

Age_mat_j = repmat(Age_vec_j,[Nyears,1]);
Age_mat_a = repmat(Age_vec_a,[Nyears,1]);

Nages = length(Age_vec);


% Combine info on size & density
for p = 1:size(Density_actual,2)
    for t = 1:size(Density_actual,1)
        X_a = squeeze(Sizeclass_actual(t,:,p));
        X_a = X_a./sum(X_a);
        X(t,:,p) = X_a.*Density_actual(t,p);
    end
end

T_obs = zeros(10,1);
for o = 1:length(Obsyears)
T_obs(Obsyears(o)) = o;
end

% generate uniform size trajectories too
Size_mat_uni = repmat(mean(Size_mat,3),[1,1,size(Size_mat,3)]);


save tomales_kalman_setup.mat