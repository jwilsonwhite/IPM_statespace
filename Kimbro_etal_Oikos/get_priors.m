function P = get_priors(Model_str,Meta)

% Define priors (and # parameters) for different model formulations.
% (updated 12/31/15)

% Deal with all abundances in log units for convenience

% Average R from dataset:
Rmean = 0.02;
Rhyper = 2; % larger hyperparameter stdev for recruitment
Roffset = 10; % 10 km - right in the middle of the bay for quadratic function

Amean = 0.2; % adult mortality estimate
Jmean = 0.2; % juvenile mortality estimate
Mhyper = 1; % mortality hyperparameter std
Moffset = 10; % same as Roffset

T = length(Meta.T);

switch Model_str.rect_time  % make these gamma....
    case 'Con' % only spatial variation in recruitment
        switch Model_str.rect_space
            case 'Con'
                 R_prior = log(Rmean);
                 R_hyper = Rhyper;
            case 'Grad'
                 R_prior = [log(Rmean), 0]; % intercept, slope
                 R_hyper = repmat(Rhyper,[1,2]);
            case 'Unim'
                 R_prior = [log(Rmean), 0, Roffset]; % intercept, slope, offset 
                 R_hyper = repmat(Rhyper,[1,3]);
        end % end switch over rect space
             
    case 'Var'
        switch Model_str.rect_space
            case 'Con'
                 R_prior = repmat(log(Rmean),[1,T]);
                 R_hyper = repmat(Rhyper,[1,T]);
            case 'Grad'
                 R_prior = [repmat(log(Rmean),[1,T]), 0]; % intercept, slope
                 R_hyper = repmat(Rhyper,[1,T+1]);
            case 'Unim'
                 R_prior = [repmat(log(Rmean),[1,T]), 0, Roffset]; % intercept, slope, offset
                 R_hyper = repmat(Rhyper,[1,T+2]);
        end % end switch over rect space

end

% Adult mortality
switch Model_str.aM_time
    case 'Con' % only spatial variation in mortality
        switch Model_str.aM_space
            case 'Con'
                 aM_prior = log(Amean);
                 aM_hyper = Mhyper;
            case 'Grad'
                 aM_prior = [log(Amean), 0]; % intercept, slope
                 aM_hyper = repmat(Mhyper,[1,2]);
            case 'Unim'
                 aM_prior = [log(Amean), 0, Moffset]; % intercept, slope, offset 
                 aM_hyper = repmat(Mhyper,[1,3]);
        end % end switch over rect space
             
    case 'Var'
        switch Model_str.aM_space
            case 'Con'
                 aM_prior = repmat(log(Amean),[1,T-1]);
                 aM_hyper = repmat(Mhyper,[1,T-1]); % T-1 because only need T-1 mortality rates
            case 'Grad'
                 aM_prior = [repmat(log(Amean),[1,T-1]), 0]; % intercept, slope
                 aM_hyper = repmat(Mhyper,[1,T]);
            case 'Unim'
                 aM_prior = [repmat(log(Amean),[1,T-1]), 0, Moffset]; % intercept, slope
                 aM_hyper = repmat(Mhyper,[1,T+1]);
        end % end switch over rect space

end    

% Juvenile mortality
switch Model_str.jM_time
    case 'Con' % only spatial variation in mortality
        switch Model_str.jM_space
            case 'Con'
                 jM_prior = log(Jmean);
                 jM_hyper = Mhyper;
            case 'Grad'
                 jM_prior = [log(Jmean), 0]; % intercept, slope
                 jM_hyper = repmat(Mhyper,[1,2]);
            case 'Unim'
                 jM_prior = [log(Jmean), 0, Moffset]; % intercept, slope, offset 
                 jM_hyper = repmat(Mhyper,[1,3]);
        end % end switch over rect space
             
    case 'Var'
        switch Model_str.jM_space
            case 'Con'
                 jM_prior = repmat(log(Jmean),[1,T-1]);
                 jM_hyper = repmat(Mhyper,[1,T-1]);
            case 'Grad'
                 jM_prior = [repmat(log(Jmean),[1,T-1]), 0]; % intercept, slope
                 jM_hyper = repmat(Mhyper,[1,T]);
            case 'Unim'
                 jM_prior = [repmat(log(Jmean),[1,T-1]), 0, Moffset]; % intercept, slope, offset
                 jM_hyper = repmat(Mhyper,[1,T+1]);
        end % end switch over rect space

end    

% Place into structure for export
P(1).R_prior = R_prior;
P(1).R_hyper = R_hyper;
P(1).aM_prior = aM_prior;
P(1).aM_hyper = aM_hyper;
P(1).jM_prior = jM_prior;
P(1).jM_hyper = jM_hyper;

P(1).pe_prior = 1e-1;
P(1).pe_hyper = 2;
P(1).oe_prior = 1e-1;
P(1).oe_hyper = 2;


