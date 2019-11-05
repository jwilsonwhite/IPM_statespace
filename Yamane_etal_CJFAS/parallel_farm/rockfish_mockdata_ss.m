function rockfish_mockdata

% Generate test dataset for MCMC IPM fitting

% Load in some metadata to extract parameters for blue rockfish
load('SMYS_Pt_Lobos_pre2007_31Aug2015_metadata.mat')

% Structure: 
% 10 replicates each for F = 0, 0.05, 0.1
% 10 replicates for F1 = 0.1, F2 = 0.0 (MPA scenario)
% 10 replicates for F = 0.05, but with only 7 or 3 years for fitting
% 10 replicates for F = 0.05, but with wider length bins
Nsim = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% values to change for each set of runs
F1 = 0.1; % F during burn-in
F2 = 0.1; % F during data collection

% Time structure of the data:
% baseline: 9 years for burnin, 9 years for fit
T1 = 9;
T2 = 9; % other values for time series duration = 3, 6, 12, 15 yrs
T = T1 + T2;

% length bin width
BinW = 3;

savename = 'SMYS_mockdata_F0.1ss200_12July2018';
savedir = 'mockdata_March2018';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean & lognormal error for recruit timeseries
Rmean = 50;
Rsd = 0.1;

% age structure for simulation
Maxage = 44;
Ages = 0:Maxage;

% age at first capture (defined as length, but easier to model as age)
% extract t0 from x0:
% t0 = (log(1 - x0/Linf))/k
t0 = (log(1 - Meta.fixparm(3)/Meta.fixparm(1)))/Meta.fixparm(2);
% La = Linf*(1 - exp(-k*(t - t0))
La = Meta.fixparm(1) .* (1 - exp( -Meta.fixparm(2) .* (Ages - t0)));
tc = min(Ages(La >= Meta.fixparm(5))); % age at first capture

Species_Names = {'SMYS'};
Site_Names = {'NULL'};
rng(2) % set seed here if want R1,R2,R3,..,RN to be different populations (in case of sample size analyses)

% Loop to create replicate datasets
for i = 1:Nsim
% Randomly drawn recruitment
R = exp( normrnd( log(Rmean), Rsd, T, 1) );

% Natural mortality rate + process error
M = Meta.fixparm(4);
sig_p = 1e-2; % Process error in mortality
Mt = normrnd(0,sig_p,T,1);
Mactual = mean(M+Mt);

% Survival vs age:
Fvec = [ones(1,T1)*F1, ones(1,T2)*F2];
Fmat = repmat(Fvec,[Maxage,1]).*repmat(Ages(1:end-1)'>=tc,[1,T]);
Mmat = ones(Maxage,T).*M;
Zmat = Fmat+Mmat;

N = nan(Maxage+1,T);
Ninit = 2*Rmean.*[1; cumprod(exp(-Mmat(:,1)))];
N(:,1) = Ninit;

% Now loop over time
for t = 2:T;
    Ztmp = max(0, Zmat(:,t) + Mt(t)); % mortality + process error
    Ztmp = diag(exp(-Ztmp));
    A = [zeros(1,Maxage); Ztmp];
    A = [A, zeros(Maxage+1,1)];
    N(:,t) = A*N(:,t-1);
    N(1,t) = R(t);
end

% convert to size distributions
Lvar = Meta.fixparm(7);
stdL = La.*Lvar;

len_vec = 1:BinW:round(Meta.fixparm(1)*1.5); % vector of length bins
len_mat = repmat(len_vec(:),[1 Maxage+1]);
len_mat2 = repmat(La(:)',[length(len_vec),1]);
len_mat3 = len_mat - len_mat2;
std_mat = repmat(stdL(:)',[length(len_vec),1]);

LM = normcdf(len_mat3+BinW/2,0,std_mat) - normcdf(len_mat3-BinW/2,0,std_mat); % prob of being in each bin
LM(isnan(LM)) = 0;

% convert to length distribution
NL = LM*N;

% Now expand out to cm-resolution to simulate the binning artifact in real
% data
LL = 1:1:max(len_vec);
Ltrans = nan(length(len_vec),1);
for j = 1:length(len_vec);
    Ltrans(j) = find(len_vec(j)==LL);
end
    
NL2 = zeros(length(LL),T);
for t = 1:T
    for j = 1:length(len_vec)
    NL2(Ltrans(j),t) = NL(j,t);
    end
end

% sample from Poisson distribution
NS = poissrnd(NL2);

% Info on the simulated dataset parameters
RF(1).F1 = F1;
RF(1).F2 = F2;
RF(1).M = M;
RF(1).Mactual = Mactual;
RF(1).Ractual = R;
RF(1).tc = tc;
RF(1).BinW = BinW;
RF(1).T1 = T1;
RF(1).T2 = T2;
RF(1).N = N;
RF(1).NL = NL;
RF(1).NL2 = NL2;
RF(1).NS = NS;

% Now make a data structure like the PISCO format so it will play nice with
% the code:
Years = 1999:(1999+T2-1);
for y = 1:T2
    Count = NS(:,y+T1);
    TL = LL;
    OK = Count>0;
    Count = Count(OK);
    TL = TL(OK);
D_str.(Species_Names{1}).(Site_Names{1})(y).data.Count = Count;
D_str.(Species_Names{1}).(Site_Names{1})(y).data.TL = TL;
D_str.(Species_Names{1}).(Site_Names{1})(y).data.numtrans = 1;
D_str.(Species_Names{1}).(Site_Names{1})(y).year = Years(y);
end

sname2 = strcat(savedir,'/',savename,'_R',num2str(i),'.mat');
save(sname2,'RF','D_str','Species_Names','Site_Names','Years');

end % end loop over replicate datasets





