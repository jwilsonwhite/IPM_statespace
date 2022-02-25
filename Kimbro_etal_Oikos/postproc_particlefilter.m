function postproc_particlefilter(modelnums)

% Post-process Tomales oyster model runs

% modelnums is vector listing which models to compare

dofig1 = 0;
dohist = 0;
dosave = 1;

DICtab = nan(length(modelnums),3);

for m = 1:length(modelnums)
    modelnum = modelnums(m);


filename = strcat(['Tomales_IPM_fits/Tomales_Model_',num2str(modelnum),'_20Jan2018_fit.mat']);
metaname = 'Tomales_20Jan2018_metadata.mat';
savename = strcat(['Tomales_IPM_fits/postproc_Tomales_Model_',num2str(modelnum),'_20Jan2018.mat']);

load(filename)
load(metaname)

M = length(mc_str(1).L); % length of chain


% convergence diagnostics (based on Gelman & Shirley 2011)

% exceptions for nonconvergence (rerun these?)
%188,

% chains that need more burn-in:
if any(modelnum == [8,10,16,17,44,64,70,83,88,100,209,118,180,170,208])
    burnin = round(M*0.75);
elseif any(modelnum == [34,52,65,66,82,102,190,172,174,210])
    burnin = round(M*0.9);
else
    burnin = round(M/2); % discard first half 
end

if any(modelnum == [10,12,13,16,18,33,66,68,69,90,125,161,195,192,156,192,174,210])
    OKchains = 1;
elseif any(modelnum == [2,15,30,34,36,59,62,71,72,84,88,96,102,108,121,132,139,143,159,...
                        175,177,203,213,215,168,208,114,126,138])
    OKchains = 2;
elseif any(modelnum == [26,48,54,64,65,70,82,83,123,179,193,197,211,188,180,206])
    OKchains = 3;
elseif any(modelnum == [14,24,41,42,47,131,167,154,156,140,190])
    OKchains = 1:2;
elseif any(modelnum == [6,23,78,103,157,204,214,168,158,170,178,196,198])
    OKchains = 2:3;
elseif any(modelnum == [5,31,49,60,79,85,100,104,106,113,141,201,120,212,176,172,150,162])
    OKchains = [1,3];
else
    OKchains = 1:3;
end

L = [];
vLwithin = nan(2,1);
for c = OKchains
    L = [L(:); mc_str(c).L(burnin:end-1)']; % due to a minor bug, the likehood has one additional entry
    vLwithin(c) = var(mc_str(c).L(burnin:end-1));
end % end loop over chains


Rhat = sqrt(var(L)/nanmean(vLwithin));

if Rhat > 1.1 % if there might be a convergence problem
    figure
    plot(L)
    modelnum
    keyboard
end


if dofig1
  
    figure
    hold on
    colors = {'k','b','r','g','c'};
    nchains = 3;
for n = 1:nchains
    plot(mc_str(n).L(burnin:end),colors{n});
end
end

% Summarize the data    
    R = []; aM = []; jM = []; pe = [];
    
    for c = OKchains
        R = [R; mc_str(c).R(burnin:end,:)];
        aM = [aM; mc_str(c).aM(burnin:end,:)];
        jM = [jM; mc_str(c).jM(burnin:end,:)];
      %  oe = [oe; mc_str(c).oe(burnin:end,:)];
        pe = [pe; mc_str(c).pe(burnin:end,:)];
    end
    
   
    % Make histograms of posteriors (does not currently work)
    if dohist
        
    end
   
        [DIC, Pd, Likelihood] = calculate_DIC_particlefilter(L,R,aM,jM,pe,...
                mc_str(1).R_vec,mc_str(1).isJuv,mc_str(1).N,...
                Meta,mc_str(1).Model_str);
    
       
       
   
        post_str(1).L = L;
        post_str(1).R = R;
        post_str(1).aM = aM;
        post_str(1).jM = jM;
      %  post_str(1).oe = oe;
        post_str(1).pe = pe;
        post_str(1).DIC = DIC;
        post_str(1).Pd = Pd;
      %  post_str(1).Pd3 = Pd3;
      %  post_str(1).DIC3 = DIC3;
        
        save(savename,'post_str')
        
        
        DICtab(m,1) = modelnum;
        DICtab(m,2:5) = [mean(L), mean(Likelihood), DIC, Pd ];

end % end loop over modelnum

DICtab(:,6) = DICtab(:,4)-min(DICtab(:,4));
if dosave
save DIC_results_1Feb2018.mat DICtab
end
%keyboard
    
    