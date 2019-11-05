function postproc_pisco_rockfish_fit_2(Meta_savename,Chains)

% This function examines the  posterior of an MCMC run, and makes a number
% of plots. It has some options to make plots for either simulated data
% runs or Pt Lobos data runs for the White et al. Ecol. Appl. paper.

% 'Chains' is a vector indicating which chains that were run, if MCMC chains
% were run manually in parallel. It can also be used to indicate only those
% chains that should be in the posterior, if some did not converge. If unset, default to the number of chains
% specified in Meta.MCMC.chains

% What plots to make?
Plotchains = true;
Plotfinal = true;
saveplots = true;

% Choose the directory to be used, depending on what type of run is
% analyzed
Dir ='runs/';
%Dir = 'mockdata_fits/'; % for mockdata runs

% Load in metadata

load(strcat(Dir,Meta_savename))
load(strcat(Dir,Meta.data_savename),'D_str','Species_Names','Years')

plot_savename = Meta.fit_savename(1:end-4); % trim off '.mat'
timeseries_plotname = strcat(plot_savename,'.eps');
F_histo_plotname = strcat(plot_savename,'_F_posterior_histo.eps');

%%%%%%%%%%%%%%%%%%%%
% Need to load in some things to plot predicted vs. actual:
meshsize = Meta.IPM.meshsize;
meshmin = Meta.IPM.meshmin;
meshmax = Meta.IPM.meshmax;
x = linspace(meshmin,meshmax,meshsize);
meshdiff = diff(x(1:2));
edges = x - meshdiff/2;

dy=x(2)-x(1);
Sy=makeSimpVec(dy,meshsize);

%timescale to run model
Tdata = Meta.Tdata; % timespan of data

Tpre = Meta.Tpre; % time prior to data that model will run
if isnan(Tpre(1)) % if it is a post2007 run with no pre-period
    T = 1:length(Tdata(:));
else
T = 1:length([Tpre(:); Tdata(:)]);
end

%get data into a histogram
Site_Names = Meta.Sites;

N = IPM_histo(D_str.(Species_Names{1}),Years,Site_Names,edges);

% strip out the years we don't need
OKyears = false(length(Years),1);
for i = 1:length(Years)
    if any(Meta.Tdata == Years(i))
        OKyears(i) = true;
    end
end
N = N(:,:,OKyears);

Ogive_b = Meta.ogive;
if isnan(Ogive_b(1)) % if there is no ogive
    Ogive = ones(size(x));
else
Ogive = 1-normcdf(x,Ogive_b(1),Ogive_b(2)); % probability of observation in the kelp forest
end


%simulate recruitment from distribution based on data
%specify max size of recruits in data (YOY)
Rsize = Meta.recruits.Rsize; 

%get sum of recruits for each year per site. 
R1 = squeeze(nansum(N(x<Rsize,:,:),1));

%find mean for each year across sites
if size(R1,1)>1 && size(R1,2) >1 % if more than one site
R1m = nanmean(R1(:,:),1);
else
R1m = R1;    % just one site, so no need to take mean
end
%log transform time series
Rfact = min(R1m(R1m>0)); % addition factor to avoid log(0)
R1lnm = log(R1m+Rfact);
%find mean and std of the time series
mu1 = mean(R1lnm);
sig1 = std(R1lnm);

R = normpdf(x,Meta.recruits.meansize,Meta.recruits.sdsize)'; %Average size of recruit is 4, with sd = 1
                         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in the MCMC results

if exist('Chains','var') && ~isnan(Chains)
    chains = length(Chains);
   % loop over runs & combine mc_str structure files
   fit_savename_tmp = Meta.fit_savename(1:end-4);
   mc_str_tmp = struct();
   for c = 1:length(Chains) %previously: c = 1:length(Chains)
       load(strcat(fit_savename_tmp,'_',num2str(c),'.mat'))
       if c == 1 
       mc_str_tmp = mc_str;
       else
       mc_str_tmp(c) = mc_str;
       end
   end
   mc_str = mc_str_tmp;
else
load(strcat(Meta.fit_savename))
chains = length(mc_str);
end

%Dir = 'NA';
if strcmp(Dir,'mockdata_fits/') % if doing mockdata
    NT = ones(1,length(T)); 
else
NT = mc_str(1).NT;
NT = NT(:,1:length(Meta.Tdata));
end


if strcmp(Dir,'mockdata_fits/')
Burn = [2.5e3 2.5e3 2.5e3 2.5e3];
else
Burn = [7.5e3 7.5e3 7.5e3 7.5e3]; % half length of chain
end
Col = {'r','b','k','g'};

if Plotchains
    figure(1)
    clf
    for i = 1:size(mc_str(1).parm_vec,2) % columns are dimension 2:  looking at chain convergence for each parameter posterior.
    subplot(7,3,i) %number of subplots should cover the ncols in parm_vec of fit
    hold on
    for c = 1:chains
    plot(mc_str(c).parm_vec(Burn(c):end,i),Col{c})
    end
    end
    
    subplot(7,3,19) %
    hold on
    for c = 1:chains
        plot(mc_str(c).L(Burn(c):end),Col{c})
    end
end % end if Plotchains

% Plot final posteriors, and compile the posteriors into a single matrix,
% Pvec.
% Also, calculate convergence statistics to ensure proper convergence.
if Plotfinal
    figure(2)
    set(gcf,'units','cent','position',[5,5,9, 15])
    clf
   Pvec = []; Lvec = [];
   for c = 1:chains 
       Ltmp = mc_str(c).L(Burn(c):end);
       Lvec = [Lvec; Ltmp(:)];
       Pvec = [Pvec; mc_str(c).parm_vec(Burn(c):end,:)];
       vLwithin(c) = var(Ltmp);
   end
   meanP=mean(Pvec(1:end,:));
   
   Rhat = sqrt(var(Lvec)/mean(vLwithin));
   if Rhat > 1.1 % there might be a convergence problem
       warning('Rhat > 1.1. Check chains for convergence.')
   end
   
if strcmp(Meta_savename(6),'N') % Naples
PlotSites = 1; % which sites to plot?

% which years to plot
if ~isnan(Meta.Tpre) % pre2012 runs
        Years = [11 12 14];%years 2000, 2001, 2003
else % post2007 runs
        Years = 2:4;
end

elseif strcmp(Meta_savename(6),'m') % simulated data
    PlotSites = 1;
    Years = [10, 17, 18];

elseif strcmp(Meta_savename(6),'A')  % Andrew Molera
PlotSites = 1

% which years to plot
    if ~isnan(Meta.Tpre) % pre2007 runs
        Years = [15,16,17]; %years 2004, 2005, 2006
    else % post2007 runs
        Years = 4:6;
    end

else strcmp(Meta_savename(6),'V')
    PlotSites = 1;
% which years to plot
    if ~isnan(Meta.Tpre) % pre2007 runs
        Years = [14,15,16]; %years 2003, 2004, 2005
    else % post2007 runs
        Years = 4:6;
    end
        end

end % end switch over site
   
 Predicted = struct();

   for s = 1:length(Site_Names)

       % switch off fishing for pre-existing MPAs
        meanP_tmp = meanP;
        meanP_tmp(1:2) = meanP_tmp(1:2).*(1-Meta.MPAstatus(s)); 
        
   [L(s), Npred, Nact] = do_IPM(meshsize,x,dy,Sy,meanP_tmp,Meta.fixparm,Rfact,R,squeeze(N(:,s,:)),T,Tpre,Ogive,Meta,NT(s,:),Site_Names{s});

   Predicted(1).(Site_Names{s}).Npred = Npred;
   Predicted(1).(Site_Names{s}).Likelihood = L;


NumY = 3; % how many years to plot

Panel = 1:(NumY*length(PlotSites));
Panel = reshape(Panel,length(PlotSites),NumY)';
FS = 8;

if any(PlotSites==s) % if we are in one of the sites to be plotted
    
for y = 1:NumY
subplot(NumY,length(PlotSites),Panel(y,find(PlotSites==s)))
hold on
sum(Nact(:,Years(y)))
sum(Npred(:,Years(y))*diff(x(1:2)).*NT(s,Years(y)-length(Tpre)))

if strcmp(Dir,'mockdata_fits/') % if doing mockdata:
% data was binned by 3 cm, but treated as continuous
bar(x,Nact(:,Years(y))./(3/diff(x(1:2))),(3/diff(x(1:2))),'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])
%bar(x,Nact(:,Years(y))./(diff(x(1:2))),(diff(x(1:2))),'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])

%keyboard
else % real data
bar(x,Nact(:,Years(y))./diff(x(1:2)),(diff(x(1:2))),'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])
end

% Aggregate IPM density to make it comparable to the scale of the data histogram
if strcmp(Dir,'mockdata_fits/') % if doing mockdata:
   % Binned by 3 cm
   xt = 2.5:3:max(x); % these are the edges that will create bins centered on the mock data bins

 
   % Translation matrix:
   xm = repmat(x(:)',[length(xt),1]);
   em = repmat(xt(:),[1,length(x)]);
   em2 = [zeros(1,length(x));em(1:end-1,:)];
   Tm = xm > em2 & xm <= em;
   
   % Do the translation:
   Npred2 = Tm*Npred(:,Years(y)).*NT(s,Years(y)-length(Tpre));
   
   % Plot the translated density:
 %  plot(xt-diff(xt(1:2)/2),Npred2,'ko')
   
end
    
% Plot the raw density:
plot(x,Npred(:,Years(y)).*NT(s,Years(y)-length(Tpre)),'k','linewidth',1) 
set(gca,'tickdir','out','ticklength',[0.02 0.02],'fontsize',FS)
set(gca,'xlim',[0 70])

if strcmp(Dir,'mockdata_fits/') 
set(gca,'ylim',[0 20])
end

if y == 1; title(Site_Names{s},'fontsize',FS); end
xlabel('Length (cm)','fontsize', FS)
ylabel('Density','fontsize',FS)
end % end loop over years for plotting
if saveplots
print(timeseries_plotname,'-depsc2','-tiff')
end
end % end if PlotSites

   end % end loop over sites
end % end if Plotfinal 

% Find the mode:
edge = linspace(0,max(Pvec(:,2)),50); 
[n, ~] = histc(Pvec(:,2),edge); 
F_mode = edge(n==max(n))+diff(edge(1:2)); % take midpoint of biggest bin
F_mode = F_mode(1); % in case there is a tie
Post(1).F_mode = F_mode;
Post(1).F_mean = mean(Pvec(:,2));
Post(1).F_median = median(Pvec(:,2));
Post(1).F_std = std(Pvec(:,2));

% Plot the histogram of historical F
figure(10)
clf
set(gcf,'units','cent','position',[14 5 8.4 4])
hold on
hist(Pvec(:,2),50)
plot([F_mode,F_mode],[0 2e3],'k-')
plot([Meta.Fprior(2,1),Meta.Fprior(2,1)],[0 2e3],'k--')
xlabel('F')
ylabel('Frequency')
ylim([0 1100])
xlim([0 4.5])
set(gca,'tickdir','out','ticklength',[0.015 0.015])
if saveplots
print(F_histo_plotname,'-depsc2','-tiff')
end
% Save post-burn chains & posterior of everything; final distribution of Npred; save F posterior:
post_savename = Meta.fit_savename(1:end-4); % trim off '.mat'
post_savename = strcat(post_savename,'_postproc.mat');
Meta.post_savename = post_savename;
save(strcat(Dir,Meta_savename),'-append','Meta')

Post(1).posterior = Pvec; % the posterior distribution of all the parameters
Post(1).predicted = Predicted; % The predicted distribution at all sites and times
Post(1).Rfact = Rfact;

save(strcat(Dir,post_savename),'Post')
