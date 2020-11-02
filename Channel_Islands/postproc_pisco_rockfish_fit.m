function postproc_pisco_rockfish_fit(Meta_savename,Chains)

% This function examines the  posterior of an MCMC run, and makes a number
% of plots. It has some options to make plots for either simulated data
% runs or Pt Lobos data runs for the White et al. Ecol. Appl. paper.

% 'Chains' is a vector indicating which chains that were run, if MCMC chains
% were run manually in parallel. It can also be used to indicate only those
% chains that should be in the posterior, if some did not converge. If unset, default to the number of chains
% specified in Meta.MCMC.chains

% What plots to make?
Plotchains = false;
Plotfinal = true;
saveplots = true;
Plothisto = false;

PPn = 1e2; % number of draws from posterior predictive

% Choose the directory to be used, depending on what type of run is
% analyzed
%Dir ='runs_for_publication_pre2007_July2015/';
%Dir = 'mockdata_fits_June2015/'; % for mockdata runs
%Dir = 'post2007_runs_July2015/'; % for post-2007 runs
Dir = '';

% Load in metadata
load(Meta_savename,'Meta')
load(Meta.data_savename,'D_str','Species_Names','Years','Site_Names')

plot_savename = Meta.fit_savename(1:end-4); % trim off '.mat'
timeseries_plotname = strcat('Timeseries_plots_Sept2020/',plot_savename,'.eps');
ppcheck_plotname = strcat('Posteriorpredictive_plots_Sept2020/',plot_savename,'.eps');
F_histo_plotname = strcat(plot_savename,'_F_posterior_histo.eps');
mcmc_plotname = strcat(plot_savename,'_MCMC.eps');

% Get a human-friendlier name
Spname_human = human_name(Species_Names{1});
Sitename_human = human_name(Site_Names{1});
plot_title = strcat([Spname_human,' : ',Sitename_human]);



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

if isempty(min(R1m(R1m>0))) % if no recruits any year
    R1m(1) = 1;
end
%log transform time series
Rfact = min(R1m(R1m>0)); % addition factor to avoid log(0)
R1lnm = log(R1m+Rfact);
%find mean and std of the time series
%mu1 = mean(R1lnm);
%sig1 = std(R1lnm);
    mu1 = 3;
    sig1 = 3;
    Rfact = 0;

R = normpdf(x,Meta.recruits.meansize,Meta.recruits.sdsize)'; %Average size of recruit is 4, with sd = 1
                         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in the MCMC results

if exist('Chains','var') && ~isnan(Chains)
    chains = length(Chains);
   % loop over runs & combine mc_str structure files
   fit_savename_tmp = Meta.fit_savename(1:end-4);
   mc_str_tmp = struct();
   for c = 1:length(Chains)
       load(strcat(Dir,fit_savename_tmp,'_',num2str(c),'.mat'),'mc_str')
       if c == 1 
       mc_str_tmp = mc_str;
       else
       mc_str_tmp(c) = mc_str;
       end
   end
   mc_str = mc_str_tmp;
else
    load(strcat(Dir,Meta.fit_savename),'mc_str')
    chains = length(mc_str);
end

if strcmp(Dir,'mockdata_fits_June2015/') % if doing mockdata
    NT = ones(1,length(T)); 
else
    NT = mc_str(1).NT;
    NT = NT(:,1:length(Meta.Tdata));
end

if strcmp(Dir,'mockdata_fits_June2015/')
    Burn = [2.5e3 2.5e3 2.5e3];
else
    
    % THIS IS WHERE YOU SET THE BURNIN LENGTH. 
    % The vector 'Burn' should have the same length as the number of
    % chains, and each entry should be the desired amount of burn-in, based
    % on the diagnostic plots and Rhat calculation.
    % A single value worked for most sites, a few required some adjustment
    if strcmp(Species_Names{1},'SMYS') && strcmp(Site_Names{1},'SCI_HAZARDS')
        Burn = [1e3 1e3 1e3];
    elseif strcmp(Species_Names{1},'SPUL') && strcmp(Site_Names{1},'ANACAPA_WEST_ISLE')
    Burn = [1e4+1 1e4+1 1e3];
    elseif strcmp(Species_Names{1},'SPUL') && strcmp(Site_Names{1},'SRI_JOHNSONS_LEE_SOUTH')
    Burn = [7e3 7e3 1e4+1];
    elseif strcmp(Species_Names{1},'PCLA') && strcmp(Site_Names{1},'ANACAPA_WEST_ISLE')
    Burn = [5e3 5e3 1e4+1];
    else
    Burn = [5e3 5e3 5e3];
    end
end

Col = {'r','b','k'};

if Plotchains
    figure(1)
    clf
    for i = 1:size(mc_str(1).parm_vec,2)
        subplot(8,3,i)
        hold on
        for c = 1:chains
            plot(mc_str(c).parm_vec(Burn(c):end,i),Col{c})
        end
    end
    
    subplot(8,3,24)
    hold on
    for c = 1:chains
        plot(mc_str(c).L(Burn(c):end),Col{c})
    end
    
    if saveplots
        print(mcmc_plotname,'-depsc2','-tiff');
    end
end % end if Plotchains

% Plot final posteriors, and compile the posteriors into a single matrix,
% Pvec.
% Also, calculate convergence statistics to ensure proper convergence.
if Plotfinal
    figure(2)
    set(gcf,'units','cent','position',[5,5,18,18])
    clf
    Pvec = []; Lvec = [];
    for c = 1:chains 
        Ltmp = mc_str(c).L(Burn(c):end);
        Lvec = [Lvec; Ltmp(:)];
        Pvec = [Pvec; mc_str(c).parm_vec(Burn(c):end,:)];
        vLwithin(c) = var(Ltmp);
   end
   meanP=mean(Pvec(1:end,:));
   
   Rhat = sqrt(var(Lvec)/mean(vLwithin))
   if Rhat > 1.1 % there might be a convergence problem
       warning('Rhat > 1.1. Check chains for convergence.')
   end
   
   % which years to plot
   Years = 13:27; % 2003-2017
   
   % Which sites to plot? (this is a legacy of the central coast version in
   % which there were multiple survey sites in a single run)
   PlotSites = 1;
   
   Predicted = struct();
 
   for s = 1:length(Site_Names)
       
       % switch off fishing for pre-existing MPAs
       meanP_tmp = meanP;
       meanP_tmp(1:2) = meanP_tmp(1:2).*(1-Meta.MPAstatus(s)); 
        
       % Calculate posterior predictive distribution
       %%%%%%%
       %%% NEEDS TO BE RUN WITHOUT DATA-FITTING
       %%% THEN COMPARE DISTRIBUTIONS TO THE DATA FITTING
       Pdist = Pvec(randsample(size(Pvec,1),PPn,'true'),:);
       Pdist = Pdist .* (1-Meta.MPAstatus(s));
       
       for i = 1:size(Pdist,1)
   [~, Npreddist(:,:,i)] = do_IPM(meshsize,x,dy,Sy,Pdist(i,:),Meta.fixparm,Rfact,R,squeeze(N(:,s,:)),T,Tpre,Ogive,Meta,NT(s,:),Site_Names{s},0,1);
       end
          [L(s), Npred, Nact] = do_IPM(meshsize,x,dy,Sy,meanP_tmp,Meta.fixparm,Rfact,R,squeeze(N(:,s,:)),T,Tpre,Ogive,Meta,NT(s,:),Site_Names{s});

   Predicted(1).(Site_Names{s}).Npred = Npred;
   Predicted(1).(Site_Names{s}).Npreddist = Npreddist;
   Predicted(1).(Site_Names{s}).Likelihood = L;
   
   NumY = length(Years); % how many years to plot
   
   Panel = 1:(NumY + mod(NumY,2));
   Panel = reshape(Panel,length(Panel)/2,2);
   PanelAddress = 1:(NumY + mod(NumY,2));
   
   PanelAddress = reshape(PanelAddress,2,length(Panel(:))/2)';
   FS = 8;
   FS_title = 14;
   YL = zeros(NumY,2); % to store ylimit values for plots
   
   if any(PlotSites==s) % if we are in one of the sites to be plotted
       for y = 1:NumY
           figure(2)
           subplot(length(Panel(:))/2,2,Panel(PanelAddress(y)))
           hold on
           
           
           if strcmp(Dir,'mockdata_fits_June2015/') % if doing mockdata:
               % data was binned by 3 cm, but treated as continuous
               bar(x,Nact(:,Years(y))./(3/diff(x(1:2))),(3/diff(x(1:2))),'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])
           
           else % real data
               bar(x,Nact(:,Years(y))./diff(x(1:2)),(diff(x(1:2))),'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])
           end
           
           % Aggregate IPM density to make it comparable to the scale of the data histogram
           if strcmp(Dir,'mockdata_fits_June2015/') % if doing mockdata:
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
               % plot(xt-diff(xt(1:2)/2),Npred2,'ko')
           end
           
           % Posterior predictive distribution:
          % keyboard
               Npred_low = quantile(Npreddist,0.25,3);
               Npred_high = quantile(Npreddist,0.75,3);
           
           % Plot the raw density:
           if ~isnan(Nact(1,Years(y)))
           plot(x,Npred(:,Years(y)).*NT(s,Years(y)-length(Tpre)),'k','linewidth',1) 
           plot(x,Npred_low(:,Years(y)).*NT(s,Years(y)-length(Tpre)),'k--','linewidth',1) 
           plot(x,Npred_high(:,Years(y)).*NT(s,Years(y)-length(Tpre)),'k--','linewidth',1) 

           YL(y,:) = get(gca,'ylim');
           else
           ylim([0 1])
           text(2,0.2,'ND','fontsize',FS)
           end
           set(gca,'tickdir','out','ticklength',[0.015 0.015],'fontsize',FS)
           set(gca,'xlim',[0 60])
           
           
           if strcmp(Dir,'mockdata_fits_June2015/') 
               set(gca,'ylim',[0 20])
           end
           
       end % end loop over years for plotting
       
     %  keyboard
       % Now set the Ylim to be the same in all panels (currently not used)
       Ylm = nanmax(YL(:,2));
       for y = 1:NumY
           subplot(length(Panel(:))/2,2,Panel(PanelAddress(y)))
          % subplot(NumY,length(PlotSites),Panel(y,find(PlotSites==s)))
       %    set(gca,'ylim',[0 Ylm])
         if ~isnan(Nact(1,Years(y)))
           text(52,YL(y,2)*0.9,num2str(Years(y)+1990),'fontsize',FS)
         else
           text(52,0.9,num2str(Years(y)+1990),'fontsize',FS)
         end
       end
       
       % Label bottom x-axis
       subplot(length(Panel(:))/2,2,Panel(PanelAddress(end-1)))
       xlabel('Length (cm)','fontsize',FS_title)
       subplot(length(Panel(:))/2,2,Panel(PanelAddress(8)))
       xlabel('Length (cm)','fontsize',FS_title)
       
       % Label y-axis
       subplot(length(Panel(:))/2,2,7)
       ylabel('Population density (fish per transect)','fontsize',FS_title)
       
       % Title top panel
       subplot(length(Panel(:))/2,2,1)
    %   title(plot_title,'fontsize',FS_title)
    title(strcat(Spname_human,' : '),'fontsize',FS_title)
    subplot(length(Panel(:))/2,2,2)
    title(Sitename_human,'fontsize',FS_title)
       
       
       
       if saveplots
           print(timeseries_plotname,'-depsc2','-tiff')
       end
       
       % Plot posterior predictive as well
       figure(3)
       set(gcf,'units','cent','position',[5,5,18,18])
       clf
       
       %keyboard
     %  Psum = squeeze(sum(Npreddist(x>=Meta.fixparm(5),:,:),1));
     %  Nsum = sum(Nact(x>=Meta.fixparm(5),:));
      
   %   Psum = squeeze(sum(Npreddist(x>=10,:,:),1));
    %  Nsum = sum(Nact(x>=10,:));
    
   % Psum = squeeze(sum(Npreddist(:,:,:).*repmat(x(:),[1,size(Npreddist,2),size(Npreddist,3)]),1))...
   %     ./squeeze(sum(Npreddist(:,:,:)));
    
   % Nsum = sum(Nact.*repmat(x(:),[1,size(Nact,2)]))./sum(Nact);
    Len = x>=10;
    Len = Len(:);
    x = x(:);
    Psum = squeeze(sum(Npreddist(Len,:,:).*repmat(x(Len),[1,size(Npreddist,2),size(Npreddist,3)]),1))...
        ./squeeze(sum(Npreddist(Len,:,:)));
    
    Nsum = sum(Nact(Len,:).*repmat(x(Len),[1,size(Nact,2)]))./sum(Nact(Len,:));
    
 %   keyboard
    
       
      for y = 1:NumY
           subplot(length(Panel(:))/2,2,Panel(PanelAddress(y)))
           
       hold on
       if ~isnan(Nact(1,Years(y)))
    %   histogram(squeeze(Psum(Years(y),:)*NT(s,Years(y)-length(Tpre))),'Normalization','probability')
       histogram(squeeze(Psum(Years(y),:)),'Normalization','probability')
       plot([Nsum(Years(y)) Nsum(Years(y))],[0 1],'r-')
       
       end
       

       set(gca,'tickdir','out','ticklength',[0.015 0.015],'fontsize',FS)
       set(gca,'xlim',[0 60])
       set(gca,'ylim',[0 0.5])
       
      end
      
      for y = 1:NumY
           subplot(length(Panel(:))/2,2,Panel(PanelAddress(y)))
         
         if ~isnan(Nact(1,Years(y)))
           text(52,0.45,num2str(Years(y)+1990),'fontsize',FS)
         else
           text(52,0.45,num2str(Years(y)+1990),'fontsize',FS)
         end
      end
       
       

      
       % Label bottom x-axis
       subplot(length(Panel(:))/2,2,Panel(PanelAddress(end-1)))
       xlabel('Mean length (cm)','fontsize',FS_title)
       subplot(length(Panel(:))/2,2,Panel(PanelAddress(8)))
       xlabel('Mean length (cm)','fontsize',FS_title)
       
       % Label y-axis
       subplot(length(Panel(:))/2,2,7)
       ylabel('Frequency','fontsize',FS_title)
       
       
       % Title top panel
       subplot(length(Panel(:))/2,2,1)
    %   title(plot_title,'fontsize',FS_title)
    title(strcat(Spname_human,' : '),'fontsize',FS_title)
    subplot(length(Panel(:))/2,2,2)
    title(Sitename_human,'fontsize',FS_title)
      
     if saveplots
           print(ppcheck_plotname,'-depsc2','-tiff')
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
if Plothisto
figure(10)
clf
set(gcf,'units','cent','position',[14 5 8.4 4])
hold on
H=histogram(Pvec(:,2),50,'Normalization','probability');
plot([F_mode,F_mode],[0 2e3],'k-')
plot([Meta.Fprior(2,1),Meta.Fprior(2,1)],[0 2e3],'k--')
xlabel('F')
ylabel('Frequency')
ylim([0 max(H.Values)*1.2])
xlim([0.0 0.5])
set(gca,'tickdir','out','ticklength',[0.015 0.015])
if saveplots
    print(F_histo_plotname,'-depsc2','-tiff')
end
end


% Save post-burn chains & posterior of everything; final distribution of Npred; save F posterior:
post_savename = Meta.fit_savename(1:end-4); % trim off '.mat'
post_savename = strcat(post_savename,'_postproc.mat');
Meta.post_savename = post_savename;
save(strcat(Dir,Meta_savename))

Post(1).posterior = Pvec; % the posterior distribution of all the parameters
Post(1).predicted = Predicted; % The predicted distribution at all sites and times
Post(1).Rfact = Rfact;




