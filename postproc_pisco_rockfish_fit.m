function postproc_pisco_rockfish_fit(Meta_savename)

% Examine posterior of an MCMC run
%Dir ='runs_for_publication_pre2007_July2015/';
%Dir = 'mockdata_fits_June2015/'; % for mockdata runs
%Dir = '';
Dir = 'post2007_runs_July2015/';

% Load in metadata
load(strcat(Dir,Meta_savename))

load(strcat(Dir,Meta.data_savename),'D_str','Species_Names','Years')
Plotchains = true;
Plotfinal = true;
saveplots = false;

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

%timescale to run model
Tdata = Meta.Tdata; % timespan of data

Tpre = Meta.Tpre; % time prior to data that model will run
if isnan(Tpre(1)) % if it is a post2007 run with no pre-period
    T = 1:length(Tdata(:));
else
T = 1:length([Tpre(:); Tdata(:)]);
end

%get data into a histogram
Site_Names = Meta.Sites

N = IPM_histo(D_str.(Species_Names{1}),Years,Site_Names,edges);
%divide total fish numbers by number of transects
%for i = 1:length(Years)
%    for j = 1:length(Site_Names)
%        NT = D_str.(Species_Names{1}).(Site_Names{j})(i).data.numtrans;
%        N(:,j,i) = N(:,j,i)./NT;
%    end
%end
% strip out the years we don't need
OKyears = false(length(Years),1);
for i = 1:length(Years)
    if any(Meta.Tdata == Years(i))
        OKyears(i) = true;
    end
end
N = N(:,:,OKyears);
% specify ogive for observations - don't see fish that are too large
% b/c of ontogenetic migration
Ogive_b = Meta.ogive;
%Ogive = 1./(1 + 1./exp(Ogive_b(1) + Ogive_b(2).*x)); % probability of observation in the kelp forest
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

%make vector of recruits to add in each year, this will be multiplied by
%the magnitude calculated above. Specify mean of pdf for size according to
%specific species, specify sd too.

R = normpdf(x,Meta.recruits.meansize,Meta.recruits.sdsize)'; %Average size of recruit is 4, with sd = 1
                         
%scale R so integrates appropriately
dy=x(2)-x(1);
Sy=makeSimpVec(dy,meshsize);
%R = R./(Sy*R).*sum(R); %rescale R so integrates correctly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(strcat(Dir,Meta.fit_savename))
%load(strcat(Meta.fit_savename))

chains = length(mc_str);
%if ~exist('mc_str(1).NT','var')
%    NT = ones(length(Site_Names),length(Tdata));
%else
if strcmp(Dir,'mockdata_fits_June2015/') % if doing mockdata
    NT = ones(1,length(T)); 
else
NT = mc_str(1).NT;
NT = NT(:,1:length(Meta.Tdata));
%NT = [zeros(size(NT,1),length(Meta.Tpre)),NT];
end


if strcmp(Dir,'mockdata_fits_June2015/')
Burn = [1e3 1e3 1e3];
else
Burn = [3e3 3e3 3e3];
end
Col = {'r','b','k'};

if Plotchains
    figure(1)
    clf
    for i = 1:size(mc_str(1).parm_vec,2)
    subplot(5,3,i)
    hold on
    for c = 1:chains
    plot(mc_str(c).parm_vec(Burn(c):end,i),Col{c})
    end
    end
    
    subplot(5,3,15)
    hold on
    for c = 1:chains
        plot(mc_str(c).L(Burn(c):end),Col{c})
    end
    
    
end

if Plotfinal
   % figure(2)
    figure(2)
    set(gcf,'units','cent','position',[5,5,17.4, 15])
    clf
   Pvec = [];
   for c = 2; %1:chains
       Pvec = [Pvec; mc_str(c).parm_vec(Burn(c):end,:)];
   end
   meanP=mean(Pvec(1:end,:));

   
if strcmp(Meta_savename(6),'P') % Point Lobos
PlotSites = [1, 3]; % which sites to plot? 

% which years to plot
if ~isnan(Meta.Tpre) % pre2007 runs
        Years = [10 17 18];%9+(1:NumY);
else % post2007 runs
        Years = 2:4;
end

elseif strcmp(Meta_savename(6),'B')  % Big Creek
  
PlotSites = [1, 3]; % which sites to plot? 
    
% which years to plot
if ~isnan(Meta.Tpre) % pre2007 runs
        Years = [12 17 18];%9+(1:NumY);
else % post2007 runs
        Years = 4:6;
end

else % White Rock
    
        if ~isnan(Meta.Tpre) % if it is a pre2007 run
    PlotSites = [4, 5]; % which sites to plot? 
        else
            if strcmp(Meta_savename(25),'f')
                PlotSites = [1,2];
            else
            PlotSites = [1,2];
            end
        end
    
% which years to plot
if ~isnan(Meta.Tpre) % pre2007 runs
    Years = [14 17 18];%9+(1:NuY);
else % post2007 runs
    Years = 4:6;
end

end % end switch over site
   
 Predicted = struct();
 %Site_Names
   for s = 1:length(Site_Names)

       % switch off fishing for pre-existing MPAs
        meanP_tmp = meanP;
        meanP_tmp(1:2) = meanP_tmp(1:2).*(1-Meta.MPAstatus(s)); 
        
   [L(s), Npred, Nact] = do_IPM(meshsize,x,dy,Sy,meanP_tmp,Meta.fixparm,Rfact,R,squeeze(N(:,s,:)),T,Tpre,Ogive,Meta,NT(s,:),Site_Names{s});


   Predicted(1).(Site_Names{s}).Npred = Npred;
   Predicted(1).(Site_Names{s}).Likelihood = L;


NumY = 3; %length(Tdata); % how many years to plot

Panel = 1:(NumY*length(PlotSites));
Panel = reshape(Panel,length(PlotSites),NumY)';
FS = 8;

if any(PlotSites==s) % if we are in one of the sites to be plotted
    
    
for y = 1:NumY
subplot(NumY,length(PlotSites),Panel(y,find(PlotSites==s)))
hold on
sum(Nact(:,Years(y)))
sum(Npred(:,Years(y))*diff(x(1:2)).*NT(s,Years(y)-length(Tpre)))

if strcmp(Dir,'mockdata_fits_June2015/') % if doing mockdata:
% data was binned by 3 cm, but treated as continuous
bar(x,Nact(:,Years(y))./3,3,'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])
else % real data
bar(x,Nact(:,Years(y))./diff(x(1:2)),'facecolor',[0.6 0.6 0.6],'edgecolor',[0.6 0.6 0.6])
end
plot(x,Npred(:,Years(y)).*NT(s,Years(y)-length(Tpre)),'k','linewidth',1) 
set(gca,'tickdir','out','ticklength',[0.02 0.02],'fontsize',FS)
%set(gca,'ylim',[0 15])
%set(gca,'ylim',[0 70])
set(gca,'xlim',[0 60])
if y == 1; title(Site_Names{s},'fontsize',FS); end
xlabel('Length (cm)','fontsize', FS)
ylabel('Relative frequency','fontsize',FS)
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
Post(1).F_std = std(Pvec(:,2))

% Plot the histogram of historical F
figure(10)
clf
set(gcf,'units','cent','position',[8 5 8.4 4])
hold on
hist(Pvec(:,2),50)
plot([F_mode,F_mode],[0 2e3],'k-')
plot([Meta.Fprior(2,1),Meta.Fprior(2,1)],[0 2e3],'k--')
%set(gca,'xtick',0:0.05:10,'tickdir','out','ticklength',[0.02 0.02])
%set(gca,'ylim',[0 1500])
xlabel('F')
ylabel('Frequency')
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



