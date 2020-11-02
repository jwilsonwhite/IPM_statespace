function plot_ChannelIslands_posteriors

% Create figure to display posterior estimates of F for Channel Islands
% analysis

do1 = false;
do2 = true;

% List of sites, in desired order:
Sites = {'SMI_TYLER_BIGHT','SMI_CROOK_POINT','SMI_HARRIS_PT_RESERVE','SMI_CUYLER',...
         'SRI_SOUTH_POINT','SRI_JOHNSONS_LEE_SOUTH','SCI_FORNEY','SCI_PAINTED_CAVE',...
         'SCI_HAZARDS','SCI_GULL_ISLE','SCI_VALLEY','SCI_YELLOWBANKS','SCI_PELICAN',...
         'SCI_CAVERN_POINT','SCI_SCORPION','ANACAPA_WEST_ISLE','ANACAPA_MIDDLE_ISLE',...
         'ANACAPA_EAST_ISLE','ANACAPA_LIGHTHOUSE_REEF'};
isMPA = [0 0 1 0  1 0 0 1 0 1 0 0 0 1 1 1 1 1 0];
       
Species = {'SATR','SMYS','SPUL','PCLA'};

Priors = [0.07, 0.03, 0.08, 0.1]; % priors on F
Priors_sd = [1.0, 1.0, 0.23, 1.0]; % sds on priors

% Which species to plot for each site?
doSpecies = [1 1 1 0 1 1 1 0 1 1 0 0 1 0 0 0 0 0 0; % SATR
             1 1 1 1 1 1 0 1 1 1 0 0 1 0 0 0 0 0 0; % SMYS
             0 0 0 0 1 1 0 1 0 1 1 0 1 1 0 1 0 1 1; % SPUL
             0 0 0 0 0 0 0 1 1 0 0 1 1 1 1 1 0 0 1];% PCLA
         
         
Colors = [0.5 0 0.5;
          0.05 0.05 1;
          0.8 0.1 0.1;
          0.1 0.8 0.1];
% Format figure (old version, not used in current ms)    
if do1
figure(1)
clf
set(gcf,'units','cent','position',[10 10 21 9])
hold on

for s = 1:length(Sites)
    
    for p = 1:length(Species)
        
                % Box to outline the MPA
        if isMPA(s)
            MPAcol = 'r';
        else
            MPAcol = 'b';
        end
        patch([s-0.4 s-0.4 s+0.4 s+0.4],[0 4 4 0],MPAcol,'FaceAlpha',0.1)
        
        % plot posterior + mean
        if doSpecies(p,s)
        fname = strcat(Species{p},'_',Sites{s},'_post2003_fit_July2019_postproc.mat');
        load(fname,'Post')
        F = Post.posterior(:,2); % Fpost
        
        xcoord = s + (p - 2.5)/8;
        
        plot(repmat(xcoord,length(F),1)+(rand(length(F),1)-0.5)/10,F,'ko','markersize',5,'markeredgecolor',Colors(p,:));
        plot(xcoord,mean(F),'kd','linewidth',2)
        end
       
        
    end % end Species
end % end Sites

% Format axes
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
set(gca,'xtick',1:length(Sites),'xticklabels',Sites,'xticklabelrotation',300)
end % end if do1


% Figure 2 - separate panels for each species, with histograms
if do2
figure(2)
clf
set(gcf,'units','cent','position',[10 10 21 9])

bins = linspace(eps,4,1e3); % bins for histogram

for p = 1:length(Species)
for s = 1:(length(Sites)+1)
   
    
    subplot(2,2,p)
    hold on
    
    if s <= length(Sites) % for plotting the posteriors
        
         Site_human{s} = human_name(Sites{s});
    % plot posterior + mean
        if doSpecies(p,s)
        % Box to outline the MPA
        if isMPA(s)
            MPAcol = 'r';
        else
            MPAcol = 'b';
        end
       patch([-0.2 4 4 -0.2],[s-0.3 s-0.3 s+0.3 s+0.3],MPAcol,'FaceAlpha',0.5)
        
        % USE ksdensity for plotting histograms
        % Add lines for prior & median
        
        
        fname = strcat(Species{p},'_',Sites{s},'_post2003_fit_July2019_postproc.mat');
        load(fname,'Post')
        F = Post.posterior(:,2); % Fpost
        
        % get ks density
        [N,Edges] = histcounts(F,'normalization','pdf');
        Midpoints = Edges(1:end-1) + diff(Edges(1:2));
        Values = N/max(N)*0.6;
        plot(Midpoints,Values + s-0.3,'k')
        fill(Midpoints,Values + s-0.3,[0.5 0.5 0.5])
        
        
        plot([median(F),median(F)],[s-0.3 s+0.3],'color',[1 1 1],'linewidth',2)
                
        else
            % Plot a ND if no data for that species/site
         %   text(0,s,'ND')
            
        end % end if doSpecies
        
    else % Plot the prior
        Site_human{s} = 'Prior';
        x = linspace(0,2,1e3);
        Pf = normpdf(log(x),log(Priors(p)),Priors_sd(p));
        plot(x,Pf,'k');
        fill(x,Pf,[0.5,0.5,0.5])
        
        
    end % end plotting posterior vs. prior
               
        
    end % end Sites
    xlim([-0.2 1.7])
    ylim([-0.4 length(Sites)+0.4])
    set(gca,'ytick',1:(length(Sites)+1),'yticklabels',Site_human)
    set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
    set(gca,'xgrid','on')

end % end Species

% Format axes

end % end if do2

