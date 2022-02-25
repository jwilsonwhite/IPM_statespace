function kimbro_plots(modelnum)

% Plot results of Kimbro model MCMC run

% Modified 27 Aug to handle latest version of MCMC results code
% Updated 7 Nov 2016 to handle updated form of IPM model

% Load file
if ~exist('modelnum','var')
modelnum = 1;
end

% load in post-processed data
filename = strcat(['Tomales_IPM_fits/postproc_Tomales_Model_',num2str(modelnum),'_20Jan2018.mat']);
load(filename)
% load in model run file, mostly to get model details
filename = strcat(['Tomales_IPM_fits/Tomales_Model_',num2str(modelnum),'_20Jan2018_fit.mat']);
load(filename)


% What years to plot:
Years = 5:10;%; % The years of the actual dataset
Distances = linspace(6,16,1e2); % range of distances from mouth of bay
Colors = 1-jet(length(Years));

% Indices to plot:
switch modelnum
    case {178, 180,210,124} % in the latest round these are all that matter

    Plot_opt = 'singlefig';
    Plot_opt2 = true;
    Hist_opt = false;
    xlab_j = {'Overall'};
    xlab_a = {'Overall'};
    xlab_r = {'Overall'};

end


% Do plots
switch Plot_opt
    
    case 'singlefig' % everything on one figure window
        figure(modelnum)
        clf 
        set(gcf,'units','cent','position',[40 10 9 16])
        
        %-----------------------------------------------------------------
        % Panel 1: Juvenile & adult mortality
        subplot(2,1,1)
        hold on
        
        % ---------------------
        % Plot juvenile mortality as a dot, plot adult mortality without
        % variation over years
        % Plot recruitment & data on different axes
        
      %  keyboard
       plot_spatial_pattern_notime(post_str.jM,Distances,mc_str(1).Model_str(1).jM_space,Colors,'Juvenile mortality rate (y^-^1)') 
        

        %-----------------------------------------------------------------
        
        %-----------------------------------------------------------------
        % Panel 2: Adult mortality
      %  subplot(3,1,2)
      %  hold on
        
        plot_spatial_pattern_notime(post_str.aM,Distances,mc_str(1).Model_str(1).aM_space,Colors,'Adult mortality rate (y^-^1)') 
        

        %-----------------------------------------------------------------
        
        %-----------------------------------------------------------------
        % Panel 3: Recruitment
        subplot(3,1,3)
        hold on
        
        plot_spatial_pattern(post_str.R,Years,Distances,mc_str(1).Model_str(1).rect_space,[Colors(2:end,:); 1 1 1],'Recruitment (no. 0.01 m^-^2)','rect') 
        ylim([0 1])
        
            % Add in observed recruitment data
            % Here are the data: (in # per 0.01 m2 to match pop
            % survey data)
            Obsrect=5:9;
          %  
          Data = importdata('Data/aggregate-recruitment-Jan2018.csv');
          Data = Data.data;
          Meanrect = Data(:,3);
          SErect = Data(:,4);
          
          
          Upper = Meanrect + 1.96*SErect;
          Lower = Meanrect - 1.96*SErect;
          Location= Data(:,2);
          Year = Data(:,1);
       %     Location = 15.75; % halfway between E4 & W4
            
         %   Bwidth = 0; % 0.2;
            
          %  Meanrect = log10(Meanrect + 1e-4);
           % Lower = log10(Lower + 1e-4);
           % Upper = log10(Upper + 1e-4);
           Pos = get(gca,'position');
           axes('position',Pos)
           set(gca,'color','none','yaxislocation','right')
           hold on
           
            
            for o = 1:length(Obsrect)
                
                YY = 1999+Obsrect(o);
                OK = Year==YY;
                
             %  patch(([0 0 Bwidth Bwidth]+Obsrect(o)),[0 Meanrect(o) Meanrect(o) 0],[0.5 0.5 0.5])
 
              % plot(Location-o/10,log10(Meanrect(o)),'d','markeredgecolor',Colors(o,:))
              % plot([Location, Location]-o/10,log10([Meanrect(o) Meanrect(o)+1.95*SErect(o)]),'-','color',Colors(Obsrect(o),:))
              % plot([Location, Location]-o/10,log10(max(1e-5,[Meanrect(o) Meanrect(o)-1.95*SErect(o)])),'-','color',Colors(Obsrect(o),:))

              
               plot(Location(OK)-o/10,(Meanrect(OK)),'o','markeredgecolor',Colors(o+1,:))
               
               Loctmp = Location(OK);
               Lowertmp = Lower(OK);
               Uppertmp = Upper(OK);
               for i = 1:4
               plot([Loctmp(i) Loctmp(i)]-o/10,([Lowertmp(i) Uppertmp(i)]),'-','color',Colors(o+1,:))
               end
               %  plot([Location, Location]-o/10,([Meanrect(o)+1.96 Meanrect(o)-1.95*SErect(o)])),'-','color',Colors(Obsrect(o),:))

               
            end
         ylim([0 30])
         set(gca,'xtick',[],'tickdir','out','ticklength',[0.015 0.015])
         set(gca,'ytick',0:10:30)
        
        % Just 2008 data:
    %    Meanrect = [0.25 1.15 3.25 17.95];
    %    SDrect = [0.33 1.00 2.07 9.68];
    %    Upper = Meanrect + 1.96*SDrect;
    %    Lower = Meanrect - 1.96*SDrect;
    %    Location = [6 8 12 16];
    %    for o = 1:length(Location)
    %    plot([Location(o), Location(o)],([Lower(o) Upper(o)]),'-','color',Colors(end-1,:))
    %    plot(Location(o),(Meanrect(o)),'d','markeredgecolor',Colors(end-1,:))
    %    end
        
        
         
end % End switch plot opt
        
 
      
            
    function plot_spatial_pattern(Post,Years,Distances,Space,Colors,Title,Type)       
            
        for y = 1:length(Years)
            YY = Years(y);
            switch Type
                case 'mort'
                    YY = YY -1; % mortality in year YY affects abundance in the *following* year
                otherwise
                 %   keyboard
            end
        
        % Determine the model type & get mortality rates 
        switch Space
            case 'Con'
                X = repmat(exp(Post(:,YY)),[1,length(Distances)]);
            case 'Grad'
                X = exp(repmat(Post(:,YY),[1,length(Distances)]) + repmat(Distances,[size(Post,1),1]).*Post(:,end));
            case 'Unim'
                X = exp(repmat(Post(:,YY),[1,length(Distances)]) + Post(:,end-1).*(repmat(Distances,[size(Post,1),1])-Post(:,end)).^2);
        end % end aM_space switch
        
      %  X = log10(X);
        
        plot(Distances,quantile(X,0.25),'color',Colors(y,:),'linestyle','--')
        plot(Distances,quantile(X,0.75),'color',Colors(y,:),'linestyle','--')
        ph(y) = plot(Distances,quantile(X,0.5),'color',Colors(y,:),'linestyle','-');
        text(max(Distances)+0.1,quantile(X(:,end),0.5),num2str(YY+1999))

        
        end % end loop over years
        
        % Maybe find way to show confidence intervals for each year?
        
        xlabel('Distance from bay mouth (km)')
        ylabel(Title)
     %   legend(ph,'2004','2005','2006','2009')
     if length(Years)==10
    % legend(ph,'2000','2001','2002','2003','2004','2005','2006','2007','2008','2009')
     else
    %  legend(ph,'2000','2001','2002','2003','2004','2005','2006','2007','2008')
     end
     
     set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
     set(gca,'xtick',[6 8 12 16],'xlim',[5 17]) % just where the sites are
    
 

    function plot_spatial_pattern_notime(Post,Distances,Space,Colors,Title)       
            
        
        % Determine the model type & get mortality rates 
        switch Space
            case 'Con'
                X = repmat(exp(Post(:,1)),[1,length(Distances)]);
            case 'Grad'
                X = exp(repmat(Post(:,1),[1,length(Distances)]) + repmat(Distances,[size(Post,1),1]).*Post(:,end));
            case 'Unim'
                X = exp(repmat(Post(:,1),[1,length(Distances)]) + Post(:,end-1).*(repmat(Distances,[size(Post,1),1])-Post(:,end)).^2);
        end % end aM_space switch
        
      %  X = log10(X);
        
        plot(Distances,quantile(X,0.25),'color','k','linestyle','--')
        plot(Distances,quantile(X,0.75),'color','k','linestyle','--')
        ph = plot(Distances,quantile(X,0.5),'color','k','linestyle','-');
        %text(max(Distances)+0.1,quantile(X(:,end),0.5),num2str(YY+1999))

       
        xlabel('Distance from bay mouth (km)')
        ylabel(Title)
     %   legend(ph,'2004','2005','2006','2009')
    
     
     set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
     set(gca,'xtick',[6 8 12 16],'xlim',[5 17]) % just where the sites are
     axis square

