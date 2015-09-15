function mockdata_plots

% Plot results of mockdata runs for state space IPM manuscript

% Fig 1a: Plot estimated values of F for 'baseline' runs

figure(1)
clf
set(gcf,'units','cent','position',[10,10,12.9 18])


% ------------------------------------------------------
% Loop over 3 'true' values then 10 runs in each one
Fs = {'0.0','0.05','0.1'};
Savename_Base1 = 'mockdata_fits_June2015/SMYS_mockdata_F';
Savename_Base2 = '_baseline_28May2015';
nchains = 2;
nruns = 10; %
burn = 2e3; % burn-in

Colors = {'k','k','k'};
Marker = {'o','o','o'};
Subplots = [1 3 5];

for f = 1:3
    
    
subplot(4,2,Subplots(f))
hold on
    
    % Plot actual value
    plot([0,11],repmat(str2num(Fs{f}),[1,2]),'k--');
   Parm_overall = []; %to hold the combined posterior
    for i = 1:nruns
        
        % Read in data
        Sname = strcat(Savename_Base1,Fs{f},Savename_Base2,'_R',num2str(i),'_fit.mat');
        load(Sname,'mc_str')
      
        % Combine values, do convergence diagnostics (based on Gelman & Shirley 2011)
        Parm = []; L = []; 
        vLwithin = nan(length(mc_str),1);
        for n = 1:nchains
        Parm = [Parm(:); mc_str(n).parm_vec(burn:end,2)];
        L = [L(:); mc_str(n).L(burn:end)'];
        vLwithin(n) = var(mc_str(n).L(burn:end));
        end

Rhat = sqrt(var(L)/mean(vLwithin));
if Rhat > 1.1 % if there might be a convergence problem
    keyboard
end
        
        plot(i,mean(Parm(:)),'marker',Marker{f},'color',Colors{f})
        plot([i i],quantile(Parm,[0.025 0.975]),'k-','color',Colors{f})
        
        Parm_overall = [Parm_overall; Parm(:)];
        
    end % end i
if f ==1
    set(gca,'ylim',[-0.005 0.03])
    set(gca,'ytick',0.0:0.01:0.3,'xtick',[])
elseif f ==2
    set(gca,'ylim',[0.03 0.08])
    set(gca,'ytick',0.0:0.01:0.3,'xtick',[])
else
    set(gca,'ylim',[0.0 0.18])
    set(gca,'ytick',0.0:0.04:0.3,'xtick',[])
end
set(gca,'xlim',[0 11])
set(gca,'tickdir','out','ticklength',[0.015 0.015])

mean(Parm_overall)
quantile(Parm_overall,[0.025 0.975])

end % end f
% ------------------------------------------------------

% ------------------------------------------------------
% Then also plot the values for split runs  - 0.05/0.0
Savename_Base = 'mockdata_fits_June2015/SMYS_mockdata_F0.05_0.0_baseline_28May2015';
nchains = 2;
nruns = 10; 
burn = 1e3; % burn-in

subplot(4,2,7)
hold on

    % Plot actual value
    plot([0,11],[0.0 0.0],'k--');
    plot([0,11],[0.05 0.05],'k--');
   
    for i = 1:nruns
        
        % Read in data
        Sname = strcat(Savename_Base,'_R',num2str(i),'_fit.mat');
        load(Sname,'mc_str')
        
        % Combine values, do convergence diagnostics (based on Gelman & Shirley 2011)
        Parm = []; L = []; 
        vLwithin = nan(length(mc_str),1);
        for n = 1:nchains
        Parm = [Parm(:); mc_str(n).parm_vec(burn:end,1)];
        L = [L(:); mc_str(n).L(burn:end)'];
        vLwithin(n) = var(mc_str(n).L(burn:end));
        end
        
        % clip out parts of the chain that did not converge:
        if i == 2 || i == 5
            L = L(1:2.9e3); Parm = Parm(1:2.9e3); vLwithin = var(L);
          %  
          elseif i == 5
            keyboard
        end
        
        Rhat = sqrt(var(L)/mean(vLwithin));
        if Rhat > 1.1 % if there might be a convergence problem
        keyboard
        end
        
        plot(i,mean(Parm(:)),'marker','d','color','b')
        plot([i i],quantile(Parm,[0.025 0.975]),'k-','color','b')
        
        % Combine values - F2
        Parm = [];
        for n = 1:nchains
        Parm = [Parm(:), mc_str(n).parm_vec(burn:end,2)];
        end
        
        if i == 2 || i == 5
            Parm = Parm(1:2.9e3); 
        
        end
        
        plot(i,mean(Parm(:)),'ko')
        plot([i i],quantile(Parm,[0.025 0.975]),'k-')
        
        
    end % end i

set(gca,'xlim',[0 11],'ylim',[-0.005 0.06])
set(gca,'ytick',0.0:0.02:0.3,'xtick',[])
set(gca,'tickdir','out','ticklength',[0.02 0.02])
    
% ------------------------------------------------------

% ------------------------------------------------------
% Fig 1b-d: Plot estimated values of F for 'worse' runs
SubP = [2 4 6];

% Loop over 3 'true' values then 10 runs in each one
Fs = {'7years','3years','10bin'};
Savename_Base1 = 'mockdata_fits_June2015/SMYS_mockdata_F0.05_';
Savename_Base2 = '_28May2015';
nchains = 2;
nruns = 10; % expand to 10 eventually
burn = 1e3; % burn-in

for f = 1:3
    
    subplot(4,2,SubP(f))
    hold on
    
    % Plot actual value
    plot([0,11],[0.05 0.05],'k--');
   
    for i = 1:nruns
        
        % Read in data
        Sname = strcat(Savename_Base1,Fs{f},Savename_Base2,'_R',num2str(i),'_fit.mat');
        load(Sname,'mc_str')
        
        % Combine values
        % Combine values, do convergence diagnostics (based on Gelman & Shirley 2011)
        Parm = []; L = []; 
        vLwithin = nan(length(mc_str),1);
        for n = 1:nchains
        Parm = [Parm(:); mc_str(n).parm_vec(burn:end,2)];
        L = [L(:); mc_str(n).L(burn:end)'];
        vLwithin(n) = var(mc_str(n).L(burn:end));
        end
        
        % clip out some chains that did not converge
        if f == 1
        elseif f == 2
         if i == 9
        L = L(1:2.9e3); Parm = Parm(1:2.9e3); vLwithin = var(L);
         end
        elseif f == 3
            if  i == 5
                L = L(7e3:end); Parm = Parm(7e3:end); vLwithin = var(L);
            end
        end
        
        Rhat = sqrt(var(L)/mean(vLwithin));
        if Rhat > 1.1 % if there might be a convergence problem
        keyboard
        end
        
        plot(i,mean(Parm),'ko')
        plot([i i],quantile(Parm,[0.025 0.975]),'k-')
        
        
    end % end i
    
set(gca,'xlim',[0 11],'ylim',[-0.01 0.08])
set(gca,'ytick',0:0.02:0.3,'xtick',[])
set(gca,'tickdir','out','ticklength',[0.015 0.015])
end % end f


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 2: plot recruitment
figure(2)
clf
set(gcf,'units','cent','position',[20,10,19 12])

Colors = {'k','b','r'};
Marker = {'o','^','d'};
Msize = 4;

Fs = {'0.0','0.05','0.1'};
Savename_Base1 = 'mockdata_fits_June2015/SMYS_mockdata_F';
Savename_Base2 = '_baseline_28May2015';
for f = 1:3
    
    % Plot actual value
    plot([0,11],repmat(0,[1,2]),'k-'); % mean of 50, lognormally distributed
    
    for i = 1:nruns
        
        % Read in data
        Sname = strcat(Savename_Base1,Fs{f},Savename_Base2,'_R',num2str(i),'_fit.mat');
        load(Sname,'mc_str')

       % Combine values
        Parm = [];
        for n = 1:nchains
        Parm = [Parm; mc_str(n).parm_vec(burn:end,3:14)];
        end
       
                % Get rid of non-converged runs:
        if f == 1
            if i == 3 || i == 10
                Parm = Parm(3.1e3:end,:); 
            elseif i == 6
                Parm = Parm(1:2.9e3,:);
            end
        elseif f == 2
            if i == 1
                Parm = Parm(3.1e3:end,:); 
            elseif i == 5
                Parm = Parm(1:2.9e3,:); 
            end
        elseif f == 3
            if i == 1
                Parm = Parm(1:2.9e3,:); 
            end
        end
        
        % Read in original data with actual values
        Sname = strcat('mockdata_May2015/SMYS_mockdata_F',Fs{f},Savename_Base2,'_R',num2str(i),'.mat');
        load(Sname)

        R = Parm(:,1:end-2) - repmat(log(RF.Ractual(9:end))',[size(Parm,1),1]); % difference between actual and estimated value
        
        subplot(2,1,1)
        hold on
        plot([0 12],[0 0],'k--')
        for j = 1:size(R,2)
        Rind = (rand-0.5)*0.4;
        plot(i+f/3+Rind,mean(R(:,j)),'marker',Marker{f},'color',Colors{f},'markersize',Msize)
        plot([i i]+f/3+Rind,quantile(R(:,j),[0.025 0.975]),'k-','color',Colors{f})
        end
        
        subplot(2,1,2)
        hold on
        plot(i+f/3,mean(Parm(:,end)),'marker',Marker{f},'color',Colors{f})
        plot([i i]+f/3,quantile(Parm(:,end),[0.025 0.975]),'k-','color',Colors{f})
        ylim([-0.01 0.3])
        
    end % end i

end % end f

for i = 1:2
    subplot(2,1,i)
set(gca,'xlim',[1 11.5])
set(gca,'tickdir','out','ticklength',[0.015 0.015])
end


% ------------------------------------------------------

    

