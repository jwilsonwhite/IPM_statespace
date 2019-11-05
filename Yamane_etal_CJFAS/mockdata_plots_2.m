function mockdata_plots_2

% Plot results of mockdata runs for state space IPM manuscript

% Fig 1a: Plot estimated values of F for 'baseline' runs

figure(1)
clf
set(gcf,'units','cent','position',[10,10,12.9 18])


% ------------------------------------------------------
% Loop over 'true' values then 10 runs in each one
Fs = {'0.15'}; %'0.0','0.05','0.1','0.2','0.25'
Savename_Base1 = 'mockdata_fits_Mar2018/PCLA_mockdata_F'
Savename_Base2 = '_8Aug2018'
nchains = 4;
nruns = 10; %
burn = 2e3; % burn-in

Colors = {'k','k','k','k'};
Marker = {'o','o','o','o'};
Subplots = [1 3 5];

for f = 1
    
subplot(4,2,Subplots(f))
hold on
    
    % Plot actual value
    plot([0,11],repmat(str2num(Fs{f}),[1,2]),'k--');
   Parm_overall = []; %to hold the combined posterior
   hold_means = [];
   hold_sds = [];
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
        last10 = Parm(990:1000);
        end

Rhat = sqrt(var(L)/mean(vLwithin));
if Rhat > 1.1 % if there might be a convergence problem
    keyboard
end
        mean(Parm(:)); %LY added
        quantile(Parm,[0.025 0.975]); %LY added
        plot(i,mean(Parm(:)),'marker',Marker{f},'color',Colors{f})
        plot([i i],quantile(Parm,[0.025 0.975]),'k-','color',Colors{f})
        
        Parm_overall = [Parm_overall; Parm(:)];
        hold_means = [hold_means; mean(Parm(:))];
        hold_sds = [hold_sds; std(Parm(:))];
    end % end i
if f ==1
    set(gca,'ylim',[0 0.3])
    set(gca,'ytick',0:0.05:0.3,'xtick',[]) 
elseif f ==2
    set(gca,'ylim',[0.0 0.3])
    set(gca,'ytick',0.0:0.05:0.3,'xtick',[])
else
    set(gca,'ylim',[0.0 0.16])
    set(gca,'ytick',0.0:0.05:0.35,'xtick',[])
end
set(gca,'xlim',[0 11])
set(gca,'tickdir','out','ticklength',[0.015 0.015])

meanF = mean(Parm_overall)
quantile(Parm_overall,[0.025 0.975])

stddevF = std(Parm_overall)
stddevmeanF = std(mean(Parm_overall)) % not used
hold_means
hold_sds %sd of each run's posterior of F, from burn-in to end
stddevholdmeans = std(hold_means)
end % end f
% ------------------------------------------------------

    

