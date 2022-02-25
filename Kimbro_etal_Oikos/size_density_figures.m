function size_density_figures

% Plot size & density of Olympia oysters over space in Tomales Bay
% For fig 1 in Kimbro et al. Tomales ms

% The sites:
Sites= {'e1','e2','e3','e4','w2','w3','w4'};
Distances = [6 8 12 15.1 8.2 12.1 16];
Years = [2004 2005 2006 2009];
Density = [0.168913399, 1.12281746, 2.309803922, 1.583123848, 0.083333333, 0.574278322, 4.487122627; ...
           0.222222222, 2.301851852, 0.640454793, 0.925925926, 0.088888889, 0.759259259, 1.460239651; ...
           0.087037037, 1.101851852, 1.295206972, 0.435185185, 0.074074074, 0.537037037, 0.259259259; ...
           0.152777778, 0.574659586, 1.561531279, 0.23733165, 0.018518519, 2.249155773, 0.18579793];
       
Density_sd = [0.114623266 0.657373328 1.46518222 1.339991662 0.076578049 0.644234363 0.681127482; ...
              0.314269681 2.611422361 0.432494797 0.641820969 0.120697561 0.686435296 0.127298663; ....
              0.035428012 0.562640015 0.648393246 0.519278544 0.045360921 0.325209651 0.159732286; ...
              0.170103454 0.44077454 1.01900055 0.261106735 0.028688766 0.578999762 0.371151455];
       
Length = [41.02222222, 40.63960396, 33.76506024, 23.80526316, 41.84444444, 40.49, 28.89020979; ....
          42.65, 43.71118881, 44.32, 26.774, 46.92222222, 44.26296296, 32.12026144;...
          50.06666667, 46.14622642, 44.72086331, 29.97234043, 48.425, 44.79137931, 32.87857143; ...
          29, 26.5862069, 25.7, 18.85714286, 25.5, 30.22321429, 14.26315789];
      
Length_sd = [9.917753 11.424903  8.733318  9.335497  9.554856  8.947146 9.309328;... 
             11.375914 10.556462  8.023521  7.411302 10.631296 8.582736  7.611940;...  
             9.702577  7.361770  8.237318  7.228558 9.342185  9.055991  8.403013;... 
             11.935725 12.571249 14.435950 9.226515  7.778175 14.864235  5.258243];
      
% Size
figure(1)
clf

Col1 = [1 3 5 7];
Col2 = [2 4 6 8];
Type1 = {'lin','lin','lin','quad'};
Type2 = {'lin','quad','lin','lin'};

for i = 1:4
    
    X = 5:17;
    
    subplot(4,2,Col1(i))
    hold on
    plot(Distances,Density(i,:),'ko')
    for j = 1:7
        plot([Distances(j) Distances(j)],Density(i,j)+[Density_sd(i,j) -Density_sd(i,j)],'k-')
    end
   % axis square
    box on
    set(gca,'ticklength',[0.02 0.03])

    
    switch Type1{i}
        case 'lin'
        [P,S] = polyfit(Distances,Density(i,:),1);
        [Y,Z] = polyval(P,X,S);
        %plot(X,Y+Z,'k--')
        %plot(X,Y-Z,'k--')
        Rinv = inv(S.R);
        V = (Rinv*Rinv')*S.normr^2/S.df;
        SE = sqrt(diag(V));
        T = abs(P(:)./SE(:));
        Pval1(i).pval = (1-tcdf(T,6)).*2;
        if Pval1(i).pval(1) < 0.05
        plot(X,Y,'k-')
        else
        plot(X,Y,'k--')
        end
        case 'quad'
        [P,S] = polyfit(Distances,Density(i,:),2);
        [Y,Z] = polyval(P,X,S);
        %plot(X,Y+Z,'k--')
        %plot(X,Y-Z,'k--')
        Rinv = inv(S.R);
        V = (Rinv*Rinv')*S.normr^2/S.df;
        SE = sqrt(diag(V));
        T = abs(P(:)./SE(:));
        Pval1(i).pval = (1-tcdf(T,6)).*2;
        if Pval1(i).pval(1) < 0.06
        plot(X,Y,'k-')
        else
        plot(X,Y,'k--')
        end
    end
    
    set(gca,'tickdir','out','ticklength',[0.03 0.03])
    set(gca,'xlim',[5 18],'xtick',5:5:15)
    set(gca,'ylim',[-0.5 5],'ytick',0:2:4)
    
    subplot(4,2,Col2(i))
    hold on
    plot(Distances,Length(i,:),'ko')
    for j = 1:7
        plot([Distances(j) Distances(j)],Length(i,j)+[Length_sd(i,j) -Length_sd(i,j)],'k-')
    end
  %  axis square
    box on
    
    switch Type2{i}
        case 'lin'
        [P,S] = polyfit(Distances,Length(i,:),1);
        [Y,Z] = polyval(P,X,S);
        %plot(X,Y+Z,'k--')
        %plot(X,Y-Z,'k--')
        Rinv = inv(S.R);
        V = (Rinv*Rinv')*S.normr^2/S.df;
        SE = sqrt(diag(V));
        T = abs(P(:)./SE(:));
        Pval2(i).pval = (1-tcdf(T,6)).*2;
        if Pval2(i).pval(1) < 0.05
        plot(X,Y,'k-')
        else
        plot(X,Y,'k--')
        end
        case 'quad'
        [P,S] = polyfit(Distances,Length(i,:),2);
        [Y,Z] = polyval(P,X,S);
        %plot(X,Y+Z,'k--')
        %plot(X,Y-Z,'k--')
        Rinv = inv(S.R);
        V = (Rinv*Rinv')*S.normr^2/S.df;
        SE = sqrt(diag(V));
        T = abs(P(:)./SE(:));
        Pval2(i).pval = (1-tcdf(T,6)).*2;
        if Pval2(i).pval(1) < 0.06
        plot(X,Y,'k-')
        else
        plot(X,Y,'k--')
        end
    end
    
    set(gca,'tickdir','out','ticklength',[0.02 0.02])
    set(gca,'xlim',[5 18],'xtick',5:5:15)
    set(gca,'ylim',[10 60],'ytick',0:20:60)
 
    box on
end
