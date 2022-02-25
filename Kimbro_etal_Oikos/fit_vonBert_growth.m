function fit_vonBert_growth

% Fit von Bertalanffy curves to Tomales oyster growth data

% Load data: 
D = importdata('Data/Tomales_site_growth_24Jan2018.csv');

% D has columns: Site, Age (months), Mean, Std
D = D.data;

% Curve:
% L(x) = Linf - (Linf - L0).*exp(-k.*x)
doFits=true;
if doFits
% Individual sites:
for s = 1:7
    Dt = D(D(:,1)==s,2:4); % collect corresponding rows
    save('Dttmp.mat','Dt')
    
    x0 = [40,4.5,0.2]; % Linf, L0, k
    
    b = fmincon(@vonBert,x0,[],[],[],[],[0,0,0],[Inf,Inf,10]); 
    % constrain all parameters to be positive, force k to be reasonable
    % (some unconstraind fits estimate k >> 100)
    
    X(s,1:3) = b;
    X(s,4) = mean(Dt(:,3)./Dt(:,2)); % CV of variability
    
end % end loop over site

% fix units:
X(:,1) = X(:,1)*10; % cm to mm
X(:,3) = X(:,3)*12; % months to years

% Save results
Meta = 'X is n x 4 matrix for n = 7 sites. Params are Linf, L0, k, CV';

save('Tomales_vB_fits.mat','X','Meta');
else
    load('Tomales_vB_fits.mat','X','Meta');
end

doPlot = true;
x = linspace(0,10,100);
figure(1)
clf
hold on
Col = {[0.8 0.2 0],[0.6 0.1 0.6],[0 0 1],...
       [1 0 0],[0.8 0.2 0],[0.6 0.1 0.6],[0 0 1]};
       
if doPlot
   for i = 1:7
       L = X(i,1).*(1-exp(-X(i,3).*x));
       ph = plot(x,L);
       set(ph,'color',Col{i},'linewidth',2)

   end
   
   set(gca,'tickdir','out','ticklength',[0.015 0.015])
   set(gca,'xcolor','k','ycolor','k','fontsize',10)
   ylabel('Length (mm)','fontsize',12)
   xlabel('Age (y)','fontsize',12)
   
end


function F = vonBert(x)

load Dttmp.mat Dt
L = x(1) - (x(1) - x(2)).*exp(-(x(3).*Dt(:,1)));

% SS difference, scaled by standard deviation:
F = sum(((L - Dt(:,2)).^2)./Dt(:,3));

    
