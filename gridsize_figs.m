function gridsize_figs
% Code to analyze effect of grid size & integration method on IPM, and make
% plots

% Get parameters
load('SMYS_Pt_Lobos_pre2007_13Dec2013_metadata') % blue rockfish
fixparm = Meta.fixparm;
fixparm_nm = fixparm;
fixparm_nm(4) = 0; % no mortality
F = 0; % no fishing

% Define parameters to vary
xSize = 10:10:200;
xPlot = [20 50 100];
Method = {'MidP','Simp'};

figure(1)
clf
set(gcf,'units','cent','position',[10,5,19,25])

Col = {'r','k','b'};
Sty = {'-',':'};

% Loop over values of gridsize (and method)
for i = 1:length(xSize)
for j = 1:length(Method)

% Create kernel
meshsize = xSize(i);
meshmin = 0; % chosen for blue RF
meshmax = fixparm(1)*2; % chosen for blue RF
x = linspace(meshmin,meshmax,meshsize);
meshdiff = diff(x(1:2));

Rvec = normpdf(x,Meta.recruits.meansize,Meta.recruits.sdsize)'; % R size vec

switch Method{j}
    case 'MidP'
        kmat = kernmatSimp(x,F,fixparm,1); % don't let the name fool you, this makes the same kernel regardless of integration method
        kmat = kmat*meshdiff;
    case 'Simp'
        kmat = kernmatSimp(x,F,fixparm,1);
        Sy = makeSimpVec(meshdiff,meshsize);
        Symat = repmat(Sy(:)',[length(Sy),1]);
        kmat = Symat.*kmat; % do the integration
end % end switch Method

% Get stable size distribution:
N = nan(length(x),100);
N(:,1) = Rvec;
for t = 2:100;
    N(:,t) = kmat*N(:,t-1) + Rvec;
end
N0 = N(:,end);


% Integrate to ensure comparing same size population
switch Method{j}
    case 'MidP'
      N00 = N0./(sum(N0*meshdiff));
    case 'Simp'
      N00 = N0./(Sy*N0);
end % end switch Method

% Plot some example sizes (color = # points, style = grid style)
if any(xSize(i) == xPlot)
    subplot(3,2,1:2)
    hold on
    plot(x,N0,'color',Col{xSize(i)==xPlot},'linestyle',Sty{j})
end

if any(xSize(i) == xPlot)
    subplot(3,2,3:4)
    hold on
    plot(x,N0,'color',Col{xSize(i)==xPlot},'linestyle',Sty{j})
end

% Get kernel w/o mortality (growth only)
switch Method{j}
    case 'MidP'
        kmat = kernmatSimp(x,F,fixparm_nm,1); % need to update kernmatMidP
        kmat = kmat*meshdiff;
    case 'Simp'
        kmat = kernmatSimp(x,F,fixparm_nm,1);
        Sy = makeSimpVec(meshdiff,meshsize);
        Symat = repmat(Sy(:)',[length(Sy),1]);
        kmat = Symat.*kmat; % do the integration
end % end switch Method

% Project using kernel
V2 = (kmat^10)*N00; % multiple generations?

% Integrate
switch Method{j}
    case 'MidP'
        Vsum = sum(V2.*meshdiff);
    case 'Simp'
        Vsum = Sy*V2;
end % end switch Method

% Store values
VS(i,j) = Vsum;
%VS2(i,j) = Vsum2;
Veig(i,j) = max(eig(kmat));
Md(i)= meshdiff;

% end loop
end % end loop over xSize
end % end loop over Method

% Plot initial part of the scale
subplot(3,2,5:6)
hold on
plot(xSize,(VS(:,1)),'k-')
plot(xSize,(VS(:,2)),'k--')

subplot(3,2,1:2)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
ylabel('Probability density')
xlabel('Size (cm)')
set(gca,'xlim',[0 60])

subplot(3,2,3:4)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
ylabel('Probability density')
xlabel('Size (cm)')
set(gca,'xlim',[5 15])

subplot(3,2,5:6)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
xlabel('Mesh size')
ylabel('Proportional loss')
ylim([0.8 1.2])
