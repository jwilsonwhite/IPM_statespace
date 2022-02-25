function plot_IPM_fit(modelnum,Year)

fit_filename = strcat(['Tomales_IPM_fits/Tomales_Model_',num2str(modelnum),'_20Jan2018_fit.mat']);
post_filename = strcat(['Tomales_IPM_fits/postproc_Tomales_Model_',num2str(modelnum),'_20Jan2018.mat']);
metaname = 'Tomales_20Jan2018_metadata.mat';

load(fit_filename)
load(post_filename)
load(metaname)

dy = Meta.IPM.meshdiff;
Sy=makeSimpVec(dy,Meta.IPM.meshsize);

%Rfact = mc_str(1).Rfact;
isJuv = mc_str(1).isJuv;
N = mc_str(1).N;
Rvec = mc_str(1).R_vec;
Model_str = mc_str(1).Model_str;

Par.R = mean(post_str.R,1);
Par.aM = mean(post_str.aM,1);
Par.jM = mean(post_str.jM,1);
Par.pe = mean(post_str.pe,1);
Par.Growth = mc_str(1).Par(1).Growth;
%Par.oe = eps; %; mean(oe,1);


y = Year;
%figure(modelnum)
figure
clf
set(gcf,'units','cent','position',[30 10 12 20])
%Par.pe = 0.00;
L = 0;
%keyboard

% Site order (in the model):
% E1 E2 E3 E4 W2 W3 W4
% Site order for plotting
% [blank] E1; W2 E2; W3 E3; W4 E4;
Plot_order = [2, 4, 6, 8, 3, 5, 7];

% TO DO: Make the histogram an actual histogram
%        Make axes black not gray, etc

for s = 1:size(N,2)

[L1, Npred,Nact] = do_Tomales_IPM(dy,Sy,Par,Rvec(:,s),isJuv(:,s),squeeze(N(:,s,:)),Meta,Model_str,s,0);
subplot(4,2,Plot_order(s))
hold on
if any(isnan(Nact(:,y))) % if no data this year:
plot(Meta.IPM.x,Npred(:,y)*dy,'k--','linewidth',1) %   
else
%plot(Meta.IPM.x,Nact(:,y)./Meta.Sample_size(s,y),'k-','linewidth',1,'color',[0.5 0.5 0.5])
bar(Meta.IPM.x,Nact(:,y)./Meta.Sample_size(s,y),dy,'facecolor',[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5])
%plot(Meta.IPM.x,Npred(:,y)*Meta.Sample_size(s,y)*dy,'k--','linewidth',1) %
plot(Meta.IPM.x,Npred(:,y)*dy,'k--','linewidth',1)
end
if y == 5
ylim([0 0.3]);
set(gca,'ytick',0:0.05:0.5)
elseif y == 6
ylim([0 0.14]);
set(gca,'ytick',0:0.02:1)
else
ylim([0 0.1]);
set(gca,'ytick',0:0.02:1)    
end
xlim([0 75])
L = L + L1;
set(gca,'tickdir','out','ticklength',[0.02 0.02],'xtick',0:20:100)
set(gca,'xcolor','k','ycolor','k')

if Plot_order(s) >= 6
    xlabel('Oyster size (mm)','fontsize',14)
end

if Plot_order(s) == 3
    ylabel('Abundance (no. 0.01 m^-^2)','fontsize',14)
end
end

subplot(4,2,1)
title(strcat('Model',num2str(modelnum),' Year ',num2str(Year+1999)))


