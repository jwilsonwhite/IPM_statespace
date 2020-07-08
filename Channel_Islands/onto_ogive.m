function onto_ogive

% create ogive describing ontogenetic migration of adult blue rockfish
% based on Rick Starr's rocky reef & PISCO kelp forest data
% for use with fitting models to PISCO data

% Load in PISCO data
load('SMYS_Pt_Lobos_data_sorted.mat')

% compile into a single histogram, focusing on post-2007 data
Site_Names = fieldnames(D_str.SMYS);
Count = [];
TL = [];
for s = 1:length(Site_Names)
    for y = 9:12
        Count = [Count; D_str.SMYS.(Site_Names{s})(y).data.Count(:)];
        TL = [Count; D_str.SMYS.(Site_Names{s})(y).data.TL(:)];
    end
end

for i = 1:length(Count)
    if Count(i) > 1
        TLtmp = repmat(TL(i),[Count(i),1]);
        TL = [TL; TLtmp(:)];
    end
end


% Load Starr data
RS = importdata('Starr_hookline_survey_data/Point_Lobos_Blues_Starr.csv');
% This file only holds data on S mystinus from Pt Lobos area sites
RS = RS.data;
RSstr(1).Y  = RS(:,end);
RSstr(1).L = RS(:,1);

% Strip out anomalously large observations
TL = TL(TL <= 53);
RSstr.L = RSstr.L(RSstr.L <= 53);

% Make into histograms:
Edges = 1:1:60;
Phist = histc(TL,Edges);
Rhist = histc(RSstr.L,Edges);

% convert into size-frequency
Phist = Phist./sum(Phist);
Rhist = Rhist./sum(Rhist);

% combine into a single dataset
N = length(TL) + length(RSstr.L);
ID = [zeros(length(TL),1); ones(length(RSstr.L),1)];
L = [TL; RSstr.L];

% Logistic regression
[b,dev,stats] = glmfit(L,ID,'binomial'); % default is logit link

% Back-transform to get ogive
LL = 0:60;
p = 1./(1+1./exp(b(1)+b(2).*LL));
%p = 1-p;

N = 1-normcdf(LL,28,7);

% Make a plot
figure
clf
axes;
ah(1) = gca;
hold on
plot(Edges,Phist,'b')
plot(Edges,Rhist,'r')
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xgrid','on','ygrid','on')

axes;
ah(2) = gca;
hold on
set(ah(2),'position',get(ah(1),'position'));
set(ah(2),'color','none','yaxislocation','right')
plot(LL,1-p,'k')
plot(LL,N,'k--')
set(gca,'tickdir','out','ticklength',[0.02 0.02])
keyboard



