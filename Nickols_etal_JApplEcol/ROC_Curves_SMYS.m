%ROC_Curves_SMYS.m
%make ROC curves for blue rockfish ms

%First step: Pick years to compare distributions
%yrs = [2 5 10 15 20];%but just run code for all years

clear

%SSD w/ recruitment variability
%output from rockfish_frwdproject_SSD.m
%pick which site loading
%load SMYS_White_Rock_07Sep16_fwdfitSSD.mat
load SMYS_Pt_Lobos_07Sep16_fwdfitSSD.mat
%load SMYS_Big_Creek_07Sep16_fwdfitSSD.mat

%add observation ogive to simulation runs with random recruitment
N_Ro=N_R.*repmat(OKlen',[1,T,RR]);
N_Fo=N_F.*repmat(OKlen',[1,T,RR]);

% % %so start with 2007 abundance
% N_Ro(:,1,:)=N_R(:,1,:);
% N_Fo(:,1,:)=N_F(:,1,:);

%sum up so get total of fish >= Lfish
N_R_sum=squeeze(sum(N_Ro(x>=Lfish,:,:),1));
N_F_sum=squeeze(sum(N_Fo(x>=Lfish,:,:),1));
SimRR_R_SSD=N_R_sum';
SimRR_F_SSD=N_F_sum';

%compare distributions for different years 
fp_SimRR_SSD=NaN(1000,T);
tp_SimRR_SSD=NaN(1000,T);
for i=1:T
    [fp_SimRR_SSD(:,i),tp_SimRR_SSD(:,i)] = ROC(SimRR_F_SSD(:,i),SimRR_R_SSD(:,i));
end

%Actual conditions (2007 monitoring data as start) w/ recruitment variability
%output from rockfish_frwdproject.m
%pick one for site that working on
%load SMYS_White_Rock_07Sep2016_frwdfit.mat
load SMYS_Pt_Lobos_07Sep2016_frwdfit.mat
%load SMYS_Big_Creek_07Sep2016_frwdfit.mat

%NOTE T goes to 20 for simulation runs

%add observation ogive to simulation runs with random recruitment
N_Ro=N_R.*repmat(OKlen',[1,T,RR]);
N_Fo=N_F.*repmat(OKlen',[1,T,RR]);

%sum up so get total of fish >= Lfish
N_R_sum=squeeze(sum(N_Ro(x>=Lfish,:,:),1));
N_F_sum=squeeze(sum(N_Fo(x>=Lfish,:,:),1));
SimRR_R=N_R_sum';
SimRR_F=N_F_sum';

%compare distributions for different years 
fp_SimRR=NaN(1000,T);
tp_SimRR=NaN(1000,T);
for i=1:T
    [fp_SimRR(:,i),tp_SimRR(:,i)] = ROC(SimRR_F(:,i),SimRR_R(:,i));
end

save SMYS_Pt_Lobos_07Sep16_ROC.mat %fp_DetRR tp_DetRR fp_SimRR tp_SimRR fp_SimRR_SSD tp_SimRR_SSD
%make file for each site