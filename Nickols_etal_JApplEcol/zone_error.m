%zone_error.m
clear
savename1={'SMYS_Big_Creek_25Oct16_fwdfitSSD.mat','SMYS_White_Rock_25Oct16_fwdfitSSD.mat','SMYS_Pt_Lobos_25Oct16_fwdfitSSD.mat'};
savename3={'SMYS_Big_Creek_zone_data_sorted.mat','SMYS_White_Rock_zone_data_sorted.mat','SMYS_Pt_Lobos_zone_data_sorted.mat'};
savename4={'SMYS_Big_Creek_zone_error.mat','SMYS_White_Rock_zone_error.mat','SMYS_Pt_Lobos_zone_error.mat'};
for ii=1:3;
    load(savename1{ii})
    load(savename3{ii})
    %Get data in useable format
    Site_Names = Meta.Sites;
    Species_Names = Meta.Species;
    %get data into a histogram for site running projection based on
    N = IPM_histo_per_zone2(D_str.(Species_Names{1}),Years,Site_Names,edges);
    N_MPA=squeeze(N(:,Meta.MPAnew,:,:));
    N_MPAtrans=nan(size(N_MPA));
    if ii==1 %Big Creek
        yrkey=[1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1];
        for jj = 1:length(Years)
            if jj==9
            NT = D_str.(Species_Names{1}).(Site_Names{Meta.MPAnew})(jj).data.transects_perzone;
            elseif jj==10
            NT = D_str.(Species_Names{1}).(Site_Names{Meta.MPAnew})(jj).data.transects_perzone;
            else
            NT=NaN;
            end
            for zn=1:length(NT)
            N_MPAtrans(:,jj,zn) = N_MPA(:,jj,zn)./NT(zn);
            end
        end
    elseif ii==2 %White Rock
        yrkey=[1 1 1 2 2 2 2 2 3 2 2 2 2 1 1 1];
        for jj = 1:length(Years)
            if yrkey(jj)==1
                NT=NaN;
            elseif yrkey(jj)==2
                NT = D_str.(Species_Names{1}).(Site_Names{Meta.MPAnew})(jj).data.transects_perzone;
            elseif yrkey(jj)==3
                NT = D_str.(Species_Names{1}).(Site_Names{Meta.MPAnew})(jj).data.total_transects_perzone;
            end
            for zn=1:length(NT)
            N_MPAtrans(:,jj,zn) = N_MPA(:,jj,zn)./NT(zn);
            end
        end
    elseif ii==3 %Pt Lobos
        yrkey=[2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
        for jj=1:length(Years)
        NT = D_str.(Species_Names{1}).(Site_Names{Meta.MPAnew})(jj).data.transects_perzone;
        
        for zn=1:length(NT)
        N_MPAtrans(:,jj,zn) = N_MPA(:,jj,zn)./NT(zn);
        end
        end
    end

%find fish open to fishery
Lfish = fixparm(5);
Lf_x=x(x>=Lfish);
Lfishv = x >= Lfish; % logical vector of fished sizes
%TOTAL ABUNDANCE > LFISH
N_pertran = squeeze(sum(N_MPAtrans(Lfishv,:,:),1)); %total fish per site per year

Nmean_MPA_data=nanmean(N_pertran');
Nstd_MPA_data=nanstd(N_pertran');

save(savename4{ii},'Nmean_MPA_data','Nstd_MPA_data')
end