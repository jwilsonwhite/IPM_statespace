%zone_error_ref.m
clear
savename1={'SMYS_Big_Creek_25Oct16_fwdfitSSD.mat','SMYS_White_Rock_25Oct16_fwdfitSSD.mat','SMYS_Pt_Lobos_25Oct16_fwdfitSSD.mat'};
savename3={'SMYS_Big_Creek_zone_data_sorted.mat','SMYS_White_Rock_zone_data_sorted.mat','SMYS_Pt_Lobos_zone_data_sorted.mat'};
savename4={'SMYS_Big_Creek_zone_error_ref.mat','SMYS_White_Rock_zone_error_ref.mat','SMYS_Pt_Lobos_zone_error_ref.mat'};
for ii=1:3
    load(savename1{ii})
    load(savename3{ii})
    %Get data in useable format
    Site_Names = Meta.Sites;
    Species_Names = Meta.Species;
    %get data into a histogram for site running projection based on
    N = IPM_histo_per_zone2(D_str.(Species_Names{1}),Years,Site_Names,edges);
    %N is in dimension size site year zone 
    if ii==1 %Big Creek - use site Esalen
        %ref_stat=logical([0 0 1 1]);
        rfst=3;
    elseif ii==2 %White Rock - use Harmony (4) or Estero (1)
        %ref_stat=logical([1 0 1 1]);
        rfst=4;
    elseif ii==3 %Pt Lobos - use Malpaso (2) or Soberanes (5)
        %ref_stat=logical([0 1 0 1 1 0]);
        rfst=2;
    end
    
    
    N_REF=squeeze(N(:,rfst,:,:));
    N_REFtrans=nan(size(N_REF));
        
        %figure out which years sampling occurred and how many months of
        %sampling occurred.
        
    if ii==1 %Big Creek
        yrkey=[1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1];
        %Esalen data 2001-2011 
        for jj = 1:length(Years)
            if yrkey(jj)==2
            NT = D_str.(Species_Names{1}).(Site_Names{rfst})(jj).data.transects_perzone;
            else
            NT=NaN;
            end
            for zn=1:length(NT)
            N_REFtrans(:,jj,zn) = N_REF(:,jj,zn)./NT(zn);
            end
        end
    elseif ii==2 %White Rock - Harmony
        yrkey=[1 1 1 1 1 1 1 1 3 2 2 1 2 1 1 1];
        %2007, sampled multiple months, 2008, 2009, 2011 normal
        for jj = 1:length(Years)
            if yrkey(jj)==1
                NT=NaN;
            elseif yrkey(jj)==2
                NT = D_str.(Species_Names{1}).(Site_Names{Meta.MPAnew})(jj).data.transects_perzone;
            elseif yrkey(jj)==3
                NT = D_str.(Species_Names{1}).(Site_Names{Meta.MPAnew})(jj).data.total_transects_perzone;
            end
            for zn=1:length(NT)
            N_REFtrans(:,jj,zn) = N_REF(:,jj,zn)./NT(zn);
            end
        end
    elseif ii==3 %Pt Lobos - Malpaso
        yrkey=[1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2];
        %data for 2006-2014
        for jj=1:length(Years)
            if yrkey(jj)==1
                NT=NaN;
            elseif yrkey(jj)==2
                NT = D_str.(Species_Names{1}).(Site_Names{Meta.MPAnew})(jj).data.transects_perzone;
            end
        
            for zn=1:length(NT)
            N_REFtrans(:,jj,zn) = N_REF(:,jj,zn)./NT(zn);
            end
        end
    end

%find fish open to fishery
Lfish = fixparm(5);
Lf_x=x(x>=Lfish);
Lfishv = x >= Lfish; % logical vector of fished sizes
%TOTAL ABUNDANCE > LFISH
N_pertran = squeeze(sum(N_REFtrans(Lfishv,:,:),1)); %total fish per site per year

Nmean_REF_data=nanmean(N_pertran');
Nstd_REF_data=std(N_pertran','omitnan');

save(savename4{ii},'Nmean_REF_data','Nstd_REF_data')
end