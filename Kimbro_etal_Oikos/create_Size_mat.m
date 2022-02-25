function create_Size_mat

load SizevsAge.csv
Data = SizevsAge;
Site = Data(:,1);
Age = Data(:,2);
Mean = Data(:,3);
Std = Data(:,4);
Sizebins = 0:4:84; % need to add one more bin to the end, since the bins are defined by their right hand edge
                  

Uni_sites = unique(Site); % the unique site IDs
Uni_ages = unique(Age); % the unique ageclasses


Size_mat = ones(length(Uni_ages),length(Sizebins)-1,length(Uni_sites))*NaN;

    for i = 1:length(Uni_sites)% this should loop over the seven sites
        for j = 1:length(Uni_ages)% this should loop over each age or year (0-10)
            
            OKrow = Site == Uni_sites(i) & Age == Uni_ages(j); % find the row of data for the correct site & age
            mean_temp = Mean(OKrow); % get the correct mean
            sd_temp = Std(OKrow); % get the correct SD
            
           
            vecbins1 = Sizebins(1:(end-1)); % size bins from the first to the second-to-last
            vecbins2 = Sizebins(2:end); % size bins from the second to the last
            
            % same step you were attempting before, but vectorized somewhat
            bin_probs = normcdf(vecbins2,mean_temp,sd_temp) - normcdf(vecbins1,mean_temp,sd_temp);
            
            Size_mat(j,:,i) = bin_probs;
            
        end % end loop over ages
        
        % Now make sure they are all nonzero proportions
    
    Size_mat(Size_mat == 0) = realmin;
    RowSum = sum(Size_mat(:,:,i),2); % the sum along each row
    Size_mat(:,:,i) = Size_mat(:,:,i)./repmat(RowSum(:),[1,length(Sizebins)-1]); % divide by rowsum

    end % end loop over sites
    save Size_mat.mat Size_mat