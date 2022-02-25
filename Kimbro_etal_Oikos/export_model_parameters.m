function export_model_parameters

% Extract parameters from a desired model

load postproc_mcmc_results_16Aug2010_modelnum_8.mat

B = post_str.B;

        % Export data
        Mort_J = struct([]);
        Mort_A = struct([]);
        Rect = struct([]);
    
        Ind_j = j_ind(1:7);   
        Ind_a = a_ind(1:7);
        Ind_r = r_ind;
        Sites = {'W2','W3','W4','E1','E2','E3','E4'};
        Siteorder = [4, 5, 1, 6, 2, 7, 3];
        Years = 1999:2008;
        
        Ind_j = Ind_j(Siteorder);
        Ind_a = Ind_a(Siteorder);
        Sites = Sites(Siteorder);
        
        % Now just get the regional values:
        
     %   Ind_j = Ind_j([1 4 6]);
      %  Ind_a = Ind_a([1 4 6]);
        
        Mort_J(1).mean = mean(B(:,Ind_j));
        Mort_J(1).std = std(B(:,Ind_j));
        Mort_J(1).quant025 = quantile(B(:,Ind_j),0.025);
        Mort_J(1).quant975 = quantile(B(:,Ind_j),0.975);
        
        Mort_A(1).mean = mean(B(:,Ind_a));
        Mort_A(1).std = std(B(:,Ind_a));
        Mort_A(1).quant025 = quantile(B(:,Ind_a),0.025);
        Mort_A(1).quant975 = quantile(B(:,Ind_a),0.975);
        
        MeanR = mean(B(:,Ind_r));
        StdR = std(B(:,Ind_r));
        Quant025R = quantile(B(:,Ind_r),0.025);
        Quant975R = quantile(B(:,Ind_r),0.975);
        
        MeanR = reshape(MeanR,[7,10]);
        StdR = reshape(StdR,[7,10]);
        Quant025R = reshape( Quant025R,[7,10]);
         Quant975R = reshape( Quant975R,[7,10]);
        
        MeanR = MeanR(Siteorder,:);
        StdR = StdR(Siteorder,:);
        Quant025R = Quant025R(Siteorder,:);
        Quant975R = Quant975R(Siteorder,:);
        
      %  MeanR = MeanR([1, 4, 6],:);
      %  StdR = StdR([1, 4, 6],:);
      %  Quant025R = Quant025R([1, 4, 6],:);
      %  Quant975R = Quant975R([1, 4, 6],:);
        
        Rect(1).mean = MeanR;
        Rect(1).std = StdR;
        Rect(1).quant025 = Quant025R;
        Rect(1).quant975 = Quant975R;
        
        %keyboard
        
        Regions = {'Outer','Mid','Inner'};
        
        save Model8_summary Mort_J Mort_A Rect Years Regions
        
        

        
        
        
        