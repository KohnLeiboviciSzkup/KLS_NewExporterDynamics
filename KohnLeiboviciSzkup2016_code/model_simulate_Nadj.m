function [r sim] = model_simulate_Nadj(m,s,r)

    s.flag_short=0;
    s.flag_extraresults=1;

%% Simulate productivity shocks

%Initalize objects
    ind_z = zeros(s.N,s.T*(1+s.burn));
    
%Initialize productivity realizations at stationary distribution
%Input -> z_pi, stationary distribution of productivity shocks, vector of size [s.z_grid_size x 1], computed by approx_shocks.m when run by model_solve.m
%Output -> ind_z(:,1), vector of indexes across productivity states for each firm, size [s.N x 1]
    k=1;
    for i=1:s.z_grid_size %length(r.log_z_grid)
        v(k:k+floor(r.z_pi(i)*s.N)-1)=i;  
        k=length(v)+1;
    end
    v(length(v):s.N)=round(mean(v));
    ind_z(:,1)=v';
    clear v;

%Compute productivity realizations after initial state
%Input -> z.P, transition matrix computed by approx_shocks.m
%Output -> ind_z(:,2:end)
    C = cumsum(r.z_P,2);
    R = rand(s.N,s.T*(1+s.burn));
    for j=2:s.T*(1+s.burn);
        f = repmat(R(:,j),1,s.z_grid_size);
        ind_z(:,j) = 1+sum( f > C(ind_z(:,j-1),:) ,2);
        clear f;
    end
    clear R C;
    ind_z = ind_z';

    
    
%Alternative way of simulating shocks (based on compecon function ddpsimul.m)
    %ind_z = ddpsimul(r.z_P,v,s.T*(1+s.burn)-1)';
        
if s.flag_short==1

        % Clearing variables soution from dynamic problem that won't be
        % essential
        

    r.y_prime=[];
    r.log_z_grid=[];
    r.z_P=[];
    r.z_pi=[];
    r.constrained=[];
    r.constrained_d=[];
    r.constrained_f=[];
    r.v=[];
    r.pi_d=[];
    r.pi_export=[];
        
        
%Initialize objects
    ind_y = zeros(s.T*(1+s.burn)+1,s.N);
    exporter = zeros(s.T*(1+s.burn)+1,s.N);   
        
    
%First period initialization
    ind_y(1,1:s.N) = round(s.y_grid_size/4); %We initialize every firm with lowest assets
    
%Compute asset and export states for t>1    
    for t=1:s.T*(1+s.burn)  
               
         ind_k = exporter(t,:)+1; %k=1 is new exporter, k=2 is continuing exporter

       
        %Construct 3-dimensional indexes from unidimensional indexes ind_y, ind_z, and ind_k     
        %The policy functions are multidimensional matrices -> we map ind_y, ind_z, and ind_k to the policy functions using ind        
            ind = s.y_grid_size*s.z_grid_size*(ind_k-1)+s.y_grid_size*(ind_z(t,:)-1)+ind_y(t,:);
        
        %Update state variables    
            exporter(t+1,:) = r.exporter(ind); % r.exporter is a decision rule given states today, sim.exporter updates simulation   
            ind_y(t+1,:)= r.y_prime_ind(ind);
            
    end

    
    %sim.y=r.y_grid(ind_y); % Labor values of the state
    exporter_decision = exporter(2:end,:); %Export decisions (as opposed to status)
    ind_y_decision = ind_y(2:end,:);
    sim.y = r.y_grid(ind_y_decision); %Labor values chosen by the firm for current production
    clear ind 
    r.y_grid=[];
       
    xhistory = zeros(s.T*(1+s.burn),s.N);
    xhistory(1,:) = exporter_decision(1,:); %First period assume all exporters are new
    xhistory(1,exporter_decision(1,:)==0) = -1;
    for t = 2:s.T*(1+s.burn)
        for j=1:s.N
           %If choose to export today
           if exporter_decision(t,j)==1 && exporter_decision(t-1,j)==0 %If new exporter
              xhistory(t,j) = 1;
           elseif exporter_decision(t,j)==1 && exporter_decision(t-1,j)==1 %If continuing exporter
                xhistory(t,j) = xhistory(t-1,j)+1;
                
             %If choose not to export today
           elseif exporter_decision(t,j)==0 && exporter_decision(t-1,j)==1 %If new non-exporter
                xhistory(t,j) = -1;
           elseif exporter_decision(t,j)==0 && exporter_decision(t-1,j)==0 %If continuing non-exporter
                xhistory(t,j) = xhistory(t-1,j) - 1;
           end
        end
    end
       
    % Burn first s.burn*s.T time periods
    ind_y = ind_y(s.T*s.burn+1:end-1,:);
    %sim.y = sim.y(s.T*s.burn+1:end-1,:);
    sim.y = sim.y(s.T*s.burn+1:end,:);
    ind_z = ind_z(s.T*s.burn+1:end,:);
    exporter = exporter(s.T*s.burn+1:end-1,:);
    exporter_decision = exporter_decision(s.T*s.burn+1:end,:);
    ind_y_decision = ind_y_decision(s.T*s.burn+1:end,:);
    xhistory = xhistory(s.T*s.burn+1:end,:);


    ind_k = exporter+1; %k=1 if new exporter, k=2 if continuing exporter

    
    
    %Initialize objects
    q_d = zeros(s.T,s.N);
    p_d = zeros(s.T,s.N);
    sim.n_f = zeros(s.T,s.N);
    sim.n_d = zeros(s.T,s.N);
    
    ind=s.y_grid_size*s.z_grid_size*(ind_k-1)+s.y_grid_size*(ind_z-1)+ind_y;
    ind_decision=s.y_grid_size*s.z_grid_size*(ind_k-1)+s.y_grid_size*(ind_z-1)+ind_y_decision;
    
    %clear ind_y ind_k 
    
    %Productivity     
    z = r.z_grid(ind_z);
    
    r.z_grid=[];
    %clear ind_z
%Quantities            
    q_f = r.q_f_export(ind_decision).*exporter_decision;
    q_d(exporter_decision==0) = r.q_d(ind_decision(exporter_decision==0));
    q_d(exporter_decision==1) = r.q_d_export(ind_decision(exporter_decision==1));
       
%Prices     
    p_f = r.p_f_export(ind_decision).*exporter_decision;
    p_d(exporter_decision==0) = r.p_d(ind_decision(exporter_decision==0));
    p_d(exporter_decision==1) = r.p_d_export(ind_decision(exporter_decision==1));

%Sales    
    sim.n_f = r.n_f_export(ind_decision).*exporter_decision;
    sim.n_d(exporter_decision==0) = r.n_d(ind_decision(exporter_decision==0));
    sim.n_d(exporter_decision==1) = r.n_d_export(ind_decision(exporter_decision==1));        
    sim.n = sim.n_f + sim.n_d;
    sim.n(xhistory==1) = sim.n(xhistory==1) + m.F + m.S;
    sim.n(xhistory>1) = sim.n(xhistory>1) + m.F;       
    r.n_d=[];
    r.n_d_export=[];
    r.n_f_export=[];
    
     
    %Sales     
    sales_f = p_f.*q_f;
    sales_d = p_d.*q_d;
    sales = sales_d+sales_f;
    
    r.p_f_export=[];
    r.p_d_export=[];
    r.p_d=[];
    r.q_f_export=[];
    r.q_d_export=[];
    r.q_d=[];
    
    clear p_f p_d 
    
    sales_growth   = [(sales_d(2:end,:)+sales_f(2:end,:)-sales_d(1:end-1,:)-sales_f(1:end-1,:))./(sales_d(1:end-1,:)+sales_f(1:end-1,:))];
    sales_d_growth = [(sales_d(2:end,:)-sales_d(1:end-1,:))./sales_d(1:end-1,:)];
    sales_f_growth = [(sales_f(2:end,:)-sales_f(1:end-1,:))./sales_f(1:end-1,:)];
    log_sales_growth = [log(sales(2:end,:))-log(sales(1:end-1,:))];    
    
    sim.sales_growth_4yravg = (1/4)*[(sales_d(4:end,:)+sales_f(4:end,:)-sales_d(1:end-3,:)-sales_f(1:end-3,:))./(sales_d(1:end-3,:)+sales_f(1:end-3,:))];  
    
    BJ_length = 3;
    for t=1:1 %Starts at 1, and can go up to s.T-BJ_length (its robust to these changes)
       sim.BJstart(t,:) = xhistory(t,:)<0 & xhistory(t+BJ_length,:)>0; 
       sim.BJstop(t,:) = xhistory(t,:)>0 & xhistory(t+BJ_length,:)<0; 
       sim.BJboth(t,:) = xhistory(t,:)>0 & xhistory(t+BJ_length,:)>0; 
       sim.BJneither(t,:) = xhistory(t,:)<0 & xhistory(t+BJ_length,:)<0; 
    end
    
    sim.sales_growth_4yravg = sim.sales_growth_4yravg(size(sim.BJstart,1),:);    
    
    sim.BJstart_sales_gr = median(sim.sales_growth_4yravg(sim.BJstart==1));
    sim.BJstop_sales_gr = median(sim.sales_growth_4yravg(sim.BJstop==1));
    sim.BJboth_sales_gr = median(sim.sales_growth_4yravg(sim.BJboth==1));
    sim.BJneither_sales_gr = median(sim.sales_growth_4yravg(sim.BJneither==1));
        
    
    %Export intensity     
    export_intensity = sales_f./(sales_d+sales_f);     
   
    %Firm-level statistics
    %share of exporters
    sim.share_exporters = mean(sum(exporter_decision,2)/s.N);    
                
     %   clear exporter_decision;
        
    %Export entry and exit rates
    sim.share_starters = mean(sum(xhistory(2:end,:)==1,2)./sum(xhistory(1:end-1,:)<=-1,2));
    sim.share_stoppers = mean(sum(xhistory(2:end,:)==-1,2)./sum(xhistory(1:end-1,:)>=1,2));    
        
    %Export intensity
    temp = export_intensity(xhistory>=1);
    sim.export_sales_ratio_med =median(temp(:));        
        
    clear temp 
    %clear export_intensity
        
    %Share of exporters that are new
    sim.sharex_new = sum(sum(xhistory==1))/sum(sum(xhistory>=1));
    
    %Exporter size premium -> exporters vs non-exporters
        %Labor            
            tempS = zeros(size(q_f));
            tempS(xhistory==1) = 1;
            temp1 = ((m.tau*q_f(xhistory>=1)+q_d(xhistory>=1))./z(xhistory>=1))+m.F+tempS(xhistory>=1)*m.S;
            temp2 = q_d(xhistory<=-1)./z(xhistory<=-1);
            sim.export_size_premium_med = median(temp1(:))/median(temp2(:));   
            clear temp1 temp2 tempS    
            
    %Continuing exporter size premium -> continuing exporters vs new exporters
        %Labor
            temp1 = ((m.tau*q_f(xhistory>1)+q_d(xhistory>1))./z(xhistory>1))+m.F;            
            temp2 = ((m.tau*q_f(xhistory==1)+q_d(xhistory==1))./z(xhistory==1))+m.F+m.S;
            sim.cont_export_size_premium_med = median(temp1(:))/median(temp2(:));  
            clear temp1 temp2   
    
% Continuing non-exporter size premium -> continuing non-exporters vs new
    % non-exporters (stoppers)
    
        %Labor
            temp1 = q_d(xhistory<-1)./z(xhistory<-1);            
            temp2 = q_d(xhistory==-1)./z(xhistory==-1);
            sim.cont_nonexport_size_premium_med = median(temp1(:))/median(temp2(:));  
            clear temp1 temp2 
            
    % New non-exporter size premium -> new non-exporters vs exporters
    
        %Labor
            tempS = zeros(size(q_f));
            tempS(xhistory==1) = 1;
            temp1 = ((m.tau*q_f(xhistory>=1)+q_d(xhistory>=1))./z(xhistory>=1))+m.F+tempS(xhistory>=1)*m.S;
            temp2 = q_d(xhistory==-1)./z(xhistory==-1);
            sim.new_nonexport_size_premium_med = median(temp2(:))/median(temp1(:));  
            clear tempS temp1 temp2 

    %New exporters / stoppers 
        tempS = zeros(size(q_f));
        tempS(xhistory==1) = 1;
        temp1 = ((m.tau*q_f(xhistory==1)+q_d(xhistory==1))./z(xhistory==1))+m.F+tempS(xhistory==1)*m.S;
        temp2 = (q_d(xhistory==-1))./z(xhistory==-1);
        sim.new_stoppers_size_premium_med = median(temp1(:))/median(temp2(:));  
        clear tempS temp1 temp2         
        
   % Statistics for cohorts of firms that didn't export for s.CD years, and exported for s.CN consecutive years 

    %Compute cohorts
    xperiods_cohorts = zeros(s.T,s.N);
    for k=1:s.T-s.CN-s.CD %Cohorts
        xperiods_exporters = xhistory(k+s.CD+s.CN-1,:)==s.CN;
        xperiods_cohorts(k:k+s.CD+s.CN-1,xperiods_exporters) = 1;
    end
    xperiods_xhistory = xperiods_cohorts.*xhistory; %First value will be any negative number
    xperiods_xhistory_gr = xperiods_xhistory(2:end,:); %Get rid of first period in the panel since cannot compute growth rates (previous obs doesn't exist)
  
    clear xperiods_exporters;
    
%Compute statistics -> Averages and medians across cohorts, by firm age within the cohort
for i = 1:s.CN    
    %Median
        %Level
            sim.xperiods_sales_f_med(i,1) = median(sales_f(xperiods_xhistory==i));
            sim.xperiods_sales_d_med(i,1) = median(sales_d(xperiods_xhistory==i));
            sim.xperiods_sales_med(i,1) = median(sales_d(xperiods_xhistory==i)+sales_f(xperiods_xhistory==i));
            sim.xperiods_export_intensity_med(i,1) = median(export_intensity(xperiods_xhistory==i));
           
        %Growth       
            sim.xperiods_sales_f_growth_med(i,1) = median(sales_f_growth(xperiods_xhistory_gr==i));
            sim.xperiods_sales_d_growth_med(i,1) = median(sales_d_growth(xperiods_xhistory_gr==i));    
            sim.xperiods_sales_growth_med(i,1) = median(sales_growth(xperiods_xhistory_gr==i));        
       
end
clear export_intensity 

sim.xperiods_sales_d_growth_med_before_x = median(sales_d_growth(xperiods_xhistory_gr<0));
sim.xperiods_sales_growth_med_before_x = median(sales_growth(xperiods_xhistory_gr<0));

sim.xperiods_sales_d_growth_med = [sim.xperiods_sales_d_growth_med_before_x ;sim.xperiods_sales_d_growth_med];
sim.xperiods_sales_growth_med = [sim.xperiods_sales_growth_med_before_x ;sim.xperiods_sales_growth_med];
   
clear xperiods_xhistory xperiods_cohorts; 

  % Feb 13th, corrected typo: we should count only firms at t if they
   % are counted at t-1, and firms at t if they can be counted at t+1.

    %Hazard rate -> computed by pooling firms from all cohorts
    xperiod_count_pool1 = zeros(1,6); % Counts only firms at t, if they are counted at t-1
    xperiod_count_pool2 = zeros(1,6); % Counts only firms at t, if they can be counted at t+1
    
    for i=1:6     
       
        xperiod_count_pool1(i) = sum(sum(xhistory(2:end,:)==i)); % Counts only firms at t, if they are counted at t-1
        xperiod_count_pool2(i) = sum(sum(xhistory(1:end-1,:)==i)); % Counts only firms at t, if they can be counted at t+1
        
    end
    
    
    for i=1:5
    
        sim.hazard_pool(i) = (xperiod_count_pool2(i)-xperiod_count_pool1(i+1))/xperiod_count_pool2(i);
    
    end
   
    clear xperiod_count_pool1  xperiod_count_pool2  
    
    
    % Statistics for cohorts of firms that exported for at least s.CD years, and did not export for s.CN consecutive years 

    %Compute cohorts
    nonxperiods_cohorts = zeros(s.T,s.N);
    for k=1:s.T-s.CN-s.CD %Cohorts
        nonxperiods_exporters = xhistory(k+s.CD+s.CN-1,:)==-s.CN;
        nonxperiods_cohorts(k:k+s.CD+s.CN-1,nonxperiods_exporters) = 1;
    end
    nonxperiods_xhistory = nonxperiods_cohorts.*xhistory; %First value will be any negative number
    nonxperiods_xhistory_gr = nonxperiods_xhistory(2:end,:); %Same as above for xperiods_xhistory_gr
    
    %Compute statistics -> Averages and medians across cohorts, by firm age within the cohort
    for i = -1:-1:-s.CN    
    %Median
        %Level
            sim.nonxperiods_sales_d_med(-i,1) = median(sales_d(nonxperiods_xhistory==i));
            sim.nonxperiods_sales_med(-i,1) = median(sales_d(nonxperiods_xhistory==i)+sales_f(nonxperiods_xhistory==i));
           
        %Growth       
            sim.nonxperiods_sales_d_growth_med(-i,1) = median(sales_d_growth(nonxperiods_xhistory_gr==i));    
            sim.nonxperiods_sales_growth_med(-i,1) = median(sales_growth(nonxperiods_xhistory_gr==i));        
       
    end

    sim.nonxperiods_sales_d_growth_med_before_x = median(sales_d_growth(nonxperiods_xhistory_gr>0));
    sim.nonxperiods_sales_growth_med_before_x = median(sales_growth(nonxperiods_xhistory_gr>0));

    
    sim.nonxperiods_sales_d_growth_med = [sim.nonxperiods_sales_d_growth_med_before_x ;sim.nonxperiods_sales_d_growth_med];
    sim.nonxperiods_sales_growth_med = [sim.nonxperiods_sales_growth_med_before_x ;sim.nonxperiods_sales_growth_med];


    % hazard for stoppers:
%     nonxperiod_count_pool = zeros(1,6);
%     for i=-1:-1:-6     
%         nonxperiod_count_pool(-i) = sum(sum(xhistory==i));
%     end
%     for i=-1:-1:-5
%         sim.nonxhazard_pool(-i) = (nonxperiod_count_pool(-i)-nonxperiod_count_pool(-i+1))/nonxperiod_count_pool(-i);
%     end    
%     
%     clear nonxperiod_count_pool;
    
     % Feb 13th, corrected typo: we should count only firms at t if they
   % are counted at t-1, and firms at t if they can be counted at t+1.

    %Hazard rate -> computed by pooling firms from all cohorts
    nonxperiod_count_pool1 = zeros(1,6); % Counts only firms at t, if they are counted at t-1
    nonxperiod_count_pool2 = zeros(1,6); % Counts only firms at t, if they can be counted at t+1
    
    for i=1:6      
       
        nonxperiod_count_pool1(i) = sum(sum(xhistory(2:end,:)==-i)); % Counts only firms at t, if they are counted at t-1
        nonxperiod_count_pool2(i) = sum(sum(xhistory(1:end-1,:)==-i)); % Counts only firms at t, if they can be counted at t+1
        
    end
    
    
    for i=1:5
    
        sim.nonxhazard_pool(i) = (nonxperiod_count_pool2(i)-nonxperiod_count_pool1(i+1))/nonxperiod_count_pool2(i);
    
    end
   
    clear nonxperiod_count_pool1  nonxperiod_count_pool2   

%Labor statistics
    wagebill = m.w*sim.n;
    sim.wagebill_acf_mat = zeros(2,s.N);
    for i=1:s.N
        sim.wagebill_acf_mat(:,i) = autocorr(log(wagebill(:,i)),1); 
    end
    sim.wagebill_acf_med = median(sim.wagebill_acf_mat(2,:));
    sim.wagebill_acf_mean = mean(sim.wagebill_acf_mat(2,:));
    
%Sales statistics   

    % calculating the median of sales statistics for each firm (within)
    sales_var = sales;
    sim.sales_acf_mat = zeros(2,s.N);
    sim.sales_std_vec = zeros(s.N,1);
    sim.sales_mean_vec = zeros(s.N,1);
    for i=1:s.N
        sim.sales_acf_mat(:,i) = autocorr(log(sales_var(:,i)),1); 
        sim.sales_std_vec(i) = std(log(sales_var(:,i)));
        sim.sales_mean_vec(i) = mean(log(sales_var(:,i)));
    end
    sim.sales_acf = median(sim.sales_acf_mat(2,:));
    sim.sales_std = median(sim.sales_std_vec);
    sim.sales_mean = median(sim.sales_mean_vec);
    sim.sales_cv = sim.sales_std/sim.sales_mean;        
    
    % calculating the median of sales statistics by pooling all firms
    sales_corr_input = zeros(s.N*(s.T-1),2);
    for i=1:s.N
        sales_corr_input((i-1)*(s.T-1)+1:i*(s.T-1),:) = [log(sales_var(1:end-1,i)) log(sales_var(2:end,i))];
    end    
    sim.sales_acf2 = corr(sales_corr_input);    
    sim.sales_std2 = std(log(sales_var(:)));
    sim.sales_mean2 = mean(log(sales_var(:)));
    sim.sales_cv2 = sim.sales_std2/sim.sales_mean2;   
        
    sales_y = sales(2:end,:);
    sales_x = sales(1:end-1,:);
    sim.sales_ar1_reg = regress(log(sales_y(:)),[ones(length(sales_y(:)),1) log(sales_x(:))]);
    
%Sales growth statistics   

    % calculating the median of sales growth statistics for each firm (within)
    sales_var = log_sales_growth;
    sim.sales_growth_acf_mat = zeros(2,s.N);
    sim.sales_growth_std_vec = zeros(s.N,1);
    for i=1:s.N
        sim.sales_growth_acf_mat(:,i) = autocorr((sales_var(:,i)),1); 
        sim.sales_growth_std_vec(i) = std(sales_var(:,i));
        sim.sales_growth_mean_vec(i) = mean(sales_var(:,i));
    end
    sim.sales_growth_acf = median(sim.sales_growth_acf_mat(2,:));
    sim.sales_growth_std = median(sim.sales_growth_std_vec);
    sim.sales_growth_mean = median(sim.sales_growth_mean_vec);
    sim.sales_growth_cv = sim.sales_growth_std/sim.sales_growth_mean; 
    
    % calculating the median of sales growth statistics by pooling all firms
    sales_growth_corr_input = zeros(s.N*(s.T-2),2);
    for i=1:s.N
        sales_growth_corr_input((i-1)*(s.T-2)+1:i*(s.T-2),:) = [(sales_var(1:end-1,i)) (sales_var(2:end,i))];
    end    
    sim.sales_growth_acf2 = corr(sales_growth_corr_input);    
    sim.sales_growth_std2 = std((sales_var(:)));   
    sim.sales_growth_mean2 = mean(sales_var(:));
    sim.sales_growth_cv2 = sim.sales_growth_std2/sim.sales_growth_mean2; 
        
%Sales size premium    
    %Exporters vs non-exporters 
        temp1 = (sales_d(xhistory>=1)+sales_f(xhistory>=1));
        temp2 = sales_d(xhistory<=-1);
        sim.export_size_premium_sales_med = median(temp1(:))/median(temp2(:));          
        clear temp1 temp2    

    %Continuing exporters vs new exporters
        temp1 = (sales_d(xhistory>1)+sales_f(xhistory>1));
        temp2 = (sales_d(xhistory==1)+sales_f(xhistory==1));
        sim.export_size_premium_cont_sales_med = median(temp1(:))/median(temp2(:));   
        clear temp1 temp2      
             
        
%Sales distribution with 3 bins    
    b = 0.33; %Size of the middle bin
    lb = (1-b)/2;
    ub = (1+b)/2;
    
    sales_sort = sort((sales(:)));
    sales_small = median(sales_sort(1:round(lb*s.T*s.N))); 
    sales_medium = median(sales_sort(round(lb*s.T*s.N):round(ub*s.T*s.N))); 
    sales_large = median(sales_sort(round(ub*s.T*s.N):end)); 
    
    sim.largemedium_size_premium_med = sales_large/sales_medium;
    sim.mediumsmall_size_premium_med = sales_medium/sales_small;
    
        
%Labor size distribution
%RECHECK VALUES USED WITH THOSE IN THE DATA
    sim.laborsize_psmall = prctile(sim.n(:),41.22);
    sim.laborsize_pmed = prctile(sim.n(:),83.34);    
    sim.laborsize = (sim.n<=sim.laborsize_psmall) + 2*(sim.n>sim.laborsize_psmall & sim.n<=sim.laborsize_pmed) + 3*(sim.n>sim.laborsize_pmed);
         
%Size-Status Distribution: for each size, distribution of export status
    %Status-Size Distribution: for each export status, distribution of size
        sim.size_distribution=zeros(3,3);
        sim.size_distribution(1,1)=sum(sum(sim.laborsize==1 & xhistory<0))/sum(sum(xhistory<0)); %Non-exporters
        sim.size_distribution(1,2)=sum(sum(sim.laborsize==2 & xhistory<0))/sum(sum(xhistory<0));
        sim.size_distribution(1,3)=sum(sum(sim.laborsize==3 & xhistory<0))/sum(sum(xhistory<0));
        sim.size_distribution(2,1)=sum(sum(sim.laborsize==1 & xhistory==1))/sum(sum(xhistory==1)); %New exporters
        sim.size_distribution(2,2)=sum(sum(sim.laborsize==2 & xhistory==1))/sum(sum(xhistory==1));
        sim.size_distribution(2,3)=sum(sum(sim.laborsize==3 & xhistory==1))/sum(sum(xhistory==1));
        sim.size_distribution(3,1)=sum(sum(sim.laborsize==1 & xhistory>1))/sum(sum(xhistory>1)); %Cont exporters
        sim.size_distribution(3,2)=sum(sum(sim.laborsize==2 & xhistory>1))/sum(sum(xhistory>1));
        sim.size_distribution(3,3)=sum(sum(sim.laborsize==3 & xhistory>1))/sum(sum(xhistory>1));
        
        sim.size_distribution2=zeros(2,3);
        sim.size_distribution2(1,1)=sum(sum(sim.laborsize==1 & xhistory<0))/sum(sum(xhistory<0)); %Non-exporters
        sim.size_distribution2(1,2)=sum(sum(sim.laborsize==2 & xhistory<0))/sum(sum(xhistory<0));
        sim.size_distribution2(1,3)=sum(sum(sim.laborsize==3 & xhistory<0))/sum(sum(xhistory<0));
        sim.size_distribution2(2,1)=sum(sum(sim.laborsize==1 & xhistory>=1))/sum(sum(xhistory>=1)); %Exporters
        sim.size_distribution2(2,2)=sum(sum(sim.laborsize==2 & xhistory>=1))/sum(sum(xhistory>=1));
        sim.size_distribution2(2,3)=sum(sum(sim.laborsize==3 & xhistory>=1))/sum(sum(xhistory>=1));      
        
    %Size-Status Distribution: for each size, distribution of export status
    %Similar to Figure 1d in Alessandria and Choi
        sim.status_distribution=zeros(3,3);
        sim.status_distribution(1,1)=sum(sum(sim.laborsize==1 & xhistory<0))/sum(sum(sim.laborsize==1)); %Small 
        sim.status_distribution(1,2)=sum(sum(sim.laborsize==1 & xhistory==1))/sum(sum(sim.laborsize==1));
        sim.status_distribution(1,3)=sum(sum(sim.laborsize==1 & xhistory>1))/sum(sum(sim.laborsize==1));
        sim.status_distribution(2,1)=sum(sum(sim.laborsize==2 & xhistory<0))/sum(sum(sim.laborsize==2)); %Medium
        sim.status_distribution(2,2)=sum(sum(sim.laborsize==2 & xhistory==1))/sum(sum(sim.laborsize==2));
        sim.status_distribution(2,3)=sum(sum(sim.laborsize==2 & xhistory>1))/sum(sum(sim.laborsize==2));
        sim.status_distribution(3,1)=sum(sum(sim.laborsize==3 & xhistory<0))/sum(sum(sim.laborsize==3)); %Large
        sim.status_distribution(3,2)=sum(sum(sim.laborsize==3 & xhistory==1))/sum(sum(sim.laborsize==3));
        sim.status_distribution(3,3)=sum(sum(sim.laborsize==3 & xhistory>1))/sum(sum(sim.laborsize==3));
    
        sim.laborsize_largemed_mean = mean(sim.n(sim.laborsize==3))/mean(sim.n(sim.laborsize==2));
        sim.laborsize_medsmall_mean = mean(sim.n(sim.laborsize==2))/mean(sim.n(sim.laborsize==1));
        sim.laborsize_largemed_med = median(sim.n(sim.laborsize==3))/median(sim.n(sim.laborsize==2));
        sim.laborsize_medsmall_med = median(sim.n(sim.laborsize==2))/median(sim.n(sim.laborsize==1));

        %This is exactly as in Alessandria and Choi Fig 1d: share of exporters by size
            sim.status_distribution2=zeros(3,2);
            sim.status_distribution2(1,1)=sum(sum(sim.laborsize==1 & xhistory<0))/sum(sum(sim.laborsize==1)); %Small 
            sim.status_distribution2(1,2)=sum(sum(sim.laborsize==1 & xhistory>=1))/sum(sum(sim.laborsize==1));
            sim.status_distribution2(2,1)=sum(sum(sim.laborsize==2 & xhistory<0))/sum(sum(sim.laborsize==2)); %Medium
            sim.status_distribution2(2,2)=sum(sum(sim.laborsize==2 & xhistory>=1))/sum(sum(sim.laborsize==2));
            sim.status_distribution2(3,1)=sum(sum(sim.laborsize==3 & xhistory<0))/sum(sum(sim.laborsize==3)); %Large
            sim.status_distribution2(3,2)=sum(sum(sim.laborsize==3 & xhistory>=1))/sum(sum(sim.laborsize==3));       

    %Alternative labor size distribution
    %Labor size distribution
        sim.ls_size1 = prctile(sim.n(:),29.74);
        sim.ls_size2 = prctile(sim.n(:),56.33);   
        sim.ls_size3 = prctile(sim.n(:),83.34);   
        sim.ls = (sim.n<=sim.ls_size1) + 2*(sim.n>sim.ls_size1 & sim.n<=sim.ls_size2) + 3*(sim.n>sim.ls_size2 & sim.n<=sim.ls_size3) + 4*(sim.n>sim.ls_size3);
         
        %Size-Status Distribution: for each size, distribution of export status
        %Status-Size Distribution: for each export status, distribution of size
        sim.size_distribution3=zeros(2,4);
        sim.size_distribution3(1,1)=sum(sum(sim.ls==1 & xhistory<0))/sum(sum(xhistory<0)); %Non-exporters
        sim.size_distribution3(1,2)=sum(sum(sim.ls==2 & xhistory<0))/sum(sum(xhistory<0));
        sim.size_distribution3(1,3)=sum(sum(sim.ls==3 & xhistory<0))/sum(sum(xhistory<0));
        sim.size_distribution3(1,4)=sum(sum(sim.ls==4 & xhistory<0))/sum(sum(xhistory<0));
        sim.size_distribution3(2,1)=sum(sum(sim.ls==1 & xhistory>=1))/sum(sum(xhistory>=1)); %Exporters
        sim.size_distribution3(2,2)=sum(sum(sim.ls==2 & xhistory>=1))/sum(sum(xhistory>=1));
        sim.size_distribution3(2,3)=sum(sum(sim.ls==3 & xhistory>=1))/sum(sum(xhistory>=1));
        sim.size_distribution3(2,4)=sum(sum(sim.ls==4 & xhistory>=1))/sum(sum(xhistory>=1));               
            
    %Sales distribution
        sim.ss_size1 = prctile(sales(:),25);
        sim.ss_size2 = prctile(sales(:),50);   
        sim.ss_size3 = prctile(sales(:),75);   
        sim.ss = (sales<=sim.ss_size1) + 2*(sales>sim.ss_size1 & sales<=sim.ss_size2) + 3*(sales>sim.ss_size2 & sales<=sim.ss_size3) + 4*(sales>sim.ss_size3);

        %Size-Status Distribution: for each size, distribution of export status
        %Status-Size Distribution: for each export status, distribution of size
        sim.Ssize_distribution=zeros(2,4);
        sim.Ssize_distribution(1,1)=sum(sum(sim.ss==1 & xhistory<0))/sum(sum(xhistory<0)); %Non-exporters
        sim.Ssize_distribution(1,2)=sum(sum(sim.ss==2 & xhistory<0))/sum(sum(xhistory<0));
        sim.Ssize_distribution(1,3)=sum(sum(sim.ss==3 & xhistory<0))/sum(sum(xhistory<0));
        sim.Ssize_distribution(1,4)=sum(sum(sim.ss==4 & xhistory<0))/sum(sum(xhistory<0));
        sim.Ssize_distribution(2,1)=sum(sum(sim.ss==1 & xhistory>=1))/sum(sum(xhistory>=1)); %Exporters
        sim.Ssize_distribution(2,2)=sum(sum(sim.ss==2 & xhistory>=1))/sum(sum(xhistory>=1));
        sim.Ssize_distribution(2,3)=sum(sum(sim.ss==3 & xhistory>=1))/sum(sum(xhistory>=1));
        sim.Ssize_distribution(2,4)=sum(sum(sim.ss==4 & xhistory>=1))/sum(sum(xhistory>=1));           
        
    clear sim.laborsize xhistory;
    
    return; %Breaks function, goes back to main.m
       
end %s.flag_short
    
    

    
%% Long simulate 

% State variables -> Labor and export status (whether the firm exported or not in the previous period)

%Initialize objects
    ind_y = zeros(s.T*(1+s.burn)+1,s.N);
    sim.exporter = zeros(s.T*(1+s.burn)+1,s.N);   
    
%First period initialization
    ind_y(1,1:s.N) = round(s.y_grid_size/4); %We initialize every firm with lowest assets
    
%Compute labor and export states for t>1    
    for t=1:s.T*(1+s.burn)  
                

        ind_k = sim.exporter(t,:)+1; %k=1 is new exporter, k=2 is continuing exporter

         
      
        ind = s.y_grid_size*s.z_grid_size*(ind_k-1)+s.y_grid_size*(ind_z(t,:)-1)+ind_y(t,:);
        
        %Update state variables    
            sim.exporter(t+1,:) = r.exporter(ind); % r.exporter is a decision rule given states today, sim.exporter updates simulation   
            ind_y(t+1,:)= r.y_prime_ind(ind);
                 
    end
    
    sim.exporter_decision = sim.exporter(2:end,:); %Export decisions (as opposed to status)
    ind_y_decision = ind_y(2:end,:);
    sim.y=r.y_grid(ind_y_decision); 
    
%Compute export history
%If x_history==k and k>0 -> Firm has exported for k consecutive periods
%If x_history==k and k<0 -> Firm has not exported for k consecutive periods
%In the first period, assumes all exporters are new
    sim.xhistory = zeros(s.T*(1+s.burn),s.N);
    sim.xhistory(1,:) = sim.exporter_decision(1,:); %First period assume all exporters are new
    sim.xhistory(1,sim.exporter_decision(1,:)==0) = -1;
    for t = 2:s.T*(1+s.burn)
        for j=1:s.N
            %If choose to export today
            if sim.exporter_decision(t,j)==1 && sim.exporter_decision(t-1,j)==0 %If new exporter
                sim.xhistory(t,j) = 1;
            elseif sim.exporter_decision(t,j)==1 && sim.exporter_decision(t-1,j)==1 %If continuing exporter
                sim.xhistory(t,j) = sim.xhistory(t-1,j)+1;
                
            %If choose not to export today
            elseif sim.exporter_decision(t,j)==0 && sim.exporter_decision(t-1,j)==1 %If new non-exporter
                sim.xhistory(t,j) = -1;
            elseif sim.exporter_decision(t,j)==0 && sim.exporter_decision(t-1,j)==0 %If continuing non-exporter
                sim.xhistory(t,j) = sim.xhistory(t-1,j) - 1;
            end
        end
    end
    
       
%% Burn first s.burn*s.T time periods

    %Exogenous state / export policy function / history-dependent statistics
        ind_z = ind_z(s.T*s.burn+1:end,:);
        sim.exporter_decision = sim.exporter_decision(s.T*s.burn+1:end,:);
        sim.xhistory = sim.xhistory(s.T*s.burn+1:end,:);                
        %sim.exporter_decision = sim.exporter_decision(s.T*s.burn+1:end,:);
        ind_y_decision = ind_y_decision(s.T*s.burn+1:end,:);
        sim.y = sim.y(s.T*s.burn+1:end,:);
        
    %Endogenous states
    %Also drop last period, which reflects state at period (1+s.burn)*s.T+1 from choice made at period (1+s.burn)*s.T
        ind_y = ind_y(s.T*s.burn+1:end-1,:);        
        sim.exporter = sim.exporter(s.T*s.burn+1:end-1,:);    
       

%% Simulate fundamental variables 

%Initialize objects
    sim.sales_f = zeros(s.T,s.N);
    sim.sales_d = zeros(s.T,s.N);
    sim.q_f = zeros(s.T,s.N);
    sim.q_d = zeros(s.T,s.N);
    sim.n_f = zeros(s.T,s.N);
    sim.n_d = zeros(s.T,s.N);
    sim.n = zeros(s.T,s.N);
    sim.p_f = zeros(s.T,s.N);
    sim.p_d = zeros(s.T,s.N);

%With sunk costs, need to take care of extra state (ie. previous export status)
%Note that ind_k is now a matrix, above it was a vector
    if m.S>0
        ind_k = sim.exporter+1; %k=1 if new exporter, k=2 if continuing exporter
     else
        ind_k = 1; %Status does not matter
    end
        
      
    ind=s.y_grid_size*s.z_grid_size*(ind_k-1)+s.y_grid_size*(ind_z-1)+ind_y;
    ind_decision=s.y_grid_size*s.z_grid_size*(ind_k-1)+s.y_grid_size*(ind_z-1)+ind_y_decision;
    
%Quantities            
     sim.q_f = r.q_f_export(ind_decision).*sim.exporter_decision;
     sim.q_d(sim.exporter_decision==0) = r.q_d(ind_decision(sim.exporter_decision==0));
     sim.q_d(sim.exporter_decision==1) = r.q_d_export(ind_decision(sim.exporter_decision==1));
       
%Prices     
     sim.p_f = r.p_f_export(ind_decision).*sim.exporter_decision;
     sim.p_d(sim.exporter_decision==0) = r.p_d(ind_decision(sim.exporter_decision==0));
     sim.p_d(sim.exporter_decision==1) = r.p_d_export(ind_decision(sim.exporter_decision==1));
        
%Labor
     sim.n_f = r.n_f_export(ind_decision).*sim.exporter_decision;
     sim.n_d(sim.exporter_decision==0) = r.n_d(ind_decision(sim.exporter_decision==0));
     sim.n_d(sim.exporter_decision==1) = r.n_d_export(ind_decision(sim.exporter_decision==1));        
        
     sim.n = sim.n_f + sim.n_d;
     sim.n(sim.xhistory==1) = sim.n(sim.xhistory==1) + m.F + m.S;
     %sim.n(sim.xhistory>0) = sim.n(sim.xhistory>0) + m.F;
     sim.n(sim.xhistory>1) = sim.n(sim.xhistory>1) + m.F;
     
%Sales     
     sim.sales_f = sim.p_f.*sim.q_f;
     sim.sales_d = sim.p_d.*sim.q_d;
     sim.sales = sim.sales_d + sim.sales_f;

%Productivity     
    sim.z = r.z_grid(ind_z); 
    
     
%% Tests

%Share of firms with assets at highest labor state
   sim.share_y_max = sum(sum(sim.y==r.y_grid(end)))/sum(sum(sim.y>0)); 
   sim.share_y_min = sum(sum(sim.y==r.y_grid(1)))/sum(sum(sim.y>0)); 
       
%Volatility across time periods of average productivity and average labor   
    sim.y_mean_std = std(mean(sim.y,2))/mean(mean(sim.y));
    sim.z_mean_std = std(mean(sim.z,2))/mean(mean(sim.z));
     
%% Variables, firm-level statistics, and aggregate statistics
     
%Variables
    %Export intensity     
        sim.export_intensity = sim.sales_f./(sim.sales_d+sim.sales_f);     

    %Growth rates     
        sim.sales_growth   = [zeros(1,s.N);(sim.sales(2:end,:)-sim.sales(1:end-1,:))./sim.sales(1:end-1,:)];
        sim.sales_d_growth = [zeros(1,s.N);(sim.sales_d(2:end,:)-sim.sales_d(1:end-1,:))./sim.sales_d(1:end-1,:)];
        sim.sales_f_growth = [zeros(1,s.N);(sim.sales_f(2:end,:)-sim.sales_f(1:end-1,:))./sim.sales_f(1:end-1,:)];
        sim.y_growth   = [zeros(1,s.N);(sim.y(2:end,:)-sim.y(1:end-1,:))./sim.y(1:end-1,:)];
        sim.y_loggrowth = [zeros(1,s.N);log(sim.y(2:end,:)./sim.y(1:end-1,:))];      
        
        %     sim.sales_growth   = [zeros(1,s.N);log(sim.sales(2:end,:)./sim.sales(1:end-1,:))];
        %     sim.sales_d_growth = [zeros(1,s.N);log(sim.sales_d(2:end,:)./sim.sales_d(1:end-1,:))];
        %     sim.sales_f_growth = [zeros(1,s.N);log(sim.sales_f(2:end,:)./sim.sales_f(1:end-1,:))];
        
%Aggregate statistics  
    %GDP
        sim.GDP = sum(sim.sales_f + sim.sales_d,2);    
    
    %Exports
        sim.exports_agg = sum(sim.sales_f,2);    
        
%Firm-level statistics
    %Number of exporters, share of exporters, average exports
        sim.exporters_num = sum(sim.exporter_decision,2); %Number of exporters per period
        sim.exports_firm = sim.exports_agg./sim.exporters_num;
        sim.share_exporters = mean(sum(sim.exporter_decision,2)/s.N);
        
    %Export entry and exit rates
        sim.share_starters = mean(sum(sim.xhistory(2:end,:)==1,2)./sum(sim.xhistory(1:end-1,:)<=-1,2));
        sim.share_stoppers = mean(sum(sim.xhistory(2:end,:)==-1,2)./sum(sim.xhistory(1:end-1,:)>=1,2));    
        
    %Median export intensity
        temp = sim.export_intensity(sim.xhistory>=1);
        sim.export_sales_ratio_med = median(temp(:));        
        clear temp
        
    %Share of exporters that are new
        sim.sharex_new = sum(sum(sim.xhistory==1))/sum(sum(sim.xhistory>=1));

    %Exporter size premium -> exporters vs non-exporters
        %Labor            
            tempS = zeros(size(sim.q_f));
            tempS(sim.xhistory==1) = 1;
            temp1 = ((m.tau*sim.q_f(sim.xhistory>=1)+sim.q_d(sim.xhistory>=1))./sim.z(sim.xhistory>=1))+m.F+tempS(sim.xhistory>=1)*m.S;
            temp2 = sim.q_d(sim.xhistory<=-1)./sim.z(sim.xhistory<=-1);
            sim.export_size_premium_med = median(temp1(:))/median(temp2(:));   
            
            sim.testF2 = median(((m.tau*sim.q_f(sim.xhistory>=1))./sim.z(sim.xhistory>=1))+m.F+tempS(sim.xhistory>=1)*m.S);
            sim.testD2 = median(((sim.q_d(sim.xhistory>=1))./sim.z(sim.xhistory>=1)));          
           
            clear temp1 temp2
            
    %Continuing exporter size premium -> continuing exporters vs new exporters
        %Labor
            temp1 = ((m.tau*sim.q_f(sim.xhistory>1)+sim.q_d(sim.xhistory>1))./sim.z(sim.xhistory>1))+m.F;            
            temp2 = ((m.tau*sim.q_f(sim.xhistory==1)+sim.q_d(sim.xhistory==1))./sim.z(sim.xhistory==1))+m.F+m.S;
            sim.cont_export_size_premium_med = median(temp1(:))/median(temp2(:));  
            clear temp1 temp2
            
    %New exporter size premium -> new exporters vs non exporters
        %Labor
            temp1 = ((m.tau*sim.q_f(sim.xhistory==1)+sim.q_d(sim.xhistory==1))./sim.z(sim.xhistory==1))+m.F+m.S;
            
            sim.testF = median(m.tau*sim.q_f(sim.xhistory==1)./sim.z(sim.xhistory==1) + m.F + m.S);
            sim.testD = median(((sim.q_d(sim.xhistory==1))./sim.z(sim.xhistory==1)));
            sim.testD3 = median((sim.n_d(sim.xhistory==1)));
            
            temp2 = sim.q_d(sim.xhistory<=-1)./sim.z(sim.xhistory<=-1);
            sim.new_export_size_premium_med = median(temp1(:))/median(temp2(:));   
            clear temp1 temp2         
            
    %New exporters / stoppers 
        tempS = zeros(size(sim.q_f));
        tempS(sim.xhistory==1) = 1;
        temp1 = ((m.tau*sim.q_f(sim.xhistory==1)+sim.q_d(sim.xhistory==1))./sim.z(sim.xhistory==1))+m.F+tempS(sim.xhistory==1)*m.S;
        temp2 = (sim.q_d(sim.xhistory==-1))./sim.z(sim.xhistory==-1);
        sim.new_stoppers_size_premium_med = median(temp1(:))/median(temp2(:));  
        clear tempS temp1 temp2           
  
            
          % New non-exporter size premium -> new non-exporters vs exporters
    
        %Labor
            tempS = zeros(size(sim.q_f));
            tempS(sim.xhistory==1) = 1;
            temp1 = ((m.tau*sim.q_f(sim.xhistory>=1)+sim.q_d(sim.xhistory>=1))./sim.z(sim.xhistory>=1))+m.F+tempS(sim.xhistory>=1)*m.S;
            temp2 = sim.n_d(sim.xhistory==-1);
            sim.new_nonexport_size_premium_med = median(temp2(:))/median(temp1(:));  
            clear tempS temp1 temp2             
            
   
    sim.sales_growth_4yravg = (1/4)*[(sim.sales_d(4:end,:)+sim.sales_f(4:end,:)-sim.sales_d(1:end-3,:)-sim.sales_f(1:end-3,:))./(sim.sales_d(1:end-3,:)+sim.sales_f(1:end-3,:))];  

            
        BJ_length = 1;
    for t=1:1 %Starts at 1, and can go up to s.T-BJ_length (its robust to these changes)
       sim.BJstart(t,:) = sim.xhistory(t,:)<0 & sim.xhistory(t+BJ_length,:)>0; 
       sim.BJstop(t,:) = sim.xhistory(t,:)>0 & sim.xhistory(t+BJ_length,:)<0; 
       sim.BJboth(t,:) = sim.xhistory(t,:)>0 & sim.xhistory(t+BJ_length,:)>0; 
       sim.BJneither(t,:) = sim.xhistory(t,:)<0 & sim.xhistory(t+BJ_length,:)<0; 
    end
    
    sim.sales_growth_4yravg = sim.sales_growth_4yravg(size(sim.BJstart,1),:);    
    
    sim.BJstart_sales_gr = median(sim.sales_growth_4yravg(sim.BJstart==1));
    sim.BJstop_sales_gr = median(sim.sales_growth_4yravg(sim.BJstop==1));
    sim.BJboth_sales_gr = median(sim.sales_growth_4yravg(sim.BJboth==1));
    sim.BJneither_sales_gr = median(sim.sales_growth_4yravg(sim.BJneither==1));
            
        
   
%% Statistics for cohorts of firms that didn't export for s.CD years, and exported for s.CN consecutive years 

%Compute cohorts
    sim.xperiods_cohorts = zeros(s.T,s.N);
    for k=1:s.T-s.CN-s.CD %Cohorts
        xperiods_exporters = sim.xhistory(k+s.CD+s.CN-1,:)==s.CN;
        sim.xperiods_cohorts(k:k+s.CD+s.CN-1,xperiods_exporters) = 1;
    end
    sim.xperiods_xhistory = sim.xperiods_cohorts.*sim.xhistory; %First value will be any negative number

%Compute statistics -> Averages and medians across cohorts, by firm age within the cohort
for i = 1:s.CN    
    %Median
        %Level
            sim.xperiods_sales_f_med(i,1) = median(sim.sales_f(sim.xperiods_xhistory==i));
            sim.xperiods_sales_d_med(i,1) = median(sim.sales_d(sim.xperiods_xhistory==i));
            sim.xperiods_sales_med(i,1) = median(sim.sales(sim.xperiods_xhistory==i));
            sim.xperiods_export_intensity_med(i,1) = median(sim.export_intensity(sim.xperiods_xhistory==i));
            sim.xperiods_z_med(i,1) = median(sim.z(sim.xperiods_xhistory==i));
            sim.xperiods_y_med(i,1) = median(sim.y(sim.xperiods_xhistory==i));  
            
        %Growth       
            sim.xperiods_sales_f_growth_med(i,1) = median(sim.sales_f_growth(sim.xperiods_xhistory==i));
            sim.xperiods_sales_d_growth_med(i,1) = median(sim.sales_d_growth(sim.xperiods_xhistory==i));    
            sim.xperiods_sales_growth_med(i,1) = median(sim.sales_growth(sim.xperiods_xhistory==i));        
            sim.xperiods_y_growth_med(i,1) = median(sim.y_growth(sim.xperiods_xhistory==i));  

    %Mean
        %Log-growth
            sim.xperiods_y_loggrowth_mean(i,1) = mean(sim.y_loggrowth(sim.xperiods_xhistory==i));   
            
end
sim.xperiods_sales_d_growth_med_before_x = median(sim.sales_d_growth(sim.xperiods_xhistory<0 & sim.xperiods_cohorts==1));
    
sim.xperiods_y_growth_med_before_x = median(sim.y_growth(sim.xperiods_xhistory<0 & sim.xperiods_cohorts==1));
sim.xperiods_y_loggrowth_mean_before_x = mean(sim.y_loggrowth(sim.xperiods_xhistory<0 & sim.xperiods_cohorts==1));

%Add pre-entry period to cohort statistics
    sim.xperiods_sales_d_growth_med = [sim.xperiods_sales_d_growth_med_before_x;sim.xperiods_sales_d_growth_med];
    sim.xperiods_y_growth_med = [sim.xperiods_y_growth_med_before_x;sim.xperiods_y_growth_med];
    sim.xperiods_y_loggrowth_mean = [sim.xperiods_y_loggrowth_mean_before_x;sim.xperiods_y_loggrowth_mean]; 

%Hazard rate -> computed by pooling firms from all cohorts
%     for i=1:6     
%         sim.xperiod_count_pool(i) = sum(sum(sim.xhistory==i));
%     end
%     for i=1:5
%         sim.hazard_pool(i) = (sim.xperiod_count_pool(i)-sim.xperiod_count_pool(i+1))/sim.xperiod_count_pool(i);
%     end

  % Feb 13th, corrected typo: we should count only firms at t if they
   % are counted at t-1, and firms at t if they can be counted at t+1.

    %Hazard rate -> computed by pooling firms from all cohorts
    sim.xperiod_count_pool1 = zeros(1,6); % Counts only firms at t, if they are counted at t-1
    sim.xperiod_count_pool2 = zeros(1,6); % Counts only firms at t, if they can be counted at t+1
    
    for i=1:6     
       
        sim.xperiod_count_pool1(i) = sum(sum(sim.xhistory(2:end,:)==i)); % Counts only firms at t, if they are counted at t-1
        sim.xperiod_count_pool2(i) = sum(sum(sim.xhistory(1:end-1,:)==i)); % Counts only firms at t, if they can be counted at t+1
        
    end
    
    
    for i=1:5
    
        sim.hazard_pool(i) = (sim.xperiod_count_pool2(i)-sim.xperiod_count_pool1(i+1))/sim.xperiod_count_pool2(i);
    
    end     
    
    
    % Hazard for stoppers -> computed by pooling firms from all cohorts


  % Feb 13th, corrected typo: we should count only firms at t if they
   % are counted at t-1, and firms at t if they can be counted at t+1.

    %Hazard rate -> computed by pooling firms from all cohorts
    sim.nonxperiod_count_pool1 = zeros(1,6); % Counts only firms at t, if they are counted at t-1
    sim.nonxperiod_count_pool2 = zeros(1,6); % Counts only firms at t, if they can be counted at t+1
    
    for i=1:6      
       
        sim.nonxperiod_count_pool1(i) = sum(sum(sim.xhistory(2:end,:)==-i)); % Counts only firms at t, if they are counted at t-1
        sim.nonxperiod_count_pool2(i) = sum(sum(sim.xhistory(1:end-1,:)==-i)); % Counts only firms at t, if they can be counted at t+1
        
    end
    
    
    for i=1:5
    
        sim.nonxhazard_pool(i) = (sim.nonxperiod_count_pool2(i)-sim.nonxperiod_count_pool1(i+1))/sim.nonxperiod_count_pool2(i);
    
    end   
%% Other variables / moments (optional)

if s.flag_extraresults==1

sim.export_sales_ratio = mean(mean(sim.export_intensity(sim.xhistory>=1)));              
sim.export_intensity_q = sim.q_f./(sim.q_d+sim.q_f);
sim.z_growth = [zeros(1,s.N);log(sim.z(2:end,:)./sim.z(1:end-1,:))];
sim.y_growth = [zeros(1,s.N);log(sim.y(2:end,:)./sim.y(1:end-1,:))];

%Exporter size premium -> exporters vs non-exporters
    %Labor         
        tempS = zeros(size(sim.q_f));
        tempS(sim.xhistory==1) = 1;
        temp1 = ((m.tau*sim.q_f(sim.xhistory>=1)+sim.q_d(sim.xhistory>=1))./sim.z(sim.xhistory>=1))+m.F+tempS(sim.xhistory>=1)*m.S;
        temp2 = sim.q_d(sim.xhistory<=-1)./sim.z(sim.xhistory<=-1);
        sim.export_size_premium2 = mean(temp1(:))/mean(temp2(:));   
        clear temp1 temp2                

    %Sales
        sim.export_size_premium_sales = mean( mean( (sim.sales_d(sim.xhistory>=1)+sim.sales_f(sim.xhistory>=1)) ) )...
            ./mean( mean( sim.sales_d(sim.xhistory<=-1) ) );

        temp1 = (sim.sales_d(sim.xhistory>=1)+sim.sales_f(sim.xhistory>=1));
        temp2 = sim.sales_d(sim.xhistory<=-1);
        sim.export_size_premium_sales_med = median(temp1(:))/median(temp2(:));          
        clear temp1 temp2

%Continuing exporter size premium -> continuing exporters vs new exporters
    %Labor
        sim.export_size_premium_cont = mean(mean(((m.tau*sim.q_f(sim.xhistory>=2)+sim.q_d(sim.xhistory>=2))./sim.z(sim.xhistory>=2))+m.F))...
            ./mean(mean(((m.tau*sim.q_f(sim.xhistory==1)+sim.q_d(sim.xhistory==1))./sim.z(sim.xhistory==1))+m.F+m.S));        

    %Sales
        sim.export_size_premium_cont_sales = mean(  mean( sim.sales_f(sim.xhistory>=2)+sim.sales_d(sim.xhistory>=2)) )...
            ./mean( mean( (sim.sales_f(sim.xhistory==1)+sim.sales_d(sim.xhistory==1)) ) );        

        temp1 = (sim.sales_d(sim.xhistory>1)+sim.sales_f(sim.xhistory>1));
        temp2 = (sim.sales_d(sim.xhistory==1)+sim.sales_f(sim.xhistory==1));
        sim.export_size_premium_cont_sales_med = median(temp1(:))/median(temp2(:));   
        clear temp1 temp2
            
%Export cohort statistics -> Averages across cohorts, by firm age within the cohort
for i = 1:s.CN 
    %Average
        %Level
            sim.xperiods_sales_f_mean(i,1) = mean(sim.sales_f(sim.xperiods_xhistory==i));
            sim.xperiods_sales_d_mean(i,1) = mean(sim.sales_d(sim.xperiods_xhistory==i));
            sim.xperiods_sales_mean(i,1) = mean(sim.sales(sim.xperiods_xhistory==i));
            sim.xperiods_export_intensity_mean(i,1) = mean(sim.export_intensity(sim.xperiods_xhistory==i));
            sim.xperiods_z_mean(i,1) = mean(sim.z(sim.xperiods_xhistory==i));
            sim.xperiods_y_mean(i,1) = mean(sim.y(sim.xperiods_xhistory==i));  
                
            
        %Growth
            sim.xperiods_sales_f_growth_mean(i,1) = mean(sim.sales_f_growth(sim.xperiods_xhistory==i));
            sim.xperiods_sales_d_growth_mean(i,1) = mean(sim.sales_d_growth(sim.xperiods_xhistory==i)); 
            sim.xperiods_sales_growth_mean(i,1) = mean(sim.sales_growth(sim.xperiods_xhistory==i));      
       
end
sim.xperiods_sales_d_growth_mean_before_x = mean(sim.sales_d_growth(sim.xperiods_xhistory<0 & sim.xperiods_cohorts==1));
               
%Sales statistics   
    for i=1:s.N
        sim.sales_acf(:,i) = autocorr(sim.sales(:,i));                
    end
    sim.sales_acf = mean(sim.sales_acf(:,i),2);
     
    sim.sales_std = std(sim.sales_growth(:));
    sim.sales_kurt = kurtosis(sim.sales_growth(:));
    sim.sales_iqr = iqr(sim.sales_growth(:))/sim.sales_std; %CHECK, NOT SURE IF CORRECT
    
    sim.sales_sort = sort(sim.sales(:));
    sim.sales_sum = sum(sim.sales(:));
    sim.sales_p1 = sum(sim.sales_sort(round(0.99*s.T*s.N):s.T*s.N))/sim.sales_sum; %Share of aggregate sales by top 1%
    sim.sales_p5 = sum(sim.sales_sort(round(0.95*s.T*s.N):s.T*s.N))/sim.sales_sum; 
    sim.sales_p10 = sum(sim.sales_sort(round(0.90*s.T*s.N):s.T*s.N))/sim.sales_sum;
    sim.sales_p20 = sum(sim.sales_sort(round(0.80*s.T*s.N):s.T*s.N))/sim.sales_sum;
    
%Labor statistics
    wagebill = m.w*sim.n;
    sim.wagebill_acf_mat = zeros(2,s.N);
    for i=1:s.N
        sim.wagebill_acf_mat(:,i) = autocorr(log(wagebill(:,i)),1); 
    end
    sim.wagebill_acf_med = median(sim.wagebill_acf_mat(2,:));
    sim.wagebill_acf_mean = mean(sim.wagebill_acf_mat(2,:));    
    
%Labor size distribution
%RECHECK VALUES USED WITH THOSE IN THE DATA
    sim.laborsize_psmall = prctile(sim.n(:),41.22);
    sim.laborsize_pmed = prctile(sim.n(:),83.34);    
    sim.laborsize = (sim.n<=sim.laborsize_psmall) + 2*(sim.n>sim.laborsize_psmall & sim.n<=sim.laborsize_pmed) + 3*(sim.n>sim.laborsize_pmed);
    
    %Status-Size Distribution: for each export status, distribution of size
        sim.size_distribution=zeros(3,3);
        sim.size_distribution(1,1)=sum(sum(sim.laborsize==1 & sim.xhistory<0))/sum(sum(sim.xhistory<0)); %Non-exporters
        sim.size_distribution(1,2)=sum(sum(sim.laborsize==2 & sim.xhistory<0))/sum(sum(sim.xhistory<0));
        sim.size_distribution(1,3)=sum(sum(sim.laborsize==3 & sim.xhistory<0))/sum(sum(sim.xhistory<0));
        sim.size_distribution(2,1)=sum(sum(sim.laborsize==1 & sim.xhistory==1))/sum(sum(sim.xhistory==1)); %New exporters
        sim.size_distribution(2,2)=sum(sum(sim.laborsize==2 & sim.xhistory==1))/sum(sum(sim.xhistory==1));
        sim.size_distribution(2,3)=sum(sum(sim.laborsize==3 & sim.xhistory==1))/sum(sum(sim.xhistory==1));
        sim.size_distribution(3,1)=sum(sum(sim.laborsize==1 & sim.xhistory>1))/sum(sum(sim.xhistory>1)); %Cont exporters
        sim.size_distribution(3,2)=sum(sum(sim.laborsize==2 & sim.xhistory>1))/sum(sum(sim.xhistory>1));
        sim.size_distribution(3,3)=sum(sum(sim.laborsize==3 & sim.xhistory>1))/sum(sum(sim.xhistory>1));
    
    %Size-Status Distribution: for each size, distribution of export status
        sim.status_distribution=zeros(3,3);
        sim.status_distribution(1,1)=sum(sum(sim.laborsize==1 & sim.xhistory<0))/sum(sum(sim.laborsize==1)); %Small 
        sim.status_distribution(1,2)=sum(sum(sim.laborsize==1 & sim.xhistory==1))/sum(sum(sim.laborsize==1));
        sim.status_distribution(1,3)=sum(sum(sim.laborsize==1 & sim.xhistory>1))/sum(sum(sim.laborsize==1));
        sim.status_distribution(2,1)=sum(sum(sim.laborsize==2 & sim.xhistory<0))/sum(sum(sim.laborsize==2)); %Medium
        sim.status_distribution(2,2)=sum(sum(sim.laborsize==2 & sim.xhistory==1))/sum(sum(sim.laborsize==2));
        sim.status_distribution(2,3)=sum(sum(sim.laborsize==2 & sim.xhistory>1))/sum(sum(sim.laborsize==2));
        sim.status_distribution(3,1)=sum(sum(sim.laborsize==3 & sim.xhistory<0))/sum(sum(sim.laborsize==3)); %Large
        sim.status_distribution(3,2)=sum(sum(sim.laborsize==3 & sim.xhistory==1))/sum(sum(sim.laborsize==3));
        sim.status_distribution(3,3)=sum(sum(sim.laborsize==3 & sim.xhistory>1))/sum(sum(sim.laborsize==3));
    
        sim.laborsize_largemed_mean = mean(sim.n(sim.laborsize==3))/mean(sim.n(sim.laborsize==2));
        sim.laborsize_medsmall_mean = mean(sim.n(sim.laborsize==2))/mean(sim.n(sim.laborsize==1));
        sim.laborsize_largemed_med = median(sim.n(sim.laborsize==3))/median(sim.n(sim.laborsize==2));
        sim.laborsize_medsmall_med = median(sim.n(sim.laborsize==2))/median(sim.n(sim.laborsize==1));
    
    %Sales distribution
        sim.ss_size1 = prctile(sim.sales(:),25);
        sim.ss_size2 = prctile(sim.sales(:),50);   
        sim.ss_size3 = prctile(sim.sales(:),75);   
        sim.ss = (sim.sales<=sim.ss_size1) + 2*(sim.sales>sim.ss_size1 & sim.sales<=sim.ss_size2) + 3*(sim.sales>sim.ss_size2 & sim.sales<=sim.ss_size3) + 4*(sim.sales>sim.ss_size3);

        %Size-Status Distribution: for each size, distribution of export status
        %Status-Size Distribution: for each export status, distribution of size
        sim.Ssize_distribution=zeros(2,4);
        sim.Ssize_distribution(1,1)=sum(sum(sim.ss==1 & sim.xhistory<0))/sum(sum(sim.xhistory<0)); %Non-exporters
        sim.Ssize_distribution(1,2)=sum(sum(sim.ss==2 & sim.xhistory<0))/sum(sum(sim.xhistory<0));
        sim.Ssize_distribution(1,3)=sum(sum(sim.ss==3 & sim.xhistory<0))/sum(sum(sim.xhistory<0));
        sim.Ssize_distribution(1,4)=sum(sum(sim.ss==4 & sim.xhistory<0))/sum(sum(sim.xhistory<0));
        sim.Ssize_distribution(2,1)=sum(sum(sim.ss==1 & sim.xhistory>=1))/sum(sum(sim.xhistory>=1)); %Exporters
        sim.Ssize_distribution(2,2)=sum(sum(sim.ss==2 & sim.xhistory>=1))/sum(sum(sim.xhistory>=1));
        sim.Ssize_distribution(2,3)=sum(sum(sim.ss==3 & sim.xhistory>=1))/sum(sum(sim.xhistory>=1));
        sim.Ssize_distribution(2,4)=sum(sum(sim.ss==4 & sim.xhistory>=1))/sum(sum(sim.xhistory>=1));      
        
%Average labor 
    sim.y_exporters_avg = sum(sim.exporter_decision.*sim.y,2)./sum(sim.exporter_decision,2);
    sim.y_nonexporters_avg = sum((1-sim.exporter_decision).*sim.y,2)./sum((1-sim.exporter_decision),2);
    sim.y_avg = sum(sim.y,2)/s.N;

end %end if, extra results
    
%% Notes
%1) We have coded an alternative way of computing the statistics for cohorts of exporters, which compute statistics by cohorts and then take averages, see old code (eg. January 2012 folder) 

end


