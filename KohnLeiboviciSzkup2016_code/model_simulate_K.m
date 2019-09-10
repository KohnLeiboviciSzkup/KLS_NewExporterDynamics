function [r sim] = model_simulate_K(m,s,r)


    s.flag_short=0;
    s.flag_extraresults=1;

%% Simulate productivity shocks

%Initalize objects
    ind_z_prod = zeros(s.N,s.T*(1+s.burn));
    ind_y = zeros(s.N,s.T*(1+s.burn));
    
%Initialize productivity realizations at stationary distribution
%Input -> z_pi, stationary distribution of productivity shocks, vector of size [s.z_grid_size x 1], computed by approx_shocks.m when run by model_solve.m
%Output -> ind_z(:,1), vector of indexes across productivity states for each firm, size [s.N x 1]
    k=1;
    for i=1:length(r.z_pi) %length(r.log_z_grid)
        v(k:k+floor(r.z_pi(i)*s.N)-1)=i;  
        k=length(v)+1;
    end
    v(length(v):s.N)=round(mean(v));
    ind_z_prod(:,1)=v';
    clear v;

%Compute productivity realizations after initial state
%Input -> z.P, transition matrix computed by approx_shocks.m
%Output -> ind_z(:,2:end)
    C = cumsum(r.z_P_original,2);
    R = rand(s.N,s.T*(1+s.burn));
    for j=2:s.T*(1+s.burn);
        f = repmat(R(:,j),1,s.z_grid_size_original);
        ind_z_prod(:,j) = 1+sum( f > C(ind_z_prod(:,j-1),:) ,2);
        clear f;
    end
    clear R C;
    ind_z_prod = ind_z_prod';
    
%Simulate F shocks
if s.y_grid_size_original>1
    
    k=1;
    for i=1:length(r.y_grid_original)
        v(k:k+floor(r.y_pi(i)*s.N)-1)=i;  
        k=length(v)+1;
    end
    v(length(v):s.N)=round(mean(v));
    ind_y(:,1)=v';
    clear v;

    C = cumsum(r.y_P,2);
    R = rand(s.N,s.T*(1+s.burn));
    for j=2:s.T*(1+s.burn);
        f = repmat(R(:,j),1,s.y_grid_size_original);
        ind_y(:,j) = 1+sum( f > C(ind_y(:,j-1),:) ,2);
        clear f;
    end
    clear R C;
    ind_y = ind_y';  
    
else
    
    ind_y = ones(s.T*(1+s.burn),s.N);
    
end
    
%Combine ind_y and ind_z_prod into ind_z
    ind_z = (ind_y-1)*s.z_grid_size_original+ind_z_prod;
    

%% Simulate y shocks, if active
    
if s.Fshocks==1 || s.Dshocks==1   

% Initalize objects
    ind_y = ones(s.N,s.T*(1+s.burn));    
    
% Same as above for productivity shocks

    k=1;
    for i=1:length(r.y_grid)
        v(k:k+floor(r.y_pi(i)*s.N)-1)=i;  
        k=length(v)+1;
    end
    v(length(v):s.N)=round(mean(v));
    ind_y(:,1)=v';
    clear v;


    C = cumsum(r.y_P,2);
    R = rand(s.N,s.T*(1+s.burn));
    for j=2:s.T*(1+s.burn);
        f = repmat(R(:,j),1,s.y_grid_size);
        ind_y(:,j) = 1+sum( f > C(ind_y(:,j-1),:) ,2);
        clear f;
    end
    clear R C;
    ind_y = ind_y';    
    

elseif s.Tshocks==1
    
    % Initalize objects
    ind_y = ones(s.N,s.T*(1+s.burn));    
    
    sim.tauRND = rand(s.N,s.T*(1+s.burn)+1)';
    ind_y = ind_y';
    
end

%% s.flag_short==1, only computes basic statistics  

if s.flag_short==1

        % Clearing variables soution from dynamic problem that won't be
        % essential
                
    if s.assets==1    
        
        r.a_prime=[];
        r.log_z_grid=[];
        r.z_P=[];
        r.z_pi=[];
        r.constrained=[];
        r.constrained_d=[];
        r.constrained_f=[];
        r.c=[];
        r.v=[];
        r.pi_d=[];
        r.pi_export=[];        
        
        %Initialize objects
        ind_a = zeros(s.T*(1+s.burn)+1,s.N);
        exporter = zeros(s.T*(1+s.burn)+1,s.N);   


        %First period initialization
        ind_a(1,1:s.N) = 1; %We initialize every firm with lowest assets

        %Compute asset and export states for t>1    
        for t=1:s.T*(1+s.burn)  

            %With sunk costs, need to take care of extra state (ie. previous export status)
                if m.S>0
                   ind_k = exporter(t,:)+1; %k=1 is new exporter, k=2 is continuing exporter
                else
                   ind_k =1; %Status does not matter
                end

            %Construct 3-dimensional indexes from unidimensional indexes ind_a, ind_z, and ind_k     
            %The policy functions are multidimensional matrices -> we map ind_a, ind_z, and ind_k to the policy functions using ind        
                ind = s.a_grid_size*s.z_grid_size*(ind_k-1)+s.a_grid_size*(ind_z(t,:)-1)+ind_a(t,:);

            %Update state variables    
                exporter(t+1,:) = r.exporter(ind); % r.exporter is a decision rule given states today, sim.exporter updates simulation   
                ind_a(t+1,:)= r.a_prime_ind(ind);

        end
    
        sim.a=r.a_grid(ind_a); %Asset values
        exporter_decision = exporter(2:end,:); %Export decisions (as opposed to status)
        clear ind 
        r.a_grid=[];
        
    end
        
    if s.assets==0
        
        r.log_z_grid=[];
        r.z_P=[];
        r.z_pi=[];
        r.v=[];
        r.pi_d=[];
        r.pi_export=[];

    % Initialize objects
        exporter = zeros(s.T*(1+s.burn)+1,s.N);   

    % Compute asset and export states for t>1    
        for t=1:s.T*(1+s.burn)  

                ind_k = exporter(t,:)+1; % k=1 is new exporter, k=2 is continuing exporter

            % Construct 3-dimensional indexes from unidimensional indexes ind_y, ind_z, and ind_k     
            % The policy functions are multidimensional matrices -> we map ind_y, ind_z, and ind_k to the policy functions using ind      

                ind = s.y_grid_size*s.z_grid_size*(ind_k-1)+s.y_grid_size*(ind_z(t,:)-1)+ind_y(t,:);

            % Update state variables    
                exporter(t+1,:) = r.exporter(ind); % r.exporter is a decision rule given states today, sim.exporter updates simulation   

                if s.Tshocks==1
                    %SEE MORE DETAILED COMMENTS BELOW, IN THE LONG VERSION OF THE CODE
                    %Non exporters continue with high cost
                        ind_y(t+1,r.exporter(ind)==0) = 1;

                    %Transition from high to low trade cost
                        temp1 = sim.tauRND(t+1,:)>=m.Trho; %Switching from high to low
                        temp2 = r.exporter(ind)==1; %Export
                        temp3 = ind_y(t,:)==1; %High cost
                        ind_y(t+1,logical(temp1.*temp2.*temp3)) = 2; %1=high, 2=low  

                    %Exporters with low cost, stay with low cost
                        temp4 = ind_y(t,:)==2;
                        ind_y(t+1,logical(temp2.*temp4)) = 2;

                    %Reset trade costs for exiting exporters
                        % we initialized at 1, so no need to reset it
                end

        end

        exporter_decision = exporter(2:end,:); %Export decisions (as opposed to status)
        clear ind 
    end
       
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
       
    if s.assets==1
        % Burn first s.burn*s.T time periods
        ind_a = ind_a(s.T*s.burn+1:end-1,:);
        sim.a = sim.a(s.T*s.burn+1:end-1,:);
        ind_z = ind_z(s.T*s.burn+1:end,:);
        exporter = exporter(s.T*s.burn+1:end-1,:);
        exporter_decision = exporter_decision(s.T*s.burn+1:end,:);
        xhistory = xhistory(s.T*s.burn+1:end,:);

        %With sunk costs, need to take care of extra state (ie. previous export status)
        %Note that ind_k is now a matrix, above it was a vector
        if m.S>0
            ind_k = exporter+1; %k=1 if new exporter, k=2 if continuing exporter
         else
            ind_k = 1; %Status does not matter
        end
        
    elseif s.assets==0
        
        % Burn first s.burn*s.T time periods
    
        ind_z = ind_z(s.T*s.burn+1:end,:);    
        if s.Tshocks==0
            ind_y = ind_y(s.T*s.burn+1:end,:);
        elseif s.Tshocks==1
            sim.tauRND = sim.tauRND(s.T*s.burn+1:end,:);
            ind_y = ind_y(s.T*s.burn+1:end-1,:);
        end
        exporter = exporter(s.T*s.burn+1:end-1,:);
        exporter_decision = exporter_decision(s.T*s.burn+1:end,:);
        xhistory = xhistory(s.T*s.burn+1:end,:);

        % Note that ind_k is now a matrix, above it was a vector
        ind_k = exporter+1; %k=1 if new exporter, k=2 if continuing exporter
        
    end                
            
    
    %Initialize objects
    q_d = zeros(s.T,s.N);
    p_d = zeros(s.T,s.N);
    sim.n_f = zeros(s.T,s.N);
    sim.n_d = zeros(s.T,s.N);

    %Construct 3-dimensional indexes from unidimensional indexes ind_a, ind_z, and ind_k     
    %The policy functions are multidimensional matrices -> we map ind_a, ind_z, and ind_k to the policy functions using ind        
    if s.assets==1
        ind=s.a_grid_size*s.z_grid_size*(ind_k-1)+s.a_grid_size*(ind_z-1)+ind_a;
    elseif s.assets==0
        ind = s.y_grid_size*s.z_grid_size*(ind_k-1)+s.y_grid_size*(ind_z-1)+ind_y;  
    end
    clear ind_a ind_k 
    
    %Productivity     
    z = r.z_grid(ind_z);
    
    r.z_grid=[];
    clear ind_z
    
    %Quantities            
    q_f = r.q_f_export(ind).*exporter_decision;
    q_d(exporter_decision==0) = r.q_d(ind(exporter_decision==0));
    q_d(exporter_decision==1) = r.q_d_export(ind(exporter_decision==1));
       
    %Prices     
    p_f = r.p_f_export(ind).*exporter_decision;
    p_d(exporter_decision==0) = r.p_d(ind(exporter_decision==0));
    p_d(exporter_decision==1) = r.p_d_export(ind(exporter_decision==1));

    %Labor    
    sim.n_f = r.n_f_export(ind).*exporter_decision;
    sim.n_d(exporter_decision==0) = r.n_d(ind(exporter_decision==0));
    sim.n_d(exporter_decision==1) = r.n_d_export(ind(exporter_decision==1));   
    sim.n = sim.n_f + sim.n_d;
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
    
    if s.assets==1
       %External finance
         sim.ext_finance = zeros(s.T,s.N);
         sim.ext_finance(exporter_decision==0) = max(m.w*(m.alpha_1*sim.n_d(exporter_decision==0))-sim.a(exporter_decision==0),0);
         sim.ext_finance(exporter_decision==1) = max(m.w*(m.alpha_1*sim.n_d(exporter_decision==1)...
                + sim.n_f(exporter_decision==1)) + m.w*m.alpha_2*(m.F+(1-exporter(exporter_decision==1))*m.S)-sim.a(exporter_decision==1),0); 
%           sim.ext_finance(exporter_decision==0) = m.w*(m.alpha_1*sim.n_d(exporter_decision==0));
%           sim.ext_finance(exporter_decision==1) = m.w*(m.alpha_1*sim.n_d(exporter_decision==1)+ sim.n_f(exporter_decision==1))...
%                                                   + m.w*m.alpha_2*(m.F+(1-exporter(exporter_decision==1))*m.S);
                                              

         
         clear exporter
         %sim.a=[];
         
         sim.GDP =  sum(sales_f + sales_d,2);
         sim.ext_finance_agg = sum(sim.ext_finance,2);
         sim.ext_finance_premium =  sim.ext_finance_agg./sim.GDP;
         sim.ext_finance_premium_mean = mean(sim.ext_finance_premium);
                  
         sim.ext_finance=[];
         sim.GDP=[];
         sim.ext_finance_agg=[];
         sim.ext_finance_premium=[];
                
        clear exporter_decision;
    end
    
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
            %temp1 = ((m.tau*q_f(xhistory>=1)+q_d(xhistory>=1))./z(xhistory>=1))+m.F+tempS(xhistory>=1)*m.S;
            %temp2 = q_d(xhistory<=-1)./z(xhistory<=-1);
            temp1 = ( m.tau*q_f(xhistory>=1)./(z(xhistory>=1).*(sim.a(xhistory>=1).^m.x)) ).^(1/(1-m.x))+( q_d(xhistory>=1)./(z(xhistory>=1).*(sim.a(xhistory>=1).^m.x)) ).^(1/(1-m.x))+m.F+tempS(xhistory>=1)*m.S;
            temp2 = ( q_d(xhistory<=-1)./(z(xhistory<=-1).*(sim.a(xhistory<=-1).^m.x)) ).^(1/(1-m.x));
            sim.export_size_premium_med = median(temp1(:))/median(temp2(:));   
            clear temp1 temp2 tempS    
            
    %Continuing exporter size premium -> continuing exporters vs new exporters
        %Labor
            %temp1 = ((m.tau*q_f(xhistory>1)+q_d(xhistory>1))./z(xhistory>1))+m.F;            
            %temp2 = ((m.tau*q_f(xhistory==1)+q_d(xhistory==1))./z(xhistory==1))+m.F+m.S;
            temp1 = ( m.tau*q_f(xhistory>1)./(z(xhistory>1).*(sim.a(xhistory>1).^m.x)) ).^(1/(1-m.x))+( q_d(xhistory>1)./(z(xhistory>1).*(sim.a(xhistory>1).^m.x)) ).^(1/(1-m.x))+m.F;
            temp2 = ( m.tau*q_f(xhistory==1)./(z(xhistory==1).*(sim.a(xhistory==1).^m.x)) ).^(1/(1-m.x))+( q_d(xhistory==1)./(z(xhistory==1).*(sim.a(xhistory==1).^m.x)) ).^(1/(1-m.x))+m.F+m.S;
            
            
            sim.cont_export_size_premium_med = median(temp1(:))/median(temp2(:));  
            clear temp1 temp2   
    
    % Continuing non-exporter size premium -> continuing non-exporters vs new
    % non-exporters (stoppers)
    
        %Labor
            %temp1 = q_d(xhistory<-1)./z(xhistory<-1);            
            %temp2 = q_d(xhistory==-1)./z(xhistory==-1);
            
            temp1 = ( q_d(xhistory<-1)./(z(xhistory<-1).*(sim.a(xhistory<-1).^m.x)) ).^(1/(1-m.x));
            temp2 = ( q_d(xhistory==-1)./(z(xhistory==-1).*(sim.a(xhistory==-1).^m.x)) ).^(1/(1-m.x));
            
            sim.cont_nonexport_size_premium_med = median(temp1(:))/median(temp2(:));  
            clear temp1 temp2 
            
    % New non-exporter size premium -> new non-exporters vs exporters
    
        %Labor
            tempS = zeros(size(q_f));
            tempS(xhistory==1) = 1;
            %temp1 = ((m.tau*q_f(xhistory>=1)+q_d(xhistory>=1))./z(xhistory>=1))+m.F+tempS(xhistory>=1)*m.S;
            %temp2 = q_d(xhistory==-1)./z(xhistory==-1);
            
            temp1 = ( m.tau*q_f(xhistory>=1)./(z(xhistory>=1).*(sim.a(xhistory>=1).^m.x)) ).^(1/(1-m.x))+( q_d(xhistory>=1)./(z(xhistory>=1).*(sim.a(xhistory>=1).^m.x)) ).^(1/(1-m.x))+m.F+tempS(xhistory>=1)*m.S;
            temp2 = ( q_d(xhistory==-1)./(z(xhistory==-1).*(sim.a(xhistory==-1).^m.x)) ).^(1/(1-m.x));
            
            sim.new_nonexport_size_premium_med = median(temp2(:))/median(temp1(:));  
            clear tempS temp1 temp2 

    %New exporters / stoppers 
        %tempS = zeros(size(q_f));
        %tempS(xhistory==1) = 1;
        %temp1 = ((m.tau*q_f(xhistory==1)+q_d(xhistory==1))./z(xhistory==1))+m.F+tempS(xhistory==1)*m.S;
        %temp2 = (q_d(xhistory==-1))./z(xhistory==-1);
        temp1 = ( m.tau*q_f(xhistory==1)./(z(xhistory==1).*(sim.a(xhistory==1).^m.x)) ).^(1/(1-m.x))+( q_d(xhistory==1)./(z(xhistory==1).*(sim.a(xhistory==1).^m.x)) ).^(1/(1-m.x))+m.F+m.S;
        temp2 = ( q_d(xhistory==-1)./(z(xhistory==-1).*(sim.a(xhistory==-1).^m.x)) ).^(1/(1-m.x));
                
        
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
    
    
    %clear export_intensity 

    sim.xperiods_sales_d_growth_med_before_x = median(sales_d_growth(xperiods_xhistory_gr<0));
    sim.xperiods_sales_growth_med_before_x = median(sales_growth(xperiods_xhistory_gr<0));

    
    sim.xperiods_sales_d_growth_med = [sim.xperiods_sales_d_growth_med_before_x ;sim.xperiods_sales_d_growth_med];
    sim.xperiods_sales_growth_med = [sim.xperiods_sales_growth_med_before_x ;sim.xperiods_sales_growth_med];

    
    clear xperiods_cohorts; 



%     %Hazard rate -> computed by pooling firms from all cohorts
%     xperiod_count_pool = zeros(1,6);
%     for i=1:6     
%         xperiod_count_pool(i) = sum(sum(xhistory==i));
%     end
%     for i=1:5
%         sim.hazard_pool(i) = (xperiod_count_pool(i)-xperiod_count_pool(i+1))/xperiod_count_pool(i);
%     end
%     sim.hazard_pool_4yr = 1 - xperiod_count_pool(4)/xperiod_count_pool(1);
%     clear xperiod_count_pool 
    
  % Feb 13th/March10 in this version, corrected typo: we should count only firms at t if they
   % are counted at t-1, and firms at t if they can be counted at t+1.

    %Hazard rate -> computed by pooling firms from all cohorts
    xperiod_count_pool1 = zeros(1,6); % Counts only firms at t, if they are observed at t-1
    xperiod_count_pool2 = zeros(1,6); % Counts only firms at t, if they can be observed at t+1
    
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


%     % hazard for stoppers:
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
    
        %Sales size distribution
            sim.sales = sales;
            sim.ss_size1 = prctile(sim.sales(:),25);
            sim.ss_size2 = prctile(sim.sales(:),50);   
            sim.ss_size3 = prctile(sim.sales(:),75);   
            sim.ss = (sim.sales<=sim.ss_size1) + 2*(sim.sales>sim.ss_size1 & sim.sales<=sim.ss_size2) + 3*(sim.sales>sim.ss_size2 & sim.sales<=sim.ss_size3) + 4*(sim.sales>sim.ss_size3);

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
            
    clear sim.laborsize;        
        
   
    
    
    
    return; %Breaks function, goes back to main.m
       
end %s.flag_short
    
    
%% State variables -> Assets and export status (whether the firm exported or not in the previous period)

%Initialize objects
    ind_a = zeros(s.T*(1+s.burn)+1,s.N);
    sim.exporter = zeros(s.T*(1+s.burn)+1,s.N);   
    
%First period initialization
    ind_a(1,1:s.N) = 1; %We initialize every firm with lowest assets
    
%Compute asset and export states for t>1    
    for t=1:s.T*(1+s.burn)  
              
        %With sunk costs, need to take care of extra state (ie. previous export status)
            if m.S>0
               ind_k = sim.exporter(t,:)+1; %k=1 is new exporter, k=2 is continuing exporter
            else
               ind_k =1; %Status does not matter
            end
         
        %Construct 3-dimensional indexes from unidimensional indexes ind_a, ind_z, and ind_k     
        %The policy functions are multidimensional matrices -> we map ind_a, ind_z, and ind_k to the policy functions using ind        
            ind = s.a_grid_size*s.z_grid_size*(ind_k-1)+s.a_grid_size*(ind_z(t,:)-1)+ind_a(t,:);
        
        %Update state variables    
            sim.exporter(t+1,:) = r.exporter(ind); % r.exporter is a decision rule given states today, sim.exporter updates simulation   
            ind_a(t+1,:)= r.a_prime_ind(ind);
        
             
    end
    sim.a=r.a_grid(ind_a); %Asset values
    sim.exporter_decision = sim.exporter(2:end,:); %Export decisions (as opposed to status)
    
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
       
  
    %Endogenous states
    %Also drop last period, which reflects state at period (1+s.burn)*s.T+1 from choice made at period (1+s.burn)*s.T
        ind_a = ind_a(s.T*s.burn+1:end-1,:);
        sim.a = sim.a(s.T*s.burn+1:end-1,:);
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
        
%Construct 3-dimensional indexes from unidimensional indexes ind_a, ind_z, and ind_k     
%The policy functions are multidimensional matrices -> we map ind_a, ind_z, and ind_k to the policy functions using ind        
        ind=s.a_grid_size*s.z_grid_size*(ind_k-1)+s.a_grid_size*(ind_z-1)+ind_a;
     
%Quantities            
     sim.q_f = r.q_f_export(ind).*sim.exporter_decision;
     sim.q_d(sim.exporter_decision==0) = r.q_d(ind(sim.exporter_decision==0));
     sim.q_d(sim.exporter_decision==1) = r.q_d_export(ind(sim.exporter_decision==1));
       
%Prices     
     sim.p_f = r.p_f_export(ind).*sim.exporter_decision;
     sim.p_d(sim.exporter_decision==0) = r.p_d(ind(sim.exporter_decision==0));
     sim.p_d(sim.exporter_decision==1) = r.p_d_export(ind(sim.exporter_decision==1));
        
%Labor
     sim.n_f = r.n_f_export(ind).*sim.exporter_decision;
     sim.n_d(sim.exporter_decision==0) = r.n_d(ind(sim.exporter_decision==0));
     sim.n_d(sim.exporter_decision==1) = r.n_d_export(ind(sim.exporter_decision==1));        
        
     sim.n = sim.n_f + sim.n_d;
     sim.n(sim.xhistory==1) = sim.n(sim.xhistory==1) + m.F + m.S;
     sim.n(sim.xhistory>1) = sim.n(sim.xhistory>1) + m.F;
        
%Sales     
     sim.sales_f = sim.p_f.*sim.q_f;
     sim.sales_d = sim.p_d.*sim.q_d;
     sim.sales = sim.sales_d + sim.sales_f;

%Productivity     
    sim.z = r.z_grid(ind_z); 
    
     
%% Tests

%Share of firms with assets at highest asset state
   sim.share_a_max = sum(sum(sim.a==r.a_grid(end)))/sum(sum(sim.a>0)); 
   sim.share_a_min = sum(sum(sim.a==r.a_grid(1)))/sum(sum(sim.a>0)); 
       
%Volatility across time periods of average productivity and average assets    
    sim.a_mean_std = std(mean(sim.a,2))/mean(mean(sim.a));
    sim.z_mean_std = std(mean(sim.z,2))/mean(mean(sim.z));
     
%% Variables, firm-level statistics, and aggregate statistics
     
%Variables
    %Export intensity     
        sim.export_intensity = sim.sales_f./(sim.sales_d+sim.sales_f);     

    %Growth rates     
        sim.sales_growth   = [(sim.sales(2:end,:)-sim.sales(1:end-1,:))./sim.sales(1:end-1,:)];
        sim.sales_d_growth = [(sim.sales_d(2:end,:)-sim.sales_d(1:end-1,:))./sim.sales_d(1:end-1,:)];
        sim.sales_f_growth = [(sim.sales_f(2:end,:)-sim.sales_f(1:end-1,:))./sim.sales_f(1:end-1,:)];
        sim.a_growth   = [(sim.a(2:end,:)-sim.a(1:end-1,:))./sim.a(1:end-1,:)];
        sim.a_loggrowth = [log(sim.a(2:end,:)./sim.a(1:end-1,:))];      
        
         sim.log_sales_growth   = [log(sim.sales(2:end,:)./sim.sales(1:end-1,:))];
         sim.log_sales_d_growth = [log(sim.sales_d(2:end,:)./sim.sales_d(1:end-1,:))];
         sim.log_sales_f_growth = [log(sim.sales_f(2:end,:)./sim.sales_f(1:end-1,:))];
    
    %External finance (total amount of resources borrowed)
         sim.ext_finance = zeros(s.T,s.N);
         sim.ext_finance(sim.exporter_decision==0) = max(m.w*(m.alpha_1*sim.n_d(sim.exporter_decision==0))-sim.a(sim.exporter_decision==0),0);
         sim.ext_finance(sim.exporter_decision==1) = max(m.w*(m.alpha_1*sim.n_d(sim.exporter_decision==1)...
                + sim.n_f(sim.exporter_decision==1)) + m.w*m.alpha_2*(m.F+(1-sim.exporter(sim.exporter_decision==1))*m.S)-sim.a(sim.exporter_decision==1),0); 
%           sim.ext_finance(sim.exporter_decision==0) = m.w*(m.alpha_1*sim.n_d(sim.exporter_decision==0));
%           sim.ext_finance(sim.exporter_decision==1) = m.w*(m.alpha_1*sim.n_d(sim.exporter_decision==1)+ sim.n_f(sim.exporter_decision==1))...
%                                                   + m.w*m.alpha_2*(m.F+(1-sim.exporter(sim.exporter_decision==1))*m.S);


         sim.ef_growth   = [(sim.ext_finance(2:end,:)-sim.ext_finance(1:end-1,:))./sim.ext_finance(1:end-1,:)];
         sim.ef_loggrowth = [log(sim.ext_finance(2:end,:)./sim.ext_finance(1:end-1,:))]; 
      
         
         
         
         
         
         
%Bernard and Jensen moments         
    sim.sales_growth_4yravg = (1/4)*[(sim.sales_d(4:end,:)+sim.sales_f(4:end,:)-sim.sales_d(1:end-3,:)-sim.sales_f(1:end-3,:))./(sim.sales_d(1:end-3,:)+sim.sales_f(1:end-3,:))];  

    BJ_length = 3;
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
         
%Aggregate statistics  
    %GDP
        sim.GDP = sum(sim.sales_f + sim.sales_d,2);
    
    %External finance
         sim.ext_finance_agg = sum(sim.ext_finance,2);  % Sum over firms (aggregate)
         sim.ext_finance_premium =  sim.ext_finance_agg./sim.GDP;
         sim.ext_finance_premium_mean = mean(sim.ext_finance_premium); % Mean over time periods
    
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
            %temp1 = ((m.tau*sim.q_f(sim.xhistory>=1)+sim.q_d(sim.xhistory>=1))./sim.z(sim.xhistory>=1))+m.F+tempS(sim.xhistory>=1)*m.S;
            %temp2 = sim.q_d(sim.xhistory<=-1)./sim.z(sim.xhistory<=-1);
            temp1 = ( m.tau*sim.q_f(sim.xhistory>=1)./(sim.z(sim.xhistory>=1).*(sim.a(sim.xhistory>=1).^m.x)) ).^(1/(1-m.x))+( sim.q_d(sim.xhistory>=1)./(sim.z(sim.xhistory>=1).*(sim.a(sim.xhistory>=1).^m.x)) ).^(1/(1-m.x))+m.F+tempS(sim.xhistory>=1)*m.S;
            temp2 = ( sim.q_d(sim.xhistory<=-1)./(sim.z(sim.xhistory<=-1).*(sim.a(sim.xhistory<=-1).^m.x)) ).^(1/(1-m.x));
            sim.export_size_premium_med = median(temp1(:))/median(temp2(:));   
            clear temp1 temp2 tempS
            
    %Continuing exporter size premium -> continuing exporters vs new exporters
        %Labor
            %temp1 = ((m.tau*sim.q_f(sim.xhistory>1)+sim.q_d(sim.xhistory>1))./sim.z(sim.xhistory>1))+m.F;            
            %temp2 = ((m.tau*sim.q_f(sim.xhistory==1)+sim.q_d(sim.xhistory==1))./sim.z(sim.xhistory==1))+m.F+m.S;
            temp1 = ( m.tau*sim.q_f(sim.xhistory>1)./(sim.z(sim.xhistory>1).*(sim.a(sim.xhistory>1).^m.x))    ).^(1/(1-m.x))+ ( sim.q_d(sim.xhistory>1) ./(sim.z(sim.xhistory>1 ).*(sim.a(sim.xhistory>1 ).^m.x)) ).^(1/(1-m.x))+m.F;
            temp2 = ( m.tau*sim.q_f(sim.xhistory==1)./(sim.z(sim.xhistory==1).*(sim.a(sim.xhistory==1).^m.x)) ).^(1/(1-m.x))+ ( sim.q_d(sim.xhistory==1)./(sim.z(sim.xhistory==1).*(sim.a(sim.xhistory==1).^m.x)) ).^(1/(1-m.x))+m.F+m.S;
     
            sim.cont_export_size_premium_med = median(temp1(:))/median(temp2(:));  
            clear temp1 temp2
        
            
            
    %New exporter size premium -> new exporters vs non exporters
        %Labor
            %temp1 = ((m.tau*sim.q_f(sim.xhistory==1)+sim.q_d(sim.xhistory==1))./sim.z(sim.xhistory==1))+m.F+m.S;
            %temp2 = sim.q_d(sim.xhistory<=-1)./sim.z(sim.xhistory<=-1);
            temp1 = ( m.tau*sim.q_f(sim.xhistory==1)./(sim.z(sim.xhistory==1).*(sim.a(sim.xhistory==1).^m.x)) ).^(1/(1-m.x))+ ( sim.q_d(sim.xhistory==1)./(sim.z(sim.xhistory==1).*(sim.a(sim.xhistory==1).^m.x)) ).^(1/(1-m.x))+m.F+m.S;
            temp2 = ( sim.q_d(sim.xhistory<=-1)./(sim.z(sim.xhistory<=-1).*(sim.a(sim.xhistory<=-1).^m.x)) ).^(1/(1-m.x));
            
            sim.new_export_size_premium_med = median(temp1(:))/median(temp2(:));   
            clear temp1 temp2
            
            
    % Continuing non-exporter size premium -> continuing non-exporters vs new
    % non-exporters (stoppers)
    
        %Labor
            %temp1 = sim.q_d(sim.xhistory<-1)./sim.z(sim.xhistory<-1);            
            %temp2 = sim.q_d(sim.xhistory==-1)./sim.z(sim.xhistory==-1);
            temp1 = ( sim.q_d(sim.xhistory<-1) ./(sim.z(sim.xhistory<-1 ).*(sim.a(sim.xhistory<-1).^m.x)) ).^(1/(1-m.x));
            temp2 = ( sim.q_d(sim.xhistory==-1)./(sim.z(sim.xhistory==-1).*(sim.a(sim.xhistory==-1).^m.x)) ).^(1/(1-m.x));
            
            sim.cont_nonexport_size_premium_med = median(temp1(:))/median(temp2(:));  
            clear temp1 temp2 
            
            
    % New non-exporter size premium -> new non-exporters vs exporters
    
        %Labor
            tempS = zeros(size(sim.q_f));
            tempS(sim.xhistory==1) = 1;
            %temp1 = ((m.tau*sim.q_f(sim.xhistory>=1)+sim.q_d(sim.xhistory>=1))./sim.z(sim.xhistory>=1))+m.F+tempS(sim.xhistory>=1)*m.S;
            %temp2 = sim.q_d(sim.xhistory==-1)./sim.z(sim.xhistory==-1);
            temp1 = ( m.tau*sim.q_f(sim.xhistory>=1)./(sim.z(sim.xhistory>=1).*(sim.a(sim.xhistory>=1).^m.x)) ).^(1/(1-m.x))+( sim.q_d(sim.xhistory>=1)./(sim.z(sim.xhistory>=1).*(sim.a(sim.xhistory>=1).^m.x)) ).^(1/(1-m.x))+m.F+tempS(sim.xhistory>=1)*m.S;
            temp2 = ( sim.q_d(sim.xhistory==-1)./(sim.z(sim.xhistory==-1).*(sim.a(sim.xhistory==-1).^m.x)) ).^(1/(1-m.x));

            sim.new_nonexport_size_premium_med = median(temp2(:))/median(temp1(:));  
            clear tempS temp1 temp2         
  
            

            
    %New exporters / stoppers 
        %tempS = zeros(size(sim.q_f));
        %tempS(sim.xhistory==1) = 1;
        %temp1 = ((m.tau*sim.q_f(sim.xhistory==1)+sim.q_d(sim.xhistory==1))./sim.z(sim.xhistory==1))+m.F+tempS(sim.xhistory==1)*m.S;
        %temp2 = (sim.q_d(sim.xhistory==-1))./sim.z(sim.xhistory==-1);
        temp1 = ( m.tau*sim.q_f(sim.xhistory==1)./(sim.z(sim.xhistory==1).*(sim.a(sim.xhistory==1).^m.x)) ).^(1/(1-m.x))+( sim.q_d(sim.xhistory==1)./(sim.z(sim.xhistory==1).*(sim.a(sim.xhistory==1).^m.x)) ).^(1/(1-m.x))+m.F+m.S;
        temp2 = ( sim.q_d(sim.xhistory==-1)./(sim.z(sim.xhistory==-1).*(sim.a(sim.xhistory==-1).^m.x)) ).^(1/(1-m.x));
                
    
        sim.new_stoppers_size_premium_med = median(temp1(:))/median(temp2(:));  
        clear temp1 temp2      
       
        
               
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
    sim.xperiods_xhistory_gr = sim.xperiods_xhistory(2:end,:);
    sim.xperiods_cohorts_gr = sim.xperiods_cohorts(2:end,:);
    
%Compute statistics -> Averages and medians across cohorts, by firm age within the cohort
for i = 1:s.CN    
    %Median
        %Level
            sim.xperiods_sales_f_med(i,1) = median(sim.sales_f(sim.xperiods_xhistory==i));
            sim.xperiods_sales_d_med(i,1) = median(sim.sales_d(sim.xperiods_xhistory==i));
            sim.xperiods_sales_med(i,1) = median(sim.sales(sim.xperiods_xhistory==i));
            sim.xperiods_export_intensity_med(i,1) = median(sim.export_intensity(sim.xperiods_xhistory==i));
            sim.xperiods_z_med(i,1) = median(sim.z(sim.xperiods_xhistory==i));
            sim.xperiods_assets_med(i,1) = median(sim.a(sim.xperiods_xhistory==i));  
            
        %Growth       
            sim.xperiods_sales_f_growth_med(i,1) = median(sim.sales_f_growth(sim.xperiods_xhistory_gr==i));
            sim.xperiods_sales_d_growth_med(i,1) = median(sim.sales_d_growth(sim.xperiods_xhistory_gr==i));    
            sim.xperiods_sales_growth_med(i,1) = median(sim.sales_growth(sim.xperiods_xhistory_gr==i));        
            sim.xperiods_a_growth_med(i,1) = median(sim.a_growth(sim.xperiods_xhistory_gr==i));  
            sim.xperiods_ef_growth_med(i,1) = median(sim.ef_growth(sim.xperiods_xhistory_gr==i & isfinite(sim.ef_growth) & sim.ef_growth~=-1)); 
            
    %Mean
        %Log-growth
            sim.xperiods_a_loggrowth_mean(i,1) = mean(sim.a_loggrowth(sim.xperiods_xhistory_gr==i));   
            
            %Average external finance growth computed only for firms that use external finance 
            %(and do not switch from not borrowing to borrowing, or from borrowing to not borrowing)            
            sim.xperiods_ef_loggrowth_mean(i,1) = mean(sim.ef_loggrowth(sim.xperiods_xhistory_gr==i & isfinite(sim.ef_loggrowth) & sim.ef_loggrowth~=-1));             
            
end
sim.xperiods_sales_d_growth_med_before_x = median(sim.sales_d_growth(sim.xperiods_xhistory_gr<0 & sim.xperiods_cohorts_gr==1));
    
sim.xperiods_a_growth_med_before_x = median(sim.a_growth(sim.xperiods_xhistory_gr<0 & sim.xperiods_cohorts_gr==1));
sim.xperiods_a_loggrowth_mean_before_x = mean(sim.a_loggrowth(sim.xperiods_xhistory_gr<0 & sim.xperiods_cohorts_gr==1));
sim.xperiods_ef_growth_med_before_x = median(sim.ef_growth(sim.xperiods_xhistory_gr<0 & sim.xperiods_cohorts_gr==1 & isfinite(sim.ef_growth) & sim.ef_growth~=-1));
sim.xperiods_ef_loggrowth_mean_before_x = mean(sim.ef_loggrowth(sim.xperiods_xhistory_gr<0 & sim.xperiods_cohorts_gr==1 & isfinite(sim.ef_loggrowth) & sim.ef_loggrowth~=-1));

%Add pre-entry period to cohort statistics
    sim.xperiods_sales_d_growth_med = [sim.xperiods_sales_d_growth_med_before_x;sim.xperiods_sales_d_growth_med];
    sim.xperiods_a_growth_med = [sim.xperiods_a_growth_med_before_x;sim.xperiods_a_growth_med];
    sim.xperiods_a_loggrowth_mean = [sim.xperiods_a_loggrowth_mean_before_x;sim.xperiods_a_loggrowth_mean];
    sim.xperiods_ef_growth_med = [sim.xperiods_ef_growth_med_before_x;sim.xperiods_ef_growth_med];
    sim.xperiods_ef_loggrowth_mean = [sim.xperiods_ef_loggrowth_mean_before_x;sim.xperiods_ef_loggrowth_mean];      
    
    
% %Hazard rate -> computed by pooling firms from all cohorts
%     
%     sim.xperiod_count_pool = zeros(1,6);
%     for i=1:6     
%         sim.xperiod_count_pool(i) = sum(sum(sim.xhistory==i));
%     end
%     for i=1:5
%         sim.hazard_pool(i) = (sim.xperiod_count_pool(i)-sim.xperiod_count_pool(i+1))/sim.xperiod_count_pool(i);
%     end
%     sim.hazard_pool_4yr = sim.xperiod_count_pool(4)/sim.xperiod_count_pool(1);
%     

  % Feb 13th/March 10 in this version, corrected typo: we should count only firms at t if they
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
    


    % Statistics for cohorts of firms that exported for at least s.CD years, and did not export for s.CN consecutive years 

    %Compute cohorts
    nonxperiods_cohorts = zeros(s.T,s.N);
    for k=1:s.T-s.CN-s.CD %Cohorts
        nonxperiods_exporters = sim.xhistory(k+s.CD+s.CN-1,:)==-s.CN;
        nonxperiods_cohorts(k:k+s.CD+s.CN-1,nonxperiods_exporters) = 1;
    end
    nonxperiods_xhistory = nonxperiods_cohorts.*sim.xhistory; %First value will be any negative number
    nonxperiods_xhistory_gr = nonxperiods_xhistory(2:end,:);
    
    %Compute statistics -> Averages and medians across cohorts, by firm age within the cohort
    for i = -1:-1:-s.CN    
    %Median
        %Level
            sim.nonxperiods_sales_d_med(-i,1) = median(sim.sales_d(nonxperiods_xhistory==i));
            sim.nonxperiods_sales_med(-i,1) = median(sim.sales_d(nonxperiods_xhistory==i)+sim.sales_f(nonxperiods_xhistory==i));
            
            sim.nonxperiods_z_med(-i,1) = median(sim.z(nonxperiods_xhistory==i));
            sim.nonxperiods_assets_med(-i,1) = median(sim.a(nonxperiods_xhistory==i)); 
           
        %Growth       
            sim.nonxperiods_sales_d_growth_med(-i,1) = median(sim.sales_d_growth(nonxperiods_xhistory_gr==i));    
            sim.nonxperiods_sales_growth_med(-i,1) = median(sim.sales_growth(nonxperiods_xhistory_gr==i)); 
            
            sim.nonxperiods_a_growth_med(-i,1) = median(sim.a_growth(nonxperiods_xhistory_gr==i));  
            sim.nonxperiods_ef_growth_med(-i,1) = median(sim.ef_growth(nonxperiods_xhistory_gr==i & isfinite(sim.ef_growth) & sim.ef_growth~=-1)); 

         %Mean
         %Log-growth
         
            sim.nonxperiods_a_loggrowth_mean(-i,1) = mean(sim.a_loggrowth(nonxperiods_xhistory_gr==i));   
             
            %Average external finance growth computed only for firms that use external finance 
            %(and do not switch from not borrowing to borrowing, or from borrowing to not borrowing)            
            sim.nonxperiods_ef_loggrowth_mean(-i,1) = mean(sim.ef_loggrowth(nonxperiods_xhistory_gr==i & isfinite(sim.ef_loggrowth) & sim.ef_loggrowth~=-1));             
    
       
    end

    sim.nonxperiods_sales_d_growth_med_before_x = median(sim.sales_d_growth(nonxperiods_xhistory_gr>0));
    sim.nonxperiods_sales_growth_med_before_x = median(sim.sales_growth(nonxperiods_xhistory_gr>0));
    sim.nonxperiods_a_growth_med_before_x = median(sim.a_growth(nonxperiods_xhistory_gr>0));
    sim.nonxperiods_a_loggrowth_mean_before_x = mean(sim.a_loggrowth(nonxperiods_xhistory_gr>0));
    sim.nonxperiods_ef_growth_med_before_x = median(sim.ef_growth(nonxperiods_xhistory_gr>0 & isfinite(sim.ef_growth) & sim.ef_growth~=-1));
    sim.nonxperiods_ef_loggrowth_mean_before_x = mean(sim.ef_loggrowth(nonxperiods_xhistory_gr>0 & isfinite(sim.ef_loggrowth) & sim.ef_loggrowth~=-1));
 
 %Add pre-entry period to cohort statistics
    
    sim.nonxperiods_sales_d_growth_med = [sim.nonxperiods_sales_d_growth_med_before_x ;sim.nonxperiods_sales_d_growth_med];
    sim.nonxperiods_sales_growth_med = [sim.nonxperiods_sales_growth_med_before_x ;sim.nonxperiods_sales_growth_med];
    sim.nonxperiods_a_growth_med = [sim.nonxperiods_a_growth_med_before_x;sim.nonxperiods_a_growth_med];
    sim.nonxperiods_a_loggrowth_mean = [sim.nonxperiods_a_loggrowth_mean_before_x;sim.nonxperiods_a_loggrowth_mean];
    sim.nonxperiods_ef_growth_med = [sim.nonxperiods_ef_growth_med_before_x;sim.nonxperiods_ef_growth_med];
    sim.nonxperiods_ef_loggrowth_mean = [sim.nonxperiods_ef_loggrowth_mean_before_x;sim.nonxperiods_ef_loggrowth_mean];      
% 
% % Hazard for stoppers -> computed by pooling firms from all cohorts
%    
%     sim.nonxperiod_count_pool = zeros(1,6);
%     for i=-1:-1:-6     
%         sim.nonxperiod_count_pool(-i) = sum(sum(sim.xhistory==i));
%     end
%     for i=-1:-1:-5
%         sim.nonxhazard_pool(-i) = (sim.nonxperiod_count_pool(-i)-sim.nonxperiod_count_pool(-i+1))/sim.nonxperiod_count_pool(-i);
%     end    
    
  % Feb 13thMarch 10th in this version, corrected typo: we should count only firms at t if they
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


%Sales statistics   

    % calculating the median of sales statistics for each firm (within)
    sales_var = sim.sales;
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
        
%Sales growth statistics   

    % calculating the median of sales growth statistics for each firm (within)
    sales_var = sim.log_sales_growth;
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
        
%Exporters vs non-exporters 
    temp1 = (sim.sales_d(sim.xhistory>=1)+sim.sales_f(sim.xhistory>=1));
    temp2 = sim.sales_d(sim.xhistory<=-1);
    sim.export_size_premium_sales_med = median(temp1(:))/median(temp2(:));          
    clear temp1 temp2    

%Continuing exporters vs new exporters
    temp1 = (sim.sales_d(sim.xhistory>1)+sim.sales_f(sim.xhistory>1));
    temp2 = (sim.sales_d(sim.xhistory==1)+sim.sales_f(sim.xhistory==1));
    sim.export_size_premium_cont_sales_med = median(temp1(:))/median(temp2(:));   
    clear temp1 temp2    

%Sales distribution with 3 bins    
    b = 0.33; %Size of the middle bin
    lb = (1-b)/2;
    ub = (1+b)/2;
    
    sim.sales_sort = sort((sim.sales(:)));
    sim.sales_small = median(sim.sales_sort(1:round(lb*s.T*s.N))); 
    sim.sales_medium = median(sim.sales_sort(round(lb*s.T*s.N):round(ub*s.T*s.N))); 
    sim.sales_large = median(sim.sales_sort(round(ub*s.T*s.N):end)); 
    
    sim.largemedium_size_premium_med = sim.sales_large/sim.sales_medium;
    sim.mediumsmall_size_premium_med = sim.sales_medium/sim.sales_small;    
    
    
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
        
    
    
    
end


