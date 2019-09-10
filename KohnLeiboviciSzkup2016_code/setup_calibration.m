%% Setup calibration

    %% Common calibration parameters

    m.sigma = 5;  
    m.gamma = 2;  

    m.r = 0.01;
    m.w = 1;   
    m.Q_d = 50; 
    m.P_d = 1; 
    m.Q_f = 50; 
    m.P_f = 1; 

    %% Initialization

    s.x_grid_size = 2; 

    %% Set parameters

    if s.model_baseline_financialfrictions==1


        
        
       % Baseline Calibration
        %Parameters  
        m.lambda = 1.64;         
        m.F = 1.70204; 
        m.S = 0; 
        m.beta = 0.83097;        
        m.alpha_1 = 0.53213; %Working capital requirement for domestic labor costs
        m.alpha_2 = 1; %Working capital requirement for fixed export cost
        m.tau = 1.51294; 

        %Productivity    
        m.log_z_sigma = 0.10559; 
        m.log_z_rho = 0.90909; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1 

        %Fixed costs
        m.log_y_sigma = 0.000000001;  
        m.log_y_rho = 0; % rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_y_mu = -(m.log_y_sigma^2)/((1-m.log_y_rho^2)*2); %Normalize average productivity to 1         
        m.F_share_low = 0.16057;
        m.F_value_low = 0.11341; %High cost = F, low cost = F_value_low*F
        s.y_grid_size = 2;
        s.y_grid_range_sd = 3;
        s.x_grid_size = 1;  
        m.TC = m.F;

% % Alternative Calibration to match export intensity growth
%         %Parameters  
%         m.lambda = 1.64;         
%         m.F = 1.4354; 
%         m.S = 0; 
%         m.beta = 0.63316;        
%         m.alpha_1 = 0.34198; %Working capital requirement for domestic labor costs
%         m.alpha_2 = 1; %Working capital requirement for fixed export cost
%         m.tau = 1.3164; 
% 
%         %Productivity    
%         m.log_z_sigma = 0.11863; 
%         m.log_z_rho = 0.91092; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
%         m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1 
% 
%         %Fixed costs
%         m.log_y_sigma = 0.000000001;  
%         m.log_y_rho = 0; % rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
%         m.log_y_mu = -(m.log_y_sigma^2)/((1-m.log_y_rho^2)*2); %Normalize average productivity to 1         
%         m.F_share_low = 0.17088;
%         m.F_value_low = 0.10999; %High cost = F, low cost = F_value_low*F
%         s.y_grid_size = 2;
%         s.y_grid_range_sd = 3;
%         s.x_grid_size = 1;  
%         m.TC = m.F;


% % Alternative Calibration setting alpha=1, calibrating lambda additionally
%         %Parameters  
%         m.lambda = 1.3185;         
%         m.F = 1.709; 
%         m.S = 0; 
%         m.beta = 0.91334;        
%         m.alpha_1 = 1; %Working capital requirement for domestic labor costs
%         m.alpha_2 = 1; %Working capital requirement for fixed export cost
%         m.tau = 1.5087; 
% 
%         %Productivity    
%         m.log_z_sigma = 0.10216; 
%         m.log_z_rho = 0.91721; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
%         m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1 
% 
%         %Fixed costs
%         m.log_y_sigma = 0.000000001;  
%         m.log_y_rho = 0; % rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
%         m.log_y_mu = -(m.log_y_sigma^2)/((1-m.log_y_rho^2)*2); %Normalize average productivity to 1         
%         m.F_share_low = 0.15712;
%         m.F_value_low = 0.12767; %High cost = F, low cost = F_value_low*F
%         s.y_grid_size = 2;
%         s.y_grid_range_sd = 3;
%         s.x_grid_size = 1;  
%         m.TC = m.F;


% % Alternative Calibration matching agg. export intensity
%         %Parameters  
%         m.lambda = 1.64;         
%         m.F = 3.4838 ; 
%         m.S = 0; 
%         m.beta = 0.84262;        
%         m.alpha_1 = 0.33743; %Working capital requirement for domestic labor costs
%         m.alpha_2 = 1; %Working capital requirement for fixed export cost
%         m.tau = 1.0959; 
% 
%         %Productivity    
%         m.log_z_sigma = 0.09036; 
%         m.log_z_rho = 0.90383; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
%         m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1 
% 
%         %Fixed costs
%         m.log_y_sigma = 0.000000001;  
%         m.log_y_rho = 0; % rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
%         m.log_y_mu = -(m.log_y_sigma^2)/((1-m.log_y_rho^2)*2); %Normalize average productivity to 1         
%         m.F_share_low = 0.11127;
%         m.F_value_low = 0.16314; %High cost = F, low cost = F_value_low*F
%         s.y_grid_size = 2;
%         s.y_grid_range_sd = 3;
%         s.x_grid_size = 1;  
%         m.TC = m.F;

    elseif s.model_baseline_sunkcosts==1

%Baseline calibration        
        %Parameters  
        m.lambda = 1e8;         
        m.F = 1.0252; 
        m.S = 4.1405; 
        m.beta = 0.97782;        
        m.alpha_1 = 0.53213; %Working capital requirement for domestic labor costs
        m.alpha_2 = 1; %Working capital requirement for fixed export cost
        m.tau = 1.5158; 

        %Productivity    
        m.log_z_sigma = 0.17847; 
        m.log_z_rho = 0.78155; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1     

        %Fixed costs
        m.log_y_sigma = 0.000000001;  
        m.log_y_rho = 0; % rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_y_mu = -(m.log_y_sigma^2)/((1-m.log_y_rho^2)*2); %Normalize average productivity to 1             
        m.F_share_low = 0.054789;
        m.F_value_low = 0.22927; %High cost = F, low cost = F_value_low*F
        s.y_grid_size = 2;
        s.y_grid_range_sd = 3;
        m.TC(1) = m.F+m.S;
        m.TC(2) = m.F;
      
% % Alternative Calibration matching agg. export intensity               
%         %Parameters  
%         m.lambda = 1e8;         
%         m.F = 4.9683; 
%         m.S = 5.9449; 
%         m.beta = 0.96790;        
%         m.alpha_1 = 0.33743; %Working capital requirement for domestic labor costs
%         m.alpha_2 = 1; %Working capital requirement for fixed export cost
%         m.tau = 1.063; 
%                               
%     %Productivity    
%         m.log_z_sigma = 0.10865; 
%         m.log_z_rho = 0.80588; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
%         m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1     
%     
%     %Fixed costs
%         m.log_y_sigma = 0.000000001;  
%         m.log_y_rho = 0; % rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
%         m.log_y_mu = -(m.log_y_sigma^2)/((1-m.log_y_rho^2)*2); %Normalize average productivity to 1             
%         m.F_share_low = 0.087069;
%         m.F_value_low = 0.096108; %High cost = F, low cost = F_value_low*F
%         s.y_grid_size = 2;
%         s.y_grid_range_sd = 3;
%         m.TC(1) = m.F+m.S;
%         m.TC(2) = m.F;        
        
    elseif s.model_extensions_Fshocks==1

        %Parameters  
        m.lambda = 1e8;         
        m.F = 1.28566; %%1.280944074970001 ; 
        m.S =13.4443; %13.176521432078587 ; 
        m.beta = 0.97782;        
        m.alpha_1 = 0.53213; %Working capital requirement for domestic labor costs
        m.alpha_2 = 1; %Working capital requirement for fixed export cost
        m.tau = 1.51772; % 1.513611411650500; 

        %Productivity    
        m.log_z_sigma = 0.199644; %0.198355361695761; 
        m.log_z_rho = 0.740256; %0.744905708254327; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1     

        %Fixed costs
        m.log_y_sigma = 0.525946; %0.523358012217615;  
        m.log_y_rho = 0; % rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_y_mu = -(m.log_y_sigma^2)/((1-m.log_y_rho^2)*2); %Normalize average productivity to 1         
        s.y_grid_size = 100;
        s.y_grid_range_sd = 3;
        
    elseif s.model_extensions_Dshocks==1

        %Parameters  
        m.lambda = 1e8;         
        m.F = 1.49730; 
        m.S = 6.88280; 
        m.beta = 0.97782;        
        m.alpha_1 = 0.53213; %Working capital requirement for domestic labor costs
        m.alpha_2 = 1; %Working capital requirement for fixed export cost
        m.tau = 1.43510; 

        %Productivity    
        m.log_z_sigma = 0.15994; 
        m.log_z_rho = 0.80540; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1     

        %Demand shocks
        m.log_y_sigma = 0.91029;  
        m.log_y_rho = 0; % rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_y_mu = -(m.log_y_sigma^2)/((1-m.log_y_rho^2)*2); %Normalize average productivity to 1         

        %Fixed costs
        m.F_share_low = 0.04727;
        m.F_value_low = 0.28273 ; %High cost = F, low cost = F_value_low*F    
        s.y_grid_size = 100;
        s.y_grid_range_sd = 3;

    elseif s.model_extensions_Tshocks==1

        %Parameters  
        m.lambda = 1e8;         
        m.F = 1.5998; 
        m.S = 3.2737; 
        m.beta = 0.97782;     %0.98   
        m.alpha_1 = 0.53213; % 0.53306 %Working capital requirement for domestic labor costs
        m.alpha_2 = 1; %Working capital requirement for fixed export cost
        m.tau =  1.5158;  %1.51

        %Productivity    
        m.log_z_sigma = 0.18004; 
        m.log_z_rho = 0.76645; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1     

        %Fixed costs
        m.F_share_low = 0.036981;
        m.F_value_low = 0.24709; %High cost = F, low cost = F_value_low*F    

        %Trade cost shocks
        s.y_grid_size = 11;
        s.y_grid_range_sd = 1;
        m.Trho = 0.69375; %Probability of transition to lower trade cost
        m.x = 0.089482;
        m.tauL = m.tau - m.x; %Value of lowest trade cost   
        m.tauH = m.tau + m.x;

       

    elseif s.model_extensions_laboradjustmentcosts==1

        %Parameters  
        m.lambda = 1e8;         
        m.F = 2.1162; 
        m.S = 1.7219; 
        m.beta = 0.97782;        
        m.alpha_1 = 0.53213; %0.53306; %Working capital requirement for domestic labor costs
        m.alpha_2 = 1; %Working capital requirement for fixed export cost
        m.tau = 1.5147; 

        %Productivity    
        m.log_z_sigma = 0.16537; 
        m.log_z_rho = 0.83582; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1     

        %Fixed costs
        m.F_share_low = 0.14613;
        m.F_value_low = 0.15777; %High cost = F, low cost = F_value_low*F    

        %Labor adjustment costs
        s.y_grid_size = 200;
        m.xi = 0.060995;

        if s.NadjFN==1
            s.adjCost = @(n,np) m.xi.*((np-n)./n).^2.*n;      
        elseif s.NadjFN==2
            s.adjCost = @(n,np) m.xi.*(np-n).^2; 
        elseif s.NadjFN==3
            s.adjCost = @(n,np) m.xi.*((np-n)./n).^2; 
        elseif s.NadjFN==4
            s.adjCost = @(n,np) m.xi.*(np-n); 
        elseif s.NadjFN==5
            s.adjCost = @(n,np) m.xi.*((np-n)./n); 
        else
            disp('Choose a type of adjustment cost');
            return; %break;               
        end        

    elseif s.model_extensions_capital==1

        %Parameters  
        m.lambda = 1.64;         
        m.F = 0.74962; 
        m.S = 0; 
        m.beta = 0.63935;        
        m.alpha_1 = 0.93233; %Working capital requirement for domestic labor costs
        m.alpha_2 = 1; %Working capital requirement for fixed export cost
        m.tau = 2.4924; 

        %Productivity    
        m.log_z_sigma = 0.23008; 
        m.log_z_rho = 0.70473; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1     

        %Capital
        m.x = 0.3; % Capital share of output
        m.delta = 0.1;
        m.r = -m.delta/(1+m.delta);  % 5 percent depreciation in a'(1+delta) -> 1+m.delta=1/(1+m.r), in dynamic problem                           

        %Fixed costs
        m.F_share_low = 0.14926;
        m.F_value_low = 0.10935; %High cost = F, low cost = F_value_low*F    
        m.TC = m.F;
        s.x_grid_size = 1; 
        s.y_grid_size = 2;

    elseif s.model_extensions_homogeneousF_ff==1

        %Parameters  
        m.lambda = 1.64;         
        m.F = 0.99512; 
        m.S = 0; 
        m.beta = 0.80165;        
        m.alpha_1 = 0.47799; %Working capital requirement for domestic labor costs
        m.alpha_2 = 1; %Working capital requirement for fixed export cost
        m.tau = 1.4547; 

        %Productivity    
        m.log_z_sigma = 0.089389; 
        m.log_z_rho = 0.93881; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1     

        %Fixed costs
        m.F_share_low = 0;
        m.F_value_low = 1; %High cost = F, low cost = F_value_low*F      
        s.x_grid_size = 1;  
        s.y_grid_size = 2;
        m.TC = m.F;

    elseif s.model_extensions_homogeneousF_sc==1
 
        %Parameters  
        m.lambda = 1e8;         
        m.F = 0.89587; %0.89713; 
        m.S = 4.169; %4.14870; 
        m.beta = 0.9235; %0.92463;        
        m.alpha_1 = 0.47799; %Working capital requirement for domestic labor costs
        m.alpha_2 = 1; %Working capital requirement for fixed export cost
        m.tau = 1.5113; 

        %Productivity    
        m.log_z_sigma = 0.16238; %0.16198; 
        m.log_z_rho = 0.77879; %0.77971; %rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_z_mu = -(m.log_z_sigma^2)/((1-m.log_z_rho^2)*2); %Normalize average productivity to 1     

        %Fixed costs
        m.log_y_sigma = 0.000000001;  
        m.log_y_rho = 0; % rho=0 is iid N(mu,sigma^2) else process AR(1), mean mu, with errors N(0,sigma^2)    
        m.log_y_mu = -(m.log_y_sigma^2)/((1-m.log_y_rho^2)*2); %Normalize average productivity to 1             
        m.F_share_low = 0;
        m.F_value_low = 1; %High cost = F, low cost = F_value_low*F
        s.y_grid_size = 2;
        s.y_grid_range_sd = 3;
        m.TC(1) = m.F+m.S;
        m.TC(2) = m.F;    

    end