function x = trade_lib_fun(tau,flag_ff,flag_sc,lambda)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the Main file to replicate the results in the paper Kohn,       %
% Leibovici, Szkup (2012).                                                %
%                                                                         % 
% The code needs the following m-files:                                   %
% main.m, model_solve.m, model_simulate.m, shocks.m, static_problem.m,    %
% dynamic_problem.m, q.m                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start


    %% Choose model

    % Paper
    s.model_baseline_financialfrictions = 0; %flag_ff;
    s.model_baseline_sunkcosts = 0; %flag_sc;

    % Appendix 
    s.model_extensions_Dshocks = 0;
    s.model_extensions_Tshocks = 0;
    s.model_extensions_Fshocks = 0;
    s.model_extensions_laboradjustmentcosts = 0;
    s.model_extensions_capital = 0;
    s.model_extensions_homogeneousF_ff = flag_ff;
    s.model_extensions_homogeneousF_sc = flag_sc;

    %% Calibrate model

    setup_flags;
    setup_calibration;
    setup_solution;

    
    %% Trade Liberalization:
    m.tau = tau;
        
    if flag_ff == 1
       m.lambda = lambda; 
    end
        
        
    
%% Solve model

[r s] = model_solve(m,s);


%% Simulate model

tic;



if s.Nadj==0
    
    if s.assets==1
        [~, sim] = model_simulate(m,s,r);
    else
        [~, sim] = model_simulate_Sy(m,s,r);
    end
    
elseif s.assets==0 && s.Nadj==1

    [~, sim] = model_simulate_Sadj(m,s,r);

end

disp('The simulation takes...');
toc;

%% GDP

sales_d = sim.tot_sales_d_avg;
sales_f = sim.tot_sales_f_avg;
GDP = mean(sim.GDP);
exports_per_firm = sim.exports_per_firm_avg;
exports_per_firm_med = sim.exports_per_firm_med;
share_exporters = sim.share_exporters;

x = [GDP sales_f sales_d exports_per_firm exports_per_firm_med share_exporters];


end
