
%% Setup flags

    %% Check

    if (s.model_baseline_financialfrictions+s.model_baseline_sunkcosts+s.model_extensions_Dshocks+s.model_extensions_Tshocks+s.model_extensions_Fshocks+s.model_extensions_laboradjustmentcosts+s.model_extensions_capital+s.model_extensions_homogeneousF_ff+s.model_extensions_homogeneousF_sc>1)
        disp('Run only one model at a time! ');
        return; %break;
    end
    if (s.model_baseline_financialfrictions+s.model_baseline_sunkcosts+s.model_extensions_Dshocks+s.model_extensions_Tshocks+s.model_extensions_Fshocks+s.model_extensions_laboradjustmentcosts+s.model_extensions_capital+s.model_extensions_homogeneousF_ff+s.model_extensions_homogeneousF_sc==0)
        disp('Choose a model! ');
        return; %break;
    end

    %% Assign flags

    s.assets = 0;                % Model with assets as a state (Sunk cost or financial frictions model)
    s.assets_Fshocks = 0;        % Shocks to fixed costs + assets (only financial friction)
    s.Fshocks = 0;               % Shocks to fixed cost. It includes the benchmark sunk cost model
    s.F_risk = 0;                %=0 if heterogeneous costs, =1 if idiosyncratic shocks to the costs -> Only for Fshocks, assets_Fshocks models
    s.Dshocks = 0;               % Shocks to foreign demand. 
    s.Tshocks = 0;               % Shocks to Iceberg costs. 
    s.Nadj = 0;                  % Labor adjustment costs.
    s.NadjFN = 0;                % Type of labor adjustment costs. 
    s.FE = 0;                    % Non-stochastic heterogeneous fixed costs for Dshocks, Tshocks, Nadj models
    s.flag_capital = 0;          % Solves problem with assets as input of production
    s.flag_unc = 0;              % Computes statistics for 'more financially constrained' firms vs 'less financially constrained' firms
    
    
    
    if s.model_baseline_financialfrictions==1
        s.assets = 1;
        s.assets_Fshocks = 1;     
        s.flag_unc = 1;
    elseif s.model_baseline_sunkcosts==1
        s.assets = 1;
        s.assets_Fshocks = 1;
    elseif s.model_extensions_Dshocks==1
        s.Dshocks = 1;
        s.FE = 1;
    elseif s.model_extensions_Tshocks==1
        s.Tshocks = 1;
        s.FE = 1;    
    elseif s.model_extensions_Fshocks==1
        s.Fshocks = 1;
        s.F_risk = 1;
    elseif s.model_extensions_laboradjustmentcosts==1
        s.Nadj = 1;
        s.NadjFN = 1; 
        s.FE = 1;
    elseif s.model_extensions_capital==1
        s.assets = 1;
        s.assets_Fshocks = 1;
        s.flag_capital = 1;
    elseif s.model_extensions_homogeneousF_ff==1
        s.assets = 1;
        s.assets_Fshocks = 1;   
    elseif s.model_extensions_homogeneousF_sc==1
        s.assets = 1;    
        s.assets_Fshocks = 1;       
    end