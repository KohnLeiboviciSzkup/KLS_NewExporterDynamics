function [r sim] = model_simulate(m,s,r)

    %Baseline model
    if s.model_baseline_financialfrictions==1 || s.model_baseline_sunkcosts==1 || s.model_extensions_homogeneousF_ff==1 || s.model_extensions_homogeneousF_sc==1
        [r sim] = model_simulate_baseline(m,s,r);
    end

    %Extensions without assets
    if s.model_extensions_Dshocks==1 || s.model_extensions_Tshocks==1 || s.model_extensions_Fshocks==1 
        [r sim] = model_simulate_Sy(m,s,r);
    end

    %Extension with labor adjustment costs
    if s.model_extensions_laboradjustmentcosts==1
        [r sim] = model_simulate_Nadj(m,s,r);
    end

    %Extension with capital        
    if s.model_extensions_capital==1
        [r sim] = model_simulate_K(m,s,r);
    end



end

