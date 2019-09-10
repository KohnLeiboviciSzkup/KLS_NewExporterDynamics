function tau_eq = trade_lib_tariff_equivalent(tau_ini,p,GDP_growth_sc,GDP_ff,lambda)

format long

% Parameters:
flag_ff = 1;
flag_sc = 0;


% Guess:
    taug = tau_ini-1.5*p;
    
% options:
    options = optimset('Display','on','TolX',1e-6);
    
% Solution:
    tau_eq = fzero(@fun1,taug,options);
    
% Auxiliary function:
    function z = fun1(tau)
        
        % liberalization
        x = trade_lib_fun(tau,flag_ff,flag_sc,lambda);
        
        % growth
        GDP_growth = (x(1)/GDP_ff-1)*100;
        
        % comparison:
        z = GDP_growth-GDP_growth_sc;
        
        display([z tau])
    end


    
end