function [q fval exit] = q_K(m,a,z,TC, exp_f,exp_d,exp1,y0)

% This functions saves for optimal quantity in the problem of the firm 
    % producing for both markets, with Capital

    options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8);
%     % Initial guess set y0 = q_unconstrained_domestic_;
   
    function y1 = fun1(q)
        
    y1 = m.alpha_1 * q ^ (1/(1-m.x)) + (m.tau^(1/(1-m.x))) * ( (   (1/m.alpha_1) * exp_d * q^(-1/exp1)  + 1 -  (1/m.alpha_1) )...
         /exp_f)^(-exp1/(1-m.x))    -  (z * a^m.x)^(1/(1-m.x)) * (m.lambda * a - m.alpha_2 * m.w * TC)/m.w;
        
    end

    [q fval exit] = fsolve(@fun1,y0,options); 


end