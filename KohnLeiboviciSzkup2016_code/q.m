function [q fval exit] = q(m,a,z,TC,num,denom,y0)

    % This functions saves for optimal quantity in the problem of the firm 
    % producing for both markets

    options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8);
%     % Initial guess set y0 = q_unconstrained_domestic_;
%     y0 = ( ((m.sigma - 1)/m.sigma) * (z/m.w) )^m.sigma * m.P_d^(m.sigma)...
%     *  m.Q_d;
   
    function y1 = fun1(x)
        
    % denom already has tau in it
    y1 = m.alpha_1 * x + m.tau * ((   (1/m.alpha_1) * (denom * x^(-1/m.sigma) - (m.tau*m.w/z) ) + m.tau * m.w/z  )...
         /num)^(-m.sigma)    -  z * (m.lambda * a - m.alpha_2 * m.w * TC)/m.w;
        
    end

    [q fval exit] = fsolve(@fun1,y0,options); 


end