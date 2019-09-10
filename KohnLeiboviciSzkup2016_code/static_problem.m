function r = static_problem(m,s,r)

% This function computes the static problem of the firm.
% It takes as input the model parameters (m.) and solution parameters (s.), 
% and returns as output the q,p,n and pi of the firm's problem in each of
% the following cases:
%
% (1) Domestic unconstrained firm
% (2) Domestic constrained firm
% (3) Exporter unconstrained firm
% (4) Exporter constrained firm
%
% If additionally s.flag_microimplications==1, it returns gains
% decomposition and extensive/intensive constraints measures.
%
% The structure of the code is the following:
% 
% First initialize matrices.
%
% Then, main code: for each state (a,z,x),
%
% (i) First solve the problem for the domestic unconstrained firm.
% (i.1) Check whether the solution satisfies the borrowing constraint. If
% it does, save the results to '*_d' matrices, and go to (iii), else go to (ii).
% (ii) Solve the problem for the domestic constrained firm. Save the
% results to '*_d' matrices. Go to (iii).
% (iii) Check if given the assets a, the firm can at least cover the fixed
% costs. If it cannot, save results to '*_export' matrices (=0 for each q_i,p_i,n_i,pi, i=d,f), if
% it can go to (iv).
% (iv) Solve the problem for the firm producing for both markets,
% unconstrained.
% (iv.1) Check whether the solution satisfies the borrowing constraint. If
% it does, save the results to '*_export' matrices, and go to (vi), else go to (v).
% (v) Solve the problem for the firm producing for both markets, and
% constrained. Save the results to '*_export' matrices. Go to (vi).
% (vi) Loop ends. Save results that we want to keep to 'r.*'.
% (vii) Solve new loop only if s.flag_microimplications==1
%
% WARNING: if we run the model with sunk costs then the indicator of
% export status and Lagrange multiplier and constraint variables
% are not correct
% In the case of sunk cost one cannot use r.constrained to determine
% whether a firm is exporter and is constrained or not. 

    
    
%% Initializing matrices
    
% p, q, n and profits (pi) for firms producing only for domestic market

p_unc_d = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);   
p_const_d = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
q_unc_d = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
q_const_d = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
n_unc_d = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
n_const_d = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
pi_unc_d = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
pi_const_d = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 

q_d  = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 
p_d  = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 
n_d  = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
pi_d = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 


% p, q, n and profits (pi) for firms producing for both markets,
% overall

q_d_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 
q_f_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 
p_d_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 
p_f_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 
n_d_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 
n_f_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 
n_export   = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 
pi_export  = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 


% p, q, n and profits (pi) for firms producing for both markets,
% unconstrained

q_unc_d_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
q_unc_f_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
p_unc_d_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
p_unc_f_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
n_unc_d_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
n_unc_f_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
n_unc_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
pi_unc_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
    
% p, q, n and profits (pi) for firms producing for both markets,
% constrained

q_const_d_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
q_const_f_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
p_const_d_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
p_const_f_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
n_const_export   = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
n_const_d_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 
n_const_f_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size); 
pi_const_export = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);    

% Variables determining whether firms are constrained or not
constrained = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
constrained_d = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
constrained_f = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);

% Variables storing values from function q.m
exit = 9*ones(s.a_grid_size,s.z_grid_size,s.x_grid_size);
fval = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
% fval is absolute difference from 0 (from fsolve), and exit is 'exit flag'

eexitflag = 0;

%% Calculating z_bar: threshold below which firm always produces domestically
% Only valid for financial frictions model (m.S=0)

    %if m.S==0 || m.S>0;
    if m.S==0

    z_bar = @(y) 1/(m.sigma-1) * (m.w*y*m.F)^(1/(m.sigma-1)) * (m.sigma/(m.P_f*m.Q_f^(1/m.sigma)))^(m.sigma/(m.sigma-1))*m.w*m.tau;        
        
    % For given z, the level of assets below which the firm that exports is constrained:
    a_ec = @(z,y) (m.w*y*m.F/m.lambda)+ (1/m.lambda)*(z/m.w)^(m.sigma-1)*((m.sigma-1)/(m.sigma))^(m.sigma)*...
        (m.alpha_1* m.P_d^m.sigma* m.Q_d+m.P_f^m.sigma* m.Q_f /(m.tau^(m.sigma-1)));

    % For given a, the level of productivity above which the firm that produces domestically is constrained:
    z_dc = @(a) ((m.lambda*a/(m.alpha_1*m.w)))^(1/(m.sigma-1)) * (  (m.sigma/(m.sigma-1))* (m.w/(m.P_d*m.Q_d^(1/m.sigma))) )^(m.sigma/(m.sigma-1));

    % For given z, the level of assets below which the firm that produces domestically is constrained:
    a_dc = @(z) (m.alpha_1/m.lambda)*(z/m.w)^(m.sigma-1) *((m.sigma-1)/m.sigma)^m.sigma *m.P_d^m.sigma* m.Q_d; 


    % Threshold levels for exporting decision

    % Threshold level of productivity to export, for firms unconstrained in both markets
    z_1 = z_bar;
    a_1 = @(y) a_ec(z_1(y),y);

    % First order condition for domestic labor, exporting constrained firm
        f_n = @(x,z,y) ((m.sigma-1)/m.sigma )*x(1)^(-1/m.sigma) * z^((m.sigma-1)/m.sigma) * m.P_d*m.Q_d^(1/m.sigma)...
         - m.alpha_1*((m.sigma-1)/m.sigma) *( ((m.lambda*x(2)) -m.alpha_1*m.w*x(1) -m.w*y*m.F)/m.w )^(-1/m.sigma) *...
         (z/m.tau)^((m.sigma-1)/m.sigma)*m.P_f*m.Q_f^(1/m.sigma) -(1-m.alpha_1)* m.w ; %A0

    % Total labor wage bill under unconstrained production
        ntot_unc = @(x,z) m.alpha_1*((z/m.w)^(m.sigma-1))* ( (((m.sigma-1)/m.sigma)^m.sigma) *m.P_d^m.sigma*m.Q_d   )  +((z/(m.w*m.tau))^(m.sigma-1))*...
            ( (((m.sigma-1)/m.sigma)^m.sigma) *m.P_f^m.sigma*m.Q_f   );

    % Total labor wage bill under constrained production
        ntot_con = @(x,z,y) m.lambda*x(2)-m.w*y*m.F;

     % Condition that has to be hold so that the firm is unconstrained in export market:
        cond = @(x,z,y) (ntot_con(x,z,y)>ntot_unc(x,z)-0.1)*1000000;


    % Export threshold for a firm that is constrained domestically
        z_bar1 = @(x,z,y) - ((m.lambda*x(2)*z /(m.alpha_1*m.w))^((m.sigma-1)/m.sigma)*m.P_d*m.Q_d^(1/m.sigma)-m.w*m.lambda*x(2)/(m.alpha_1*m.w))...
                 +(x(1)*z)^((m.sigma-1)/m.sigma)*m.P_d*m.Q_d^(1/m.sigma)+((m.lambda*x(2) -m.alpha_1* m.w*x(1)-m.w*y*m.F)/m.w)^((m.sigma-1)/m.sigma)...
                 *(z/m.tau)^((m.sigma-1)/m.sigma)*m.P_f*m.Q_f^(1/m.sigma) -m.lambda*x(2) -(1-m.alpha_1)* m.w*x(1); %A1

    % Export threshold for a firm that is unconstrained domestically
        z_bar2 = @(x,z,y) -( ( ((m.sigma-1)^(m.sigma-1))/(m.sigma^m.sigma) )*m.Q_d*m.P_d^m.sigma*( (z/m.w)^(m.sigma-1) ) )...
                 +(x(1)*z)^((m.sigma-1)/m.sigma)*m.P_d*m.Q_d^(1/m.sigma)+((m.lambda*x(2) -m.alpha_1* m.w*x(1)-m.w*y*m.F)/m.w)^((m.sigma-1)/m.sigma)...
                 *(z/m.tau)^((m.sigma-1)/m.sigma)*m.P_f*m.Q_f^(1/m.sigma) -m.lambda*x(2) -(1-m.alpha_1)* m.w*x(1); %A2

    % Export threshold:    
        z_bar_fun = @(x,z,y) z_bar2(x,z,y)*(z<z_dc(x(2)))+z_bar1(x,z,y)*(z>=z_dc(x(2)));

    % Function that gives the threshold:    
        fun = @(x,z,y) [f_n(x,z,y) z_bar_fun(x,z,y) cond(x,z,y)];

    end


%% Main loop of static problem

for k = 1:s.x_grid_size  % loop over export status
    
    count=1;
    a_z=max(r.a_grid);
    
    for j = 1:s.z_grid_size  % loop over productivity       
       
        
        z=r.z_grid(j);
        y=r.y_grid(j);       
        
        if mod(j,s.z_grid_size_original+1)==0
            count = 1;
        end
        
        if m.S==0
            
            if z>=z_bar(y) 

                % initial condition
                a0=a_1(y); 
                n0=( ( (m.sigma - 1)/m.sigma  * r.z_grid(j)...
                               * m.P_d / m.w  ) ^ m.sigma) * m.Q_d; %unconstrained domestic labor
                x0=[n0, a0];

                % initial conditions using ealier solution
                if count>1 
                    x0=sol;
                end

                % solves for a and labor
                options = optimset('Display','off','Algorithm','levenberg-marquardt','TolX',1e-7);
                [x ffval eexitflag] = fsolve(@(x) fun(x,z,y),x0,options);

                % saving a
                sol = x ;
                a_z=sol(2);
                count=count+1;

            end
            
        end
       
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % 1. Domestic Market (unconstrained)  -Step (i)-               
          % See 'Case 1: Domestic unconstrained firm' in paper
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
          % quantity if unconstrained
          
          q_unc_d(:,j,k) = ( ( (m.sigma - 1)/m.sigma  * r.z_grid(j)...
                           * m.P_d / m.w  ) ^ m.sigma) * m.Q_d;
                    
          % price if unconstrained
          p_unc_d(:,j,k) = q_unc_d(:,j,k).^(-1/m.sigma) ...
                           .* (m.Q_d )^(1/m.sigma) .* m.P_d;
                    
          % labour if unconstrained
          n_unc_d(:,j,k) = q_unc_d(:,j,k) ./ r.z_grid(j);
                    
          % profits if unconstrained
          pi_unc_d(:,j,k) = p_unc_d(:,j,k) .* q_unc_d(:,j,k) - m.w .* n_unc_d(:,j,k);

          
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % 3. Domestic and Foreign (unconstrained) -Step (iv)-  
          % See 'Case 3: Exporting unconstrained firm' in paper
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
              
             % Quantities
             q_unc_d_export(:,j,k) = ( ( (m.sigma - 1)/m.sigma...
                                  * r.z_grid(j)* m.P_d / m.w  )...
                                  ^ m.sigma) * m.Q_d;
                              
             q_unc_f_export(:,j,k) = ( ( (m.sigma - 1)/m.sigma...
                                  * r.z_grid(j)* m.P_f / (m.tau * m.w))...
                                  ^ m.sigma) * m.Q_f;
                    
             % Domestic and export prices
          
             p_unc_d_export(:,j,k) = q_unc_d_export(:,j,k).^(-1/m.sigma)...
                                  .*m.Q_d^(1/m.sigma) .* m.P_d;
                    
             p_unc_f_export(:,j,k) = q_unc_f_export(:,j,k).^(-1/m.sigma)...
                                  .*m.Q_f^(1/m.sigma) .* m.P_f;
                    
             % Domestic, export and total labor
             n_unc_d_export(:,j,k) = q_unc_d_export(:,j,k)./r.z_grid(j);
             n_unc_f_export(:,j,k) = m.tau .* q_unc_f_export(:,j,k)./r.z_grid(j);
             n_unc_export(:,j,k) = n_unc_d_export(:,j,k) + n_unc_f_export(:,j,k);
                    
             % Profits
             pi_unc_export(:,j,k) = p_unc_d_export(:,j,k) .* q_unc_d_export(:,j,k)...
                               + p_unc_f_export(:,j,k) .* q_unc_f_export(:,j,k)...
                               - m.w.*n_unc_export(:,j,k) - m.w .* y .* m.TC(k);                   

       for i = 1:s.a_grid_size  % loop over assets
                
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
          % Checking if the constraint is not binding - Step (i.1)-
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          if m.alpha_1 * m.w * n_unc_d(i,j,k) < m.lambda * r.a_grid(i)  
                        
             q_d(i,j,k) = q_unc_d(i,j,k); % quantity
             p_d(i,j,k) = p_unc_d(i,j,k); % price
             n_d(i,j,k) = n_unc_d(i,j,k); % labour
             pi_d(i,j,k) = pi_unc_d(i,j,k); % profits   
             constrained(i,j,k) = 0; % indicator   (overall)
             constrained_d(i,j,k) = 0; % Indicator for specific market
                        
          else  % If the constraint is binding
                    
              
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % 2. Domestic Market (constrained)  -Step (ii)-               
          % See 'Case 2: Domestic constrained firm' in paper
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
             % labour if constrained
             n_const_d(i,j,k) = (m.lambda * r.a_grid(i))/(m.w * m.alpha_1);
                    
             % quantity if constrained
             q_const_d(i,j,k) = r.z_grid(j) * (n_const_d(i,j,k));
                    
             % price if constrained
             p_const_d(i,j,k) = (q_const_d(i,j,k))^(-1/m.sigma)...
                                 * (m.Q_d )^(1/m.sigma) * m.P_d;
                    
             % profits if constrained
             pi_const_d(i,j,k) = p_const_d(i,j,k) * q_const_d(i,j,k)...
                 - m.w * n_const_d(i,j,k);
                    
             % Saving results    
              q_d(i,j,k) = q_const_d(i,j,k); % quantity
              p_d(i,j,k) = p_const_d(i,j,k); % price
              n_d(i,j,k) = n_const_d(i,j,k); % labour
              pi_d(i,j,k) = pi_const_d(i,j,k); % profits   
              constrained_d(i,j,k) = 1; % indicator (for specific market)
                    
          end
                             
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Domestic and Foreign Sales               
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Check if given the assets a, the firm can at least 
          % cover the fixed costs.   -Step (iii)-            
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
          
          %Updated this part of the code on March 13th
          %We now check that the asset thresholds obtained above are actual
          %solutions to the system of equations that determine them
          %(eexitflag>0)
          if m.lambda * r.a_grid(i) - m.alpha_2 * m.w * y * m.TC(k) <= 0 ...
                   || (m.S==0 && r.z_grid(j)<z_bar(y))  || ( (m.S==0 && r.a_grid(i)<a_z)  ...
                   && eexitflag>0 )
                    
                    
             q_d_export(i,j,k) = 0; % domestic quantity
             q_f_export(i,j,k) = 0; % foreign quantity
             q_f_export(i,j,k) = 0; % foreign quantity
             p_d_export(i,j,k) = 0; % domestic price
             p_f_export(i,j,k) = 0; % foreign price
             n_d_export(i,j,k) = 0; % labour domestic
             n_f_export(i,j,k) = 0; % labour export
             pi_export(i,j,k)  = 0; % profits
             
             if m.S==0 && r.z_grid(j)>=z_bar(y)
                constrained_f(i,j,k) = 1; % indicator (for specific market)
             end
                
          else 
                    
 
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
          % Checking if the constraint is not binding - Step (iv.1)-
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
             if m.w * (m.alpha_1 * n_unc_d_export(i,j,k) + n_unc_f_export(i,j,k) )...
                     <= m.lambda * r.a_grid(i) - m.alpha_2 * m.w * y * m.TC(k) 
                        
                q_d_export(i,j,k) = q_unc_d_export(i,j,k); % domestic quantity
                q_f_export(i,j,k) = q_unc_f_export(i,j,k); % foreign quantity
                p_d_export(i,j,k) = p_unc_d_export(i,j,k); % domestic price
                p_f_export(i,j,k) = p_unc_f_export(i,j,k); % foreign price
                n_d_export(i,j,k) = n_unc_d_export(i,j,k); % labour domestic
                n_f_export(i,j,k) = n_unc_f_export(i,j,k); % labour foreign
                pi_export(i,j,k) = pi_unc_export(i,j,k); % profits
                constrained_f(i,j,k) = 0;
            
             else           

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % 4. Domestic and Export (constrained) -Step (v)-  
          % See 'Case 4: Exporting constrained firm' in paper
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful expressions (q.m, below):
num = (1-1/m.sigma) * m.P_f * (m.Q_f) ^ (1/m.sigma);
denom = m.tau * (1-1/m.sigma) * m.P_d * ( m.Q_d) ^ (1/m.sigma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
             
             if j>1 && q_const_d_export(i,j-1,k)~=0 
                 
                 y0= q_const_d_export(i,j-1,k);       
                 
             else
                 
                 % Initial guess set y0 = q_unconstrained_domestic;
                 y0 = ( ((m.sigma - 1)/m.sigma) * (r.z_grid(j)/m.w) )^m.sigma * m.P_d^(m.sigma)...
                 *  m.Q_d;
                              
             end

             % Domestic output
                [q_const_d_export(i,j,k) fval(i,j,k) exit(i,j,k)] =...
                    q(m,r.a_grid(i),r.z_grid(j),y*m.TC(k), num, denom,y0); % fsolve 
                   %WARNING: WE PASS y*TC INTO q FUNCTION, SO NO NEED TO ADJUST FIXED COST INSIDE q
                               
             % Export output
                q_const_f_export(i,j,k) = ( (  (1/m.alpha_1) * (denom * ...
                    q_const_d_export(i,j,k)^(-1/m.sigma) - m.tau * m.w/ r.z_grid(j))...
                    + m.tau * m.w/r.z_grid(j) ) /num)^(-m.sigma);
                    
             % Domestic and export prices
                p_const_d_export(i,j,k) = q_const_d_export(i,j,k)^(-1/m.sigma)...
                    *m.Q_d^(1/m.sigma) * m.P_d;
                            
                p_const_f_export(i,j,k) = q_const_f_export(i,j,k)^(-1/m.sigma)...
                    *m.Q_f^(1/m.sigma) * m.P_f;
                            
             % Domestic, foreign and total labor
                n_const_d_export(i,j,k) = q_const_d_export(i,j,k)/r.z_grid(j);
                            
                n_const_f_export(i,j,k) = m.tau * q_const_f_export(i,j,k)/r.z_grid(j);
                            
                n_const_export(i,j,k) =n_const_d_export(i,j,k) + n_const_f_export(i,j,k);
                
             % Profits
                pi_const_export(i,j,k) = p_const_d_export(i,j,k) * q_const_d_export(i,j,k)...
                            + p_const_f_export(i,j,k) * q_const_f_export(i,j,k)...
                            - m.w * n_const_export(i,j,k) - m.w * y * m.TC(k);
                                                      
             % Saving the results:                            
                q_d_export(i,j,k) = q_const_d_export(i,j,k); % domestic quantity
                q_f_export(i,j,k) = q_const_f_export(i,j,k); % foreign quantity
                p_d_export(i,j,k) = p_const_d_export(i,j,k); % domestic price
                p_f_export(i,j,k) = p_const_f_export(i,j,k); % foreign price
                n_export(i,j,k) = n_const_export(i,j,k); % labour
                n_d_export(i,j,k) = n_const_d_export(i,j,k); % labour domestic
                n_f_export(i,j,k) = n_const_f_export(i,j,k); % labour foreign
                pi_export(i,j,k) = pi_const_export(i,j,k); % profits                            
                constrained_f(i,j,k) = 1; % indicator (for specific market)
                     
                            
                        
                        
             end 
                        
                        
                        
          end
                            
                    
                    
      end % end of the loop over assets
   end % loop over productivity
end % loop over export status

 
%% Saving results that want to keep outside the function

% p, q, n and profits (pi) for firms producing only for domestic market

r.q_d  = q_d; 
r.p_d  = p_d; 
r.n_d  = n_d;
r.pi_d = pi_d; 

% p, q, n and profits (pi) for firms producing for both markets,
% overall

r.q_d_export = q_d_export; 
r.q_f_export = q_f_export; 
r.p_d_export = p_d_export; 
r.p_f_export = p_f_export; 
r.n_d_export = n_d_export; 
r.n_f_export = n_f_export; 
r.pi_export  = pi_export;      
        
% Variables determining whether firms are constrained or not

r.constrained = constrained;
r.constrained_d = constrained_d;
r.constrained_f = constrained_f;
        
        
end %end function

