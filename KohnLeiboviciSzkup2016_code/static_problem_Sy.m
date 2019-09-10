function r = static_problem_Sy(m,s,r)
%% Initializing matrices
    
% p, q, n and profits (pi) for firms producing only for domestic market
q_d  = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size); 
p_d  = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size); 
n_d  = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);
pi_d = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size); 

% p, q, n and profits (pi) for firms producing for both markets,
% overall
q_d_export = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);
q_f_export = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);
p_d_export = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);
p_f_export = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);
n_d_export = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);
n_f_export = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);
n_export = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);
n_f_export_TC = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);
pi_export = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);
    

    
%% Main loop of static problem

for i = 1:s.y_grid_size
   for j = 1:s.z_grid_size  % loop over productivity
      for k = 1:s.x_grid_size  % loop over export status
                    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % 1. Domestic Market (unconstrained)  -Step (i)-               
          % See 'Case 1: Domestic unconstrained firm' in paper
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
          % quantity if unconstrained
          
          q_d(i,j,k) = ( ( (m.sigma - 1)/m.sigma  * r.z_grid(j)...
                           * m.P_d / m.w  ) ^ m.sigma) * m.Q_d;
                    
          % price if unconstrained
          p_d(i,j,k) = q_d(i,j,k)^(-1/m.sigma) ...
                           * (m.Q_d )^(1/m.sigma) * m.P_d;
                    
          % labour if unconstrained
          n_d(i,j,k) = q_d(i,j,k) / r.z_grid(j);
                    
          % profits if unconstrained
          pi_d(i,j,k) = p_d(i,j,k) * q_d(i,j,k) - m.w * n_d(i,j,k);
                                                                         
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Domestic and Foreign Sales               
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
      
             % Domestic quantities and prices
             q_d_export(i,j,k) = ( ( (m.sigma - 1)/m.sigma...
                                  * r.z_grid(j)* m.P_d / m.w  )...
                                  ^ m.sigma) * m.Q_d;
                              
             p_d_export(i,j,k) = q_d_export(i,j,k)^(-1/m.sigma)...
                                  *m.Q_d^(1/m.sigma) * m.P_d;                              
                              
             % Export quantities and prices
                                        
             if s.Dshocks==1 
                    
                    q_f_export(i,j,k) = ( ( (m.sigma - 1)/m.sigma...
                                      * r.z_grid(j)* m.P_f / (m.tau * m.w))...
                                      ^ m.sigma) * m.Q_f*r.y_grid(i) ;

                    p_f_export(i,j,k) = q_f_export(i,j,k)^(-1/m.sigma)...
                                      *(m.Q_f*r.y_grid(i))^(1/m.sigma) * m.P_f;       

      
             elseif s.Fshocks==1

                q_f_export(i,j,k) = ( ( (m.sigma - 1)/m.sigma...
                                  * r.z_grid(j)* m.P_f / (m.tau * m.w))...
                                  ^ m.sigma) * m.Q_f ;
            
                p_f_export(i,j,k) = q_f_export(i,j,k)^(-1/m.sigma)...
                                  *m.Q_f^(1/m.sigma) * m.P_f;              
                              
                              
             elseif s.Tshocks==1
                 
                q_f_export(i,j,k) = ( ( (m.sigma - 1)/m.sigma...
                                  * r.z_grid(j)* m.P_f / (r.y_grid(i) * m.w))...
                                  ^ m.sigma) * m.Q_f ;
            
                p_f_export(i,j,k) = q_f_export(i,j,k)^(-1/m.sigma)...
                                  *m.Q_f^(1/m.sigma) * m.P_f;                    
                 
             end
                              
                
                   
             % Domestic, export and total labor
             
             n_d_export(i,j,k) = q_d_export(i,j,k)/r.z_grid(j);
             
             if s.Tshocks==1
                n_f_export(i,j,k) = r.y_grid(i) * q_f_export(i,j,k)/r.z_grid(j); 
             elseif s.Tshocks==0
                n_f_export(i,j,k) = m.tau * q_f_export(i,j,k)/r.z_grid(j);
             end 
                          
             n_export(i,j,k) = n_d_export(i,j,k) + n_f_export(i,j,k);
             
             
             
             if s.Fshocks==0 && s.FE==0
                
                n_f_export_TC(i,j,k) = n_f_export(i,j,k) + m.TC(k);  
             
                % Profits
                pi_export(i,j,k) = p_d_export(i,j,k) * q_d_export(i,j,k)...
                                   + p_f_export(i,j,k) * q_f_export(i,j,k)...
                                   - m.w*n_export(i,j,k) - m.w * m.TC(k); 
             
             elseif s.Fshocks==1 || s.FE==0
                 
                n_f_export_TC(i,j,k) = n_f_export(i,j,k) + m.TC(k,i);       
        
                % Profits
                pi_export(i,j,k) = p_d_export(i,j,k) * q_d_export(i,j,k)...
                                    + p_f_export(i,j,k) * q_f_export(i,j,k)...
                                   - m.w*n_export(i,j,k) - m.w * m.TC(k,i); 
             
             elseif s.FE==1
                 
                n_f_export_TC(i,j,k) = n_f_export(i,j,k) + m.TC(k,j);       
        
                % Profits
                pi_export(i,j,k) = p_d_export(i,j,k) * q_d_export(i,j,k)...
                                    + p_f_export(i,j,k) * q_f_export(i,j,k)...
                                   - m.w*n_export(i,j,k) - m.w * m.TC(k,j);                 
                 
                               
             end
             
             
                  
             
      end % end of the loop over export status
   end % loop over productivity
end %loop over additional shock

 
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
r.n_f_export_TC = n_f_export_TC;


        
        
end %end function

