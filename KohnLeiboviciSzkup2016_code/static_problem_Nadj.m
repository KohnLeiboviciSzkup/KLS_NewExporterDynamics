function r = static_problem_Nadj(m,s,r)
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
                    
      % Only domestic production (e'=0):
      
      q_d(i,j,k)  = r.z_grid(j)*r.y_grid(i);
      p_d(i,j,k)  = m.P_d*m.Q_d^(1/m.sigma)*q_d(i,j,k)^(-1/m.sigma);
      n_d(i,j,k)  = r.y_grid(i);
      pi_d(i,j,k) = p_d(i,j,k)*q_d(i,j,k) - m.w*r.y_grid(i);
      
      
      % Both domestic and foreign production (e'=1):
      
      ratio = ((m.P_f*m.Q_f^(1/m.sigma))/(m.tau*m.P_d*m.Q_d^(1/m.sigma)))^m.sigma;

      if r.y_grid(i)>m.TC(k,j) 

          q_d_export(i,j,k) = r.z_grid(j)*(r.y_grid(i)-m.TC(k,j)) * (1/ (1+m.tau*ratio));
          q_f_export(i,j,k) = r.z_grid(j)*(r.y_grid(i)-m.TC(k,j)) * (ratio/(1+m.tau*ratio));  
          p_d_export(i,j,k) = m.P_d*m.Q_d^(1/m.sigma) * q_d_export(i,j,k)^(-1/m.sigma);
          p_f_export(i,j,k) = m.P_f*m.Q_f^(1/m.sigma) * q_f_export(i,j,k)^(-1/m.sigma);
          n_d_export(i,j,k) = q_d_export(i,j,k)/r.z_grid(j);
          n_f_export(i,j,k) = m.tau*q_f_export(i,j,k)/r.z_grid(j);
          n_f_export_TC(i,j,k) = m.tau*q_f_export(i,j,k)/r.z_grid(j) + m.TC(k,j);

      end        
        
      pi_export(i,j,k) = p_d_export(i,j,k)*q_d_export(i,j,k) + p_f_export(i,j,k)*q_f_export(i,j,k) - m.w*r.y_grid(i);
      
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

