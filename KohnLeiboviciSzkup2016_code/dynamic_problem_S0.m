function r = dynamic_problem_S0(m,s,r,k,v_initial,ac_iter,dispn)

% Comments (24.02.2014):

% This function solves problem for S=0 for a given export status and initial value function and
% returns policies and value function

   z_P = r.z_P(1:s.z_grid_size/2,1:s.z_grid_size/2);
   exporter_ff = r.pi_export>=r.pi_d;
    
    
        
   l=1;
   v_new_L = v_initial(:,1:s.z_grid_size/2,1);
   v_old = v_new_L;
   v_diff = 1;
   u_export = zeros(s.a_grid_size,s.a_grid_size);
      
   while v_diff>s.eps      
            
         for j = 1:s.z_grid_size/2
                        
             c = (exporter_ff(:,j,k).*r.pi_export(:,j,k) + (1-exporter_ff(:,j,k)).*r.pi_d(:,j,k) + r.a_grid')*ones(1,s.a_grid_size) - ones(s.a_grid_size,1)*(r.a_grid/(1+m.r));  
             u_export(c>0) = ((c(c>0)).^(1-m.gamma))./(1-m.gamma); 
             u_export(c<=0) = -1e10; % Only allow for positive consumption
             [v_export,index_export] = max(u_export + ones(s.a_grid_size,1)*(m.beta*z_P(j,:)*v_new_L(:,:,1)'),[],2); 
                    
             v_new_L(:,j,1) = v_export;                 
             r.exporter_L(:,j,k) = exporter_ff(:,j,k);   
             r.a_prime_L(:,j,k)  = r.a_grid(index_export);
             r.a_prime_ind_L(:,j,k)=index_export;
             r.c_L(:,j,k) = exporter_ff(:,j,k).*r.pi_export(:,j,k) + (1-exporter_ff(:,j,k)).*r.pi_d(:,j,k) + r.a_grid' - r.a_prime_L(:,j,k)/(1+m.r);% recovering policy function for consumption                                                                   
                    
         end
            



         aind=r.a_prime_ind_L(:,:,k);
         c=r.pi_export(:,1:s.z_grid_size/2,k).*r.exporter_L(:,:,k) + (1-r.exporter_L(:,:,k)).*r.pi_d(:,1:s.z_grid_size/2,k) + r.a_grid'*ones(1,s.z_grid_size/2) - r.a_grid(aind)./(1+m.r);
         u = ((c).^(1-m.gamma))./(1-m.gamma); 

         v_n= v_new_L;
         for iter=1:ac_iter;
             for j = 1:s.z_grid_size/2
                 v_new_L(:,j,1) = u(:,j) + (m.beta*z_P(j,:)*v_n(aind(:,j),:,1)')';
             end
             v_n= v_new_L;
         end    

        
         v_diff = max(max(max(abs(v_new_L-v_old))));
         v_old = v_new_L;         %s.adj*v_new + (1-s.adj)*v_old; 
            
         if mod(l,dispn)==0
            disp(['Delta V: ' num2str(v_diff)]);
         end    
         l=l+1;            
   end %while v_diff>s.eps 
   
   
   
   % Second half of the problem
   
   v_diff  = 1;
   
   if k==1
      
       v_new_H =v_new_L;
   
   else
       
       v_new_H =v_initial(:,s.z_grid_size/2+1:end,1);
   
   end
   
   v_old   = v_new_H;
   l=1;
    
   while v_diff>s.eps
        
         for j = s.z_grid_size/2+1:s.z_grid_size
                   
             c = (exporter_ff(:,j,k).*r.pi_export(:,j,k) + (1-exporter_ff(:,j,k)).*r.pi_d(:,j,k) + r.a_grid')*ones(1,s.a_grid_size) - ones(s.a_grid_size,1)*(r.a_grid/(1+m.r));  
             u_export(c>0) = ((c(c>0)).^(1-m.gamma))./(1-m.gamma); 
             u_export(c<=0) = -1e10; % Only allow for positive consumption
             [v_export,index_export] = max(u_export + ones(s.a_grid_size,1)*(m.beta*z_P(j-s.z_grid_size/2,:)*v_new_H(:,:,1)'),[],2); 
                    
              v_new_H(:,j-(s.z_grid_size/2),1) = v_export;                 
              r.exporter_H(:,j-(s.z_grid_size/2),k) = exporter_ff(:,j,k);   
              r.a_prime_H(:,j-(s.z_grid_size/2),k)  = r.a_grid(index_export);
              r.a_prime_ind_H(:,j-(s.z_grid_size/2),k)=index_export;
              r.c_H(:,j-(s.z_grid_size/2),k) = exporter_ff(:,j,k).*r.pi_export(:,j,k) + (1-exporter_ff(:,j,k)).*r.pi_d(:,j,k) + r.a_grid' - r.a_prime_H(:,j-(s.z_grid_size/2),k)/(1+m.r);% recovering policy function for consumption                                                                   
                    
         end      
        
 
         aind=r.a_prime_ind_H(:,:,k);
         c=r.pi_export(:,(s.z_grid_size/2)+1:end,k).*r.exporter_H(:,:,k) + (1-r.exporter_H(:,:,k)).*r.pi_d(:,(s.z_grid_size/2)+1:end,k) + r.a_grid'*ones(1,s.z_grid_size/2) - r.a_grid(aind)./(1+m.r);
         u = ((c).^(1-m.gamma))./(1-m.gamma); 

         v_n= v_new_H;
         for iter=1:ac_iter;
             for j = 1:s.z_grid_size/2
                 v_new_H(:,j,1) = u(:,j) + (m.beta*z_P(j,:)*v_n(aind(:,j),:,1)')';
             end
             v_n= v_new_H;
         end        
        


        
         v_diff = max(max(max(abs(v_new_H-v_old))));
         v_old = v_new_H;         %s.adj*v_new + (1-s.adj)*v_old; 
        
         if mod(l,dispn)==0
            disp(['Delta V: ' num2str(v_diff)]);
         end    
         l=l+1;
        
   end
   
   r.v_new_L=v_new_L;
   r.v_new_H=v_new_H;
   
   
end