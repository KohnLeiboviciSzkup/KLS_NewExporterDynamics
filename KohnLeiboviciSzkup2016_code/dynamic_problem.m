function r = dynamic_problem(m,s,r)


%% Dynamic problem: value functions and policy functions.    

r.exporter = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
r.a_prime = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);  
r.a_prime_ind = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);  
r.c = zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);  

r.exporter_L = zeros(s.a_grid_size,s.z_grid_size/2,s.x_grid_size);
r.a_prime_L = zeros(s.a_grid_size,s.z_grid_size/2,s.x_grid_size);  
r.a_prime_ind_L = zeros(s.a_grid_size,s.z_grid_size/2,s.x_grid_size);  
r.c_L = zeros(s.a_grid_size,s.z_grid_size/2,s.x_grid_size);      

r.exporter_H = zeros(s.a_grid_size,s.z_grid_size/2,s.x_grid_size);
r.a_prime_H = zeros(s.a_grid_size,s.z_grid_size/2,s.x_grid_size);  
r.a_prime_ind_H = zeros(s.a_grid_size,s.z_grid_size/2,s.x_grid_size);  
r.c_H = zeros(s.a_grid_size,s.z_grid_size/2,s.x_grid_size);      
    

    
%% Solving dynamic problem for m.S=0, k=1

ac_iter=30;  %=0 -> no accelerator
dispn=1;
k=1;
% Initial guess
v_initial = ((r.pi_d(:,:,1).^(1-m.gamma))./(1-m.gamma))./(1-m.beta);

r = dynamic_problem_S0(m,s,r,k,v_initial,ac_iter,dispn);

v_new_L = r.v_new_L;
v_new_H = r.v_new_H;
    
   
    %% if m.S>0
    
if m.S>0  % m.S>0
   
    
   % Solving dynamic problem for m.S=0, k=2 
   
   k=2;
   % initial guess for k=2
   v_initial(:,1:s.z_grid_size/2,1) = v_new_L;
   v_initial(:,s.z_grid_size/2+1:end,1) = v_new_H;
      
   r = dynamic_problem_S0(m,s,r,k,v_initial,ac_iter,dispn); 
   
   v_new_L2 = r.v_new_L;
   v_new_H2 = r.v_new_H;
   
   
   % Loop for m.S>0 problem
   dispn=500;
   
   % Initial guesses:
   v_new=zeros(s.a_grid_size,s.z_grid_size,s.x_grid_size);
   v_new(:,1:s.z_grid_size/2,1)=v_new_L;
   v_new(:,1:s.z_grid_size/2,2)=v_new_L2;
   v_new(:,s.z_grid_size/2+1:end,1)=v_new_H;
   v_new(:,s.z_grid_size/2+1:end,2)=v_new_H2;
   
   v_old   = v_new;
   l=1;
   v_diff = 1;
   
   
   
   % Useful matrices to create u_export and u_noexport outside the loop
   pi_export=r.pi_export(:);
   pi_d=r.pi_d(:);
   az_grid(:,:,1)=r.a_grid'*ones(1,s.z_grid_size);
   az_grid(:,:,2)=az_grid(:,:,1);
   az_grid=az_grid(:);
   
   % if exports:
   c=(pi_export + az_grid)*ones(1,s.a_grid_size)-ones(s.a_grid_size*s.z_grid_size*s.x_grid_size,1)*(r.a_grid/(1+m.r));
   u_export = (c.^(1-m.gamma))./(1-m.gamma); 
   u_export(c<=0) = -1e10; % Only allow for positive consumption
   
   % If it does not export:
   c = (pi_d + az_grid)*ones(1,s.a_grid_size) -  ones(s.a_grid_size*s.z_grid_size*s.x_grid_size,1)*(r.a_grid/(1+m.r));                                         
   u_noexport = (c.^(1-m.gamma))./(1-m.gamma); 
   u_noexport(c<=0) = -1e10;   
   
   clear pi_export pi_d az_grid c
   
   
   while v_diff>s.eps
       
       for k = 1:s.x_grid_size     
           for j = 1:s.z_grid_size
           
               
               rinit=(k-1)*s.z_grid_size*s.a_grid_size+ (j-1)*s.a_grid_size+1;
               rend=(k-1)*s.a_grid_size*s.z_grid_size+j*s.a_grid_size;
               
               % if exports:
               [v_export,index_export] = max(u_export(rinit:rend,:) + ones(s.a_grid_size,1)*(m.beta*r.z_P(j,:)*v_new(:,:,2)'),[],2);                         
               
               % If does not export:
               [v_noexport,index_noexport] = max(u_noexport(rinit:rend,:) + ones(s.a_grid_size,1)*(m.beta*r.z_P(j,:)*v_new(:,:,1)'),[],2);               
                        
               
               % This condition checkes if firm exports and also if the
               % firm can afford entry costs with its assets
               cond = v_export >= v_noexport & (r.a_grid'-(m.alpha_2*m.w*m.TC(k,j))/m.lambda > 0);
            
               % The following condition imposes monotonicity in export
               % decision, on assets (if firm exports with assets i, it
               % should export for all assets j>i. 
               ind = find(cond==1,1,'first');
               if isempty(ind)==0
                cond(ind+1:end) = 1;
               end
               
               v_new(cond,j,k)   = v_export(cond);
               v_new(~cond,j,k)  = v_noexport(~cond);
               r.exporter(cond,j,k)  = 1;
               r.exporter(~cond,j,k) = 0;
               r.a_prime(cond,j,k)  = r.a_grid(index_export(cond));
               r.a_prime(~cond,j,k) = r.a_grid(index_noexport(~cond));
               r.a_prime_ind(cond,j,k)  = index_export(cond);
               r.a_prime_ind(~cond,j,k) = index_noexport(~cond);
               r.c(cond,j,k)  = r.pi_export(cond,j,k)+ r.a_grid(cond)'  - r.a_prime(index_export(cond),j,k)/(1+m.r);
               r.c(~cond,j,k) = r.pi_d(~cond,j,k)    + r.a_grid(~cond)' - r.a_prime(index_noexport(~cond),j,k)/(1+m.r);
               
          end
       end

       v_diff = max(abs(v_new(:)-v_old(:))); 
       v_old = v_new;        
        
       if mod(l,dispn)==0
          disp(['Delta V: ' num2str(v_diff)]);
       end
       
       l=l+1;
        
    end   
   
end  % if m.S>0 
  

    
if m.S==0    
    for xstate = 1:s.x_grid_size
        r.v(:,:,xstate) = [v_new_L(:,:,1) v_new_H(:,:,1)];
        r.exporter(:,:,xstate) = [r.exporter_L(:,:,1) r.exporter_H(:,:,1)];
        r.a_prime(:,:,xstate) = [r.a_prime_L(:,:,1) r.a_prime_H(:,:,1)];  
        r.a_prime_ind(:,:,xstate) = [r.a_prime_ind_L(:,:,1)  r.a_prime_ind_H(:,:,1)];  
        r.c(:,:,xstate) = [r.c_L(:,:,1) r.c_H(:,:,1)] ;  
    end
end
    
    
    
    
    
    
end