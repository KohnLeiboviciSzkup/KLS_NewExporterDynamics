function r = dynamic_problem_Nadj(m,s,r)

% Dynamic problem: value functions and policy functions.    

    r.exporter = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);
    r.y_prime = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);  
    r.y_prime_ind = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);   
   
   
    % Initial guess
   v_initial = ((r.pi_d(:,:,1).^(1-m.gamma))./(1-m.gamma))./(1-m.beta);
    
  
   v_new(:,:,1)=v_initial;
   v_new(:,:,2)=v_initial;
   
   v_old  = v_new;
   v_diff = 1;    
    
    
    yadj = zeros(s.y_grid_size,s.y_grid_size);


    for i = 1:s.y_grid_size

        yadj(:,i) = m.xi.*((r.y_grid-r.y_grid(i))./r.y_grid(i)).^2.*r.y_grid(i);

    end


    pi_export=r.pi_export(:);
    pi_d=r.pi_d(:);
    u_export = zeros(length(pi_export),s.y_grid_size);
    u_noexport = zeros(length(pi_export),s.y_grid_size);

    for k = 1:s.x_grid_size
       for j = 1:s.z_grid_size 

           rinit=(k-1)*s.z_grid_size*s.y_grid_size+ (j-1)*s.y_grid_size+1;
           rend=(k-1)*s.y_grid_size*s.z_grid_size+j*s.y_grid_size;         

           u_export(rinit:rend,:) =  ( ((r.y_grid> m.TC(k)).*pi_export(rinit:rend) + (r.y_grid<=m.TC(k)).* -1e8 )*ones(1,s.y_grid_size) -  yadj )';
           u_noexport(rinit:rend,:) = (pi_d(rinit:rend)*ones(1,s.y_grid_size)-  yadj)';

       end
    end


    iter=1;
    
    while v_diff>s.eps
        
        for k = 1:s.x_grid_size
            for j = 1:s.z_grid_size

                   rinit=(k-1)*s.z_grid_size*s.y_grid_size+ (j-1)*s.y_grid_size+1;
                   rend=(k-1)*s.z_grid_size*s.y_grid_size+j*s.y_grid_size;

                   % if exports:
                   [v_export,index_export] = max(u_export(rinit:rend,:) + ones(s.y_grid_size,1)*(m.beta*r.z_P(j,:)*v_new(:,:,2)'),[],2);                         

                   % If does not export:
                   %if k==1 || ( r.exporter(1,j,1) == 0 && r.exporter(end,j,1) == 0) 
                   [v_noexport,index_noexport] = max(u_noexport(rinit:rend,:) + ones(s.y_grid_size,1)*(m.beta*r.z_P(j,:)*v_new(:,:,1)'),[],2); 

                   cond = v_export >= v_noexport; % & (r.y_grid> m.TC(k) );
                        
                   v_new(cond,j,k)   = v_export(cond);
                   v_new(~cond,j,k)  = v_noexport(~cond);
                   r.exporter(cond,j,k)  = 1;
                   r.exporter(~cond,j,k) = 0;
                   r.y_prime(cond,j,k)  = r.y_grid(index_export(cond));
                   r.y_prime(~cond,j,k) = r.y_grid(index_noexport(~cond));
                   r.y_prime_ind(cond,j,k)  = index_export(cond);
                   r.y_prime_ind(~cond,j,k) = index_noexport(~cond);
                  

            end
        end
        
        
        
        
        v_diff = max(abs(v_new(:)-v_old(:)));
        v_old = v_new;        
        
        if mod(iter,100)==0
            disp(['Delta V: ' num2str(v_diff)]);
        end
        iter=iter+1;
    end
    r.v = v_new;
    
    
end