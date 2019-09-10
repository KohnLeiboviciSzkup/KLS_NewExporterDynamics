function r = dynamic_problem_Sy(m,s,r)

% Dynamic problem: value functions and policy functions.    

    r.exporter = zeros(s.y_grid_size,s.z_grid_size,s.x_grid_size);   
    r.exp = zeros(s.y_grid_size,s.z_grid_size);
    v_new = (r.pi_d)./(1-m.beta);
    v_old = v_new; 
    v_diff = 1;
    
    while v_diff>s.eps


        for k = 1:s.x_grid_size               

            if s.Tshocks==0

                % If firm exports:
                v_export = r.pi_export(:,:,k) + m.beta*r.y_P*v_new(:,:,2)*r.z_P';

                % If firm does not export:                        
                v_noexport = r.pi_d(:,:,k) + m.beta*r.y_P*v_new(:,:,1)*r.z_P';

                % Export decision:
                r.exp(v_export>=v_noexport) = 1;
                r.exp(v_export<v_noexport)  = 0;
                r.exporter(:,:,k) = r.exp;

                % Update:
                v_new(:,:,k) = max(v_export , v_noexport); 

            elseif s.Tshocks==1

                % If firm exports:
                v_export = r.pi_export(:,:,k) + m.beta*r.y_P*v_new(:,:,2)*r.z_P';

                % If firm does not export:        
                v_noexport= r.pi_d(:,:,k) + m.beta*r.y_PNX*v_new(:,:,1)*r.z_P';

                % Export decision:
                r.exp(v_export>=v_noexport) = 1;
                r.exp(v_export<v_noexport)  = 0;
                r.exporter(:,:,k) = r.exp;

                % Update:
                v_new(:,:,k) = max(v_export , v_noexport); 

            end


       end
     
        
        v_diff = max(max(max(abs(v_new-v_old))));
        v_old = v_new; 
        

    end
    r.v = v_new;
    
    
end