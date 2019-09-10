function [r s] = model_solve(m,s)

    %% Setup asset and productivity grids

    %Asset grid
    if s.assets==1
        r.a_grid = logspace(log10(s.a_grid_lb),log10(s.a_grid_ub),s.a_grid_size);
    end

    %Productivity shocks    
    [r.log_z_grid,r.z_P,r.z_pi] = approx_shocks(s.z_grid_size,s.z_grid_range_sd,m.log_z_rho,m.log_z_sigma,m.log_z_mu);  % r.z_P is the transition, r.z_pi is the limiting distribution.
    r.z_grid = exp(r.log_z_grid); % recovering levels from lognormal distribution

    r.z_P_original = r.z_P;
    s.z_grid_size_original = s.z_grid_size;

    %% Setup y-grids

    if s.model_baseline_financialfrictions==1 || s.model_baseline_sunkcosts==1 || s.model_extensions_homogeneousF_ff==1 || s.model_extensions_homogeneousF_sc==1 || s.model_extensions_capital==1
        r.y_pi = [1-m.F_share_low;m.F_share_low];
        r.y_P = [1 0;0 1];
        r.y_grid = [1;m.F_value_low];        
        s.y_grid_size = 2;    

        r.y_grid_original = r.y_grid;
        s.y_grid_size_original = s.y_grid_size;        

        %Combine grids into one
        P = zeros(s.y_grid_size*s.z_grid_size);
        for yrow=1:s.y_grid_size 
            for zrow=1:s.z_grid_size
                for ycol=1:s.y_grid_size 
                    for zcol=1:s.z_grid_size
                        irow = s.z_grid_size*(yrow-1)+zrow;
                        jcol = s.z_grid_size*(ycol-1)+zcol;
                        P(irow,jcol) = r.z_P(zrow,zcol)*r.y_P(yrow,ycol);
                    end
                end
            end
        end

        r.z_P = P;

        r.z_grid = r.z_grid(:,ones(1,s.y_grid_size));
        r.z_grid = r.z_grid(:);

        r.y_grid = r.y_grid(:,ones(1,s.z_grid_size))';
        r.y_grid = r.y_grid(:);

        s.z_grid_size = length(r.z_grid); %s.z_grid_size*s.y_grid_size;        
        s.y_grid_size = length(r.y_grid); %s.z_grid_size*s.y_grid_size;     


        m.TC(1,1:s.z_grid_size_original) = (m.F + m.S);
        m.TC(2,1:s.z_grid_size_original) = m.F;

        m.TC(1,s.z_grid_size_original+1:2*s.z_grid_size) = (m.F + m.S)*m.F_value_low;
        m.TC(2,s.z_grid_size_original+1:2*s.z_grid_size) = m.F*m.F_value_low;         


    end

    if s.model_extensions_Fshocks==1
        [r.log_y_grid,r.y_P,r.y_pi] = approx_shocks(s.y_grid_size,s.y_grid_range_sd,m.log_y_rho,m.log_y_sigma,m.log_y_mu);  % r.y_P is the transition, r.y_pi is the limiting distribution.
        r.y_grid = exp(r.log_y_grid); % recovering levels from lognormal distribution 

        r.y_grid_original = r.y_grid;
        s.y_grid_size_original = s.y_grid_size;      

        for i = 1:s.y_grid_size
           m.TC(1,i) = (m.F + m.S)*r.y_grid(i);
           m.TC(2,i) = m.F*r.y_grid(i);
        end    
    end

    if s.model_extensions_Dshocks==1
        [r.log_y_grid,r.y_P,r.y_pi] = approx_shocks(s.y_grid_size,s.y_grid_range_sd,m.log_y_rho,m.log_y_sigma,m.log_y_mu);  % r.y_P is the transition, r.y_pi is the limiting distribution.
        r.y_grid = exp(r.log_y_grid); % recovering levels from lognormal distribution    

        r.y_grid_original = r.y_grid;
        s.y_grid_size_original = s.y_grid_size;        

        r.FE_pi = [1-m.F_share_low;m.F_share_low];
        r.FE_P = [1 0;0 1];
        r.FE_grid = [1;m.F_value_low];         

        %Combine grids into one
        P = zeros(s.z_grid_size);
        for FErow=1:2 
            for zrow=1:s.z_grid_size
                for FEcol=1:2
                    for zcol=1:s.z_grid_size
                        irow = s.z_grid_size*(FErow-1)+zrow;
                        jcol = s.z_grid_size*(FEcol-1)+zcol;
                        P(irow,jcol) = r.z_P(zrow,zcol)*r.FE_P(FErow,FEcol);
                    end
                end
            end
        end

        r.z_P = P;

        r.z_grid = r.z_grid(:,ones(1,2));
        r.z_grid = r.z_grid(:);

        r.z_pi = [(1-m.F_share_low)*r.z_pi;m.F_share_low*r.z_pi];

        r.FE_grid = r.FE_grid(:,ones(1,s.z_grid_size))';
        r.FE_grid = r.FE_grid(:);


        s.z_grid_size = s.z_grid_size*2;    

        m.TC(1,1:s.z_grid_size_original) = (m.F + m.S);
        m.TC(2,1:s.z_grid_size_original) = m.F;

        m.TC(1,s.z_grid_size_original+1:2*s.z_grid_size) = (m.F + m.S)*m.F_value_low;
        m.TC(2,s.z_grid_size_original+1:2*s.z_grid_size) = m.F*m.F_value_low;            
    end            

    if s.model_extensions_Tshocks==1
        r.y_grid = linspace(m.tauH,m.tauL,s.y_grid_size);

        r.y_P = zeros(s.y_grid_size,s.y_grid_size);
        for n=1:s.y_grid_size-1
            r.y_P(n,n) = 1-m.Trho;
            r.y_P(n,n+1) = m.Trho;
        end
        r.y_P(s.y_grid_size,s.y_grid_size) = 1;

        r.y_PNX = zeros(s.y_grid_size,s.y_grid_size);
        for n=1:s.y_grid_size
            r.y_PNX(n,1) = 1;
        end    

        r.FE_pi = [1-m.F_share_low;m.F_share_low];
        r.FE_P = [1 0;0 1];
        r.FE_grid = [1;m.F_value_low];         

        %Combine grids into one
        P = zeros(s.z_grid_size);
        for FErow=1:2 
            for zrow=1:s.z_grid_size
                for FEcol=1:2
                    for zcol=1:s.z_grid_size
                        irow = s.z_grid_size*(FErow-1)+zrow;
                        jcol = s.z_grid_size*(FEcol-1)+zcol;
                        P(irow,jcol) = r.z_P(zrow,zcol)*r.FE_P(FErow,FEcol);
                    end
                end
            end
        end

        r.z_P = P;

        r.z_grid = r.z_grid(:,ones(1,2));
        r.z_grid = r.z_grid(:);

        r.z_pi = [(1-m.F_share_low)*r.z_pi;m.F_share_low*r.z_pi];

        r.FE_grid = r.FE_grid(:,ones(1,s.z_grid_size))';
        r.FE_grid = r.FE_grid(:);


        s.z_grid_size = s.z_grid_size*2;      

       m.TC(1,1:s.z_grid_size_original) = (m.F + m.S);
       m.TC(2,1:s.z_grid_size_original) = m.F;

       m.TC(1,s.z_grid_size_original+1:2*s.z_grid_size) = (m.F + m.S)*m.F_value_low;
       m.TC(2,s.z_grid_size_original+1:2*s.z_grid_size) = m.F*m.F_value_low;                
    end

    if s.model_extensions_laboradjustmentcosts==1
        s.y_grid_lb = 0.01; 

        n_d_unc = ((m.sigma/(m.sigma-1))*m.w/(m.P_d*m.Q_d^(1/m.sigma)))^(-m.sigma) * max(r.z_grid)^(m.sigma-1);
        n_f_unc = ((m.sigma/(m.sigma-1))*m.w/(m.P_f*m.Q_f^(1/m.sigma)))^(-m.sigma) * m.tau^(1-m.sigma) * max(r.z_grid)^(m.sigma-1);
        s.y_grid_ub = 2*(n_d_unc+n_f_unc);                

        r.y_grid = linspace(0,1,s.y_grid_size)';
        r.y_grid = r.y_grid.^2;
        r.y_grid = r.y_grid.*(s.y_grid_ub-s.y_grid_lb) + s.y_grid_lb;     

        r.FE_pi = [1-m.F_share_low;m.F_share_low];
        r.FE_P = [1 0;0 1];
        r.FE_grid = [1;m.F_value_low];         

        %Combine grids into one
        P = zeros(s.z_grid_size);
        for FErow=1:2 
            for zrow=1:s.z_grid_size
                for FEcol=1:2
                    for zcol=1:s.z_grid_size
                        irow = s.z_grid_size*(FErow-1)+zrow;
                        jcol = s.z_grid_size*(FEcol-1)+zcol;
                        P(irow,jcol) = r.z_P(zrow,zcol)*r.FE_P(FErow,FEcol);
                    end
                end
            end
        end

        r.z_P = P;

        r.z_grid = r.z_grid(:,ones(1,2));
        r.z_grid = r.z_grid(:);

        r.z_pi = [(1-m.F_share_low)*r.z_pi;m.F_share_low*r.z_pi];

        r.FE_grid = r.FE_grid(:,ones(1,s.z_grid_size))';
        r.FE_grid = r.FE_grid(:);

        s.z_grid_size = s.z_grid_size*2;    

       m.TC(1,1:s.z_grid_size_original) = (m.F + m.S);
       m.TC(2,1:s.z_grid_size_original) = m.F;

       m.TC(1,s.z_grid_size_original+1:2*s.z_grid_size) = (m.F + m.S)*m.F_value_low;
       m.TC(2,s.z_grid_size_original+1:2*s.z_grid_size) = m.F*m.F_value_low;                
    end
    
    
    

    
    
    

    %% Solve static and dynamic problems

    if s.model_baseline_financialfrictions==1 || s.model_baseline_sunkcosts==1 || s.model_extensions_homogeneousF_ff==1 || s.model_extensions_homogeneousF_sc==1          
        r = static_problem(m,s,r);   
        r = dynamic_problem(m,s,r);   
    end
    
    if s.model_extensions_Dshocks==1 || s.model_extensions_Tshocks==1 || s.model_extensions_Fshocks==1
        r = static_problem_Sy(m,s,r);     
        r = dynamic_problem_Sy(m,s,r);
    end

    if s.model_extensions_laboradjustmentcosts==1
        r = static_problem_Nadj(m,s,r); 
        r = dynamic_problem_Nadj(m,s,r);
    end

    if s.model_extensions_capital==1
        r = static_problem_K(m,s,r);   
        r = dynamic_problem_K(m,s,r);      
    end
    
 
    
    if s.flag_unc==1 %only if s.model_baseline_financialfrictions==1 
        
       m_unc = m;
       r_unc = r;

       m_unc.lambda = 1e8;
       r_unc = static_problem(m_unc,s,r_unc);         

       r.pi_d_unc = r_unc.pi_d;
       r.pi_export_unc = r_unc.pi_export;

       r.exporter_ff_const = r.pi_export>=r.pi_d;
       r.pi_const = r.exporter_ff_const.*r.pi_export + (1-r.exporter_ff_const).*r.pi_d;

       r.exporter_ff_unc = r.pi_export_unc>=r.pi_d_unc;
       r.pi_unc = r.exporter_ff_unc.*r.pi_export_unc + (1-r.exporter_ff_unc).*r.pi_d_unc;   
       
    end

end
