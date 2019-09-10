
%% Solution parameters

    %Productivity     
    s.z_grid_size = 100; % number of grid points for productivity
    s.z_grid_range_sd = 3; % range of values for productivity states (in number of stdev)

    %Assets    
    if s.assets==1
        s.a_grid_size = 200; % number of grid points for assets
        s.a_grid_lb = 0.01; %%0.01; % lowest value for assets  
        s.a_grid_ub = max(m.Q_d,m.Q_f); %Range between lowest and highest asset
    end

    %Simulations    
    s.T = 13;   % Number of periods
    s.N =50000;  % Number of firms         
    s.CN = 4;
    s.CD = 1;     % Only works for 0 and 1 for now (xhistory negative)
    s.burn = 50; % Number of time periods burnt = s.burn*s.T
    s.seed = 88;
%    RandStream.setDefaultStream(RandStream('mt19937ar','seed',s.seed)); % Sets seed for randomizations       
     RandStream.setGlobalStream(RandStream('mt19937ar','seed',s.seed)); % Sets seed for randomizations       
    
    %Convergence parameters     
    s.eps = 1e-8;

