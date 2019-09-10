%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the Main file to replicate the results in the paper Kohn,       %
% Leibovici, Szkup (2016).                                                %
% Please cite as:                                                         %
% Kohn, D., Leibovici, F., & Szkup, M. (2016). Financial frictions and    % 
%   new exporter dynamics. International economic review, 57(2), 453-486. %
% © Copyright 2016 by Kohn, D., Leibovici, F. and M. Szkup.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Start

    clear all;
    close all;  
    clc;
    tic;

    %% Choose model

    % Paper
    s.model_baseline_financialfrictions = 1;
    s.model_baseline_sunkcosts = 0;

    % Appendix 
    s.model_extensions_Fshocks = 0;
    s.model_extensions_Dshocks = 0;
    s.model_extensions_Tshocks = 0;
    s.model_extensions_laboradjustmentcosts = 0;
    s.model_extensions_capital = 0;
    s.model_extensions_homogeneousF_ff = 0;
    s.model_extensions_homogeneousF_sc = 0;

    %% Calibrate model

    setup_flags;
    setup_calibration;
    setup_solution;

    %% Solve model

    [r s] = model_solve(m,s);

    %% Simulate model

    [r sim] = model_simulate(m,s,r);

    %% Store results 

    store_results;

    %% End

    disp('Code takes...');
    toc;


