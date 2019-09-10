March, 2016

This code replicates the results in the paper Kohn, Leibovici, and Szkup (2016)                                                %

Please cite as:                                                         
  Kohn, D., Leibovici, F., and Szkup, M. (2016). Financial frictions and     
   new exporter dynamics. International economic review, 57(2), 453-486. 
 © Copyright 2016 by Kohn, D., Leibovici, F. and M. Szkup.               

In this document, we provide instructions to use the code, although the code itself is documented and should be straightforward to use.
Email us if you have any questions or comments: davidkohn@uc.cl, fleibovici@gmail.com, michal.szkup@ubc.ca

The code is composed of the following files:
(1) main: This is the file you need to run. Select the model (i.e. baseline_financialfrictions, baseline_sunkcosts, etc) and run.
(2) model_solve: Calls discretization function, static and dynamic problems given model selection
(3) model_simulate: Calls simulate function given model selection
(4) approx_shocks: discretizes AR(1) process for productivity shocks
(5) dynamic_problem: Solves baseline dynamic problem
(6) dynamic_problem_S0: Solves dynamic problem for S=0 for a given export status and initial value function. Used in dynamic_problem.m
(7) static_problem: Solves baseline static problem
(8) q: Saves for optimal quantity in the problem of the firm producing for both markets
(9) model_simulate_baseline: Simulates baseline model and computes statistics
(10) setup_flags: updates flags given model selection 
(11) setup_calibration: updates calibration given model selection
(12) setup_solution: solution parameters
(13) store_results: saves results 

Additional files required to reproduce results in Section 6 of the paper (Policy Analysis):
(14) trade_lib_main: Reproduces all the numbers needed for Table 7 and Table 8 in the paper. Uses trade_lib_fun.m and trade_lib_tariff_equivalent.m and 
(15) trade_lib_fun: Computes results of trade liberalization
(16) trade_lib_tariff_equivalent: Computes results of the tariff equivalent to financial frictions

Additional files for the extensions with shocks to export costs, demand shocks and iceberg costs shocks:
(17) dynamic_problem_Sy: Solves dynamic problem
(18) static_problem_Sy: Solves static problem
(19) model_simulate_Sy: Simulates model and computes statistics

Additional files for the extension with labor adjustment costs:
(20) dynamic_problem_Nadj: Solves dynamic problem
(21) static_problem_Nadj: Solves static problem
(22) model_simulate_Nadj: Simulates model and computes statistics

Additional files for the extension where assets are used as an input of production (i.e. physical capital):
(23) dynamic_problem_K: Solves dynamic problem
(24) static_problem_K: Solves static problem
(25) q_K: Saves for optimal quantity in the problem of the firm producing for both markets
(26) model_simulate_K: Simulates model and computes statistics

Instructions:
1- Unzip
2- Open main.m
3- Select model you want to run
4- Run


 



