% This code reproduces all the numbers needed for Table 7 and Table 8.
%
% Functions:
%  (1) fun_trade_lib(tau,flag_ff,flag_sc)
%      This function computes GDP, sales_d, sales_f, share of exporters and
%      exports per firm for given tau. 
%      
%      -> flag_ff = 1 and flag_sc = 0 then results correspond to financial friction model
%      -> flag_ff = 0 and flag_sc = 1 then results correspond to sunk cost model
%
%  (2) tau_eq(tau_ff,delta_tau,GDP_growth_sc,GDP_ff)
%      This function computes the tariff equivalent of financial frictions
%      when in the sunk cost model variable trade cost, tau, decreased by
%      delta_tau percent.

clc
clear all
close all

tau_ff = 1.4547; % benchmark tau
tau_sc = 1.5113; % benchmark tau

% size of liberalization:
prc = [0 0.05,0.1];
x_ff = zeros(numel(prc),6);
x_sc = zeros(numel(prc),6);
lambda = 1.64;


%% TRADE LIBERALIZATION:

for i = 1:numel(prc)
    
    % Trade Liberalization (Financial Frictions Model)    
    tau_ff_lib = tau_ff - prc(i);
    flag_ff = 1;
    flag_sc = 0;
    x_ff(i,:) = trade_lib_fun(tau_ff_lib,flag_ff,flag_sc,lambda);

    % Trade Liberalization (Sunk Cost Model)    
    tau_sc_lib = tau_sc - prc(i);
    flag_ff = 0;
    flag_sc = 1;
    x_sc(i,:) = trade_lib_fun(tau_sc_lib,flag_ff,flag_sc,lambda);

end


%% TARIFF EQUIVALENT:

GDP_growth_sc = zeros(numel(prc)-1,1); % GDP growth in sunk cost model following trade liberalization
GDP_ff = x_ff(1,1);                    % GDP in the financial frictions model before liberalization
tau_eq = zeros(numel(prc)-1,1);        % tariff equivalent vector


for i = 2:numel(prc)
    
    GDP_growth_sc(i) = (x_sc(i,1)/x_sc(1,1)-1)*100;
    tau_eq(i) = trade_lib_tariff_equivalent(tau_ff,prc(i),GDP_growth_sc(i),GDP_ff,lambda);
    
end
    
    
%% Lambda 1.77

lambda = 1.77;

for i = 1:numel(prc)
    
    % Trade Liberalization (Financial Frictions Model)    
    tau_ff_lib = tau_ff - prc(i);
    flag_ff = 1;
    flag_sc = 0;
    x_ff(i,:) = trade_lib_fun(tau_ff_lib,flag_ff,flag_sc,lambda);


end


    

