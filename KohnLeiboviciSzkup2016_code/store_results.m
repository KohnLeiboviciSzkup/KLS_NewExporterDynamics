
% Calibration moments
results_moments_cal.starters = sim.share_starters;
results_moments_cal.stoppers = sim.share_stoppers;
results_moments_cal.x_size_premium = sim.export_size_premium_med;
results_moments_cal.x_sales_ratio = sim.export_sales_ratio_med;
if s.assets==1
    results_moments_cal.extfinancepremium = sim.ext_finance_premium_mean;
end

% Sales distribution
results_distribution_cal = sim.Ssize_distribution';

% dynamics of new exporters:
results_dynamics.newx_sales_d_growth_med = sim.xperiods_sales_d_growth_med;
results_dynamics.newx_sales_f_growth_med = sim.xperiods_sales_f_growth_med;
results_dynamics.newx_export_intensity_med = sim.xperiods_export_intensity_med;
results_dynamics.newx_hazard = sim.hazard_pool';
if s.assets==1
    results_dynamics.newx_ef_log_growth = sim.xperiods_ef_loggrowth_mean;
end
results_dynamics.nonx_hazard = sim.nonxhazard_pool';

% Labor size premiums:
results_moments.newnonx_exporters_sprem=sim.new_nonexport_size_premium_med;
results_moments.newnonx_nonexporters_sprem=sim.new_export_size_premium_med/sim.new_stoppers_size_premium_med ;
results_moments.newnonx_newexporters_sprem=1/sim.new_stoppers_size_premium_med;


% Extensions
if s.model_extensions_Dshocks==1
   results_moments_cal.salesFD_SDratio = sim.salesFD_SDratio;
   tmoments.salesFD_SDratio = 1.406;
end

if s.model_extensions_laboradjustmentcosts ==1
   results_moments_cal.wagebill_acf_med = sim.wagebill_acf_med;
   tmoments.wagebill_acf_med = 0.62;
end




if s.model_baseline_financialfrictions == 1 || s.model_baseline_sunkcosts == 1;

    % Online appendix
    % Bernard and Jensen Stats
    appendix_BJ(1,1) = sim.BJstart_sales_gr1;
    appendix_BJ(1,2) = sim.BJstop_sales_gr1;
    appendix_BJ(1,3) = sim.BJboth_sales_gr1;
    appendix_BJ(1,4) = sim.BJneither_sales_gr1;
    
    appendix_BJ(2,1) = sim.BJstart_sales_gr4;
    appendix_BJ(2,2) = sim.BJstop_sales_gr4;
    appendix_BJ(2,3) = sim.BJboth_sales_gr4;
    appendix_BJ(2,4) = sim.BJneither_sales_gr4;
    
    appendix_BJ(3,1) = sim.BJstart_sales_gr8;
    appendix_BJ(3,2) = sim.BJstop_sales_gr8;
    appendix_BJ(3,3) = sim.BJboth_sales_gr8;
    appendix_BJ(3,4) = sim.BJneither_sales_gr8;

end    
    
if s.model_baseline_financialfrictions == 1
    
    % Working capital asymmetry
    appendix_workcapital_int=sim.work_capital_int_mean;

    % dynamics of financial constrained and unconstrained firms:
    appendix_dynamics_con.hazard=sim.hazard_pool_const;
    appendix_dynamics_con.export_intensity=sim.xperiods_export_intensity_med_const;
    appendix_dynamics_con.sales_f_growth=sim.xperiods_sales_f_growth_med_const;
    appendix_dynamics_con.sales_d_growth=sim.xperiods_sales_d_growth_med_const;
    appendix_dynamics_unc.hazard =sim.hazard_pool_unc;
    appendix_dynamics_unc.export_intensity=sim.xperiods_export_intensity_med_unc;
    appendix_dynamics_unc.sales_f_growth =sim.xperiods_sales_f_growth_med_unc;
    appendix_dynamics_unc.sales_d_growth =sim.xperiods_sales_d_growth_med_unc;

end


