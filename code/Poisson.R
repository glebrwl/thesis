library(fixest)
# Get the estimates:
IV_Pois_est = function(fml, df) {
  # Run 2SLS:
  fe_ols = feols(fml, df)
  # Get the residuals from the 1st stage:
  resid_1st_stage = resid(summary(fe_ols, stage=1))
  # Change 2nd stage OLS with Poisson, add residuals as regressors:
  y = model.matrix(fe_ols, type = 'lhs')
  RHS = model.matrix(fe_ols, type=c('iv.endo', 'iv.exo'))
  RHS_control_function = cbind(RHS, resid_1st_stage)
  # Introduce FEs for 2nd stage:
  FE = if(!is.null(fe_ols$fixef_vars)) model.matrix(fe_ols, type = "fixef") else NULL
  # Run 2nd stage
  feglm.fit(y, RHS_control_function, FE, family='poisson')
}

# Get the VCOV matrix:
IV_Pois_VarCov = function(fml, df, rep) {
  n = nrow(df)
  # Prepare 
  all_coef = vector('list', rep)
  # Sample half of data rep times
  for (i in 1:rep) {
    all_coef[[i]] = coef(IVPois_est(fml, df[sample(n, n/2), ]))
  }
  all_coef_mrx = do.call('rbind', all_coef)
  var(all_coef_mrx)
}

# Function to get estimates and SEs:
IV_Pois = function(fml, df, rep=200) {
  res = IV_Pois_est(fml, df)
  vcov_mrx = IV_Pois_VarCov(fml, df, rep)
  summary(res, .vcov=vcov_mrx)
}

z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn + z_zgjij

# Basin and year fixed effects
est_noIV_byfe = fepois(btl_p_y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | HYBAS_ID + year, df1)
est_IV_byfe = IV_Pois(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | HYBAS_ID + year | dams_pby ~ RGxD_hat, df1)
est_IV_noboot_byfe = IV_Pois_est(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | HYBAS_ID + year | dams_pby ~ RGxD_hat, df1)

# Country and year fixed effects
est_noIV_cyfe = fepois(btl_p_y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | Country + year, df1)
est_IV_cyfe  = IV_Pois(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | Country + year | dams_pby ~ RGxD_hat, df1)
est_IV_noboot_cyfe  = IV_Pois_est(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | Country + year | dams_pby ~ RGxD_hat, df1)

# Country^year fixed effects
est_noIV_cpyfe = fepois(btl_p_y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | Country^year, df1)
est_IV_cpyfe = IV_Pois(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | Country^year | dams_pby ~ RGxD_hat, df1)
est_IV_noboot_cpyfe = IV_Pois_est(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | Country^year | dams_pby ~ RGxD_hat, df1)

# Basins + Country^year fixed effects
est_noIV_allfe = fepois(btl_p_y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | HYBAS_ID + Country^year, df1)
est_IV_allfe = IV_Pois(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df1)
est_IV_noboot_allfe = IV_Pois_est(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df1)


# Basin and year fixed effects
etable(est_noIV_byfe, est_IV_byfe, est_IV_noboot_byfe)
# Country and year fixed effects
etable(est_noIV_cyfe, est_IV_cyfe, est_IV_noboot_cyfe)
# Country^year fixed effects
etable(est_noIV_cpyfe, est_IV_cpyfe, est_IV_noboot_cpyfe)
# Basins + Country^year fixed effects
etable(est_noIV_allfe, est_IV_allfe, est_IV_noboot_allfe)

etable(est_noIV_byfe, est_noIV_cyfe, est_noIV_cpyfe, est_noIV_allfe)
etable(est_IV_noboot_byfe, est_IV_noboot_cyfe, est_IV_noboot_cpyfe, est_IV_noboot_allfe)
etable(est_IV_byfe, est_IV_cyfe, est_IV_cpyfe, est_IV_allfe)






# Basin and year fixed effects
est_noIV_byfe_oth = fepois(btl_p_y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + year, df1)
est_IV_byfe_oth = IV_Pois(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + year | dams_pby ~ RGxD_hat, df1)
est_IV_noboot_byfe_oth = IV_Pois_est(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + year | dams_pby ~ RGxD_hat, df1)

# Country and year fixed effects
est_noIV_cyfe_oth = fepois(btl_p_y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year, df1)
est_IV_cyfe_oth  = IV_Pois(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df1)
est_IV_noboot_cyfe_oth  = IV_Pois_est(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df1)

# Country^year fixed effects
est_noIV_cpyfe_oth = fepois(btl_p_y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country^year, df1)
est_IV_cpyfe_oth = IV_Pois(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country^year | dams_pby ~ RGxD_hat, df1)
est_IV_noboot_cpyfe_oth = IV_Pois_est(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country^year | dams_pby ~ RGxD_hat, df1)

# Basins + Country^year fixed effects
est_noIV_allfe_oth = fepois(btl_p_y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + Country^year, df1)
est_IV_allfe_oth = IV_Pois(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df1)
est_IV_noboot_allfe_oth = IV_Pois_est(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df1)

etable(est_IV_byfe_oth, est_IV_cyfe_oth, est_IV_cpyfe_oth, est_IV_allfe_oth, tex=TRUE)
