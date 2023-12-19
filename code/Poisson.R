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
    all_coef[[i]] = coef(IV_Pois_est(fml, df[sample(n, n/2), ]))
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


### All conflicts:
# Country and year fixed effects
est_noIV_cyfe_oth_conf = fepois(confs_py ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year, df)
est_IV_cyfe_oth_conf  = IV_Pois(confs_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df)
est_IV_noboot_cyfe_oth_conf  = IV_Pois_est(confs_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df)
etable(est_noIV_cyfe_oth_conf, est_IV_cyfe_oth_conf, est_IV_noboot_cyfe_oth_conf)
# Country, year and Country^year fixed effects
est_noIV_cpyfe_oth_conf = fepois(confs_py ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year, df)
est_IV_cpyfe_oth_conf = IV_Pois(confs_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year | dams_pby ~ RGxD_hat, df)
est_IV_noboot_cpyfe_oth_conf = IV_Pois_est(confs_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year | dams_pby ~ RGxD_hat, df)
etable(est_noIV_cpyfe_oth_conf, est_IV_cpyfe_oth_conf, est_IV_noboot_cpyfe_oth_conf)

etable(est_IV_noboot_cyfe_oth_conf, est_IV_noboot_cpyfe_oth_conf, est_IV_Con_cyfe_conf, est_IV_Con_allfe_conf, est_IV_cyfe_oth_conf, est_IV_cpyfe_oth_conf, tex=TRUE)


### Battles:
# Country and year fixed effects
est_noIV_cyfe_oth_btl = fepois(battles_py ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year, df)
est_IV_cyfe_oth_btl  = IV_Pois(battles_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df)
est_IV_noboot_cyfe_oth_btl  = IV_Pois_est(battles_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df)
etable(est_noIV_cyfe_oth_btl, est_IV_cyfe_oth_btl, est_IV_noboot_cyfe_oth_btl)
# Country, year and Country^year fixed effects
est_noIV_cpyfe_oth_btl = fepois(battles_py ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year, df)
est_IV_cpyfe_oth_btl = IV_Pois(battles_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year | dams_pby ~ RGxD_hat, df)
est_IV_noboot_cpyfe_oth_btl = IV_Pois_est(battles_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year | dams_pby ~ RGxD_hat, df)
etable(est_noIV_cpyfe_oth_btl, est_IV_cpyfe_oth_btl, est_IV_noboot_cpyfe_oth_btl)

etable(est_IV_noboot_cyfe_oth_btl, est_IV_noboot_cpyfe_oth_btl, est_IV_Con_cyfe_btl, est_IV_Con_allfe_btl, est_IV_cyfe_oth_btl, est_IV_cpyfe_oth_btl, tex=TRUE)


### All riots:
# Country and year fixed effects
est_noIV_cyfe_oth_riot = fepois(riots_py ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year, df)
est_IV_cyfe_oth_riot  = IV_Pois(riots_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df)
est_IV_noboot_cyfe_oth_riot  = IV_Pois_est(riots_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df)
etable(est_noIV_cyfe_oth_riot, est_IV_cyfe_oth_riot, est_IV_noboot_cyfe_oth_riot)
# Country, year and Country^year fixed effects
est_noIV_cpyfe_oth_riot = fepois(riots_py ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year, df)
est_IV_cpyfe_oth_riot = IV_Pois(riots_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year | dams_pby ~ RGxD_hat, df)
est_IV_noboot_cpyfe_oth_riot = IV_Pois_est(riots_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year | dams_pby ~ RGxD_hat, df)
etable(est_noIV_cpyfe_oth_riot, est_IV_cpyfe_oth_riot, est_IV_noboot_cpyfe_oth_riot)

etable(est_IV_noboot_cyfe_oth_riot, est_IV_noboot_cpyfe_oth_riot, est_IV_Con_cyfe_riot, est_IV_Con_allfe_riot, est_IV_cyfe_oth_riot, est_IV_cpyfe_oth_riot, tex=TRUE)