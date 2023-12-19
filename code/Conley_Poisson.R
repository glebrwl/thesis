library(fixest)
library(stringr)

# Get the estimates:
IV_Pois_Conley_est = function(fml, df) {
  # Run 2SLS:
  fe_ols = feols(fml, df)
  # Extract the residuals from the 1st stage:
  resid_1st_stage = resid(summary(fe_ols, stage=1))
  # Change 2nd stage OLS with Poisson, add residuals as regressors:
  y = model.matrix(fe_ols, type = 'lhs')
  RHS = model.matrix(fe_ols, type=c('iv.endo', 'iv.exo'))
  RHS_control_function = cbind(RHS, resid_1st_stage)
  # Introduce FEs for 2nd stage:
  FE = if(!is.null(fe_ols$fixef_vars)) model.matrix(fe_ols, type = "fixef") else NULL
  # Add coordinates:
  lat = df$lat
  lon = df$lon
  # Produce final data:
  new_data = cbind(y, RHS_control_function, FE, lat, lon)
  # Parse formula for feglm:
  fml_str = deparse(fml, width.cutoff = 500)
  new_str = strsplit(fml_str, '\\|')[[1]][2]
  new_fml = as.formula(paste('y ~', paste(colnames(RHS_control_function), collapse = " + "), '|', new_str))
  # Run 2nd stage:
  feglm(new_fml, new_data, family='poisson')
}

# Function to get estimates and SEs:
IV_Pois_Conley = function(fml, df, cutoff=100) {
  res = IV_Pois_Conley_est(fml, df)
  vcov_mrx = vcov_conley(IV_Pois_Conley_est(fml, df), cutoff = cutoff, lat = "lat", lon = "lon")
  summary(res, .vcov=vcov_mrx)
}

### All conflicts:
# Country + year Fixed Effects:
est_IV_Con_cyfe_conf = IV_Pois_Conley(confs_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df)
# Country + year + Country^year Fixed Effects:
fe_ols = feols(confs_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year | dams_pby ~ RGxD_hat, df)
resid_1st_stage = resid(summary(fe_ols, stage=1))
y = model.matrix(fe_ols, type = 'lhs')
RHS = model.matrix(fe_ols, type=c('iv.endo', 'iv.exo'))
RHS_control_function = cbind(RHS, resid_1st_stage)
FE = data.frame(Country = df$Country, year = df$year)
lat = df$lat
lon = df$lon
new_data = cbind(y, RHS_control_function, FE, lat, lon)
res = feglm(y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn + resid_1st_stage| Country + year + Country^year , new_data, family='poisson')
vcov_mrx = vcov_conley(res, cutoff = cutoff, lat = "lat", lon = "lon")
est_IV_Con_allfe_conf = summary(res, .vcov=vcov_mrx)
etable(est_IV_Con_cyfe_conf, est_IV_Con_allfe_conf)


### Battles:
# Country + year Fixed Effects:
est_IV_Con_cyfe_btl = IV_Pois_Conley(battles_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df)
# Country + year + Country^year Fixed Effects:
#### Manually calculate functions for interaction of FEs:
#cutoff = 100
## Country^year FE:
#fe_ols = feols(battles_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country^year | dams_pby ~ RGxD_hat, df)
#resid_1st_stage = resid(summary(fe_ols, stage=1))
#y = model.matrix(fe_ols, type = 'lhs')
#RHS = model.matrix(fe_ols, type=c('iv.endo', 'iv.exo'))
#RHS_control_function = cbind(RHS, resid_1st_stage)
#FE = data.frame(Country = df$Country, year = df$year)
#lat = df$lat
#lon = df$lon
#new_data = cbind(y, RHS_control_function, FE, lat, lon)
#res = feglm(y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn + resid_1st_stage| Country^year , new_data, family='poisson')
#vcov_mrx = vcov_conley(res, cutoff = cutoff, lat = "lat", lon = "lon")
#est_IV_Con_cpyfe = summary(res, .vcov=vcov_mrx)

# Country^year FE:
fe_ols = feols(battles_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year | dams_pby ~ RGxD_hat, df)
resid_1st_stage = resid(summary(fe_ols, stage=1))
y = model.matrix(fe_ols, type = 'lhs')
RHS = model.matrix(fe_ols, type=c('iv.endo', 'iv.exo'))
RHS_control_function = cbind(RHS, resid_1st_stage)
FE = data.frame(Country = df$Country, year = df$year)
lat = df$lat
lon = df$lon
new_data = cbind(y, RHS_control_function, FE, lat, lon)
res = feglm(y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn + resid_1st_stage| Country + year + Country^year , new_data, family='poisson')
vcov_mrx = vcov_conley(res, cutoff = cutoff, lat = "lat", lon = "lon")
est_IV_Con_allfe_btl = summary(res, .vcov=vcov_mrx)
etable(est_IV_Con_cyfe_btl, est_IV_Con_allfe_btl)


### All riots:
# Country + year Fixed Effects:
est_IV_Con_cyfe_riot = IV_Pois_Conley(riots_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df)
# Country + year + Country^year Fixed Effects:
fe_ols = feols(riots_py ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year + Country^year | dams_pby ~ RGxD_hat, df)
resid_1st_stage = resid(summary(fe_ols, stage=1))
y = model.matrix(fe_ols, type = 'lhs')
RHS = model.matrix(fe_ols, type=c('iv.endo', 'iv.exo'))
RHS_control_function = cbind(RHS, resid_1st_stage)
FE = data.frame(Country = df$Country, year = df$year)
lat = df$lat
lon = df$lon
new_data = cbind(y, RHS_control_function, FE, lat, lon)
res = feglm(y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn + resid_1st_stage| Country + year + Country^year , new_data, family='poisson')
vcov_mrx = vcov_conley(res, cutoff = cutoff, lat = "lat", lon = "lon")
est_IV_Con_allfe_riot = summary(res, .vcov=vcov_mrx)
etable(est_IV_Con_cyfe_riot, est_IV_Con_allfe_riot)
