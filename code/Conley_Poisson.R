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

# Get the VCOV matrix:
IV_Pois_Conley_VarCov = function(fml, df, cutoff=100) {
  vcov_conley(IV_Pois_Conley_est(fml, df), cutoff = cutoff, lat = "lat", lon = "lon")
}

# Function to get estimates and SEs:
IV_Pois_Conley = function(fml, df, cutoff=100) {
  res = IV_Pois_Conley_est(fml, df)
  vcov_mrx = IV_Pois_Conley_VarCov(fml, df, cutoff)
  summary(res, .vcov=vcov_mrx)
}

est_IV_Con_byfe = IV_Pois_Conley(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + year | dams_pby ~ RGxD_hat, df1)
est_IV_Con_cyfe = IV_Pois_Conley(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country + year | dams_pby ~ RGxD_hat, df1)
#est_IV_Con_cpyfe = IV_Pois_Conley(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country^year | dams_pby ~ RGxD_hat, df1)             - does not operate due to interaction of FEs 
#est_IV_Con_allfe = IV_Pois_Conley(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df1)  - does not operate due to interaction of FEs 

### Manually calculate functions for interaction of FEs:
cutoff = 100
# Country^year FE:
fe_ols = feols(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country^year | dams_pby ~ RGxD_hat, df1)
resid_1st_stage = resid(summary(fe_ols, stage=1))
y = model.matrix(fe_ols, type = 'lhs')
RHS = model.matrix(fe_ols, type=c('iv.endo', 'iv.exo'))
RHS_control_function = cbind(RHS, resid_1st_stage)
FE = data.frame(Country = df1$Country, year = df1$year)
lat = df1$lat
lon = df1$lon
new_data = cbind(y, RHS_control_function, FE, lat, lon)
res = feglm(y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn + resid_1st_stage| Country^year , new_data, family='poisson')
vcov_mrx = vcov_conley(res, cutoff = cutoff, lat = "lat", lon = "lon")
est_IV_Con_cpyfe = summary(res, .vcov=vcov_mrx)

# HYBAS_ID + Country^year FE:
fe_ols = feols(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df1)
resid_1st_stage = resid(summary(fe_ols, stage=1))
y = model.matrix(fe_ols, type = 'lhs')
RHS = model.matrix(fe_ols, type=c('iv.endo', 'iv.exo'))
RHS_control_function = cbind(RHS, resid_1st_stage)
FE = data.frame(HYBAS_ID = df1$HYBAS_ID, Country = df1$Country, year = df1$year)
lat = df1$lat
lon = df1$lon
new_data = cbind(y, RHS_control_function, FE, lat, lon)
res = feglm(y ~ dams_pby + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn + resid_1st_stage| HYBAS_ID + Country^year , new_data, family='poisson')
vcov_mrx = vcov_conley(res, cutoff = cutoff, lat = "lat", lon = "lon")
est_IV_Con_allfe = summary(res, .vcov=vcov_mrx)

etable(est_IV_Con_byfe, est_IV_Con_cyfe, est_IV_Con_cpyfe, est_IV_Con_allfe, tex=TRUE)
