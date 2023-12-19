library(sf)
#library(lfe)
#library(dplyr)
library(fixest)

df = st_read("5rivers_fully_prepared_data.shp")


##################################################################
##################################################################
################# Run all possible combinations: ################# 

dependent_vars <- c("had_fight", "btl_p_y")
independent_vars <- c("z_fobki", "z_hyeah", "z_vcjei", "z_nlvsk", "z_ahjvn", "z_zgjij")
conley_values <- c(100, 200, 400)

results <- data.frame(dependent_var = character(), 
                      independent_vars = character(), 
                      conley = numeric(),
                      estimate = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)
results1 <- data.frame(dependent_var = character(), 
                      independent_vars = character(), 
                      conley = numeric(),
                      estimate = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

get_combinations <- function(vars, max_length) {
  unlist(lapply(1:max_length, function(n) combn(vars, n, simplify = FALSE)), recursive = FALSE)
}
restricted_combinations <- get_combinations(independent_vars[3:6], 3)
i = 0

for (dep_var in dependent_vars) {
  for (indep_comb in restricted_combinations) {
    i <- i+1
    print(i/length(restricted_combinations))
    for (conley_value in conley_values) {
      if (dep_var == "had_fight" && identical(indep_comb, c("z_vcjei", "z_nlvsk", "z_ahjvn", "z_zgjij"))  && conley_value == 100) {
        next
      }
      if (dep_var == "had_fight" && identical(indep_comb, c("z_fobki", "z_vcjei", "z_nlvsk", "z_ahjvn", "z_zgjij")) && conley_value == 100) {
        next
      }
      formula <- as.formula(paste(dep_var, "~", paste(c("z_fobki", "z_hyeah", indep_comb), collapse = " + "), 
                                  "| HYBAS_ID + Country^year | dams_pby ~ RGxD_hat"))
      try <- feols(formula, df, conley(conley_value, dist="spherical"))
      results <- rbind(results, data.frame(dependent_var = dep_var, 
                                           independent_comb = paste(c("z_fobki", "z_hyeah", indep_comb), collapse = ", "), 
                                           conley = conley_value,
                                           estimate = try$coeftable[1, 1],
                                           p_value = try$coeftable[1, 4]))
    }
  }
}

write.csv(results, "2023_all_models_feols,csv",row.names = FALSE)
write.csv(results, "2017_all_models_feols,csv",row.names = FALSE)

#################### Simulations are finished ####################
##################################################################
##################################################################

# Simple OLS:
ols1 = feols(had_conf ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + year, df, vcov = conley(100, distance = "spherical"))
ols2 = feols(had_conf ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country, df, vcov = conley(100, distance = "spherical"))
ols3 = feols(had_conf ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country + Country^year, df, vcov = conley(100, distance = "spherical"))
ols4 = feols(had_conf ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + Country^year, df, vcov = conley(100, distance = "spherical"))
etable(ols1, ols2, ols3, ols4)

ols5 = feols(had_fight ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + year, df, vcov = conley(100, distance = "spherical"))
ols6 = feols(had_fight ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country, df, vcov = conley(100, distance = "spherical"))
ols7 = feols(had_fight ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country + Country^year, df, vcov = conley(100, distance = "spherical"))
ols8 = feols(had_fight ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + Country^year, df, vcov = conley(100, distance = "spherical"))
etable(ols5, ols6, ols7, ols8)

ols9 = feols(had_riot ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + year, df, vcov = conley(100, distance = "spherical"))
ols10 = feols(had_riot ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country, df, vcov = conley(100, distance = "spherical"))
ols11 = feols(had_riot ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country + Country^year, df, vcov = conley(100, distance = "spherical"))
ols12 = feols(had_riot ~ dams_pby + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + Country^year, df, vcov = conley(100, distance = "spherical"))
etable(ols9, ols10, ols11, ols12)

# Reduced forms:
rf1 = feols(had_conf ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + year, df, vcov = conley(100, distance = "spherical"))
rf2 = feols(had_conf ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country, df, vcov = conley(100, distance = "spherical"))
rf3 = feols(had_conf ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country + Country^year, df, vcov = conley(100, distance = "spherical"))
rf4 = feols(had_conf ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + Country^year, df, vcov = conley(100, distance = "spherical"))
etable(rf1, rf2, rf3, rf4)

rf5 = feols(had_fight ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + year, df, vcov = conley(100, distance = "spherical"))
rf6 = feols(had_fight ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country, df, vcov = conley(100, distance = "spherical"))
rf7 = feols(had_fight ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country + Country^year, df, vcov = conley(100, distance = "spherical"))
rf8 = feols(had_fight ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + Country^year, df, vcov = conley(100, distance = "spherical"))
etable(rf5, rf6, rf7, rf8)

rf9 = feols(had_riot ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + year, df, vcov = conley(100, distance = "spherical"))
rf10 = feols(had_riot ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country, df, vcov = conley(100, distance = "spherical"))
rf11 = feols(had_riot ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| year + Country + Country^year, df, vcov = conley(100, distance = "spherical"))
rf12 = feols(had_riot ~ RGxD_hat + z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + Country^year, df, vcov = conley(100, distance = "spherical"))
etable(rf9, rf10, rf11, rf12)

# 2SLS with FE:
tsls_1 = feols(had_conf ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + year | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
tsls_2 = feols(had_conf ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | year + Country | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
tsls_3 = feols(had_conf ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | year + Country + Country^year | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
tsls_4 = feols(had_conf ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
etable(tsls_1, tsls_2, tsls_3, tsls_4)

tsls_5 = feols(had_fight ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + year | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
tsls_6 = feols(had_fight ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | year + Country | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
tsls_7 = feols(had_fight ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | year + Country + Country^year | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
tsls_8 = feols(had_fight ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
etable(tsls_1, tsls_2, tsls_3, tsls_4)

tsls_9 = feols(had_riot ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + year | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
tsls_10 = feols(had_riot ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | year + Country | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
tsls_11 = feols(had_riot ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | year + Country + Country^year | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
tsls_12 = feols(had_riot ~ z_oxjpe + z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df, vcov = conley(100, distance = "spherical"))
etable(tsls_9, tsls_10, tsls_11, tsls_12)


# All conflicts:
etable(ols2, ols3, rf2, rf3, tsls_2, tsls_3, tex=TRUE)
# All fights:
etable(ols6, ols7, rf6, rf7, tsls_6, tsls_7, tex=TRUE)
# All riots:
etable(ols10, ols11, rf10, rf11, tsls_10, tsls_11, tex=TRUE)


# Trying the same estimations as in try_11, try_12, try_13 and try_14 but with btl_p_y as dependent var hasn't produced any significant results for dams_pby
# Trying the same estimations as in try_11, try_12, try_13 and try_14 but with SUB_AREA + tot_riv_di + av_elev + shr_l_25 + shr_25_50 + shr_50_1k as independent has resulted in similar estimates