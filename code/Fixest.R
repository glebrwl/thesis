library(sf)
library(lfe)
library(dplyr)
library(fixest)

df <- st_read("5rivers_fully_prepared_data_all_conflicts_centroid.shp")
df1<-df[!(df$year==2018 | df$year==2019| df$year==2020| df$year==2021| df$year==2022 | df$year==2023),]

try_1 = feols(had_fight ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df, conley(400, distance = "spherical"))
try_1 = feols(btl_p_y ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df, NW(5) ~ HYBAS_ID + Country^year)
try_1 = feols(had_fight ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df)
z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn + z_zgjij

try_1 = feols(btl_p_y ~ z_fobki + z_hyeah + z_vcjei | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df1, conley(200, distance = "spherical"))
summary(try_1)
try_1$coeftable[1, 1]
summary(try, stage = 1)
try$coeftable[1,4]
summary(try_11, "conley")


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


#try_11 = feols(had_fight ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij| HYBAS_ID + year | dams_pby ~ RGxD_hat, df1, vcov = conley(990, distance = "spherical"))
#try_12 = feols(had_fight ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | year + Country | dams_pby ~ RGxD_hat, df1, vcov = conley(100, distance = "spherical"))
#try_13 = feols(had_fight ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | Country^year | dams_pby ~ RGxD_hat, df1, vcov = conley(100, distance = "spherical"))
#try_14 = feols(had_fight ~ z_fobki + z_hyeah + z_vcjei + z_ahjvn + z_zgjij | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df1, vcov = conley(100, distance = "spherical"))
#etable(try_11, try_12, try_13, try_14)

try_11 = feols(had_fight ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn| HYBAS_ID + year | dams_pby ~ RGxD_hat, df1, vcov = conley(100, distance = "spherical"))
try_12 = feols(had_fight ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | year + Country | dams_pby ~ RGxD_hat, df1, vcov = conley(100, distance = "spherical"))
try_13 = feols(had_fight ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | Country^year | dams_pby ~ RGxD_hat, df1, vcov = conley(100, distance = "spherical"))
try_14 = feols(had_fight ~ z_fobki + z_hyeah + z_vcjei + z_nlvsk + z_ahjvn | HYBAS_ID + Country^year | dams_pby ~ RGxD_hat, df1, vcov = conley(100, distance = "spherical"))
etable(try_11, try_12, try_13, try_14, tex=TRUE)

# Trying the same estimations as in try_11, try_12, try_13 and try_14 but with btl_p_y as dependent var hasn't produced any significant results for dams_pby
