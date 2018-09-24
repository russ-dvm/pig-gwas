
## this was used in round 2 models.

covariates <- data.frame("sample" = as.character(
    c(
      "sal_shed_v2",
      "sal_shed_v3",
      "sal_shed_v4",
      "sal_shed_v5",
      "sal_shed_v2-5",
      "sal_shed_v2-6",
      
      "sal_sero_v2",
      "sal_sero_v3",
      "sal_sero_v4",
      "sal_sero_v5",
      
      "sal_titre_v2",
      "sal_titre_v3",
      "sal_titre_v4",
      "sal_titre_v5"
    )
  ), "covars" = as.character(
    c(
      "sal2,seas2_int", #shed_v2
      "sal3", #shed_v3
      "sal4,age4,trial_int,sero4_int", #shed_v4
      "sal5,age5,seas5_int", #shed_v5
      "diet_int,trial_int,sal25", #shed_v2-5
      "diet_int,trial_int,sal26", #shed_v2-6
      
      "seas2_int,sero2", #sero_v2
      "farm_int,sero3", #sero_v3
      "farm_int,sero4,sal4", #sero_v4
      "farm_int,sero5", #sero_v5
      
      "seas2_int", #titre_v2
      "farm_int,sal3,diet_int", #titre_v3
      "sal4,trial_int", #titre_v4
      "farm_int" #titre_v5
      
    )
  )
)



assign_covars <- function(x) {
  
  sample_row <- which(covariates[,1] == x)
  
  covars <- strsplit(as.character(covariates[sample_row, 2]), ",")[[1]]
  
  return(covars)
  
}
