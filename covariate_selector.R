
## this was used in round 2 models.

covariates <- data.frame("sample" = as.character(
    c(
      "sal_shed_v2",
      "sal_shed_v3",
      "sal_shed_v4",
      "sal_shed_v5",
      "sal_shed_v6",
      "sal_shed_v2-5",
      "sal_shed_v2-6",
      "sal_shed_v5-6",
      
      "sal_sero_v2",
      "sal_sero_v3",
      "sal_sero_v4",
      "sal_sero_v5",
      "sal_sero_v3-5",
      
      "sal_titre_v2",
      "sal_titre_v3",
      "sal_titre_v4",
      "sal_titre_v5"
    )
  ), "covars" = as.character(
    c(
      "sal2,seas2_int", #shed_v2
      "sal3", #shed_v3
      "sal4,age4,trial_int,sero4", #shed_v4
      "sal5,age5,farm,seas5_int", #shed_v5
      "sal6,trial,age6,seas6_int", #shed_v6
      "sal25,diet_int,trial_int", #shed_v2-5
      "sal26", #shed_v2-6
      "sal56,trial", #shed_v5-6
      
      "sero2,seas2_int", #sero_v2
      "sero3,farm", #sero_v3
      "sero4,farm,sal4", #sero_v4
      "sero5,farm", #sero_v5
      "sero35,farm", #sero_v3-5
      
      "seas2_int", #titre_v2
      "farm,sal3,diet_int", #titre_v3
      "sal4,trial_int", #titre_v4
      "farm" #titre_v5
      
    )
  )
)



assign_covars <- function(x) {
  
  sample_row <- which(covariates[,1] == x)
  
  covars <- strsplit(as.character(covariates[sample_row, 2]), ",")[[1]]
  
  return(covars)
  
}
