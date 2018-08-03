options_theta <- c("univariate_theta","multivariate_theta")
options_lambda <- c("yes","none")
options_Z <- c("naive_symmetric","gibbs","simple_flip","years","from_lambda")



all_options <- expand.grid("theta"=options_theta,"lambda"=options_lambda,"Z"=options_Z)

scale_pars <- c("indiv_propn","year_propn","adaptive","swap")

indiv_propn_values <- c(1,0.5,0.1)0
year_propn_values <- c(1,0.5,0.1)
adaptive <- c(TRUE,FALSE)
swap <- c(TRUE,FALSE)

all_tuning <- expand.grid("theta"=options_theta,"lambda"=options_lambda,"Z"=options_Z,
  "indiv"=indiv_propn_values, "year"=year_propn_values,"adaptive"=adaptive,"swap"=swap)


## Run 10 chains for each
n_indivs <- c(10,100,500)
n_coins <- c(10, 100, 200)