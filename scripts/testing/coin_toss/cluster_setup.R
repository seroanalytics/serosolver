setwd("~/net/home/coinflip")


## Submit credentials to didewin cluster
options(didehpc.credentials = "~/.smbcredentials",
        didehpc.cluster = "fi--didemrchnb")
                                        # didehpc.cluster = "fi--dideclusthn")
## Local location of the antibodyKinetics and lazymcmc packages
src <- provisionr::package_sources()
## Full file path to run submission function
sources <- c("~/net/home/coinflip/scripts/run_test.R",
             "~/net/home/coinflip/scripts/model_funcs.R",
             "~/net/home/coinflip/scripts/probability_funcs.R",
             "~/net/home/coinflip/scripts/proposal_funcs.R",
             "~/net/home/coinflip/scripts/proposal_theta.R",
             "~/net/home/coinflip/scripts/mcmc_funcs.R",
             "~/net/home/coinflip/scripts/mcmc_funcs_gibbs.R")

## Setup contexts
context::context_log_start()
root <- "contexts"
packages <- c("coda","plyr","reshape2")

ctx <- context::context_save(packages=packages,path=root, sources=sources,package_sources=src)

## Submit setup to cluster
obj <- didehpc::queue_didehpc(ctx)
