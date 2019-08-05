Introduction
------------

In disease systems where a single infection is expected in a given
period of time, serological data can be combined with models of antibody
kinetics to infer past infection times and the force of infection on a
population \[1–4\]. However, simple models of antibody kinetics are
often incomplete for viruses, such as influenza outside of a single
season, due to the antigenically variable nature of the virus. Whereas
raised antibody titres against a pathogen with a single circulating
phenotype are indicative of infection, cross-reactive antibodies from
antigenically similar influenza viruses results in a measurable antibody
response against strains that an individual may not have seen \[5–7\].
Understanding the immunological processes which build an individual’s
antibody repertoire is further complicated by the interference of a
pre-existing memory B cell response on the production of new B cell
populations against the infecting strain \[8,9\].

Progress has been made in understanding the immunological mechanisms
which link observed immune responses with disease incidence through
experimental systems and human population studies \[10–14\]. Taking into
account phenomena such as original antigenic sin (OAS), antigenic
seniority, back-boosting and age-specific antibody responses have been
used to significantly improve our interpretation of serological data
\[15–18\]. The full repertoire of humoral responses by one individual to
all previous and future strains can be described by an antibody
landscape, which quantifies the expected antibody titre against viruses
on an antigenic map \[17\]. Quantifying these mechanisms, along with
antibody boosting and cross reactivity, can therefore be used to better
predict antibody landscapes for humans \[19\].

### Model Description

`serosolver` uses a dynamic model of antibody dynamics, previously
described by Kucharski *et al.* \[20\] to infer infection histories and
attack rates from cross-sectional or longitudinal serological data. The
model infers individual infection histories, historical attack rates,
and patterns of antibody dynamics by accounting for short-term, broadly
reactive antibody boosting following exposure; long-term,
narrowly-reactive antibody boosting; boosting suppression through
antigenic seniority; and measurement error.

### Vignette Outline

This quick start guide is designed to provide the basic instructions to
use `serosolver` on serological data. For more in depth instructions
about the different capabilities of `serosolver` see the case study
vignettes.The following steps are covered in this guide:

1.  Formatting and inputting a serological data set
2.  Specifying inputs
    -   Parameter values
    -   Antigenic map
    -   Age mask
    -   Starting infection histories
3.  Running MCMC
4.  Processing outputs and generating plots

Install serosolver
------------------

The easiest way to install the development version of `serosolver` is to
use the `devtools` package:

    #devtools::install_github("seroanalytics/serosolver")
    library(serosolver)
    library(plyr)
    library(data.table)
    library(ggplot2)

1. Data Format
--------------

`serosolver` expects a data frame as input in long format with the
following columns:

1.  individual: consecutive integer ID of individuals from 1, …, *N*,
    where *N* is the total number of individuals
2.  samples: integer value of the time each sample was collected
    -   For annual samples, the sample time is simply the year each
        sample was collected in. For example, data collected between
        2009 and 2011 would have values of 2009, 2010, or 2011.
    -   For semi-annual samples, the sample time is the year the sample
        was collected multiplied by 2 (since there are two semesters in
        every year). For example, data collected between 2009 and 2011
        would have values of (2009 × 2) + 1 = 4019,
        (2009.5 × 2) + 1 = 4020, (2010 × 2) + 1 = 4021,
        (2010.5 × 2) + 1 = 4022, (2011 × 2) + 1 = 4023,
        (2011.5 × 2) + 1 = 4024.
    -   Similarly, for quarterly or monthly samples, the sample time
        follows the same pattern. For example for a sample collected in
        July of 2009, the quarterly sample time would be
        (2009.75 × 4) + 1 = 8040 (because July is in the third quarter
        of the year) and the monthly sample time would be
        (2009.5833 × 12) + 1 = 24116 (because July is the seventh month
        and 7/12 = 0.5833).  
3.  virus: numeric time of when the virus was circulating (should be in
    the same format as samples)
4.  titre: integer of titre value against the given virus at that
    sampling time (needs to be on log scale)
5.  run: integer value of the repeat number of each sample (if there are
    no repeat samples, then run should be 1 for every sample)
6.  DOB: integer value of date of birth (should be in the same format as
    samples)
7.  group: integer value giving the population membership of the
    individual

<!-- -->

      titre_dat <- data.frame(individual=c(rep(1,4),rep(2,4)),
                              samples=c(8039,8040,8044,8047,8039,8041,8045,8048),
                              virus=c(rep(8036,8)),
                              titre=c(0,0,7,7,0,5,6,5),
                              run=c(rep(1,8)),
                              DOB=c(rep(8036,8)),
                              group=c(rep(1,8))
                             )
      knitr::kable(head(titre_dat))

For the remainder of this quick start guide we refer to the input data
set as `titre_dat`.

2. Specifying inputs
--------------------

Along with a data file, serosolver expects several additional inputs
namely, a parameter input file for starting values of antibody kinetics
parameters, an antigenic map, an age mask, and starting infection
histories. The following section will walk the users through creating
each of these inputs.

### Parameter values

Several input files are stored within the `serosolver` package in the
‘inputs’ directory, however users can create their own. The default
parameter file ‘parTab\_base.csv’ is shown below:

<table>
<thead>
<tr class="header">
<th></th>
<th style="text-align: left;">names</th>
<th style="text-align: right;">values</th>
<th style="text-align: right;">fixed</th>
<th style="text-align: right;">steps</th>
<th style="text-align: right;">lower_bound</th>
<th style="text-align: right;">upper_bound</th>
<th style="text-align: right;">lower_start</th>
<th style="text-align: right;">upper_start</th>
<th style="text-align: right;">type</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1</td>
<td style="text-align: left;">mu</td>
<td style="text-align: right;">1.80</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">1.00</td>
<td style="text-align: right;">3.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td>2</td>
<td style="text-align: left;">mu_short</td>
<td style="text-align: right;">2.70</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">2.00</td>
<td style="text-align: right;">3.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="odd">
<td>3</td>
<td style="text-align: left;">tau</td>
<td style="text-align: right;">0.05</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.01</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td>4</td>
<td style="text-align: left;">wane</td>
<td style="text-align: right;">0.20</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.01</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="odd">
<td>5</td>
<td style="text-align: left;">sigma1</td>
<td style="text-align: right;">0.10</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.10</td>
<td style="text-align: right;">0.2</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td>6</td>
<td style="text-align: left;">sigma2</td>
<td style="text-align: right;">0.03</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.01</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="odd">
<td>7</td>
<td style="text-align: left;">MAX_TITRE</td>
<td style="text-align: right;">8.00</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">8.00</td>
<td style="text-align: right;">8.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td>8</td>
<td style="text-align: left;">error</td>
<td style="text-align: right;">1.00</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">25</td>
<td style="text-align: right;">0.50</td>
<td style="text-align: right;">2.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="odd">
<td>9</td>
<td style="text-align: left;">alpha</td>
<td style="text-align: right;">1.00</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">0.00</td>
<td style="text-align: right;">1000.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td>10</td>
<td style="text-align: left;">beta</td>
<td style="text-align: right;">1.00</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">0.00</td>
<td style="text-align: right;">1000.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="odd">
<td>11</td>
<td style="text-align: left;">kappa</td>
<td style="text-align: right;">0.90</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.80</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td>12</td>
<td style="text-align: left;">t_change</td>
<td style="text-align: right;">6.00</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">0.40</td>
<td style="text-align: right;">12.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="odd">
<td>13</td>
<td style="text-align: left;">boost_limit</td>
<td style="text-align: right;">6.00</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">0.00</td>
<td style="text-align: right;">8.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td>14</td>
<td style="text-align: left;">gradient</td>
<td style="text-align: right;">0.20</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.00</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="odd">
<td>15</td>
<td style="text-align: left;">titre_dependent</td>
<td style="text-align: right;">0.00</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.00</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td>16</td>
<td style="text-align: left;">wane_type</td>
<td style="text-align: right;">0.00</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1.00</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td>18</td>
<td style="text-align: left;">mu_mean</td>
<td style="text-align: right;">2.00</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">-100</td>
<td style="text-align: right;">100</td>
<td style="text-align: right;">1.00</td>
<td style="text-align: right;">3.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td>19</td>
<td style="text-align: left;">mu_sd</td>
<td style="text-align: right;">0.50</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">0.10</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="odd">
<td>20</td>
<td style="text-align: left;">rho_mean</td>
<td style="text-align: right;">0.00</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">-10</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">-2.00</td>
<td style="text-align: right;">2.0</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td>21</td>
<td style="text-align: left;">rho_sd</td>
<td style="text-align: right;">1.00</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">1.00</td>
<td style="text-align: right;">2.0</td>
<td style="text-align: right;">1</td>
</tr>
</tbody>
</table>

In the parameter file, the column ‘fixed’ indicates if each parameter
value is fixed or should be estimated. For values that are not fixed,
starting values need to be generated. The below code generates a random
starting value for unfixed parameter values.

    # generate starting parameter values
      start_tab <- par_tab
        for(i in 1:nrow(start_tab)){
          if(start_tab[i,"fixed"] == 0){
            start_tab[i,"values"] <- runif(1,start_tab[i,"lower_start"], 
                                      start_tab[i,"upper_start"])
          }
        }

### Antigenic Map

The antigenic map specifies the coordinates of different circulating
viruses in antigenic space \[21\], informing the model about how
antigenically similar these viruses are. Viruses that are further apart
in this 2-dimensional space are expected to cross react less than
viruses that are close together.

<img src="/Library/Frameworks/R.framework/Versions/3.6/Resources/library/serosolver/extdata/antigenic_map.tiff" alt="Assumed antigenic locations of historical strains in model between 1968 and 2012 (from Kucharski *et al.* @kucharski2018)." width="75%" />
<p class="caption">
Assumed antigenic locations of historical strains in model between 1968
and 2012 (from Kucharski *et al.* \[20\]).
</p>

An antigenic map can be created using `generate_antigenic_map()` and
input by the user depending on the source of their data. We use the
antigenic map created by Fonville *et al.* \[17\] in our analyses.

    ## Read in raw coordinates
      antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
      antigenic_coords <- read.csv(antigenic_coords_path, stringsAsFactors=FALSE)

    ## Convert to form expected by serosolver
      antigenic_map <- generate_antigenic_map(antigenic_coords, buckets = 4)

    # unique strain circulation times
      strains_isolation_times <- unique(antigenic_map$inf_years)

### Age Mask

An age mask is required to indicate when each individual could have had
their first infection. For example an individual born in 1997 could not
have been infected with a strain that circulated in 1990. The function
`create_age_mask()` outputs a vector of indices of the first
`strains_isolation_times` that an individual could have been infected
with.

      unique_indiv <- titre_dat[!duplicated(titre_dat$individual),]
      ageMask <- create_age_mask(unique_indiv$DOB, strains_isolation_times)

The final inputs required by `serosolver` are starting infection
histories for each individual. The function
`setup_infection_histories_new_2()` creates a matrix of infection
histories given a matrix of titre data by looking at an individual’s
titre against each strain. Where titres are raised, it suggests an
infection. In this way, it can propose plausible initial infection
histories from which to begin MCMC sampling.

      start_inf <- setup_infection_histories_new_2(titre_dat, strains_isolation_times)

3. Run MCMC
-----------

`seroslver` uses an Adaptive Metropolis-within-Gibbs algorithm. Given a
starting point and the necessary MCMC parameters, `run_MCMC()` performs
a random-walk of the posterior space to produce an MCMC chain that can
be used to generate MCMC density and iteration plots (see section 4).

Using the default options, the MCMC can be run using the following.

    # run MCMC  
      res <- run_MCMC(par_tab = start_tab, titre_dat = titre_dat, antigenic_map = antigenic_map, 
                      start_inf_hist = start_inf, CREATE_POSTERIOR_FUNC = create_posterior_func,
                      version = 2)

If you want to change any default options, such as increasing or
decreasing the number of iterations and the adaptive period you can pass
a vector of function argument values to `run_MCMC()`.

      mcmc_pars <- c("iterations"=50000,"popt"=0.44,"popt_hist"=0.44,
                    "opt_freq"=1000,"thin"=1,
                    "adaptive_period"=10000, "save_block"=1000, 
                    "thin_hist"=10, "hist_sample_prob"=0.5,
                    "switch_sample"=2, "inf_propn"=0.5, 
                    "move_size"=5, "hist_opt"=1,
                    "swap_propn"=0.5,"hist_switch_prob"=0.5,
                    "year_swap_propn"=0.5)
      res <- run_MCMC(par_tab = start_tab, titre_dat = titre_dat, antigenic_map = antigenic_map, 
                    mcmc_pars = mcmc_pars, start_inf_hist = start_inf,
                    CREATE_POSTERIOR_FUNC = create_posterior_func, version = 2)

4. Processing Outputs
---------------------

`serosolver` contains numerous plotting functions for visualising
results and evaluating model performance. Below are examples of several
useful plotting functions.

### Diagnostic plots

To evaluate model performance, trace and density plots can be produced
by plotting the MCMC chains. In the below example the adaptive period
and burn-in from the MCMC are excluded in the diagnostic plots.

    # Density/trace plots
    chain1 <- read.csv(res$chain_file)
    chain1 <- chain1[chain1$sampno >= (mcmc_pars["adaptive_period"]+mcmc_pars["burnin"]),]
    plot(coda::as.mcmc(chain1[,c("mu","mu_short","wane")]))

### Results plots

Plot inferred attack rates.

    # Read in infection history file from MCMC output
      inf_chain <- data.table::fread(res$history_file,data.table=FALSE)
    # Remove adaptive period and burn-in
      inf_chain <- inf_chain[inf_chain$sampno >= (mcmc_pars["adaptive_period"]+mcmc_pars["burnin"]),]
      inf_chain1 <- setDT(inf_chain)
    # Define year range
      xs <- min(strains_isolation_times):max(strains_isolation_times)
    # Plot inferred attack rates
      arP <- plot_attack_rates(infection_histories = inf_chain1, titre_dat = titre_dat, strain_isolation_times=xs)

Plot inferred infection histories.

    # Plot infection histories
      IH_plot <- plot_infection_histories(chain = chain1, infection_histories = inf_chain, 
                                          titre_dat = titre_dat, individuals = c(1:5),
                                          antigenic_map = antigenic_map, par_tab = start_tab1)

Plot inferred antibody titres.

    # Plot inferred antibody titres
      titre_preds <- get_titre_predictions(chain = chain1, infection_histories = inf_chain, 
                                           titre_dat = titre_dat, individuals = c(1:5),
                                           antigenic_map = antigenic_map, par_tab = start_tab1)
      to_use <- titre_preds$predictions
      
      titre_pred_p <- ggplot(to_use[to_use$individual %in% 1:5,])+
                      geom_line(aes(x=samples, y=median))+
                      geom_point(aes(x=samples, y=titre))+
                      geom_ribbon(aes(x=samples,ymin=lower, ymax=upper),alpha=0.2,col='red')+
                      facet_wrap(~individual)

Generate cumulative incidence plots.

      y <- generate_cumulative_inf_plots(inf_chain,burnin = 0,1:5,nsamp=100,
                                       strain_isolation_times = strains_isolation_times,
                                       pad_chain=FALSE,number_col = 2,subset_years = NULL)

### Parameter estimates

The following code produces a table of the parameter estimates.

    # Table of antibody kinetics parameters
      myresults <- matrix(c(rep(0,3*7)),nrow=3)
      rownames(myresults) <- c("mu_short","wane","error")
      colnames(myresults) <- c("mean","sd","2.5%","25%","50%","75%","97.5%")
      
      myresults[,"mean"] <- round(apply(chain1[,c("mu_short","wane","error")],2,mean),3)
      myresults[,"sd"] <- round(apply(chain1[,c("mu_short","wane","error")],2,sd),3)  
      myresults[,3:7] <- t(round(apply(chain1[,c("mu_short","wane","error")],2,
                                       quantile,probs=c(0.025,0.25,0.5,0.75,0.975)),3))  

References
----------

1. White MT, Griffin JT, Akpogheneta O, Conway DJ, Koram KA, Riley EM,
et al. Dynamics of the antibody response to Plasmodium falciparum
infection in African children. J Infect Dis. 2014;210: 1115–1122.

2. Metcalf CJE, Farrar J, Cutts FT, Basta NE, Graham AL, Lessler J, et
al. Use of serological surveys to generate key insights into the
changing global landscape of infectious disease. Lancet. 2016;388:
728–730.
doi:[10.1016/S0140-6736(16)30164-7](https://doi.org/10.1016/S0140-6736(16)30164-7)

3. Borremans B, Hens N, Beutels P, Leirs H, Reijniers J, Khan A, et al.
Estimating Time of Infection Using Prior Serological and Individual
Information Can Greatly Improve Incidence Estimation of Human and
Wildlife Infections. Salathé M, editor. PLOS Comput Biol. Public Library
of Science; 2016;12: e1004882.
doi:[10.1371/journal.pcbi.1004882](https://doi.org/10.1371/journal.pcbi.1004882)

4. Pepin KM, Kay SL, Golas BD, Shriner SS, Gilbert AT, Miller RS, et al.
Inferring infection hazard in wildlife populations by linking data
across individual and population scales. Ecol Lett. 2017;
doi:[10.1111/ele.12732](https://doi.org/10.1111/ele.12732)

5. Beyer WEP, Palache AM, Lüchters G, Nauta J, Osterhaus ADME.
Seroprotection rate, mean fold increase, seroconversion rate: which
parameter adequately expresses seroresponse to influenza vaccination?
Virus Res. 2004;103: 125–32.
doi:[10.1016/j.virusres.2004.02.024](https://doi.org/10.1016/j.virusres.2004.02.024)

6. Beest D te, Bruin E de, Imholz S, Wallinga J, Teunis P, Koopmans M,
et al. Discrimination of Influenza Infection (A/2009 H1N1) from Prior
Exposure by Antibody Protein Microarray Analysis. Tang JW, editor. PLoS
One. Cambridge University Press; 2014;9: e113021.
doi:[10.1371/journal.pone.0113021](https://doi.org/10.1371/journal.pone.0113021)

7. FREEMAN G, PERERA RAPM, NGAN E, FANG VJ, CAUCHEMEZ S, IP DKM, et al.
Quantifying homologous and heterologous antibody titre rises after
influenza virus infection. Epidemiol Infect. 2016;144: 2306–2316.
doi:[10.1017/S0950268816000583](https://doi.org/10.1017/S0950268816000583)

8. Andrews SF, Huang Y, Kaur K, Popova LI, Ho IY, Pauli NT, et al.
Immune history profoundly affects broadly protective B cell responses to
influenza. Sci Transl Med. 2015;7: 316ra192–316ra192.
doi:[10.1126/scitranslmed.aad0522](https://doi.org/10.1126/scitranslmed.aad0522)

9. Zarnitsyna VI, Lavine J, Ellebedy A, Ahmed R, Antia R. Multi-epitope
Models Explain How Pre-existing Antibodies Affect the Generation of
Broadly Protective Responses to Influenza. Lauring AS, editor. PLOS
Pathog. Public Library of Science; 2016;12: e1005692.
doi:[10.1371/journal.ppat.1005692](https://doi.org/10.1371/journal.ppat.1005692)

10. Li Y, Myers JL, Bostick DL, Sullivan CB, Madara J, Linderman SL, et
al. Immune history shapes specificity of pandemic H1N1 influenza
antibody responses. J Exp Med. 2013;210.

11. Linderman SL, Chambers BS, Zost SJ, Parkhouse K, Li Y, Herrmann C,
et al. Potential antigenic explanation for atypical H1N1 infections
among middle-aged adults during the 2013-2014 influenza season. Proc
Natl Acad Sci U S A. 2014;111: 15798–15803.

12. Gostic KM, Ambrose M, Worobey M, Lloyd-Smith JO. Potent protection
against H5N1 and H7N9 influenza via childhood hemagglutinin imprinting.
Science (80- ). 2016;354.

13. Cobey S, Hensley SE. Immune history and influenza virus
susceptibility. Curr Opin Virol. 2017;
doi:[10.1016/j.coviro.2016.12.004](https://doi.org/10.1016/j.coviro.2016.12.004)

14. Monto AS, Malosh RE, Petrie JG, Martin ET. The Doctrine of Original
Antigenic Sin: Separating Good From Evil. J Infect Dis. 2017;215:
1782–1788.

15. Kim JH, Skountzou I, Compans R, Jacob J. Original Antigenic Sin
Responses to Influenza Viruses. J Immunol. 2009;183: 3294–3301.
doi:[10.4049/jimmunol.0900398](https://doi.org/10.4049/jimmunol.0900398)

16. Lessler J, Riley S, Read JM, Wang S, Zhu H, Smith GJD, et al.
Evidence for antigenic seniority in influenza A (H3N2) antibody
responses in southern China. PLoS Pathog. Public Library of Science;
2012;8: e1002802.
doi:[10.1371/journal.ppat.1002802](https://doi.org/10.1371/journal.ppat.1002802)

17. Fonville JM, Wilks SH, James SL, Fox A, Ventresca M, Aban M, et al.
Antibody landscapes after influenza virus infection or vaccination.
Science (80- ). 2014;346: 7–9.

18. Mosterín Höpping A, McElhaney J, Fonville JM, Powers DC, Beyer WEP,
Smith DJ, et al. The confounded effects of age and exposure history in
response to influenza vaccination. Vaccine. NIH Public Access; 2016;34:
540–546.
doi:[10.1016/j.vaccine.2015.11.058](https://doi.org/10.1016/j.vaccine.2015.11.058)

19. Kucharski AJ, Lessler J, Read JM, Zhu H, Jiang CQ, Guan Y, et al.
Estimating the Life Course of Influenza A(H3N2) Antibody Responses from
Cross-Sectional Data. PLoS Biol. Public Library of Science; 2015;13:
1–16.
doi:[10.1371/journal.pbio.1002082](https://doi.org/10.1371/journal.pbio.1002082)

20. Kucharski AJ, Lessler J, Cummings DAT, Riley S. Timescales of
influenza a/h3n2 antibody dynamics. PLOS Biology. 2018;16: 1–19.
doi:[10.1371/journal.pbio.2004974](https://doi.org/10.1371/journal.pbio.2004974)

21. Smith DJ. Mapping the Antigenic and Genetic. 2004;305: 371–376.
doi:[10.1126/science.1097211](https://doi.org/10.1126/science.1097211)
