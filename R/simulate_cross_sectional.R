#' Simulates multiple cohorts in a cross-sectional framework.
#' 
#' @param parTab the full parameter table controlling parameter ranges and values
#' @param group which group index to give this simulated data
#' @param n_indiv number of individuals to simulate
#' @param strainIsolationTimes vector of strain circulation times
#' @param samplingTimes possible sampling times for the individuals
#' @param antigenicMap the raw antigenic map with colnames x_coord, y_coord and inf_years
#' @param titreSensoring what proportion of titres are randomly missing?
#' @param ageMin minimum age to simulate
#' @param ageMax maximum age to simulate
#' @param attackRates a vector of attackRates to be used in the simulation
#' @return a list with: 1) the data frame of titre data as returned by \code{\link{simulate_group}}; 2) a matrix of infection histories as returned by \code{\link{simulate_infection_histories}}; 3) a vector of ages
#' @export

simulate_cross_sectional <- function(parTab, n_indiv, buckets = 12,
                                     strainIsolationTimes, samplingTimes,
                                     antigenicMap, ageMin = 5, ageMax = 80,
                                     attackRates, group = 1,
                                     sampleSensoring = 0, titreSensoring = 0) {

  ## What is the range of strain isolation times
  xs <- min(strainIsolationTimes):max(strainIsolationTimes)

  ## How many cohorts to simulate
  N.cohorts <- length(samplingTimes)

  ## Order sampling times in ascending order
  samplingTimes <- sort(samplingTimes, decreasing = TRUE)

  if (N.cohorts == 1) { # if there is only one cohort then we just need to call this line once
    if (max(strainIsolationTimes) != samplingTimes) warning("Cross-sectional study of one year is not testing strains cirulating in that year")
    dat <- simulate_data(parTab, 1, n_indiv, buckets, strainIsolationTimes,
      samplingTimes, length(samplingTimes),
      antigenicMap = fit_dat, sampleSensoring, titreSensoring, ageMin * buckets, ageMax * buckets, attackRates
    )
    titreDat.final <- dat[[1]]
    # titreDat.final <- titreDat[titreDat$virus %in% viruses,]

    infectionHistories.final <- infHist.final <- infectionHistories <- infHist <- dat[[2]]
    ages.final <- ages <- dat[[3]]
    AR.final <- AR <- dat[[4]]

    n_alive <- sapply(xs, function(x) nrow(ages[ages$DOB <= x, ]))
    numbers.final <- AR.final[, 2] * n_alive
  } else {

    # first cohort
    dat <- simulate_data(
      parTab, 1, n_indiv, buckets, strainIsolationTimes,
      samplingTimes[1], length(samplingTimes[1]), antigenicMap, sampleSensoring, titreSensoring, ageMin * buckets, ageMax * buckets, attackRates
    )
    titreDat.final <- dat[[1]]
    # titreDat.final <- titreDat <- titreDat[titreDat$virus %in% viruses,]

    infectionHistories.final <- infHist.final <- infectionHistories <- infHist <- dat[[2]]
    ages.final <- ages <- dat[[3]]
    AR.final <- AR <- dat[[4]]

    n_alive <- sapply(xs, function(x) nrow(ages.final[ages.final$DOB <= x, ]))
    numbers.final <- AR.final[, 2] * n_alive


    for (i in 2:N.cohorts) {
      tmp.dat <- simulate_data(
        parTab, 1, n_indiv, buckets, strainIsolationTimes,
        samplingTimes[i], length(samplingTimes[i]), antigenicMap, sampleSensoring, titreSensoring, ageMin * buckets, ageMax * buckets, attackRates
      )
      tmp.titreDat <- tmp.dat[[1]]
      # tmp.titreDat <- tmp.titreDat[tmp.titreDat$virus %in% viruses,]

      tmp.infectionHistories <- tmp.infHist <- tmp.dat[[2]]
      tmp.ages <- tmp.dat[[3]]
      tmp.AR <- tmp.dat[[4]]

      n_alive <- sapply(xs, function(x) nrow(tmp.ages[tmp.ages$DOB <= x, ]))
      tmp.numbers <- tmp.AR[, 2] * n_alive
      # if there are observed samples, then append files, if not, skip
      if (dim(tmp.titreDat)[1] != 0) {
        # adjust ages
        tmp.titreDat$individual <- tmp.titreDat$individual + max(titreDat.final$individual)
        tmp.ages$individual <- tmp.ages$individual + max(ages.final$individual)

        # tmp.infectionHistories$individual


        # append everything
        titreDat.final <- rbind(titreDat.final, tmp.titreDat)
        infectionHistories.final <- rbind(infectionHistories.final, tmp.infectionHistories)
        infHist.final <- rbind(infHist.final, tmp.infHist)
        ages.final <- rbind(ages.final, tmp.ages)
        AR.final <- rbind(AR.final, tmp.AR)
        numbers.final <- rbind(numbers.final, tmp.numbers)
      }
    }
  }

  return(list(titreDat.final, infHist.final, ages.final, AR.final, numbers.final))
}
