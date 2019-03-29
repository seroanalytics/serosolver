#' Swap infection history years
#'
#' Swaps the entire contents of two columns of the infection history matrix, adhering to age and strain limitations.
#' @param infection_history matrix of 1s and 0s to swap around
#' @param age_mask the first index in infection_history that each individual (row) could be infected in
#' @param strain_mask the last index in infection_history that each individual (row) could be infected in ie. the time of the latest blood sample
#' @param swap_propn what proportion of infections should be swapped?
#' @param move_size How many time points away should be chosen as candidate swaps?
#' @return the same infection_history matrix, but with two columns swapped
#' @export
inf_hist_swap <- function(infection_history, age_mask, strain_mask, swap_propn, move_size) {
  ## Choose a column
  y1 <- sample(1:ncol(infection_history), 1)

  ## Propose another column some random distance, but not 0, away
  move <- 0
  while (move == 0) move <- sample((-move_size):move_size, 1)

  ## Need to adjust if we've proposed too far away
  y2 <- y1 + move
  while (y2 < 1) y2 <- y2 + ncol(infection_history)
  if (y2 > ncol(infection_history)) y2 <- y2 - floor(y2 / ncol(infection_history)) * ncol(infection_history)

  ## Get the first and last year chronologically
  small_year <- min(y1, y2)
  big_year <- max(y1, y2)

  ## Find individuals that are alive/sampled in both years and choose the lesser of swap_propn*n_indivs and
  ## the number that are actually able to be infected in both years
  indivs <- 1:nrow(infection_history)
  alive_indivs <- indivs[intersect(which(age_mask <= small_year), which(strain_mask >= big_year))]
  samp_indivs <- sample(alive_indivs, floor(length(alive_indivs) * swap_propn))

  ## Swap contents
  tmp <- infection_history[samp_indivs, y1]
  infection_history[samp_indivs, y1] <- infection_history[samp_indivs, y2]
  infection_history[samp_indivs, y2] <- tmp
  return(infection_history)
}
#' Swap infection history years with phi term
#'
#' Swaps the entire contents of two columns of the infection history matrix, adhering to age and strain limitations. Also swaps the values of phi that correspond to these years
#' @inheritParams inf_hist_swap
#' @param phis vector of force of infection parameters for each column
#' @param n_alive number of individuals alive in each entry of phis
#' @return a list: the same infection_history matrix, but with two columns swapped; also the swapped phis
#' @seealso \code{\link{inf_hist_swap}}
#' @export
inf_hist_swap_phi <- function(infection_history, phis, age_mask, strain_mask, swap_propn, move_size, n_alive) {
  ## This first bit of code is the same as inf_hist_swap
  y1 <- sample(1:ncol(infection_history), 1)
  move <- 0
  while (move == 0) move <- sample((-move_size):move_size, 1)
  y2 <- y1 + move

  while (y2 < 1) y2 <- y2 + ncol(infection_history)
  if (y2 > ncol(infection_history)) y2 <- y2 - floor(y2 / ncol(infection_history)) * ncol(infection_history)
  small_year <- min(y1, y2)
  big_year <- max(y1, y2)

  indivs <- 1:nrow(infection_history)

  alive_indivs <- indivs[intersect(which(age_mask <= small_year), which(strain_mask >= big_year))]
  samp_indivs <- sample(alive_indivs, floor(length(alive_indivs) * swap_propn))

  tmp <- infection_history[samp_indivs, y1]

  ## Get number of infections amongst selected individuals in the two years
  no_infs_y1 <- sum(tmp)
  no_infs_y2 <- sum(infection_history[samp_indivs, y2])

  ## Total number of infections in these years
  total_infs_y1 <- sum(infection_history[, y1])
  total_infs_y2 <- sum(infection_history[, y2])

  ## How many infections will there be in these two years after the swap?
  new_total_infs_y1 <- total_infs_y1 - no_infs_y1 + no_infs_y2
  new_total_infs_y2 <- total_infs_y2 - no_infs_y2 + no_infs_y1

  ## Adjust the FOI parameters proportional to the number of infections gained/lost
  lost_foi_y1 <- gained_foi_y1 <- 0
  if (total_infs_y1 > 0) {
    lost_foi_y1 <- no_infs_y1 / total_infs_y1
  }
  if (new_total_infs_y1 > 0) {
    gained_foi_y1 <- no_infs_y2 / new_total_infs_y1
  }
  phis[y1] <- phis[y1] - phis[y1] * lost_foi_y1 + phis[y1] * gained_foi_y1

  lost_foi_y2 <- gained_foi_y2 <- 0
  if (total_infs_y2 > 0) {
    lost_foi_y2 <- no_infs_y2 / total_infs_y2
  }
  if (new_total_infs_y2 > 0) {
    gained_foi_y2 <- no_infs_y1 / new_total_infs_y2
  }
  phis[y2] <- phis[y2] - phis[y2] * lost_foi_y2 + phis[y2] * gained_foi_y2

  ## And finish the swap
  infection_history[samp_indivs, y1] <- infection_history[samp_indivs, y2]
  infection_history[samp_indivs, y2] <- tmp

  return(list(infection_history, phis))
}


#' Brute force infection history proposal
#'
#' Performs a flipping/swapping infection history update step for a matrix of infection histories. 50/50 chance of performing a flip or a swap
#' @param new_inf_hist a matrix of infection histories - rows for individuals, columns for infection epochs. Contents should be 1s and 0s
#' @param sampled_indivs a vector of indices describing rows in the infection history matrix that should be updated
#' @param age_mask a vector (one value for each individual) giving the first infection epoch that an individual could have been exposed in. That is, if an individual was born in the 7th epoch, their entry in age_mask would be 7.
#' @param strain_mask a vector (one value for each individual) giving the last infection epoch that an individual could have been exposed in.
#' @param move_sizes when performing a move step, how far should two epochs be swapped?
#' @param n_infs number of infection epochs to flip
#' @param rand_ns pre-computed random numbers (0-1) for each individual, deciding whether to do a flip or swap
#' @return a matrix of infection histories matching the input new_inf_hist
#' @export
infection_history_symmetric <- function(new_inf_hist, sampled_indivs, age_mask, strain_mask, move_sizes, n_infs, rand_ns) {
  new_inf <- new_inf_hist
  ks <- rpois(length(sampled_indivs), n_infs)
  ## For each individual
  for (i in 1:length(sampled_indivs)) {
    indiv <- sampled_indivs[i]
    ## Isolate infection history
    x <- new_inf_hist[indiv, age_mask[indiv]:strain_mask[indiv]]
    max_i <- length(x)

    ## Flip or swap with prob 50%
    if (rand_ns[i] < 1 / 2) {
      ## Choose a location and turn 0 -> 1 or 1 -> 0
      ## Poisson number of changes?
      k <- min(max_i, max(ks[i], 1))
      locs <- sample(length(x), k)
      x[locs] <- !x[locs]
    } else {
      ## Choose a location
      id1 <- sample(length(x), 1)
      move_max <- move_sizes[indiv]

      ## Choose another location up to move_max epochs away
      move <- sample(-move_max:move_max, 1)
      id2 <- id1 + move

      ## Control boundaries
      if (id2 < 1) id2 <- (move %% max_i) + id1
      if (id2 > max_i) id2 <- (id2 - 1) %% max_i + 1

      ## Swap the contents of these locations
      tmp <- x[id1]
      x[id1] <- x[id2]
      x[id2] <- tmp
    }
    new_inf[indiv, age_mask[indiv]:strain_mask[indiv]] <- x
  }
  return(new_inf)
}



#' Beta binomial infection history update
#'
#' Samples a new infection history from a beta binomial distribution for the specified number of individuals. Note that only one epoch is updated each time with either a flip or a swap step.
#' @param new_inf_hist a matrix of infection histories - rows for individuals, columns for infection epochs. Contents should be 1s and 0s
#' @param sampled_indivs a vector of indices describing rows in the infection history matrix that should be updated
#' @param age_mask a vector (one value for each individual) giving the first infection epoch that an individual could have been exposed in. That is, if an individual was born in the 7th epoch, their entry in age_mask would be 7.
#' @param strain_mask a vector (one value for each individual) giving the last infection epoch that an individual could have been exposed in.
#' @param move_sizes when performing a move step, how far should two epochs be swapped?
#' @param alpha alpha parameter of beta binomial
#' @param beta beta parameter of beta binomial
#' @return a matrix of infection histories matching the input new_inf_hist
#' @export
infection_history_betabinom <- function(new_inf_hist, sampled_indivs, age_mask, strain_mask, move_sizes, alpha, beta) {
  new_inf <- new_inf_hist
  prob_ratio <- rep(1, nrow(new_inf_hist))
  ## For each individual
  for (indiv in sampled_indivs) {
    # x <- new_inf_hist[indiv, age_mask[indiv]:ncol(new_inf_hist)]
    x <- new_inf_hist[indiv, age_mask[indiv]:strain_mask[indiv]]
    max_i <- length(x)
    rand1 <- runif(1)
    ## With prob 0.5 swap or move/add
    if (rand1 < 1 / 2) {
      ## Choose a location
      loc <- sample(length(x), 1)
      x_new <- x_old <- x
      old <- x[loc]
      ## This location can either be a 1 or a 0
      x_new[loc] <- 1
      x_old[loc] <- 0
      ## probA <- dbb(sum(x_new), length(x), alpha, beta)/choose(length(x), sum(x_new))
      ## probB <- dbb(sum(x_old), length(x), alpha, beta)/choose(length(x), sum(x_old))
      ## ratio <- probA/(probA + probB)

      ## Add a 1 some proportion of the time depending on how many 1s are already present
      ## (the sum(x[-loc]) gives the number of 1s in the vector less the location to be
      ## updated
      prob1 <- (alpha + sum(x[-loc])) / (alpha + beta + (length(x) - 1))

      if (runif(1) < prob1) {
        x[loc] <- 1
        ## Otherwise, make a 0
      } else {
        x[loc] <- 0
      }
      if (x[loc] == old) {
        prob_ratio[indiv] <- 1
      } else if (x[loc] == 1) {
        prob_ratio[indiv] <- (beta + length(x) - sum(x[-loc]) - 1) / (alpha + beta + length(x) - 1)
      } else {
        prob_ratio[indiv] <- (alpha + sum(x[-loc])) / (alpha + beta + length(x) - 1)
      }
    } else {
      ## Choose a location
      id1 <- sample(length(x), 1)
      move_max <- move_sizes[indiv]

      ## Propose a location up to move_max epochs away
      move <- sample(-move_max:move_max, 1)
      id2 <- id1 + move

      ## Control boundary conditions
      if (id2 < 1) id2 <- (move %% max_i) + id1
      if (id2 > max_i) id2 <- (id2 - 1) %% max_i + 1

      ## Swap contents
      tmp <- x[id1]
      x[id1] <- x[id2]
      x[id2] <- tmp
    }
    new_inf[indiv, age_mask[indiv]:strain_mask[indiv]] <- x
    # new_inf[indiv,age_mask[indiv]:ncol(new_inf_hist)]=x
  }
  return(list(new_inf, prob_ratio))
}


#' DEPRECATED - implemented in Cpp for speed
#' @export
infection_history_betabinom_group <- function(new_inf_hist, sampled_indivs, age_mask, strain_mask, move_sizes, n_infs, alpha, beta) {
  new_inf <- new_inf_hist
  for (indiv in sampled_indivs) {
    # x <- new_inf[indiv,age_mask[indiv]:ncol(new_inf_hist)]
    x <- new_inf_hist[indiv, age_mask[indiv]:strain_mask[indiv]]

    max_i <- length(x)
    if (runif(1) < 1 / 2) {
      k <- n_infs[indiv]
      locs <- sample(length(x), k)
      number_1s <- sum(x[-locs])
      n <- length(x[-locs])

      for (i in 1:k) {
        ratio <- (alpha + number_1s) / (alpha + beta + n)
        if (runif(1) < ratio) {
          x[locs[i]] <- 1
          number_1s <- number_1s + 1
        } else {
          x[locs[i]] <- 0
        }
        n <- n + 1
      }
    } else {
      for (i in 1:n_infs[indiv]) {
        id1 <- sample(length(x), 1)
        move_max <- move_sizes[indiv]
        move <- sample(-move_max:move_max, 1)
        id2 <- id1 + move
        if (id2 < 1) id2 <- (move %% max_i) + id1
        if (id2 > max_i) id2 <- (id2 - 1) %% max_i + 1

        tmp <- x[id1]
        x[id1] <- x[id2]
        x[id2] <- tmp
      }
    }
    # new_inf[indiv,age_mask[indiv]:ncol(new_inf_hist)]=x
    new_inf[indiv, age_mask[indiv]:strain_mask[indiv]] <- x
  }
  return(new_inf)
}



#' DEPRECATED - propose inf hist from phi
#' @export
inf_hist_prob_phi <- function(new_inf, sampled_indivs, age_mask, strain_mask, n_infs, phis) {
  # ks <- rpois(length(sampled_indivs),n_infs)
  for (i in 1:length(sampled_indivs)) {
    indiv <- sampled_indivs[i]
    x <- new_inf[indiv, age_mask[indiv]:strain_mask[indiv]]
    probs <- phis[age_mask[indiv]:length(phis)]
    max_i <- length(x)
    # k <- min(max_i, max(ks[i],1))
    # locs <- sample(length(x), k)
    # x[locs] <- rbinom(rep(1, k),1,probs[locs])
    x <- rbinom(rep(1, max_i), 1, probs)
    new_inf[indiv, age_mask[indiv]:strain_mask[indiv]] <- x
  }
  return(new_inf)
}

#' DEPRECATED - proposal for phi based on number of infections
#' @export
phi_proposal <- function(current_pars, infection_history, years, js, alpha, beta, n_alive) {
  proposed <- current_pars
  if (length(years) > 1) {
    infs <- colSums(infection_history[, years])
    proposed[js] <- rbeta(length(years), alpha + infs, beta + (n_alive[years] - infs))
  } else {
    infs <- sum(infection_history[, years])
    proposed[js] <- rbeta(1, alpha + infs, beta + (n_alive[years] - infs))
  }
  forward <- sum(dbeta(proposed[js], alpha, beta, log = TRUE))
  back <- sum(dbeta(current_pars[js], alpha, beta, log = TRUE))

  ratio <- back - forward
  return(list(proposed, ratio))
}

#' Infection history proposal - original version
#'
#' Proposes new infection histories for a vector of infection histories, where rows represent individuals and columns represent years. Proposals are either removal, addition or switching of infections.
#' Also requires the indices of sampled individuals, the vector of strain isolation times, and a vector of age masks (ie. which index of the strain_isolation_times vector is the first year in which
#' an individual *could* be infected).
#' NOTE - MIGHT NEED TO UPDATE THIS FOR GROUPS
#' @param new_infection_histories an n*m matrix of 1s & 0s indicating infection histories, where n is individuals and m i strains
#' @param sampled_indivs the indices of sampled individuals to receive proposals
#' @param strain_isolation_times the vector of strain isolation times in real time
#' @param age_mask the vector of indices for each individual specifiying which index of strain_isolation_times is the first strain each individual could have seen
#' @return a new matrix matching new_infection_histories in dimensions with proposed moves
#' @export
infection_history_proposal <- function(new_infection_histories, sampled_indivs, strain_isolation_times, age_mask, strain_mask, n_infs) {
  new_inf <- new_infection_histories
  # ks <- rpois(length(sampled_indivs), n_infs)
  for (indiv in sampled_indivs) { # Resample subset of individuals
    rand1 <- runif(1)
    # x=new_infection_histories[indiv,age_mask[indiv]:length(strain_isolation_times)] # Only resample years individual was alive
    x <- new_infection_histories[indiv, age_mask[indiv]:strain_mask[indiv]] # Only resample years individual was alive

    ## Remove infection
    if (rand1 < 1 / 3) {
      infect_id <- which(x > 0)
      # Number of 1s in first place
      n_1 <- length(infect_id)
      #k_f <- min(n_1, max(n_infs[indiv], 1))
      if (n_1 > 0) {
        # x[sample(infect_id,k_f)]=0 # Why double? DEBUG
        x[sample(c(infect_id), 1)] <- 0
      }
    }
    ## Add infection
    if (rand1 > 1 / 3 & rand1 < 2 / 3) {
      n_infect_id <- which(x == 0)
      n_0 <- length(n_infect_id)
      #k_f <- min(n_0, max(n_infs[indiv], 1))
      if (n_0 > 0) {
        x[sample(c(n_infect_id), 1)] <- 1
        # x[sample(n_infect_id,k_f)]=1
      }
    }
    ## Move infection position
    if (rand1 > 2 / 3) {
      infect_id <- which(x > 0)
      n_infect_id <- which(x == 0)
      n_1 <- length(infect_id)
      n_0 <- length(n_infect_id)
      if (n_1 > 0 & n_0 > 0) {
        x[sample(c(infect_id), 1)] <- 0
        x[sample(c(n_infect_id), 1)] <- 1
        # x[sample(infect_id,1)]=0
        # x[sample(n_infect_id,1)]=1
      }
    }
    new_inf[indiv, age_mask[indiv]:strain_mask[indiv]] <- x # Only =1 if individual was alive
  } # end loop over individuals
  return(new_inf)
}
