#' Swap infection history years
#'
#' Swaps the entire contents of two columns of the infection history matrix, adhering to age and strain limitations.
#' @param infection_history matrix of 1s and 0s to swap around representing the infection history
#' @param age_mask the first index in infection_history that each individual (row) could be infected in
#' @param strain_mask the last index in infection_history that each individual (row) could be infected in ie. the time of the latest blood sample
#' @param swap_propn what proportion of infections should be swapped?
#' @param move_size How many time points away should be chosen as candidate swaps?
#' @param proposal_ratios optional NULL. Can set the relative sampling weights of the infection state times. Should be an integer vector of length matching nrow(antigenic_map). Otherwise, leave as NULL for uniform sampling.
#' @return the same infection_history matrix, but with two columns swapped
#' @family proposals
#' @examples
#' data(example_inf_hist)
#' data(example_titre_dat)
#' data(example_antigenic_map)
#' times <- example_antigenic_map$inf_times
#' ages <- unique(example_titre_dat[,c("DOB","individual")])
#' age_mask <- create_age_mask(ages$DOB, times)
#' strain_mask <- create_strain_mask(example_titre_dat, times)
#' new_inf_hist <- inf_hist_swap(example_inf_hist, age_mask,strain_mask, 1,3)[[1]]
#' @export
inf_hist_swap <- function(infection_history, vaccination_histories_mat, age_mask, strain_mask, swap_propn, move_size, proposal_ratios=NULL, output_list = TRUE) {
    use_ratios <- NULL
    if(!is.null(proposal_ratios)) {
        use_ratios <- proposal_ratios / sum(proposal_ratios)
    }
        
    ## Choose a column
    y1 <- sample(1:ncol(infection_history), 1, prob = use_ratios) # this is a time point j

    ## Propose another column some random distance, but not 0, away
    move <- 0
    while (move == 0) move <- sample((-move_size):move_size, 1) # propose a new point move away

    ## Need to adjust if we've proposed too far away
    y2 <- y1 + move
    if(y2 < 1) y2 = -y2 + 2
    if(y2 > ncol(infection_history)) y2 = ncol(infection_history) - y2 + ncol(infection_history) # adjustements
    
    ##while (y2 < 1) y2 <- y2 + ncol(infection_history)
    ##while(y2 > ncol(infection_history)) y2 <- y2 - ncol(infection_history)

    ## Get the first and last year chronologically
    small_year <- min(y1, y2)
    big_year <- max(y1, y2)

    ## Find individuals that are alive/sampled in both years and choose the lesser of swap_propn*n_indivs and
    ## the number that are actually able to be infected in both years
    indivs <- 1:nrow(infection_history)
    alive_indivs <- indivs[intersect(which(age_mask <= small_year), which(strain_mask >= big_year))]
    # if either y1 or y2 for a person is a vaccination year then need to stop the swap happening, of those alive individuals
    if (is.null(vaccination_histories_mat)) {
      vac_times <- NULL
    } else {
      vac_times <- union(which(vaccination_histories_mat[, y1] > 0), which(vaccination_histories_mat[, y2] > 0))
    }
    alive_indivs <- setdiff(alive_indivs, vac_times)
    samp_indivs <- sample(alive_indivs, floor(length(alive_indivs) * swap_propn))

    ## Swap contents
    tmp <- infection_history[samp_indivs, y1]
    infection_history[samp_indivs, y1] <- infection_history[samp_indivs, y2]
    infection_history[samp_indivs, y2] <- tmp

    if(output_list)
      return(list(infection_history))
    else 
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
#' @family proposals
#' @examples
#' data(example_inf_hist)
#' data(example_titre_dat)
#' data(example_antigenic_map)
#' times <- example_antigenic_map$inf_times
#' ages <- unique(example_titre_dat[,c("DOB","individual")])
#' age_mask <- create_age_mask(ages$DOB, times)
#' strain_mask <- create_strain_mask(example_titre_dat, times)
#' phis <- runif(length(times))
#' n_alive <- get_n_alive(example_titre_dat,times)
#' new_inf_hist <- inf_hist_swap_phi(example_inf_hist, phis, age_mask,strain_mask, 1,3, n_alive)
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
#' @param swap_propn threshold for deciding if swap or add/remove step
#' @return a matrix of infection histories matching the input new_inf_hist
#' @family proposals
#' @examples
#' data(example_inf_hist)
#' data(example_titre_dat)
#' data(example_antigenic_map)
#'
#' times <- example_antigenic_map$inf_times
#' ages <- unique(example_titre_dat[,c("DOB","individual")])
#'
#' ## Create age and strain mask
#' age_mask <- create_age_mask(ages$DOB, times)
#' strain_mask <- create_strain_mask(example_titre_dat, times)
#'
#' ## Index of individuals to be resampled
#' indivs <- 1:nrow(example_inf_hist)
#' n_indiv <- length(indivs)
#' 
#' ## Parameters controlling proposal sizes for each individual
#' move_sizes <- rep(3, n_indiv)
#' n_infs <- rep(10, n_indiv)
#' 
#' ## Pre-compute random numbers
#' rand_ns <- runif(n_indiv)
#'
#' new_inf_hist <- infection_history_symmetric(example_inf_hist, indivs,age_mask ,strain_mask, move_sizes, n_infs, rand_ns, 0.5)
#' @export
infection_history_symmetric <- function(new_inf_hist, sampled_indivs, age_mask, strain_mask, move_sizes, n_infs, rand_ns, swap_propn = 0.5) {
  new_inf <- new_inf_hist
  ks <- rpois(length(sampled_indivs), n_infs)
  # ks <- n_infs
  ## For each individual
  for (i in 1:length(sampled_indivs)) {
    indiv <- sampled_indivs[i]
    ## Isolate infection history
    x <- new_inf_hist[indiv, age_mask[indiv]:strain_mask[indiv]]
    max_i <- length(x)

    ## Flip or swap with prob 50%
    if (rand_ns[i] > swap_propn) {
      ## Choose a location and turn 0 -> 1 or 1 -> 0
      ## Poisson number of changes?
      k <- min(max_i, max(ks[i], 1))
      locs <- sample(length(x), k)
      x[locs] <- !x[locs]
    } else {
      ## Choose a location
      infs <- which(x == 1)
      if (length(infs > 0)) {
        id1 <- sample(infs, 1)
        # id1 <- sample(length(x),1)
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
    }
    new_inf[indiv, age_mask[indiv]:strain_mask[indiv]] <- x
  }
  return(new_inf)
}
