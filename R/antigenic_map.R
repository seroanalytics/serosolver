#' Setup antigenic map for serosolver
#'
#' Cleans up an antigenic_map data frame based on provided inputs. Two checks are carried out. First, it will check if an antigenic map is provided, and if so, it will align its entries with possible_exposure_times. If no antigenic map is provided, it will create a dummy map where all pathogens have the same position on the map. Second, it will enumerate the antigenic map for each unique biomarker group, unless the antigenic map has already been enumerated. 
#' @param antigenic_map the antigenic map data frame
#' @param possible_exposure_times a vector of possible exposure times, which will be used to align the antigenic map
#' @param n_biomarker_groups the number of biomarker groups in the antigenic map
#' @param unique_biomarker_groups a vector of unique biomarker groups in the antigenic map
#' @param verbose if TRUE, prints messages about the process
#' @return list with three entries: 1) the updated antigenic map, 2) the updated possible_exposure_times vector, 3) a set of indices matching possible_exposure_times to entries in the antigenic map
#' @export
#' @family antigenic_maps
setup_antigenic_map <- function(antigenic_map=NULL, possible_exposure_times=NULL, n_biomarker_groups=1,unique_biomarker_groups=c(1), verbose=TRUE){
  ## Check if an antigenic map is provided. If not, then create a dummy map where all pathogens have the same position on the map
  if (!is.null(antigenic_map)) {
    possible_exposure_times_tmp <- unique(antigenic_map$inf_times) # How many strains are we testing against and what time did they circulate
    if(!is.null(possible_exposure_times) & !identical(possible_exposure_times, possible_exposure_times_tmp)){
      if(verbose) message(cat("Warning: provided possible_exposure_times argument does not match entries in the antigenic map. Please make sure that there is an entry in the antigenic map for each possible circulation time. Using the antigenic map times.\n"))
    }
    ## If possible exposure times was not specified, use antigenic map times instead
    if(is.null(possible_exposure_times)) {
      infection_history_mat_indices <- match(possible_exposure_times_tmp, possible_exposure_times_tmp)-1
    } else {
      infection_history_mat_indices <- match(possible_exposure_times, possible_exposure_times_tmp)-1
    }
    possible_exposure_times <- possible_exposure_times_tmp
    
    ## If no observation types assumed, set all to 1.
    if (!("biomarker_group" %in% colnames(antigenic_map))) {
      if(verbose) message(cat("Note: no biomarker_group detection in antigenic_map. Setting automatically.\n"))
      antigenic_map_tmp <- replicate(n_biomarker_groups,antigenic_map,simplify=FALSE)
      for(biomarker_group in unique_biomarker_groups){
        antigenic_map_tmp[[biomarker_group]]$biomarker_group <- biomarker_group
      }
      antigenic_map <- do.call(rbind,antigenic_map_tmp)
    }
    
  } else {
    ## Create a dummy map with entries for each observation type
    antigenic_map <- data.frame("x_coord"=1,"y_coord"=1,
                                "inf_times"=rep(possible_exposure_times, n_biomarker_groups), 
                                "biomarker_group"=rep(unique_biomarker_groups,each=length(possible_exposure_times)))
    infection_history_mat_indices <- match(possible_exposure_times, possible_exposure_times)-1
  }
  return(list(antigenic_map=antigenic_map, possible_exposure_times=possible_exposure_times, infection_history_mat_indices=infection_history_mat_indices))
}



#' @export
euc_distance <- function(i1, i2, fit_data) {
  return(sqrt((fit_data[i1, "x_coord"] - fit_data[i2, "x_coord"])^2 + (fit_data[i1, "y_coord"] - fit_data[i2, "y_coord"])^2))
}


#' Create useable antigenic map
#'
#' Creates an antigenic map from an input data frame that can be used to calculate cross reactivity. This will end up being an NxN matrix, where there are N strains circulating.
#' @param anti.map.in can either be a 1D antigenic line to calculate distance from, or a two dimensional matrix with x and y coordinates on an antigenic map
#' @return the euclidean antigenic distance between each pair of viruses in anti.map.in
#' @export
melt_antigenic_coords <- function(anti.map.in) { # anti.map.in can be vector or matrix - rows give inf_times, columns give location
  # Calculate antigenic distances
  if (is.null(dim(anti.map.in))) { # check if input map is one or 2 dimensions
    # If 1D antigenic 'line' defined, calculate distances directory from input
    (dmatrix <- sapply(anti.map.in, function(x) {
      y <- abs(anti.map.in - x)
      y
    }))
  } else { # If 2D antigenic map defined, calculate distances directory from input
    (dmatrix <- apply(
      anti.map.in, 1,
      function(x) {
        y <- sqrt(colSums(apply(anti.map.in, 1, function(y) {
          (y - x)^2
        })))
        y
      }
    ))
  }
}

#' Generate antigenic map, flexible
#'
#' Fits a smoothing spline through a set of antigenic coordinates, and uses this to predict antigenic coordinates for all potential infection time points. This version is more flexible than \code{\link{generate_antigenic_map}}, and allows the user to specify "clusters" to assume that strains circulating in a given period are all identical, rather than on a continuous path through space as a function of time.
#' @param antigenic_distances a data frame of antigenic coordinates, with columns labelled X, Y and Strain for x coord, y coord and Strain label respectively. "Strain" should be a single number giving the year of circulation of that strain. See \code{\link{example_antigenic_map}}
#' @param buckets = 1 the number of epochs per year. 1 means that each year has 1 strain; 12 means that each year has 12 strains (monthly resolution)
#' @param clusters = NULL a data frame of cluster labels, indicating which cluster each circulation year belongs to. Note that each row (year) gets repeated by the number of buckets. Column names "year" and "cluster_used"
#' @param use_clusters = FALSE if TRUE, uses the clusters data frame above, otherwise just returns as normal
#' @param spar = 0.3 to be passed to smooth.spline
#' @param year_min = 1968 first year in the antigenic map (usually 1968)
#' @param year_max = 2016 last year in the antigenic map
#' @return a fitted antigenic map
#' @family antigenic_maps
#' @examples
#' \dontrun{
#' antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
#' antigenic_coords <- read.csv(antigenic_coords_path, stringsAsFactors=FALSE)
#' antigenic_coords$Strain <- c(68,72,75,77,79,87,89,92,95,97,102,104,105,106) + 1900
#' antigenic_map <- generate_antigenic_map_flexible(antigenic_coords, buckets=1, year_min=1968, year_max=2015,spar=0.3)
#' 
#' times <- 1968:2010
#' n_times <- length(times)
#' clusters <- rep(1:5, each=10)
#' clusters <- clusters[1:n_times]
#' clusters <- data.frame(year=times, cluster_used=clusters)
#' antigenic_map <- generate_antigenic_map_flexible(antigenic_coords, buckets=1, 
#'                                                 clusters=clusters,use_clusters=TRUE,
#'                                                 year_min=1968, year_max=2010,spar=0.5)
#' }
#' @seealso \code{\link{generate_antigenic_map}}
#' @export
generate_antigenic_map_flexible <- function(antigenic_distances, buckets = 1, clusters = NULL,
                                            use_clusters = FALSE, spar = 0.3, year_min = 1968, year_max = 2016) {
  ## Convert strains to correct time dimensions
  antigenic_distances$Strain <- antigenic_distances$Strain * buckets
  ## Fit spline through antigenic coordinates
  fit <- smooth.spline(antigenic_distances$X, antigenic_distances$Y, spar = spar)
  
  ## Work out relationship between strain circulation time and x coordinate
  x_line <- lm(data = antigenic_distances, X ~ Strain)
  
  ## Enumerate all strains that could circulate
  Strain <- seq(year_min * buckets, year_max * buckets - 1, by = 1)
  
  ## Predict x and y coordinates for each possible strain from this spline
  x_predict <- predict(x_line, data.frame(Strain))
  y_predict <- predict(fit, x = x_predict)
  
  fit_data <- data.frame(x = y_predict$x, y = y_predict$y)
  fit_data$strain <- Strain
  colnames(fit_data) <- c("x_coord", "y_coord", "inf_times")
  
  ## If using clusters
  if (use_clusters) {
    ## Repeat each row by the number of buckets per year
    clusters <- clusters[rep(seq_len(nrow(clusters)), each = buckets), ]
    
    ## Enumerate out such that each row has a unique time
    clusters$year <- seq(year_min * buckets, length.out = nrow(clusters))
    ## Which time point does each cluster start?
    cluster_starts <- clusters %>% group_by(cluster_used) %>% filter(year == min(year)) 
    cluster_starts <- cluster_starts %>% rename(first_cluster_year=year)
    clusters1 <- merge(clusters, cluster_starts, by = c("cluster_used"))
    clusters1 <- clusters1[, c("cluster_used", "first_cluster_year", "year")]
    colnames(fit_data)[3] <- "first_cluster_year"
    
    ## Merge on "inf_times" with first_cluster_year, such that all viruses
    ## in a cluster have the same location as the first virus in that cluster
    fit_data <- fit_data[fit_data$first_cluster_year %in% clusters1$first_cluster_year, ]
    fit_data <- merge(clusters1, fit_data, by = "first_cluster_year")
    fit_data <- fit_data[, c("x_coord", "y_coord", "year")]
    colnames(fit_data)[3] <- "inf_times"
    fit_data <- fit_data[order(fit_data$inf_times), ]
  }
  
  return(fit_data)
}
