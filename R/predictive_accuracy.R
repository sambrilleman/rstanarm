
#' @export
predictive_accuracy <- function(object, ...) {
  UseMethod("predictive_accuracy")
}

#' @export
#' 
predictive_accuracy.stanjm <- function(object, newdataLong = NULL, newdataEvent = NULL,
                                       t, delta_t, type = "auc", draws = NULL, seed = NULL) {
  validate_stanjm_object(object)
  if (missing(t)) 
    t <- NULL
  if (missing(delta_t))
    delta_t <- NULL
  M <- get_M(object)
  
  if (is.null(t) || is.null(delta_t))
    stop("'t' and 'delta_t' must be specified when calculating the ",
         "predictive accuracy for the event submodel.")
  if (delta_t <= 0)
    stop("'delta_t' must be positive.")
  u <- t + delta_t
  
  # Construct prediction data
  # ndL: dataLong to be used in predictions
  # ndE: dataEvent to be used in predictions
  if (!identical(is.null(newdataLong), is.null(newdataEvent)))
    stop("Both newdataLong and newdataEvent must be supplied together.")
  if (is.null(newdataLong)) { # user did not specify newdata
    dats <- get_model_data(object)
    ndL <- dats[1:M]
    ndE <- dats[["Event"]]
  } else { # user specified newdata
    newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
    ndL <- newdatas[1:M]
    ndE <- newdatas[["Event"]]   
  }
  
  # Obtain event time and status variable, from event submodel formula
  fm_LHS <- formula(object, m = "Event")[[2L]]
  len <- length(fm_LHS)
  event_tvar <- as.character(fm_LHS[[len - 1L]])
  event_dvar <- as.character(fm_LHS[[len]])
  id_var <- object$id_var
  
  # Subset prediction data to only include individuals surviving up
  # to time t and their longitudinal data observed prior to time t
  subdats <- get_conditional_prediction_data(object, ndL, ndE, t = t)
  ndL <- subdats[1:M]
  ndE <- subdats[["Event"]]
  
  # Observed y: event status at time u
  y <- ndE[, c(id_var, event_tvar, event_dvar), drop = FALSE]
  
  # Predicted y: conditional survival probability at time u
  ytilde <- posterior_survfit(
    object, 
    newdataLong = ndL, 
    newdataEvent = ndE,
    times = u,
    last_time = t,
    condition = TRUE,
    extrapolate = FALSE,
    draws = draws,
    seed = seed)
  ytilde <- ytilde[, c(id_var, "survpred"), drop = FALSE]
  
  y <- merge(y, ytilde, by = id_var)
  
  ids <- y[[id_var]]
  if (any(duplicated(ids)))
    stop2("Bug found: y should only contain one row per individual.")
  eventtime <- y[[event_tvar]]
  status    <- y[[event_dvar]]
  survpred  <- y[["survpred"]]
  names(ids) <- names(eventtime) <- names(status) <- names(survpred) <- ids
  
  pairs <- combn(as.character(ids), 2)
  ids_i <- pairs[1L,]
  ids_j <- pairs[2L,]
  times_i    <- eventtime[ids_i]
  times_j    <- eventtime[ids_j]
  status_i   <- status[ids_i]
  status_j   <- status[ids_j]
  survpred_i <- survpred[ids_i]
  survpred_j <- survpred[ids_j]
  
  ind1 <- (times_i <= u & status_i == 1) & (times_j > u)
  ind2 <- (times_i <= u & status_i == 0) & (times_j > u)
  ind3 <- (times_i <= u & status_i == 1) & (times_j <= u & status_j == 0)
  ind4 <- (times_i <= u & status_i == 0) & (times_j <= u & status_j == 0)
  names(ind1) <- names(ind2) <- names(ind3) <- names(ind4) <- 
    paste(ids_i, ids_j, sep = "_")
  ind <- ind1 | ind2 | ind3 | ind4
  
  if (any(ind2)) {
    nams <- strsplit(names(ind2[ind2]), "_")
    nams_i <- sapply(nams, `[`, 1L)
    unq_nams_i <- unique(nams_i)
    ndE_2i <- ndE[ndE[[id_var]] %in% unq_nams_i, , drop = FALSE]
    ndL_2i <- lapply(ndL, function(x) x[x[[id_var]] %in% unq_nams_i, , drop = FALSE])
    last_time_2i <- ndE_2i[[event_tvar]]
    names(last_time_2i) <- ndE_2i[[id_var]]
    survpred_2i <- posterior_survfit(object, newdataLong = ndL_2i, newdataEvent = ndE_2i, 
                                     times = u, last_time = last_time_2i, 
                                     draws = draws, seed = seed)
    ind[ind2] <- ind[ind2] * (1 - survpred_2i[nams_i, "survpred"])
  }
  
  if (any(ind3)) {
    nams <- strsplit(names(ind3[ind3]), "_")
    nams_j <- sapply(nams, "[", 2)
    unq_nams_j <- unique(nams_j)
    ndE_3j <- ndE[ndE[[id_var]] %in% unq_nams_j, , drop = FALSE]
    ndL_3j <- lapply(ndL, function(x) x[x[[id_var]] %in% unq_nams_j, , drop = FALSE])
    last_time_3j <- ndE_3j[[event_tvar]]
    names(last_time_3j) <- ndE_3j[[id_var]]
    survpred_3j <- posterior_survfit(object, newdataLong = ndL_3j, newdataEvent = ndE_3j, 
                                     times = u, last_time = last_time_3j,
                                     draws = draws, seed = seed)
    ind[ind3] <- ind[ind3] * survpred_3j[nams_j, "survpred"]
  }
  
  if (any(ind4)) {
    nams <- strsplit(names(ind4[ind4]), "_")
    nams_i <- sapply(nams, "[", 1)
    nams_j <- sapply(nams, "[", 2)
    unq_nams_i <- unique(nams_i)
    unq_nams_j <- unique(nams_j)
    ndE_4i <- ndE[ndE[[id_var]] %in% unq_nams_i, , drop = FALSE]
    ndE_4j <- ndE[ndE[[id_var]] %in% unq_nams_j, , drop = FALSE]
    ndL_4i <- lapply(ndL, function(x) x[x[[id_var]] %in% unq_nams_i, , drop = FALSE])
    ndL_4j <- lapply(ndL, function(x) x[x[[id_var]] %in% unq_nams_j, , drop = FALSE])
    last_time_4i <- ndE_4i[[event_tvar]]
    last_time_4j <- ndE_4j[[event_tvar]]
    names(last_time_4i) <- ndE_4i[[id_var]]
    names(last_time_4j) <- ndE_4j[[id_var]]
    survpred_4i <- posterior_survfit(object, newdataLong = ndL_4i, newdataEvent = ndE_4i, 
                                     times = u, last_time = ndE_4i[[event_tvar]], 
                                     draws = draws, seed = seed)
    survpred_4j <- posterior_survfit(object, newdataLong = ndL_4j, newdataEvent = ndE_4j, 
                                     times = u, last_time = ndE_4j[[event_tvar]], 
                                     draws = draws, seed = seed)
    ind[ind4] <- ind[ind4] * (1 - survpred_4i[nams_i, "survpred"]) * (survpred_4j[nams_j, "survpred"])
  }
  
  auc <- sum((survpred_i < survpred_j) * c(ind), na.rm = TRUE) / sum(ind, na.rm = TRUE) 
  auc
}
  


# Subset the prediction data to only include:
# - patients who survived up to a specified time t
# - their longitudinal observations that were observed prior to time t
#
# @param object A stanjm object.
# @param dataLong A list of data frames, one for each of the longitudinal submodels.
# @param dataEvent The data frame for the event submodel.
# @param t Scalar specifying the time patients must have survived up to, as well as
#   being the time up to which longitudinal measurements should be used. Patients
#   who died or were censored before time t are discarded. Longitudinal measurements
#   observed in the data after time t are also discarded.
get_conditional_prediction_data <- function(object, dataLong, dataEvent, t) {
  if (!is.stanjm(object))
    stop2("This function is only for stanjm objects.")
  if (!is(dataLong, "list"))
    dataLong <- list(dataLong)
  fm_lhs <- formula(object, m = "Event")[[2L]]
  y_time_var <- object$time_var
  e_time_var <- as.character(fm_lhs[[length(fm_lhs) - 1L]])
  id_var <- object$id_var
  
  # patients who survived up to time t
  sel <- which(dataEvent[[e_time_var]] > t)
  dataEvent <- dataEvent[sel, , drop = FALSE]
  ids <- dataEvent[[id_var]]
  
  # patients with longitudinal observations before time t
  dataLong <- lapply(dataLong, function(x) {
    sel <- which(x[[y_time_var]] <= t)
    x[sel, , drop = FALSE]
  })  
  
  for (i in 1:length(dataLong))
    ids <- intersect(dataLong[[i]][[id_var]], ids)
  if (!length(ids))
    stop2("No individuals still at risk at time 't' and ",
          "with longitudinal measurements prior to 't'.")
  
  dataEvent <- dataEvent[dataEvent[[id_var]] %in% ids, , drop = FALSE]
  dataLong <- lapply(dataLong, function(x) {
    x[x[[id_var]] %in% ids, , drop = FALSE]
  })
  ret <- c(dataLong, list(dataEvent))
  list_nms(ret, M = get_M(object), stub = get_stub(object))
}

