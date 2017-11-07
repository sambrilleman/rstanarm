# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2017 Sam Brilleman
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Predictive error and accuracy measures for an event outcome
#' 
#' Calculate predictive error and accuracy measures (discrimination, calibration)
#' for the event outcome (i.e. the event submodel) based on a fitted joint model 
#' for longitudinal and time-to-event data estimated using \code{\link{stan_jm}}.
#' 
#' @export
#' @templateVar stanjmArg object
#' @template args-stanjm-object
#' @param ... Ignored for the \code{stanjm} method.
#' 
#' @details
#' The predictive error or accuracy measure is calculated only using individuals
#' who are still at risk of the event at time \code{t} (i.e. they have not 
#' experienced the event, or been censored, prior to time {t}). If there are 
#' individuals in the data who have experienced the event, or been censored, 
#' prior to time {t}, then these individuals are discarded before calculating the 
#' predictive error or accuracy measure. Moreover, only longitudinal data
#' observed between time 0 (i.e. baseline) and time \code{t} is used in the
#' calculation. That is, any longitudinal data observed after time {t} is also
#' discarded before calculating the predictive error or accuracy measure.
#' 
#' The \code{predictive_accuracy} function can be used to calculate the 
#' following types of measures:
#' 
#' \subsection{Time-dependent AUC measure (\code{type = "auc"})}{
#' To be completed. Describe the time-dependent area under the receiver 
#' operating characteristic curve (AUC) measure.
#' }
#' \subsection{Prediction error (\code{type = "error"})}{
#' To be completed. Describe the prediction error calculation for event submodel.
#' }
#' 
#' @seealso \code{\link{posterior_survfit}} 
#' 
#' @examples
#' \donttest{
#'   # Run example model if not already loaded
#'   if (!exists("example_jm")) example(example_jm)
#'   
#'   # Time-dependent AUC
#'   predictive_accuracy(example_jm, t = 5, u = 8)
#'   
#'   # Prediction error
#'   predictive_accuracy(example_jm, t = 5, u = 8, type = "error")
#' }
#' 
predictive_accuracy <- function(object, ...) {
  UseMethod("predictive_accuracy")
}

#' @rdname predictive_accuracy
#' @export
#' @param newdataLong,newdataEvent Optionally, a data frame (or in the case of 
#'   \code{newdataLong} this can be a list of data frames) in which to look 
#'   for variables with which to predict. If omitted, the model matrices are used. 
#'   If new data is provided, then it should also contain the longitudinal 
#'   outcome data on which to condition when drawing the new group-specific 
#'   coefficients for individuals in the new data. Note that there is only
#'   allowed to be one row of data for each individual in \code{newdataEvent}, 
#'   that is, time-varying covariates are not allowed in the prediction data for
#'   the event submodel.
#' @param t,u The argument \code{t} specifies the time up to which individuals must 
#'   have survived as well as being the time up to which the longitudinal data 
#'   in \code{newdata} is available. The argument \code{u} specifies the time 
#'   horizon up to which the prediction error or accuracy measure should be calculated.
#' @param type The type of predictive error or accuracy measure to calculate. 
#'   Can currently be any of the following:
#'   \describe{
#'     \item{\code{"auc"}}{a time-dependent AUC measure evaluating the 
#'     discriminatory ability of the model between times \code{t} and
#'     \code{u}.}
#'     \item{\code{"error"}}{the prediction error for the event outcome evaluated
#'     between times \code{t} and \code{u} and allowing for censoring.}
#'   }
#' @param loss_function The loss function to use when \code{type = "error"}.
#'   Can be \code{"square"} (the default), \code{"absolute"}, or a user-defined
#'   function.
#' @param draws Currently ignored.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#'
predictive_accuracy.stanjm <- function(object, newdataLong = NULL, newdataEvent = NULL,
                                       t, u, type = c("auc", "error"), 
                                       loss_function = "square",
                                       draws = NULL, seed = NULL, ...) {
  
  validate_stanjm_object(object)
  if (missing(t)) 
    t <- NULL
  if (missing(u))
    u <- NULL
  M <- get_M(object)
  type <- match.arg(type)
  
  if (is.null(t) || is.null(u))
    stop("'t' and 'u' must be specified when calculating the ",
         "predictive error for the event submodel.")
  if (u <= t)
    stop("'u' must be greater than 't'.")
  
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
    last_time2 = event_tvar,
    condition = TRUE,
    extrapolate = FALSE,
    draws = draws,
    seed = seed)
  ytilde <- ytilde[, c(id_var, "survpred", "survpred_eventtime"), drop = FALSE]
  
  # Merge observed event status and predicted surv prob
  y <- merge(y, ytilde, by = id_var)
  
  if (type == "auc") {
    
    # Extract necessary components as named vectors
    ids <- y[[id_var]]
    if (any(duplicated(ids)))
      stop2("Bug found: y should only contain one row per individual.")
    eventtime <- y[[event_tvar]]
    status    <- y[[event_dvar]]
    survpred  <- y[["survpred"]]           # last_time = t
    survpred2 <- y[["survpred_eventtime"]] # last_time = known censoring time
    names(ids) <- names(eventtime) <- names(status) <- 
      names(survpred) <- names(survpred2) <- ids
    
    # Ensure eventtime_i <= eventtime_j (for subjects i and j)
    ord <- order(eventtime)
    ids       <- ids[ord]
    eventtime <- eventtime[ord]
    status    <- status[ord]
    survpred  <- survpred[ord]
    survpred2 <- survpred2[ord]
    
    # Extract necessary components for each pair of subject (i,j)
    pairs <- combn(as.character(ids), 2)
    ids_i <- pairs[1L,]
    ids_j <- pairs[2L,]
    times_i    <- eventtime[ids_i]
    times_j    <- eventtime[ids_j]
    status_i   <- status[ids_i]
    status_j   <- status[ids_j]
    survpred_i <- survpred[ids_i]
    survpred_j <- survpred[ids_j]
    
    # Calculate AUC
    ind1 <- (times_i <= u & status_i == 1) & (times_j > u)
    ind2 <- (times_i <= u & status_i == 0) & (times_j > u)
    ind3 <- (times_i <= u & status_i == 1) & (times_j <= u & status_j == 0)
    ind4 <- (times_i <= u & status_i == 0) & (times_j <= u & status_j == 0)
    nms_ind <- paste0(ids_i, "_", ids_j)
    names(ind1) <- names(ind2) <- names(ind3) <- names(ind4) <- nms_ind
    ind <- ind1 | ind2 | ind3 | ind4
    
    if (any(ind2)) {
      # if subject i was censored, then weight the contribution to AUC
      # based on the (failure) probability that subject i would have
      # died in interval between their known censoring time and time u
      nms <- strsplit(names(ind2[ind2]), "_")
      nms_i <- sapply(nms, `[`, 1L)
      ind[ind2] <- ind[ind2] * (1 - survpred2[nms_i])
    }
    
    if (any(ind3)) {
      # if subject j was censored, then weight the contribution to AUC
      # based on the (survival) probability that subject j would have
      # survived from their known censoring time up until time u
      nms <- strsplit(names(ind3[ind3]), "_")
      nms_j <- sapply(nms, "[", 2)
      ind[ind3] <- ind[ind3] * survpred2[nms_j]
    }
    
    if (any(ind4)) {
      # if subjects i and j were both censored, then weight the contribution 
      # to AUC based on the (failure) probability for subject i and the
      # (survival) probability for subject j, taken between each of their 
      # known censoring times and the time horizon u
      nms <- strsplit(names(ind4[ind4]), "_")
      nms_i <- sapply(nms, "[", 1)
      nms_j <- sapply(nms, "[", 2)
      ind[ind4] <- ind[ind4] * (1 - survpred2[nms_i]) * (survpred2[nms_j])
    }
    
    # use same method for ties as survival pkg: see ?survival::survConcordance
    agree    <- sum((survpred_i < survpred_j) * c(ind), na.rm = TRUE)
    disagree <- sum((survpred_i > survpred_j) * c(ind), na.rm = TRUE)
    tied     <- sum((survpred_i == survpred_j) * c(ind), na.rm = TRUE)
    total <- agree + disagree + tied
    if (sum(ind, na.rm = TRUE) != total)
      stop2("Bug found: total and sum(ind) should be equal (= number of comparable pairs).")
    val <- (agree + (0.5 * tied)) / total 
    
  } else if (type == "error") {
    
    loss <- switch(loss_function,
                   square = function(x) {x*x},
                   absolute = function(x) {abs(x)})
    
    # Calculate mean prediction error
    y$status <- as.integer(y[[event_dvar]])
    y$ind    <- as.integer(y[[event_tvar]] > u)
    y$ind1 <- y$ind                        # died/censored after u
    y$ind2 <- y$status * (1 - y$ind)       # died before u
    y$ind3 <- (1 - y$status) * (1 - y$ind) # censored before u
    y$val <- 
      y$ind1 * loss(1 - y$survpred) +
      y$ind2 * loss(0 - y$survpred) +
      y$ind3 * (y$survpred_eventtime * loss(1- y$survpred) + 
                  (1 - y$survpred_eventtime) * loss(0- y$survpred))
    val <- mean(y$val)
  }
  
  ret <- structure(val, 
                   n_subjects_at_risk = nrow(y),
                   t = t, u = u, type = type,
                   loss_function = if (type == "error") loss_function,
                   stanreg_name = deparse(substitute(object)),
                   class = "predictive_accuracy.stanjm")
  return(ret)
}
  

# ------------------ exported but doc kept internal

#' Generic print method for \code{predictive_accuracy.stanjm} objects
#' 
#' @rdname print.predictive_accuracy.stanjm
#' @method print predictive_accuracy.stanjm
#' @keywords internal
#' @export
#' @param x An object of class \code{predictive_accuracy.stanjm}, returned by a call to 
#'   \code{\link[predictive_accuracy.stanjm]{predictive_accuracy}}.
#' @param digits Number of digits to use for formatting.
#' @param ... Ignored.
#' 
print.predictive_accuracy.stanjm <- function(x, digits = 4, ...) {
  val <- format(round(x, digits), nsmall = digits)
  type <- attr(x, "type")
  stanreg_nm <- attr(x, "stanreg_name")
  if (type == "auc") {
    msg1 <- paste0("Event prediction accuracy (discrimination) for model '", stanreg_nm, "'")
    msg2 <- paste0("\nEstimated time-dependent AUC:           ", val)
  } else if (type == "error") {
    msg1 <- paste0("Event prediction error for model '", stanreg_nm, "'")
    msg2 <- paste0("\nEstimated prediction error:             ", val)
  } else {
    stop2("Bug found: unknown type of predictive_accuracy measure.")
  }
  cat(msg1)
  cat("\n------")  
  cat(msg2)
  cat("\nPrediction is calculated at time (t):  ", attr(x, "u"))
  cat("\nUsing longitudinal data up to time (u):", attr(x, "t"))
  cat("\nNum. subjects still at risk at time t: ", attr(x, "n_subjects_at_risk"))
  if (type == "error")
    cat("\nLoss function used:                    ", attr(x, "loss_function"))
  cat("\n------")
  invisible(x)
}

# ------------------ internal

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

