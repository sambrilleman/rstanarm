# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2016, 2017 Sam Brilleman
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

# Return the design matrices required for evaluating the linear predictor or
# log-likelihood in post-estimation functions for a \code{stan_jm} model
#
# @param object A stanmvreg object
# @param newdataLong A data frame or list of data frames with the new 
#   covariate data for the longitudinal submodel
# @param newdataEvent A data frame with the new covariate data for the
#   event submodel
# @param ids An optional vector of subject IDs specifying which individuals
#   should be included in the returned design matrices.
# @param etimes An optional vector of times at which the event submodel
#   design matrices should be evaluated (also used to determine the 
#   quadrature times). If NULL then times are taken to be the eventimes in
#   the fitted object (if newdataEvent is NULL) or in newdataEvent.
# @param long_parts,event_parts A logical specifying whether to return the
#   design matrices for the longitudinal and/or event submodels.
# @return A named list (with components M, Npat, ndL, ndE, yX, tZt, 
#   yZnames, eXq, assoc_parts) 
.pp_data_jm <- function(object, newdataLong = NULL, newdataEvent = NULL, 
                        ids = NULL, etimes = NULL, long_parts = TRUE, 
                        event_parts = TRUE) {
  M <- get_M(object)
  id_var   <- object$id_var
  time_var <- object$time_var
  
  # validate newdata - return value may be NULL
  newdatas <- .validate_newdata_for_jm(object, newdataLong, newdataEvent)
  
  # drop variables from pp_data that are not required
  datas <- .clean_pp_data_jm(object, newdatas = newdatas)
  
  # prediction data
  ndL <- datas[1:M]        # glmer (longitudinal) submodels
  ndE <- datas[["Event"]]  # event submodel
  
  # possibly subset
  if (!is.null(ids)) {
    ndL <- subset_ids(object, ndL, ids)
    ndE <- subset_ids(object, ndE, ids)
  }
  id_list <- unique(ndE[[id_var]]) # unique subject id list
  
  # evaluate the last known survival time and status
  if (!is.null(newdataEvent) && is.null(etimes)) {
    # prediction data for the event submodel was provided but  
    # no event times were explicitly specified by the user, so
    # they must be evaluated using the data frame
    surv <- eval(formula(object, m = "Event")[[2L]], ndE)
    etimes  <- unclass(surv)[,"time"]
    estatus <- unclass(surv)[,"status"]    
  } else if (is.null(etimes)) {
    # if no prediction data was provided then event times are 
    # taken from the fitted model
    etimes  <- object$eventtime[as.character(id_list)]
    estatus <- object$status[as.character(id_list)]
  } else { 
    # otherwise, event times ('etimes') are only directly specified for dynamic   
    # predictions via posterior_survfit in which case the 'etimes' correspond 
    # to the last known survival time and therefore we assume everyone has survived
    # up to that point (ie, set estatus = 0 for all individuals), this is true 
    # even if there is an event indicated in the data supplied by the user.
    estatus <- rep(0, length(etimes))
  }
  res <- nlist(M, Npat = length(id_list), ndL, ndE)
  
  if (long_parts && event_parts)
    lapply(ndL, function(x) {
      if (!time_var %in% colnames(x)) 
        STOP_no_var(time_var)
      if (!id_var %in% colnames(x)) 
        STOP_no_var(id_var)
      if (any(x[[time_var]] < 0))
        stop2("Values for the time variable (", time_var, ") should not be negative.")
      mt <- tapply(x[[time_var]], factor(x[[id_var]]), max)
      if (any(mt > etimes))
        stop2("There appears to be observation times in the longitudinal data that ",
              "are later than the event time specified in the 'etimes' argument.")      
    }) 
  
  # response and design matrices for longitudinal submodels
  if (long_parts) {
    y <- lapply(1:M, function(m) eval(formula(object, m = m)[[2L]], ndL[[m]]))
    y_pp_data <- lapply(1:M, function(m) pp_data(object, ndL[[m]], m = m))
    yX  <- fetch(y_pp_data, "x")
    yZt <- fetch(y_pp_data, "Zt")
    yZ_names <- fetch(y_pp_data, "Z_names")
    flist <- lapply(ndL, function(x) factor(x[[id_var]]))
    res <- c(res, nlist(y, yX, yZt, yZ_names, flist))
  }
  
  # design matrices for event submodel and association structure
  if (event_parts) {
    qnodes <- object$qnodes
    qq     <- get_quadpoints(qnodes)
    qtimes <- uapply(qq$points,  unstandardise_qpts, 0, etimes)
    qwts   <- uapply(qq$weights, unstandardise_qwts, 0, etimes)
    starttime <- deparse(formula(object, m = "Event")[[2L]][[2L]], 500L)
    e_dat <- prepare_data_table(ndE, id_var = id_var, time_var = starttime)
    times <- c(etimes, qtimes) # times used to design event submodel matrices
    ids2 <- rep(id_list, qnodes + 1)
    e_pp_data <- rolling_merge(e_dat, ids = ids2, times = times)
    eXq <- .pp_data_mer_x(object, newdata = e_pp_data, m = "Event")
    assoc_parts <- lapply(1:M, function(m) {
      ymf <- ndL[[m]]
      grp_stuff <- object$grp_stuff[[m]]
      if (grp_stuff$has_grp) {
        grp_stuff <- get_extra_grp_info( # update grp_info with new data
          grp_stuff, flist = ymf, id_var = id_var, 
          qnodes = qnodes, grp_assoc = object$grp_assoc)
      }
      ymf <- prepare_data_table(ymf, id_var = id_var, time_var = time_var,
                                grp_var = grp_stuff$grp_var) # NB grp_var may be NULL
      make_assoc_parts(
        ymf, assoc = object$assoc[,m], id_var = id_var, time_var = time_var, 
        ids = ids2, times = times, grp_stuff = grp_stuff,
        use_function = pp_data, object = object, m = m)
    })
    assoc_attr <- nlist(.Data = assoc_parts, qnodes, qtimes, qwts, etimes, estatus)
    assoc_parts <- do.call("structure", assoc_attr)
    res <- c(res, nlist(eXq, assoc_parts))
  }
  
  return(res)
}

# Return a data frame for each submodel that:
# (1) only includes variables used in the model formula.
# (2) ensures that additional variables that are required
#     such as the ID variable or variables used in the 
#     interaction-type association structures, are included.
# (3) if new data is not provided by the user, then the cleaned
#     pp_data should only include rows that were contained in the 
#     glmod/coxmod model frames (ie. after any subsetting and/or
#     rows with NAs the model variables were removed).
#
# It is necessary to drop unneeded variables though so that 
# errors are not encountered if the original data contained 
# NA values for variables unrelated to the model formula.
# We generate a data frame here for in-sample predictions 
# rather than using a model frame, since some quantities will
# need to be recalculated at quadrature points etc, for example
# in posterior_survfit.
#
# @param object A stanmvreg object.
# @param newdatas A data frame, or a list of data frames. If NULL
#   then the user did not pass any new data to the prediction function.
#   If the user did provide new data then this will have been passed
#   through the .validate_newdata_for_jm function.
# @param m Integer specifying which submodel to get the
#   prediction data frame for.
# @param keep_response Whether to include the response variable(s)
#   in the returned data frame.
# @return A data frame or list of data frames with all the
#   (unevaluated) variables required for predictions.
.clean_jm_pp_data <- function(object, newdatas = NULL, m = NULL,
                              keep_response = TRUE,
                              clean_newdata = FALSE) {
  validate_stanmvreg_object(object)
  M <- get_M(object)
  terms <- terms(object, fixed.only = FALSE)
  
  if (is.jm(object) && is(newdatas, "list")) {
    # if the new data is a list, then it is being used to construct 
    # pp_data for all submodels as required by posterior_survfit so we
    # must identify any extra variables that are related to the 
    # association structure for the joint model and therefore need 
    # to be added to each submodels terms objects
    extra_vars <- lapply(1:M, function(m) {
      # for each submodel loop over the four possible assoc  
      # interaction formulas and collect any variables used
      forms_m <- object$assoc["which_formulas",][[m]]
      uapply(forms_m, function(x) {
        if (length(x)) {
          rownames(attr(terms.formula(x), "factors")) 
        } else NULL
      })
    })
    # also ensure that id_var is included in the event data
    extra_vars$Event <- object$id_var
    
    if (!identical(length(terms), length(extra_vars)))
      stop2("Bug found: terms and extra_vars should be same length.")
    
    # add the extra variables to the terms formula for each submodel
    terms <- xapply(terms, extra_vars, FUN = function(x, y) {
      lhs <- x[[2L]]
      rhs <- deparse(x[[3L]], 500L)
      if (!is.null(y))
        rhs <- c(rhs, y)
      reformulate(rhs, response = if (keep_response) lhs else NULL)
    })
  }
  
  # collect the data frames
  if (is.null(newdatas)) {
    # extract data from model object    
    datas <- if (is.jm(object)) 
      c(object$dataLong, list(object$dataEvent)) else object$data
    # identify rows that were in the model frame
    row_nms <- lapply(model.frame(object), rownames)
    # only keep rows that were in the model frame, and 
    # the variables that were used in the model formulas
    mfs <- xapply(w = terms, x = datas, y = row_nms,
                  FUN = function(w, x, y) 
                    get_all_vars(w, x)[y, , drop = FALSE])
    na_check <- uapply(mfs, is.na)
    if (any(na_check))
      stop2("Bug found: model data should not have contained NAs?")
  } else if (!clean_newdata) {
    return(newdata)
  } else if (is(newdatas, "list")) {
    # newdatas is a list and therefore is being used to construct pp_data
    # for all submodels as required by posterior_survfit, so validate
    # all new data frames against the model variables for each submodel
    data <- newdatas
    # only keep the variables that were used in the model formulas
    mfs <- xapply(terms, datas, FUN = get_all_vars(w, x))
    # ensure there are no NAs in the prediction data (which now
    # only contains the model variables required for prediction)
    na_check <- uapply(mfs, is.na)
    if (any(na_check))
      stop2("NAs are not allowed in the prediction data.")
  } else {
    # newdatas is not a list and therefore is a data frame for just one submodel 
    # as required by posterior_predict or posterior_traj so we only need to 
    # validate the new data against the model variables for one submodel (m)
    if (is.null(m))
      stop2("Bug found: new data is not a list, but m is NULL, so not sure ",
            "which submodel the new data corresponds to.")
    mf <- get_all_vars(terms[[m]], datas)
    na_check <- is.na(mf)
    if (any(na_check))
      stop2("NAs are not allowed in the prediction data.")
    return(mf) # one user supplied data frame returned for posterior_{predict,traj}
  }
  
  mfs <- list_nms(mfs, M, stub = get_stub(object))
  if (is.null(m)) mfs else mfs[[m]]
}


# Validate newdata, newdataLong and newdataEvent arguments
#
# @param object A stanmvreg object.
# @param newdata A data frame.
# @param newdataLong A data frame, or a list of data frames.
# @param newdataEvent A data frame.
# @param duplicate_ok A logical. If FALSE then only one row per individual is
#   allowed in the newdataEvent data frame.
# @return A (list of) validated data frames.
.validate_newdata_jm <- function(object, newdata = NULL, newdataLong = NULL, 
                                 newdataEvent = NULL, duplicate_ok = FALSE) {
  validate_stanmvreg_object(object)
  id_var   <- object$id_var
  time_var <- object$time_var
  
  # no new data specified by the user
  if (is.null(newdata) && is.null(newdataLong) && is.null(newdataEvent)) {
    return(NULL)
  }
  
  # a single user-supplied data frame for posterior_{predict,traj}
  if (!is.null(newdata)) {
    if (!is.data.frame(newdata))
      stop2("If 'newdata' is specified it must be a data frame.")
    return(newdata)
  }
  
  # otherwise validate new data for posterior_survfit or log_lik
  newdatas <- list()
  if (!is.null(newdataLong)) {
    if (!is(newdataLong, "list"))
      newdataLong <- rep(list(newdataLong), get_M(object))
    df_check <- sapply(newdataLong, is.data.frame)
    if (!all(df_check))
      stop2("'newdataLong' must be a data frame or list of data frames.")
    newdatas <- c(newdatas, newdataLong)
  }
  if (!is.null(newdataEvent)) {
    if (!is.data.frame(newdataEvent))
      stop2("'newdataEvent' must be a data frame.")
    if (!id_var %in% colnames(newdata))
      STOP_no_var(id_var)
    if (!duplicate_ok && any(duplicated(newdataEvent[[id_var]])))
      stop("'newdataEvent' should only contain one row per individual, since ",
           "time varying covariates are not allowed in the prediction data.")
    newdatas <- c(newdatas, list(Event = newdataEvent))
  }
  if (is.jm(object)) {
    id_check <- sapply(newdatas, function(x) id_var %in% colnames(x)) 
    if (!all(id_check)) 
      STOP_no_var(id_var)
    ids <- lapply(newdatas, function(x) unique(x[[id_var]]))
    sorted_ids <- lapply(ids, sort)
    if (!length(unique(sorted_ids)) == 1L) 
      stop("The same subject ids should appear in each new data frame.")
    if (!length(unique(ids)) == 1L) 
      stop("The subject ids should be ordered the same in each new data frame.")  
  }
  return(newdatas)
}

