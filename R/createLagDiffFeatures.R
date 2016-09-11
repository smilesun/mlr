# dates has to be defined to avoid a warning in R CMD CHECK
globalVariables("dates")
#' @title Generate lags and differences for feature variables
#'
#' @description Replace all variables with their generated lagged and differenced variables.
#'
#' @template arg_taskdf
#' @template arg_taskdf_target
#' @param lag [\code{integer}]\cr
#' An integer vector of lag lengths.
#' @param difference [\code{integer}]\cr
#' An integer of the order of differencing
#' @param cols [\code{character}]\cr
#' A character vector of columns to create lag features for.
#' Default is to use all columns. NOTE: For forecast regression tasks, it is not
#' a good idea to make lags of your target variable. So if cols are not specied by
#' the user, createLagDiffFeatures will return a regr task.
#' @param seasonal.lag [\code{integer}]\cr
#' An integer vector of seasonal lag lengths, made as \code{seasonal.lag * frequency}
#' @param seasonal.difference [\code{integer}]\cr
#' An integer of the seasonal order of difference, made as \code{seasonal.difference * frequency}
#' @param frequency [\code{integer}]\cr
#' An integer representing the periodicity in the time series. If frequency is declared in the task,
#' the task frequency will be used.
#' @param na.pad [\code{logical}]\cr
#' A logical to denote whether the data should be padded to the original size with NAs
#' @param difference.lag [\code{integer}]\cr
#' An integer denoting the period to difference over
#' @param seasonal.difference.lag [\code{integer}]\cr
#' An integer denoting the period to seasonaly difference over
#' @param return.nonlag [\code{logical}]\cr
#' A logical to denote whether the original unlagged features should be returned
#' @param grouping [\code{character}]\cr
#' The name of the column to be passed to data.table's \code{by} function. This will take lags and differences wrt the groups.
#' @param date.col [code{data.frame}]
#' The dates for each observation. In the case of a forecasting task, these will be taken from the task description.
#' @export
#' @family eda_and_preprocess
#' @examples
#' set.seed(1234)
#' dat = data.frame(arima.sim(model = list(ar = c(.5,.2), ma = c(.4), order = c(2,0,1)), n = 200))
#' times = as.POSIXct("1992-01-14") + 0:199
#' colnames(dat) = c("arma_test")
#' regr.task = makeRegrTask(id = "Lagged ML model", data = dat, target = "arma_test")
#' regr.task.lag = createLagDiffFeatures(regr.task, lag = 1L:10L, difference = 0L, date.col = times)
createLagDiffFeatures = function(obj, target = character(0L), lag = 0L, difference = 0L, difference.lag = 0L,
  cols = NULL, seasonal.lag = 0L, seasonal.difference = 0L,
  seasonal.difference.lag = 0L, frequency = 1L,
  na.pad = FALSE, return.nonlag = FALSE, grouping = NULL, date.col) {

  assertInteger(lag, lower = 0L, upper = 1000L)
  assertInteger(difference, lower = 0L, upper = 1000L, len = 1L)
  assertInteger(difference.lag, lower = 0L, upper = 1000L, len = 1L)
  assertInteger(seasonal.lag, lower = 0L, upper = 1000L)
  assertInteger(seasonal.difference, lower = 0L, upper = 1000L, len = 1L)
  assertInteger(seasonal.difference.lag, lower = 0L, upper = 1000L, len = 1L)
  assertLogical(na.pad)
  assert(checkClass(obj, "data.frame"), checkClass(obj, "Task"))
  assertCharacter(target, any.missing = FALSE)
  if (!is.null(cols))
    assertCharacter(cols, any.missing = FALSE)

  UseMethod("createLagDiffFeatures")
}

#' @export
createLagDiffFeatures.data.frame = function(obj, target = character(0L), lag = 0L, difference = 0L, difference.lag = 0L,
  cols = NULL, seasonal.lag = 0L, seasonal.difference = 0L,
  seasonal.difference.lag = 0L, frequency = 1L,
  na.pad = FALSE, return.nonlag = FALSE, grouping = NULL, date.col) {

  work.cols = colnames(obj)
  if (missing(date.col)) {
    stop("Dates must be given")
  } else {
    data = as.data.table(obj)
    data[, dates := date.col]
    suppressWarnings(setkeyv(data, c(grouping, "dates")))
  }

  if (!is.null(cols)) {
    if (!(target %in% cols))
      cols = c(cols, target)
    assertSubset(cols, work.cols)
    x = data[, c(cols, grouping), with = FALSE]
  } else {
    cols = work.cols
    x = data[, cols, with = FALSE]
  }
  cols = cols[!(cols %in% grouping)]
  lag.diff.full.names = vector(mode = "character")

  if (any(lag > 0)) {
    lag.vars = cols
    lag.value = lag
    lag.levels = as.vector(vapply(lag.value, function(levels) rep(levels, length(lag.vars)), c(rep(1.0, length(lag.vars)))))
    lag.names = paste0(lag.vars, ".lag.", lag.levels)
    x[, c(lag.names) := shift(.SD, lag.value), by = eval(c(grouping)), .SDcols = lag.vars]
    lag.diff.full.names = c(lag.diff.full.names, lag.names)
  }

  pad  = function(x, n) {
    len.diff = n - length(x)
    c(rep(NA, len.diff), x)
  }

  if (any(difference > 0) | any(difference.lag > 0)) {

    diff.vars = vlapply(x[, c(cols), with = FALSE], is.numeric)
    diff.vars = cols[diff.vars]
    # Since the default for both is zero, which would throw an error, set the zero one to 1
    if (any(difference > 0)) {
      difference.value = difference
    } else {
      difference.value = 1
    }
    if (any(difference.lag > 0)) {
      difference.lag.value = difference.lag
    } else {
      difference.lag.value = 1
    }
    diff.table = expand.grid(difference.value, difference.lag.value)
    diff.full.names = list()[seq_len(nrow(diff.table))]
    for (i in seq_len(nrow(diff.table))) {
      diff.iter = as.numeric(diff.table[i, ])
      diff.lag.names = as.numeric(rep(diff.iter[2], length(diff.vars)))
      diff.diff.names = as.numeric(rep(diff.iter[1], length(diff.vars)))
      diff.names = paste0(diff.vars, ".diff.", diff.diff.names, ".lag.", diff.lag.names)
      x[, c(diff.names) := lapply(.SD,
        function(xx) pad(diff(xx, lag = diff.iter[2], differences = diff.iter[1]), length(xx))),
        by = eval(c(grouping)), .SDcols = diff.vars]
      diff.full.names[[i]] = diff.names
    }
    lag.diff.full.names = c(lag.diff.full.names, unlist(diff.full.names))
  }


  if (frequency > 1L) {
    if (any(seasonal.lag > 0)) {

      seasonal.lag.vars = cols
      seasonal.lag.value = seasonal.lag * frequency
      seasonal.lag.levels = as.vector(vapply(seasonal.lag.value, function(levels) rep(levels, length(seasonal.lag.vars)), c(rep(1.0, length(seasonal.lag.vars)))))
      seasonal.lag.names = paste0(seasonal.lag.vars, ".lag.", seasonal.lag.levels)
      x[, c(seasonal.lag.names) := shift(.SD, seasonal.lag.value), by = eval(c(grouping)), .SDcols = seasonal.lag.vars]
      lag.diff.full.names = c(lag.diff.full.names, seasonal.lag.names)
    }

    if (any(seasonal.difference > 0) | any(seasonal.difference.lag > 0)) {

      diff.vars = vlapply(x[, c(cols), with = FALSE], is.numeric)
      diff.vars = cols[diff.vars]
      # Since the default for both is zero, which would throw an error, set the zero one to 1
      if (any(seasonal.difference > 0)) {
        seasonal.difference.value = seasonal.difference * frequency
      } else {
        seasonal.difference.value = 1
      }
      if (any(seasonal.difference.lag > 0)) {
        seasonal.difference.lag.value = seasonal.difference.lag * frequency
      } else {
        seasonal.difference.lag.value = 1
      }
      seasonal.diff.vars = vlapply(x[, c(cols), with = FALSE], is.numeric)
      seasonal.diff.vars = cols[seasonal.diff.vars]

      seasonal.diff.table = expand.grid(seasonal.difference.value, seasonal.difference.lag.value)
      seasonal.diff.full.names = list()[seq_len(nrow(seasonal.diff.table))]
      for (i in seq_len(nrow(seasonal.diff.table))) {
        seasonal.diff.iter = as.numeric(seasonal.diff.table[i, ])
        seasonal.diff.lag.names = as.numeric(rep(seasonal.diff.iter[2], length(seasonal.diff.vars)))
        seasonal.diff.diff.names = as.numeric(rep(seasonal.diff.iter[1], length(seasonal.diff.vars)))
        seasonal.diff.names = paste0(seasonal.diff.vars, ".diff.", seasonal.diff.diff.names, ".lag.", seasonal.diff.lag.names)
        x[, c(seasonal.diff.names) := lapply(.SD,
          function(xx) pad(diff(xx, lag = seasonal.diff.iter[2], differences = seasonal.diff.iter[1]), length(xx))),
          by = eval(c(grouping)), .SDcols = seasonal.diff.vars]
        seasonal.diff.full.names[[i]] = seasonal.diff.names
      }
      lag.diff.full.names = c(lag.diff.full.names, unlist(seasonal.diff.full.names))
    }
  }

  max.shift = 1:(max(lag, seasonal.lag * frequency,
    max(difference, 1) * max(difference.lag, 1),
    max(seasonal.difference * frequency) * max(seasonal.difference.lag * frequency)))


  if (return.nonlag) {
    data = data[, c(setdiff(work.cols, cols), "dates"), drop = FALSE, with = FALSE]
    if (ncol(data) != 0) {
      data = cbind(data, x)
    } else {
      data = x
    }
  } else {
    data = cbind(data[, unique(c(setdiff(work.cols, cols), "dates", target)), with = FALSE], x[, c(lag.diff.full.names), with = FALSE])
  }
  if (!na.pad) {
    data = data[, .SD[-max.shift, ], by = eval(c(grouping))]
  }
  setkey(data, "dates")
  data$dates = NULL
  return(as.data.frame(data))
}

#' @export
createLagDiffFeatures.Task = function(obj, target = character(0L), lag = 0L, difference = 0L, difference.lag = 0L,
  cols = NULL, seasonal.lag = 0L, seasonal.difference = 0L,
  seasonal.difference.lag = 0L, frequency = 1L,
  na.pad = FALSE, return.nonlag = FALSE, grouping = NULL, date.col) {

  target = getTaskTargetNames(obj)
  data = getTaskData(obj)

  td = getTaskDesc(obj)
  target = getTaskTargetNames(obj)
  if (!is.null(td$frequency) && frequency == 1L)
    frequency = td$frequency
  # We store the original columns as we need them for forecasting
  data.original = data

  if (missing(date.col)) {
    if (is.null(td$dates)) {
      stop("Dates must be supplied")
    } else {
      date.col = td$dates
    }
  }
  data = createLagDiffFeatures( obj = data, target = target, lag = lag, difference = difference,
    difference.lag = difference.lag,
    cols = cols,
    seasonal.lag = seasonal.lag,
    seasonal.difference = seasonal.difference,
    seasonal.difference.lag = seasonal.difference.lag,
    frequency = frequency, na.pad = na.pad,
    return.nonlag = return.nonlag, grouping = grouping, date.col = date.col)

  obj = changeData(obj, data = data)

  max.shift = max(lag, seasonal.lag * frequency,
    max(difference) * max(difference.lag),
    max(seasonal.difference * frequency) * max(seasonal.difference.lag * frequency))
  data.original = data.table(data.original)
  data.original = data.original[,.SD[ (.N - max.shift):.N,], by = eval(c(grouping))]

  obj$task.desc$pre.proc$data.original = data.original
  obj$task.desc$pre.proc$par.vals = list(lag = lag, difference = difference,
    difference.lag = difference.lag,
    cols = cols, target = target,
    seasonal.lag = seasonal.lag,
    seasonal.difference = seasonal.difference,
    seasonal.difference.lag = seasonal.difference.lag,
    frequency = frequency, na.pad = na.pad,
    return.nonlag = return.nonlag, grouping = grouping, date.col = date.col)
  obj
}

