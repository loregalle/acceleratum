#' Mine stationary accelerometer windows via reservoir sampling
#'
#' Streams a burst-format accelerometer CSV in chunks, reconstructs
#' per-sample timestamps (deriving each row's sampling rate from its
#' timestamp and the next row's), and evaluates a time-based sliding window
#' spanning burst-row boundaries for stationarity (VeDBA, mean per-axis
#' variance, and minimum magnitude). Passing windows are kept via uniform
#' (Vitter) reservoir sampling: only mean vectors stay in memory, while raw
#' window samples are written to a fixed-record scratch binary file for
#' later retrieval by slot offset.
#'
#' @param object burst object or path to a burst-format CSV file.
#' @param data_col Name of the column holding burst data cells.
#' @param ts_col Name of the column holding per-row timestamps.
#' @param axes Axis specification.
#' @param window_sec Sliding window length in seconds.
#' @param sr Nominal sampling rate (Hz). If \code{NULL}, estimated as the
#'   median per-row rate over the first chunk.
#' @param sr_tol Maximum allowed deviation of a row's derived rate from
#'   \code{sr} before its samples are discarded and the window is reset.
#' @param vedba_thresh,var_thresh,mag_thresh Stationarity thresholds (VeDBA,
#'   mean per-axis variance, minimum window mean-vector magnitude).
#' @param min_samples_per_window Minimum samples required before a window
#'   is evaluated.
#' @param reservoir_size Maximum number of windows retained.
#' @param write_to_file If TRUE, writes the reservoir to an external
#' scratch file.
#' @param slot_size Fixed sample capacity per scratch-file record; estimated
#'   from \code{window_sec} and \code{sr} if \code{NULL}.
#' @param scratch_path Path for the scratch binary file.
#' @param overwrite if TRUE, the scratch file is overwritten when present.
#' @param rng_seed rng seed for repeatability purposes.
#' @param chunk_size Rows read per chunk.
#' @param delim CSV delimiter.
#' @param sep Separator used within burst data cells.
#' @param tz Timezone used when parsing non-numeric timestamps.
#'
#' @return A list with reservoir mean vectors (\code{means}), counts, the
#'   scratch file path/slot size needed for reconstruction, and mining
#'   metadata (\code{nominal_sr}, \code{axes}, etc).
#' @export
mine_reservoir <- function(object,
                           data_col = NULL,
                           ts_col = NULL,
                           axes = NULL,
                           window_sec = 5,
                           sr = NULL,
                           sr_tol = 0.5,
                           vedba_thresh = 0.05,
                           var_thresh = 2e-4,
                           mag_thresh = 0.8,
                           min_samples_per_window = 3,
                           reservoir_size = 4000,
                           write_to_file = FALSE,
                           slot_size = NULL,
                           scratch_path = "./tmp",
                           overwrite = FALSE,
                           rng_seed = 42,
                           chunk_size = 10000,
                           delim = ",",
                           sep = " ",
                           tz = "UTC") {

  if(!is.logical(write_to_file) && length(write_to_file) != 1) {
    stop("`write_to_file` must be a logical vector of length 1.")
  }

  if (write_to_file) {
    if (file.exists(scratch_path)) {
      if (overwrite) {
        file.remove(scratch_path)
      } else {
        stop(sprintf(
          paste(
            "A file named '%s' already exists in the working directory.",
            "Remove it or choose a different `scratch_path` before proceeding."
          ),
          scratch_path
        ), call. = FALSE)
      }
    }
  }

  if (!inherits(object, "aclrtm_burst") &&
      !(is.character(object) && length(object) == 1)) {
    stop("`object` must be either an `aclrtm_burst` class object ",
         "or a valid path to a file containing burst-like formatted data.")
  }

  # if burst object, extract attributes
  if (inherits(object, "aclrtm_burst")) {
    if (is.null(data_col)) {
      data_col <- attr(object, "data_col")
    }

    if (is.null(ts_col)) {
      ts_col <- attr(object, "ts_col")
    }

    if (is.null(axes)) {
      axes <- attr(object, "axes")
    }
  }

  # check for data_col and ts_col
  if (is.null(data_col)) {
    stop("`data_col` not defined.")
  }

  if (is.null(ts_col)) {
    stop("`ts_col` not defined.")
  }

  if (is.null(axes)) {
    stop("`axes` not defined.")
  }

  axes_chr <- .parse_axes(axes)
  n_axes   <- length(axes_chr)

  set.seed(rng_seed)

  state <- new.env()
  state$row_counter   <- 0L
  state$pending_row   <- NULL
  state$nominal_sr     <- sr
  state$sr_tol         <- sr_tol
  state$window_sec     <- window_sec
  state$min_samples_per_window <- min_samples_per_window
  state$vedba_thresh   <- vedba_thresh
  state$var_thresh     <- var_thresh
  state$mag_thresh     <- mag_thresh
  state$reservoir_size <- as.integer(reservoir_size)
  state$reservoir_means  <- matrix(NA_real_, nrow = state$reservoir_size,
                                   ncol = n_axes)
  state$reservoir_filled <- 0L
  state$passing_seen     <- 0L
  state$windows_total    <- 0L
  state$write_to_file    <- write_to_file

  init_cap <- max(64L, min_samples_per_window * 4L)
  state$buf_cap   <- init_cap
  state$buf_ts    <- numeric(init_cap)
  state$buf_raw   <- matrix(0, nrow = init_cap, ncol = n_axes)
  state$buf_row   <- integer(init_cap)
  state$buf_start <- 1L
  state$buf_end   <- 0L

  state$slot_size         <- slot_size
  state$scratch_con       <- NULL
  state$scratch_open      <- FALSE
  state$record_size_bytes <- NULL
  state$scratch_path      <- scratch_path
  state$reservoir_raw     <- NULL   # only used when write_to_file = FALSE
  state$n_axes            <- n_axes
  state$bootstrap_needed  <- is.null(sr) # should nominal_sr be estimated

  finalize_setup <- function() {
    if (state$write_to_file) {
      if (is.null(state$slot_size)) {
        state$slot_size <- as.integer(
          floor(
            state$window_sec * (state$nominal_sr + state$sr_tol)
          ) + 1
        )
      }
      state$record_size_bytes <- .record_size_bytes(state$slot_size, n_axes)
      state$scratch_con  <- file(state$scratch_path, open = "wb")
      state$scratch_open <- TRUE
    } else {
      state$reservoir_raw <- vector("list", state$reservoir_size)
    }
  }

  # Ensure the scratch connection (if opened) is always closed, even on error.
  # Tracked via an explicit flag rather than isOpen(), since re-querying an
  # already-closed connection object can itself error ("invalid connection").
  on.exit({
    if (state$scratch_open) {
      close(state$scratch_con)
      state$scratch_open <- FALSE
    }
  }, add = TRUE)

  if (!state$bootstrap_needed) finalize_setup()

  process_chunk <- function(chunk, pos) {

    n_rows <- nrow(chunk)
    if (n_rows == 0L) return(invisible(NULL))

    # burst data to matrix
    if (inherits(object, "aclrtm_burst")) {
      mats <- chunk[[data_col]]
    } else {
      dat_chr <- chunk[[data_col]]
      mats <- lapply(dat_chr, .parse_burst_string, n_axes = n_axes, sep = sep)
    }

    n_rows  <- length(mats)
    ts_chr  <- chunk[[ts_col]]

    # Timestamp column to numeric seconds. Tries "already numeric" first
    # (e.g. unix seconds), falls back to POSIXct string parsing.
    num_try <- suppressWarnings(as.numeric(ts_chr))
    if (!anyNA(num_try)) {
      ts_num <- num_try
    } else {
      ts_num <- as.numeric(as.POSIXct(ts_chr, tz = tz))
    }

    # if sr is not defined, estimate nominal_sr from first chunk
    if (state$bootstrap_needed) {
      sr_vals <- numeric(0)
      if (n_rows >= 2L) {
        for (i in seq_len(n_rows - 1L)) {
          n_i <- nrow(mats[[i]])
          if (n_i > 0L) sr_vals <- c(sr_vals, n_i / (ts_num[i + 1L] - ts_num[i]))
        }
      }
      if (length(sr_vals) == 0L) {
        stop(
          "Not enough rows in the first chunk to estimate a nominal sampling ",
          "rate. Supply `sr` explicitly or increase `chunk_size`.",
          call. = FALSE
        )
      }
      state$nominal_sr <- stats::median(sr_vals)
      state$bootstrap_needed <- FALSE

      if (n_rows < 40L) {
        message(
          sprintf(
            paste(
              "Note: few samples in the first chunk,",
              "the estimated sampling rate (%s Hz) "
            ),
            round(state$nominal_sr, 2)
          ),
          "might be unreliable. Providing a nominal sampling rate ",
          "with argument `sr` is recommended."
        )
      }

      finalize_setup()
    }

    # process chunk rows
    for (i in seq_len(n_rows)) {
      state$row_counter <- state$row_counter + 1L
      current <- list(ts = ts_num[i],
                      mat = mats[[i]],
                      row_id = state$row_counter)
      if (is.null(state$pending_row)) {
        state$pending_row <- current
      } else {
        .resolve_and_process_row(state, state$pending_row, current$ts)
        state$pending_row <- current
      }
    }

    invisible(NULL)
  }

  if (inherits(object, "aclrtm_burst")) {
    process_chunk(object)
  } else if (is.character(object) && length(object == 1)) {

    col_types <- do.call(
      readr::cols_only,
      stats::setNames(list(readr::col_character(),
                           readr::col_character()),
                      c(ts_col, data_col))
    )

    readr::read_delim_chunked(
      file      = object,
      delim     = delim,
      col_types = col_types,
      callback  = readr::SideEffectChunkCallback$new(process_chunk),
      chunk_size = chunk_size
    )
  }

  # The very last row of the file has no successor; timestamp it using the
  # file's nominal sampling rate directly.
  if (!is.null(state$pending_row)) {
    .resolve_and_process_row(state, state$pending_row, sr_override = TRUE)
  }

  if (isTRUE(state$scratch_open)) {
    close(state$scratch_con)
    state$scratch_open <- FALSE
  }

  if (state$reservoir_filled == 0L) {
    if (state$write_to_file &&  file.exists(state$scratch_path)) {
      file.remove(state$scratch_path)
    }
    stop("No stationary windows found with the current thresholds. ",
         "Adjust thresholds and retry.",
         call. = FALSE)
  }

  means_out <- state$reservoir_means[seq_len(state$reservoir_filled), , drop = FALSE]
  colnames(means_out) <- axes_chr

  out <- list(
    means          = means_out,
    n_filled       = state$reservoir_filled,
    windows_total  = state$windows_total,
    passing_seen   = state$passing_seen,
    nominal_sr     = state$nominal_sr,
    axes           = axes_chr,
    n_axes         = n_axes,
    reservoir_size = state$reservoir_size,
    write_to_file  = state$write_to_file
  )

  if (state$write_to_file) {
    out$scratch_path <- state$scratch_path
    out$slot_size    <- state$slot_size
  } else {
    out$raw <- state$reservoir_raw[seq_len(state$reservoir_filled)]
  }

  out
}

#' Select a diverse subset of orientations via Farthest Point Sampling
#'
#' Normalises reservoir mean vectors to the unit sphere and greedily selects
#' \code{k} of them by angular distance (Farthest Point Sampling), repeating
#' with random restarts and keeping the restart maximising mean pairwise
#' angular distance among the selection.
#'
#' @param reservoir A reservoir object as returned by \code{mine_reservoir()}.
#' @param k Number of orientations to select.
#' @param restarts Number of random FPS restarts.
#' @param rng_seed RNG seed for restart selection.
#'
#' @return A list with the selected reservoir indices (\code{selected_idx}),
#'   their unit orientation vectors (\code{orientations}), and the winning
#'   restart's mean pairwise angular distance (\code{score}).
#' @export
select_fps <- function(reservoir, k, restarts = 8, rng_seed = 42) {
  means <- reservoir$means
  n <- nrow(means)
  if (k > n) {
    stop(
      sprintf(
        "Requested %d selections but only %d candidates are available.",
        k, n
      ),
      call. = FALSE
    )
  }

  norms <- sqrt(rowSums(means^2))
  norms[norms < 1e-12] <- 1e-12
  unit_vecs <- means / norms

  run_once <- function(start_idx) {
    sel <- integer(k)
    sel[1] <- start_idx
    dists <- as.vector(
      .angular_distance(unit_vecs,
                        unit_vecs[start_idx, , drop = FALSE])
    )
    if (k > 1) {
      for (m in 2:k) {
        idx <- which.max(dists)
        sel[m] <- idx
        new_dists <- as.vector(
          .angular_distance(unit_vecs,
                            unit_vecs[idx, , drop = FALSE])
        )
        dists <- pmin(dists, new_dists)
      }
    }
    sel
  }

  set.seed(rng_seed)
  best_score <- -1
  best_sel   <- NULL
  for (r in seq_len(restarts)) {
    start <- sample.int(n, 1L)
    sel   <- run_once(start)

    if (length(sel) < 2) {
      sc <- 0
    } else {
      X <- unit_vecs[sel, , drop = FALSE]
      D <- .angular_distance(X, X)
      sc <- mean(D[upper.tri(D)])
    }

    if (sc > best_score) {
      best_score <- sc
      best_sel   <- sel
    }
  }

  orientations <- unit_vecs[best_sel, , drop = FALSE]
  colnames(orientations) <- reservoir$axes

  list(selected_idx = best_sel, orientations = orientations, score = best_score)
}

#' Reconstruct raw sample data for selected reservoir windows
#'
#' Reads raw sample data for the given reservoir slots directly from the
#' scratch binary file by known byte offset (slot index -> offset), with no
#' second pass over the original CSV.
#'
#' @param reservoir A reservoir object as returned by \code{mine_reservoir()}.
#' @param selected_idx Reservoir slot indices to reconstruct (e.g. from
#'   \code{select_fps()$selected_idx}).
#' @param delete_scratch If \code{TRUE}, delete the scratch file after
#'   reading. Default \code{FALSE}, since a reservoir's scratch file may be
#'   read from more than once.
#'
#' @return A list, one element per selected index, each with
#'   \code{timestamp}, \code{data} (raw sample matrix), and \code{source_row}.
#' @export
reconstruct_selected <- function(reservoir, selected_idx,
                                 delete_scratch = FALSE) {
  n_axes    <- reservoir$n_axes

  if (reservoir$write_to_file) {
    slot_size <- reservoir$slot_size
    rec_bytes <- .record_size_bytes(slot_size, n_axes)

    con <- file(reservoir$scratch_path, open = "rb")
    windows <- tryCatch({
      out <- vector("list", length(selected_idx))
      for (k in seq_along(selected_idx)) {
        slot   <- selected_idx[k]
        offset <- (slot - 1) * rec_bytes
        seek(con, where = offset, origin = "start")
        n        <- readBin(con, integer(), n = 1, size = 4)
        ts       <- readBin(con, double(), n = slot_size, size = 8)
        raw_flat <- readBin(con, double(), n = slot_size * n_axes, size = 8)
        raw      <- matrix(raw_flat, nrow = slot_size,
                           ncol = n_axes, byrow = TRUE)
        row_id   <- readBin(con, integer(), n = slot_size, size = 4)

        raw_mat <- raw[seq_len(n), , drop = FALSE]
        colnames(raw_mat) <- reservoir$axes

        out[[k]] <- list(timestamp = ts[seq_len(n)],
                         data = raw_mat,
                         source_row = row_id[seq_len(n)],
                         mean = apply(raw_mat, 2, mean))
      }
      out
    }, finally = close(con))

    if (delete_scratch && file.exists(reservoir$scratch_path)) {
      file.remove(reservoir$scratch_path)
    }
  } else {
    windows <- reservoir$raw[selected_idx]

    for (i in seq_along(selected_idx)) {
      windows[[i]]$mean <- reservoir$means[selected_idx[i],]
      colnames(windows[[i]]$data) <- reservoir$axes
    }
  }

  windows
}

#' Mine, select, and reconstruct diverse stationary windows
#'
#' Convenience wrapper chaining \code{mine_reservoir()},
#' \code{select_fps()}, and \code{reconstruct_selected()}.
#'
#' @param object,data_col,ts_col,axes,window_sec,... Passed to
#'   [mine_reservoir()].
#' @param k Number of orientations to select.
#' @param restarts,rng_seed Passed to [select_fps()].
#' @param delete_scratch Passed to [reconstruct_selected()]; defaults
#'   to \code{TRUE} here, unlike \code{reconstruct_selected()}'s own default.
#'
#' @return A list with FPS-selected \code{orientations}, their raw mean
#'   vectors (\code{means_raw}), and reconstructed raw \code{windows}.
#' @export
mine_and_select <- function(object,
                            k,
                            data_col = NULL,
                            ts_col = NULL,
                            axes = NULL,
                            window_sec = 5,
                            restarts = 8, rng_seed = 42,
                            delete_scratch = TRUE, ...) {
  reservoir <- mine_reservoir(
    object = object, data_col = data_col, ts_col = ts_col, axes = axes,
    window_sec = window_sec, rng_seed = rng_seed, ...
  )
  fps     <- select_fps(reservoir, k = k, restarts = restarts,
                        rng_seed = rng_seed)
  windows <- reconstruct_selected(reservoir, fps$selected_idx,
                                  delete_scratch = delete_scratch)

  windows
}
