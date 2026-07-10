# ============================================================
# mine_fps_burst.R
#
# Streams a deployment-scale burst-format accelerometer CSV in
# chunks (via readr::read_delim_chunked), mines stationary
# windows using a true time-based sliding window (which may
# span multiple burst rows), and selects a geometrically
# diverse subset of orientations via Farthest Point Sampling
# (FPS) -- an R port of the design discussed for
# AccelerationCalibrationMiner10.py, adapted to the
# aclrtm_burst representation and to streaming/chunked reading.
#
# Design notes / departures from the Python original
# ----------------------------------------------------
# * Sliding window is sample-level and time-based, spans burst
#   row boundaries freely -- mirrors _sliding_windows() in the
#   Python script.
# * Per-row sampling rate (sr) is NOT read from a column. It is
#   derived the way aclrtm_burst / burst_to_accelerometry()
#   already does elsewhere in this package:
#       sr_i = nrow(data_i) / (ts_{i+1} - ts_i)
#   This requires a one-row lookahead, which must be carried
#   across chunk boundaries (the last row of each chunk is held
#   back as `pending_row` until the next chunk's first row
#   arrives).
# * A row whose derived sr_i deviates from the file's nominal
#   sr by more than `sr_tol` cannot be reliably timestamped
#   internally, so its samples are discarded entirely and the
#   sliding window is cleared -- this plays the same role as
#   the Python script's gap_thresh reset, just detected via the
#   rate check instead of a raw elapsed-time check. No separate
#   gap_thresh is needed as a result.
# * The nominal sampling rate (if not supplied by the user) is
#   estimated as the median of sr_i computed from the pairs of
#   consecutive rows fully contained within the FIRST CHUNK
#   only (a one-time "bootstrap" pass over that chunk before
#   any window-mining happens). A message is issued if the
#   first chunk is small (<100 rows).
# * Reservoir sampling (Vitter-style, uniform) keeps only the
#   MEAN VECTOR of each passing window in memory. The raw
#   samples of every window admitted to the reservoir are
#   written to a scratch binary file using FIXED-SIZE records,
#   so that "reservoir slot index" directly maps to a byte
#   offset (slot i -> offset (i-1)*record_size) with no
#   separate offset table and no second pass over the original
#   CSV. Only the k windows selected by FPS are ever read back,
#   directly from their known offsets.
# * Stationarity test uses VeDBA (Euclidean norm of the dynamic
#   acceleration, per sample, averaged over the window) instead
#   of ODBA (L1 norm), plus mean-of-per-axis-variance, plus a
#   minimum window mean-vector magnitude (to avoid unstable
#   direction estimates from near-zero vectors) -- all three
#   must pass for a window to be admitted to the reservoir.
# * Identity calibration only -- no calibration hook.
# * No diagnostics CSV, no plots. Output is purely: selected
#   unit orientations, their raw mean vectors, and the raw
#   sample data for each of the k selected windows.
#
# Composable functions
# ---------------------
#   mine_reservoir()       -- Pass 1: stream + mine + reservoir-sample
#   select_fps()           -- FPS selection on the reservoir means
#   reconstruct_selected() -- read the k selected windows' raw data
#   mine_and_select()       -- convenience wrapper chaining all three
#
# Dependency: readr (Imports). Add to DESCRIPTION:
#     Imports:
#         readr
# ============================================================


# ============================================================
# mine_reservoir()
#
# Pass 1: stream the CSV in chunks, reconstruct per-sample
# timestamps, run the burst-spanning sliding window, test each
# window for stationarity (VeDBA + variance + magnitude), and
# keep a uniform reservoir sample of admitted windows. Only
# mean vectors are kept in memory; raw window samples are
# written to a scratch binary file (fixed-size records) as
# they are admitted/replaced.
# ============================================================
mine_reservoir <- function(path,
                           data_col,
                           ts_col,
                           axes,
                           window_sec = 5,
                           sr = NULL,
                           sr_tol = 0.5,
                           vedba_thresh = 0.005,
                           var_thresh = 2e-4,
                           mag_thresh = 0.3,
                           min_samples_per_window = 3,
                           reservoir_size = 80000,
                           slot_size = NULL,
                           scratch_path = "tmp",
                           rng_seed = 42,
                           chunk_size = 10000,
                           delim = ",",
                           sep = " ",
                           tz = "UTC") {

  if (file.exists(scratch_path)) {
    stop(sprintf(
      paste(
        "A file named '%s' already exists in the working directory.",
        "Remove it or choose a different `scratch_path` before proceeding."
      ),
      scratch_path
    ), call. = FALSE)
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
  state$n_axes            <- n_axes
  state$bootstrap_needed  <- is.null(sr)

  finalize_setup <- function() {
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
  }

  # Ensure the scratch connection (if opened) is always closed, even on error.
  # Tracked via an explicit flag rather than isOpen(), since re-querying an
  # already-closed connection object can itself error ("invalid connection").
  on.exit({
    if (isTRUE(state$scratch_open)) {
      close(state$scratch_con)
      state$scratch_open <- FALSE
    }
  }, add = TRUE)

  if (!state$bootstrap_needed) finalize_setup()

  process_chunk <- function(chunk, pos) {
    ts_chr  <- chunk[[ts_col]]
    dat_chr <- chunk[[data_col]]
    n_rows  <- length(ts_chr)
    if (n_rows == 0L) return(invisible(NULL))

    # Timestamp column to numeric seconds. Tries "already numeric" first
    # (e.g. unix seconds), falls back to POSIXct string parsing.
    num_try <- suppressWarnings(as.numeric(ts_chr))
    if (!anyNA(num_try)) {
      ts_num <- num_try
    } else {
      ts_num <- as.numeric(as.POSIXct(ts_chr, tz = tz))
    }

    mats <- lapply(dat_chr, .parse_burst_string, n_axes = n_axes, sep = sep)

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

  col_types <- do.call(
    readr::cols_only,
    stats::setNames(list(readr::col_character(),
                         readr::col_character()),
                    c(ts_col, data_col))
  )

  readr::read_delim_chunked(
    file      = path,
    delim     = delim,
    col_types = col_types,
    callback  = readr::SideEffectChunkCallback$new(process_chunk),
    chunk_size = chunk_size
  )

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
    if (file.exists(state$scratch_path)) file.remove(state$scratch_path)
    stop("No stationary windows found with the current thresholds. ",
         "Adjust thresholds and retry.",
         call. = FALSE)
  }

  means_out <- state$reservoir_means[seq_len(state$reservoir_filled), , drop = FALSE]
  colnames(means_out) <- axes_chr

  list(
    means          = means_out,
    n_filled       = state$reservoir_filled,
    windows_total  = state$windows_total,
    passing_seen   = state$passing_seen,
    scratch_path   = state$scratch_path,
    slot_size      = state$slot_size,
    nominal_sr     = state$nominal_sr,
    axes           = axes_chr,
    n_axes         = n_axes,
    reservoir_size = state$reservoir_size
  )
}

# ============================================================
# select_fps()
#
# Normalises reservoir mean vectors to the unit sphere and runs
# greedy Farthest Point Sampling (angular distance) with random
# restarts, keeping the restart whose selection maximises mean
# pairwise angular distance. Vectors with magnitude below
# mag_thresh were already excluded from the reservoir in Pass 1
# (mine_reservoir()), so no further filtering happens here.
# ============================================================
select_fps <- function(reservoir, k, restarts = 8, rng_seed = 42) {
  means <- reservoir$means
  n <- nrow(means)
  if (k > n) {
    stop(sprintf("Requested %d selections but only %d candidates are available.", k, n),
         call. = FALSE)
  }

  norms <- sqrt(rowSums(means^2))
  norms[norms < 1e-12] <- 1e-12
  unit_vecs <- means / norms

  run_once <- function(start_idx) {
    sel <- integer(k)
    sel[1] <- start_idx
    dists <- as.vector(.angular_distance(unit_vecs, unit_vecs[start_idx, , drop = FALSE]))
    if (k > 1) {
      for (m in 2:k) {
        idx <- which.max(dists)
        sel[m] <- idx
        new_dists <- as.vector(.angular_distance(unit_vecs, unit_vecs[idx, , drop = FALSE]))
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

# ============================================================
# reconstruct_selected()
#
# Reads raw sample data for the FPS-selected reservoir slots
# directly from the scratch file (slot index -> known byte
# offset), with no second pass over the original CSV. Deletes
# the scratch file afterwards by default.
# ============================================================
reconstruct_selected <- function(reservoir, selected_idx,
                                 delete_scratch = FALSE) {
  n_axes    <- reservoir$n_axes
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
      raw      <- matrix(raw_flat, nrow = slot_size, ncol = n_axes, byrow = TRUE)
      row_id   <- readBin(con, integer(), n = slot_size, size = 4)

      raw_mat <- raw[seq_len(n), , drop = FALSE]
      colnames(raw_mat) <- reservoir$axes

      out[[k]] <- list(timestamp = ts[seq_len(n)],
                       data = raw_mat,
                       source_row = row_id[seq_len(n)])
    }
    out
  }, finally = close(con))

  if (delete_scratch && file.exists(reservoir$scratch_path)) {
    file.remove(reservoir$scratch_path)
  }

  windows
}

# ============================================================
# mine_and_select()
#
# Convenience wrapper chaining mine_reservoir() -> select_fps()
# -> reconstruct_selected(). Returns a list with the FPS-
# selected unit orientations, their raw mean vectors, and the
# raw sample data for each selected window.
# ============================================================
mine_and_select <- function(path, data_col, ts_col, axes, window_sec, k,
                            restarts = 8, rng_seed = 42,
                            delete_scratch = TRUE, ...) {
  reservoir <- mine_reservoir(
    path = path, data_col = data_col, ts_col = ts_col, axes = axes,
    window_sec = window_sec, rng_seed = rng_seed, ...
  )
  fps     <- select_fps(reservoir, k = k, restarts = restarts, rng_seed = rng_seed)
  windows <- reconstruct_selected(reservoir, fps$selected_idx,
                                  delete_scratch = delete_scratch)

  list(
    orientations = fps$orientations,
    means_raw    = reservoir$means[fps$selected_idx, , drop = FALSE],
    windows      = windows
  )
}
