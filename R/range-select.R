#' Select and annotate accelerometry bursts
#'
#' A function that opens an interactive Shiny app, lets the
#' user select named x-axis ranges, and returns them as a
#' list when the app is closed.
#'
#' @param acc An object of class \code{accelerometry}, \code{matrix}, or
#'   \code{data.frame}
#' @param axes A string composed of any non-repeating combination of "x", "y",
#'   "z" (e.g. "xyz", "xz", "y") or \code{NULL}. Columns must be present in
#'   \code{acc}.
#' @param ... passed on to methods
#' @export
range_select <- function(acc, axes = NULL, ...) {
  UseMethod("range_select")
}

#' @export
#' @rdname range_select
range_select.matrix <- function(acc, axes = NULL, ...) {

  # Check for shiny and plotly
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop(
      "Package 'shiny' is required for this function. ",
      "Install it with: install.packages('shiny')"
    )
  }
  if (!requireNamespace("shinyjs", quietly = TRUE)) {
    stop(
      "Package 'shinyjs' is required for this function. ",
      "Install it with: install.packages('shinyjs')"
    )
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "Package 'plotly' is required for this function. ",
      "Install it with: install.packages('plotly')"
    )
  }

  # Input validation
  if (length(dim(acc)) != 2L) {
    stop("`acc` must be a matrix-like object", call. = FALSE)
  }
  if (dim(acc)[2L] < 1L | dim(acc)[2L] > 3L) {
    stop("`acc` must have between 1 and 3 columns. Found ",
         dim(acc)[2L],
         ".",
         call. = FALSE)
  }
  if (dim(acc)[1L] < 1L) {
    stop("`acc` must have at least 1 row. Found ",
         dim(acc)[1L],
         ".",
         call. = FALSE)
  }

  # Get attributes
  sr <- attr(acc, "sampling_rate")
  st <- attr(acc, "start_time")

  # Check column names
  if (is.null(colnames(acc))) {
    axnms <- c("x", "y", "z")
    colnames(acc) <- axnms[1L:dim(acc)[2L]]

    tobe <- ifelse(dim(acc)[2L] > 1L, "s are", " is")
    message(
      "`acc` does not have column names.\n",
      sprintf(
        "The %i column%s assumed to be %s.",
        dim(acc)[2],
        tobe,
        paste(colnames(acc), collapse = ", ")
      )
    )
  }

  # resolve axes to plot
  if (!is.null(axes)) {
    requested  <- .parse_axes(axes)
    available  <- colnames(acc)
    to_plot    <- intersect(requested, available)

    if (length(to_plot) == 0L) {
      stop(
        "None of the requested axes (\"", paste(requested, collapse = ""),
        "\") are present in `acc` (available: \"",
        paste(available, collapse = ""), "\")."
      )
    }
  } else {
    to_plot <- intersect(c("x", "y", "z"), colnames(acc))

    if (length(to_plot) == 0L) {
      stop(
        "No column named \"x\", \"y\", or \"z\" is present in `acc`."
      )
    }
  }

  # plot axes labels
  if (is.null(sr)) {
    x_label <- "Sample"
    x <- seq_len(nrow(acc))
  } else {
    duration <- (nrow(acc) - 1L) / sr

    breaks <- data.frame(
      limit    = c(600, 36000, 864000, 10368000, Inf),
      unit     = c("s", "min", "h", "days", "months"),
      divisor  = c(1, 60, 3600, 86400, 2592000)
    )

    unit_entry <- breaks[which(breaks$limit > duration)[1L], ]

    x_label   <- paste0("Time (", unit_entry$unit, ")")
    x <- seq(0,
             by = (1 / sr) / unit_entry$divisor,
             length.out = nrow(acc))
  }

  y_label <- "acc. (g)"

  clrs <- c(x = "#66C2A5",
            y = "#FC8D62",
            z = "#8DA0CB")

  # UI
  ui <- shiny::fluidPage(

    shinyjs::useShinyjs(),

    shiny::tags$head(shiny::tags$style(shiny::HTML("
      body { font-family: 'Helvetica Neue', sans-serif; padding: 20px; }
      .toolbar { display: flex; align-items: center; gap: 10px;
                 margin-bottom: 12px; flex-wrap: wrap; }
      .mode-pill { font-size: 12px; color: #555; padding: 5px 12px;
                   background: #f0f4f8; border-radius: 20px; }
      .range-box { background: #f8f9fa; border: 1px solid #dee2e6;
                   border-radius: 6px; padding: 16px; margin-top: 16px; }
      .range-box h4 { margin-top: 0; color: #333; }
      pre { background: #1e1e1e; color: #dcdcdc; padding: 14px;
            border-radius: 4px; font-size: 12px; }
    "))),

    # Toolbar
    shiny::div(class = "toolbar",
               shiny::actionButton("clear_btn",  "Clear All",
                                   class = "btn btn-outline-danger"),
               shiny::actionButton("done_btn",   "Done - Return to R",
                                   class = "btn btn-success",
                                   onclick = "setTimeout(function(){window.close();},500);")
    ),

    # Plot
    plotly::plotlyOutput("plot", height = "460px"),

    # Live preview of captured ranges
    shiny::div(class = "range-box",
               shiny::h4("Captured Ranges"),
               shiny::verbatimTextOutput("ranges_display")
    )
  )

  # Server
  server <- function(input, output, session) {

    ranges      <- shiny::reactiveVal(data.frame(label = character(),
                                                 xmin = double(),
                                                 xmax = double()))
    pending_sel <- shiny::reactiveVal(NULL)

    # Plot
    output$plot <- plotly::renderPlotly({
      rng <- ranges()

      shapes <- lapply(seq_len(nrow(rng)), function(i) {
        r <- rng[i,]
        list(type = "rect",
             x0 = r$xmin,
             x1 = r$xmax,
             y0 = 0,
             y1 = 1,
             yref = "paper",
             fillcolor = "rgba(255, 165, 0, 0.25)",
             line  = list(color = "darkorange", width = 1.5),
             layer = "below")
      })

      annots <- lapply(seq_len(nrow(rng)), function(i) {
        r <- rng[i,]
        list(x = (r$xmin + r$xmax) / 2,
             y = 1,
             yref = "paper",
             text = r$label,
             showarrow = FALSE,
             font = list(size = 11, color = "darkorange"),
             xanchor = "center", yanchor = "bottom")
      })

      p <- plotly::plot_ly(type = "scatter",
                           mode = "lines+markers",
                           # line = list(),
                           hovertemplate = paste0(
                             x_label, ": %{x:.3f}<br>",
                             y_label, ": %{y:.3f}<extra></extra>"
                           ))

      # dummy trace to enable select tool - single invisible marker
      p <- plotly::add_trace(p,
                             x = unname(x[1L]),
                             y = unname(acc[1L,1L]),
                             mode = "markers",
                             marker = list(opacity = 0, size = 0),
                             showlegend = FALSE,
                             hoverinfo = "skip",
                             name = "_dummy")

      # add line traces
      for (nm in colnames(acc)) {
        p <- plotly::add_trace(p,
                               x = unname(x),
                               y = unname(acc[, nm]),
                               mode = "lines",
                               line = list(color = unname(clrs[nm]),
                                           width = .7),
                               name = nm)
      }

      p |>
        plotly::layout(
          dragmode        = "select",
          selectdirection = "h",
          xaxis           = list(title = x_label),
          yaxis           = list(title = y_label),
          shapes          = if (length(shapes) > 0) shapes else list(),
          annotations     = if (length(annots)  > 0) annots  else list(),
          showlegend      = TRUE,
          legend          = list(x = 1, y = 1,
                                 xanchor = "right",
                                 yanchor = "top"),
          margin          = list(t = 50)
        ) |>
        plotly::config(scrollZoom = TRUE)
    })

    # Box-select event - modal prompt
    shiny::observeEvent(plotly::event_data("plotly_brushed"), {
      d <- plotly::event_data("plotly_brushed")
      if (is.null(d) || is.null(d$x)) return()

      xmin <- min(d$x)
      xmax <- max(d$x)
      pending_sel(list(xmin = xmin, xmax = xmax))

      shiny::showModal(shiny::modalDialog(
        title = "Capture this range?",
        shiny::p(shiny::strong("X range detected:")),
        shiny::p(sprintf("%.4f  -  %.4f", xmin, xmax),
                 style = "font-family:monospace; font-size:15px; color:#0066cc;"),

        shiny::hr(),

        shiny::textInput("range_label", "Label for this range:",
                         placeholder = "Write here..."),

        footer = shiny::tagList(
          shiny::actionButton("cancel_btn", "Cancel"),
          shiny::actionButton("confirm_btn", "Capture", class = "btn-primary")
        )
      ))
    })

    # Cancel button
    shiny::observeEvent(input$cancel_btn, {
      shinyjs::runjs("Plotly.relayout('plot', {'selections': []})")
      shiny::removeModal()
    })

    # Confirm label - store range
    shiny::observeEvent(input$confirm_btn, {
      sel <- pending_sel()
      if (is.null(sel)) { shiny::removeModal(); return() }

      label <- trimws(input$range_label)

      new_range <- data.frame(label = label,
                              xmin = sel$xmin,
                              xmax = sel$xmax)
      if (is.null(sr)) {
        new_range$imin <- ceiling(new_range$xmin)
        new_range$imax <- floor(new_range$xmax)
      } else {
        new_range$imin <- max(
          1L,
          ceiling(new_range$xmin * unit_entry$divisor * sr) + 1L
        )
        new_range$imax <- min(
          nrow(acc), floor(new_range$xmax * unit_entry$divisor * sr) + 1L
        )
      }
      current   <- ranges()
      current   <- rbind(current, new_range)
      ranges(current)

      shiny::removeModal()
    })

    # Clear
    shiny::observeEvent(input$clear_btn, { ranges(data.frame(label = character(),
                                                             xmin = double(),
                                                             xmax = double())) })

    # Done: stop the app and return ranges to R
    shiny::observeEvent(input$done_btn, {
      shiny::stopApp(returnValue = ranges())
    })

    # Live display
    output$ranges_display <- shiny::renderPrint({
      rng <- ranges()
      if (nrow(rng) == 0) {
        cat("No ranges captured yet.\n")
        cat("Switch to 'Select Range' mode and draw a box on the plot.\n")
        return(invisible(NULL))
      }
      cat("\n")
      for (i in seq_len(nrow(rng))) {
        r   <- rng[i,]
        sep1 <- ifelse(r$label == "", "", " - ")
        sep2 <- ifelse(i < nrow(rng), ",", "")
        cat(sprintf("range %i: from %.4f to %.4f%s%s%s\n",
                    i, r$xmin, r$xmax, sep1, r$label, sep2))
      }
      cat("\n")
    })
  }

  # Launch (blocking) and return result
  out <- shiny::runApp(
    shiny::shinyApp(ui, server)
  )

  print(out)
  invisible(out)
}
