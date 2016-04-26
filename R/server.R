#' Shiny server code for interactive visualization of MMD distances,
#' peak plots, and MMD distances by sample.
#'
#' @inheritParams compDists
#' @inheritParams plotPeak
#'
#' @examples
#'  if(interactive()){
#'     data("MMD")
#'     runShinyMMDiff2(MMD)
#' }
#'
#' @import RColorBrewer shiny
#
# David Kuo
# March 2016
server.MMDiff2 <- function(MD, whichContrast=1) {
  function(input, output) {
    #out_df <- out_df[complete.cases(out_df), ]
    out_df <- plotDists(MD, shiny_df_opt = TRUE)
    out_df$peakPos <- rownames(out_df)
    #rownames(out_df) <- 1:nrow(out_df)

    # -------------------------------------------------------------------
    # Linked plots (middle and right)
    # -------------------------------------------------------------------

    #set up reactive values for ranges
    ranges <- reactiveValues(x = NULL, y = NULL)
    if (!is.null(out_df)) {
      output$plot_static <- renderPlot({
        condition <- NULL
        sig <- NULL
        means <- NULL
        rm(list=c(condition, sig, means))
        ggplot(out_df,
               aes(
                 means,
                 distance,
                 shape = factor(condition),
                 color = factor(sig)
               )) +
          scale_shape_manual(values = c(20, 4)) +
          theme_bw() +
          geom_point()
      })

      output$plot_dynamic <- renderPlot({
        p <-
          ggplot(out_df,
                 aes(
                   means,
                   distance,
                   shape = factor(condition),
                   color = factor(sig)
                 )) +
          scale_shape_manual(values = c(20, 4)) +
          geom_point() +
          theme_bw() +
          coord_cartesian(xlim = ranges$x, ylim = ranges$y)
        if (is.null(click_saved$singleclick)) {
          p
        } else {
          p + geom_vline(aes(xintercept = click_saved$singleclick$x),
                         linetype = 3) +
            geom_hline(aes(yintercept = click_saved$singleclick$y),
                       linetype = 3)
        }
      })
    }

    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observe({
      brush <- input$plot_static_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }

    })

    click_saved <- reactiveValues(singleclick = NULL)
    observeEvent(eventExpr = input$plot_click, handlerExpr = {
      click_saved$singleclick <- input$plot_click
    })

    # -------------------------------------------------------------------
    # Plot chosen peak
    # -------------------------------------------------------------------
    output$peak_plot <- renderPlot({
      if (is.null(click_saved$singleclick)) {
        selected_peak_id <- 'W-27'
      } else {
        selected_peak_id <-
          row.names(
            nearPoints(
              out_df,
              click_saved$singleclick,
              xvar = "means",
              yvar = "distance",
              maxpoints = 1

          ))
      }

      #if (selected_peak_id > dim(as.data.frame(Dists(MD)))[1]) {
      #  selected_peak_id <- selected_peak_id - dim(as.data.frame(Dists(MD)))[1]
      #}
      Peak.id = strsplit(selected_peak_id,'-')[[1]][2]

      plotPeak(
        MD = MD,
        Peak.id = Peak.id,
        NormMethod = NULL,
        plot.input = FALSE,
        whichContrast = whichContrast
      )

    })


    # -------------------------------------------------------------------
    # Plot distances for chosen peak
    # -------------------------------------------------------------------

    output$plot_dist_for_peak <- renderPlot({
      if (is.null(click_saved$singleclick)) {
        selected_peak_id <- 'W-27'
      } else {
        selected_peak_id <-
          row.names(
            nearPoints(
              out_df,
              click_saved$singleclick,
              xvar = "means",
              yvar = "distance",
              maxpoints = 1
          ))
      }

      #if (selected_peak_id > dim(as.data.frame(Dists(MD)))[1]) {
      #  selected_peak_id <- selected_peak_id-dim(as.data.frame(Dists(MD)))[1]
      #}
      # Log
      Peak.id = strsplit(selected_peak_id,'-')[[1]][2]
      #Peak.id = names(Regions(MD))[selected_peak_id]
      plotDISTS4Peak(
        MD,
        Peak.id = Peak.id,
        dist.method = 'MMD',
        whichContrast = whichContrast,
        Zoom = TRUE
      )
    })

    # -------------------------------------------------------------------
    # UCSC Genome Browser link generation
    # -------------------------------------------------------------------
    output$ucsc_link <- renderUI({
      if (is.null(click_saved$singleclick)) {
        selected_peak_id <- 'W-27'
      } else {
        selected_peak_id <-
          row.names(
            nearPoints(
              out_df,
              click_saved$singleclick,
              xvar = "means",
              yvar = "distance",
              maxpoints = 1
          ))
      }
      Peak.id = strsplit(selected_peak_id,'-')[[1]][2]
      ii <- match(Peak.id,names(Regions(MD)))
      ucsc_coord <- Regions(MD)[ii]
      chrom <- seqnames(ucsc_coord)
      left_coord <- start(ucsc_coord)
      right_coord <- end(ucsc_coord)
      if (grepl("hg18", Genome(MD))) {
        org_str <- "human"
        db_str <- "hg18"
      } else if (grepl("hg19", Genome(MD))) {
        db_str <- "hg19"
        org_str <- "human"
      } else if (grepl("hg38", Genome(MD))) {
        db_str <- "hg38"
        org_str <- "human"
      } else if (grepl("mm9", Genome(MD))) {
        db_str <- "mm9"
        org_str <- "mouse"
      } else if (grepl("mm10", Genome(MD))) {
        db_str <- "mm10"
        org_str <- "mouse"
      } else {
        db_str <- NULL
      }
      if (!is.null(db_str)) {
        peak_link <-
          sprintf(
            "http://genome.ucsc.edu/cgi-bin/hgTracks?org=%s&db=%s&position=%s:%s-%s",
            org_str,
            db_str,
            chrom,
            left_coord,
            right_coord
          )
      } else {
        peak_link <- "http://genome.ucsc.edu"
      }
      tags$a(href = peak_link, peak_link, target="_blank")
    })
  }
}
