#' ui component for interactive visualization of MMD,GMD and Pearson Difference
#' as well as plotting peaks
#'
#' @param MD DBAmmd object
#'
#' @examples
#'  if(interactive()){
#' load(system.file("data/MMD.RData", package="MMDiff2"))
#' runShinyMMDiff2(MMD)
#'}
#'
#' @import RColorBrewer
#' @import shiny
#
# David Kuo
# March 2016

ui.MMDiff2 <- function(MD) {
  fluidPage(
    fluidRow(column(
      width = 12,
      class = "well",
      selectInput(
        "select",
        label = h3("Select Distance"),
        choices = list(
          "MMD Distance" = "plot_MMD"
        ),
        selected = "plot_MMD"
      )
    )),
    fluidRow(column(
      width = 12,
      class = "well",
      fluidRow(
        column(
          width = 6,
          h4("Zoom with box and Double-click to reset"),
          plotOutput(
            "plot_static",
            height = 250,
            brush = brushOpts(id = "plot_static_brush",
                              resetOnNew = TRUE)
          )
        ),
        h4("Select Peak Profiles"),
        column(
          width = 6,
          plotOutput(
            "plot_dynamic",
            height = 250,
            click = clickOpts(id = "plot_click")
          )
        )
      )
    )),
    #row
    fluidRow(column(
      width = 12, class = "well",
      fluidRow(column(
        width = 6,
        h4("Region Plot for all Samples"),
        plotOutput("peak_plot",
                   height=250)
      ),
      h4("Distances for Selected Peak"),
      column(
        width = 6,
        plotOutput(
          "plot_dist_for_peak",
          height = 250)
      )
      )
    )),
    fluidRow(column(
      width = 12,
      class = "well",
      h5(
        "UCSC Link:    ",
        htmlOutput("ucsc_link")
      )
    ))#row
  ) #page
} #function
