#' Shiny Application for interactive visualization of MMD,GMD and
#' Pearson Difference as well as plotting peaks
#'
#' @inheritParams DBAmmd-Accessors
#' @inheritParams getPeakReads
#'
#' @examples
#'  if(interactive()){
#'   data("MMD")
#' runShinyMMDiff2(MMD)
#'}
#'
#' @import RColorBrewer
#' @import shiny
#' @import ggplot2
#'
#' @export
#'
#
# David Kuo
# March 2016
runShinyMMDiff2 <- function(MD, whichContrast=1) {
  shinyMMDiff2.app <- list(ui = ui.MMDiff2(MD),
                           server = server.MMDiff2(MD, whichContrast))
  runApp(shinyMMDiff2.app)
}
