#' Shiny Application for interactive visualization of MMD,GMD and
#' Pearson Difference as well as plotting peaks
#'
#' @param MD DBAmmd object
#' @param con (DEFAULT: NULL)
#'
#' @examples
#'  if(interactive()){
#' load(system.file("data/MMD.RData", package="MMDiff2"))
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
runShinyMMDiff2 <- function(MD, con=NULL) {
  shinyMMDiff2.app <- list(ui = ui.MMDiff2(MD),
                           server = server.MMDiff2(MD, con))
  runApp(shinyMMDiff2.app)
}
