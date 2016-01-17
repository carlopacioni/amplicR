#' Set up environment to use \code{amplicR}
#' 
#' Install/update additional packages needed that are not installed the 
#' conventional way.
#' 
#' \code{Bioconductor} recommends installing its packages using 
#' 
#' \code{source("https://bioconductor.org/biocLite.R")} 
#' 
#' and 
#' 
#' \code{biocLite()}.
#' 
#' This is what I do, but this also implies that \code{amplicR}'s dependencies 
#' are not installed when \code{amplicR} is installed. To circumvent this issue,
#' this function will check whether \code{Bioconductor} is installed and whether
#' its packages need to be updated. Also, because \code{dada2} is not currently 
#' available from CRAN nor from \code{Bioconductor}, it is installed by this 
#' function via \code{devtools}. If you have installed or want to install these 
#' packages manually, you don't have to use this function. If you do run this 
#' function, make sure that you do so when you first load \code{amplicR} (i.e. a
#' fresh R session) so that no \code{Bioconductor}'s packages are in use and it
#' is possible to update them if needed.
#' 
#' @export
setup <- function() {
  pkgchk <- function(pkg) {
    if (require(pkg, character.only=TRUE) == T) {
      message(paste("Package ", pkg, " is installed"))
    } else {
      message(paste("Package ", pkg, " is not installed", sep=""))
      message(paste("Installing", pkg), sep="")
      biocLite(pkg)
    }
  }
  source("https://bioconductor.org/biocLite.R")
  biocLite()

  pkgs<-c("ShortRead")
  
  for (pkg in pkgs) {
  pkgchk(pkg)
  }
  if (require("dada2", character.only=TRUE) == T) {
    message(paste("Package ", "dada2", " is installed"))
  } else {
    message(paste("Package ", "dada2", " is not installed", sep=""))
    message(paste("Installing", "dada2"), sep=" ")
    devtools::install_github("benjjneb/dada2")
  }
}
