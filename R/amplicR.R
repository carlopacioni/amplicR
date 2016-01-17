#------------------------------------------------------------------------------#
#Package info 
#------------------------------------------------------------------------------#
#' amplicR 
#' 
#' An R package to process amplicon data.
#' 
#' This package has a number of functions to filter, dereplicate and error and 
#' chimera check NGS
#' reads stored in fastq files. This funcitonality is currently limted to
#' single-read analysis. The retained reads after the data processing described 
#' above (or your own reads if you have done
#' this already in another way) can then be compared against reference sequences
#' and the number of mismatch is reported. This may be useful when, for example,
#' screening samples for particular taxa (e.g. a pathogen).
#' 
#' \code{amplicR} is mainly a wrapper of several functions provided by other R 
#' packages in order to automate some common analyses, so if you use it, please 
#' make sure to also cite the relevant packages (these are generally indicated 
#' in the documentation of \code{amplicR}'s functions). To find out the correct 
#' citation for a package, you can use the function:
#' \code{citation("package_name")} where you have to replace package_name with
#' the actual name of the package you are interested into.
#' 
#' 
#' @section Documentations: Use \code{help(package = "amplicR")} for a list of 
#'   \code{amplicR} functions and their specific documentations.
#'   
#'   A more detailed description of the package and functions can be opened 
#'   with: 
#'   
#'   \code{vignette(package="amplicR", topic="amplicR-package")} 
#'   \emph{ - This is pending... hang in there}
#'      
#' @section Citation: If you use \code{amplicR}, please cite: 
#' 
#' \emph{This is also pending}
#'   
#'   
#' @section Get in touch: Please, use 
#'   \url{https://github.org/carlopacioni/amplicR/issues} to report any issues 
#'   with \code{amplicR}. If unsure, or for feedback, contact me at: 
#'   carlo.pacioni 'at' gmail.com.
#'   
#' @docType package
#'   
#' @name amplicR
NULL
