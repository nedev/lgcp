##' .onAttach function
##'
##' A function to print a welcome message on loading package  
##'
##' @param libname libname argument
##' @param pkgname pkgname argument
##' @return ...
##' @export

.onAttach <- function(libname, pkgname)
{
	packageStartupMessage("\n Welcome to 'lgcp': Log-Gaussian Cox Process\n B. Taylor & T.M. Davies & B. Rowlingson & P. Diggle\ntype '?lgcp' for details, or 'vignette(\"lgcp\")' to view the package vignette. Type 'citation(\"lgcp\")' to view the citation for this package.)", appendLF=T)
}

