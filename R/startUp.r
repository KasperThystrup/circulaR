.onAttach <- function(libname, pkgname){
  packageStartupMessage("Welcome to circulaR package")
  # Windows mclapply
  if (Sys.info()["sysname"] == "Windows"){
    warning("Inferior operating system detected, redefining mclapply to support Windows!" )
    setMethod(f = mclapply,
              definition = function(X, FUN, ..., mc.cores){
      tryCatch(expr = {
        cl <- parallel::makeCluster(getOption("cl.cores", mc.cores))
        return(parallel::parLapply(cl = cl, X = X, fun = FUN))
      }, finally = parallel::stopCluster(cl))
    })
  }
}
