system.time(expr = mclapply(X = 1:20, FUN = function(x) Sys.sleep(10), mc.cores = 20))
