### Restart file outputs -------
extract_treelist <-
  function(sam.start,
           sam.end,
           outdir,
           runroot,
           filebase,
           var.vec.re,
           filterFile,
           finalyear) {
    library(ncdf4)
    filter.arr <-
      read.table(file.path(outdir, filterFile), header = F) # on server
    sam.vec <-
      c(sam.start:sam.end)[filter.arr$V1]
    nsam <- length(sam.vec)
    # base lists
    all.sam.list <- rep(list(), nsam)
    for (i in 1:nsam) {
      sample <- sam.vec[i]
      casename <- paste(filebase, ".", sample, sep = "")
      filetag <- paste("elm.r.", finalyear, "-", sep = "")
      filename <-
        paste0(runroot,
               "/",
               casename,
               "/run/",
               casename,
               ".",
               filetag,
               "01-01-00000.nc")
      nc <- nc_open(filename, write = F)
allvar.vec <- names(nc$var)
    write.table(
      allvar.vec,
      file = file.path(outdir, "allvar.vec.txt"),
      row.names = FALSE
    ) 
      res.arr <- vector("list", length = length(var.vec.re))  
      for (v in 1:length(var.vec.re)) {
        var.name <- var.vec.re[v]
        names(res.arr)[v] <- var.name
        res.arr[[v]] <- ncvar_get(nc, var.name)
      }
      nc_close(nc)
      res.all.df <- do.call(cbind.data.frame, res.arr[-1])
      res.all.df$nsam <- i
      all.sam.list[[i]] <- res.all.df
    }
    all.sam.var <- do.call(rbind, all.sam.list)
    all.sam.var$finalyear <- finalyear
    write.table(
      all.sam.var,
      file = file.path(outdir, paste0("restart_var_outputs.txt")),
      row.names = FALSE
    )
    return(TRUE)
  }
