### Restart file outputs -------
extract_moisture <-
  function(sam.start,
           sam.end,
           outdir,
           runroot,
           filebase,
           var.vec.re,
           filterFile,
           finalyear,
           fates_res,
           fates_c2b,
           fates_levscls) {
    library(ncdf4)
    library(tidyverse)
    filter.arr <-
      read.table(file.path(outdir, filterFile), header = F) # on server
    sam.vec <-
      c(sam.start:sam.end)[filter.arr$V1]
    nsam <- length(sam.vec)
    # base lists
    all.sam.list <- rep(list(), nsam)
    for (i in 1:nsam) {
      sample <- sam.vec[i]
      casename <- paste0(filebase, ".", sample)
      filetag <- paste0("elm.h0.", finalyear, "-12.nc")
      filename <-
        paste0(runroot,
               "/",
               casename,
               "/run/",
               casename,
               ".",
               filetag)
      nc <- nc_open(filename, write = F)
      res.arr <- vector("list", length = length(var.vec.re))  
      for (v in 1:length(var.vec.re)) {
        var.name <- var.vec.re[v]
        names(res.arr)[v] <- var.name
        res.arr[[v]] <- ncvar_get(nc, var.name)
      }
      res.all.df <- do.call(cbind.data.frame, res.arr)
      res.all.df$nsam <- i
      all.sam.list[[i]] <- res.all.df
      nc_close(nc)
    }

    all.sam.var <- do.call(rbind, all.sam.list)
    write.table(
      all.sam.var,
      file = file.path(outdir, paste0("livefuel.moisturei.raw.txt")),
      row.names = FALSE
    )
    nzmean <- function(x) {
    	zvals <- x==0
    	if (all(zvals)) 0 else mean(x[!zvals])
    }
    moisture <- all.sam.var %>% 
      group_by(nsam) %>%
      # for each nsam, fates_levscls values are for each PFT
      mutate(fates_pft = rep(1:c(n()/fates_levscls), each = fates_levscls)) %>%
      ungroup() %>%
      group_by(nsam, fates_pft) %>%
      #Values are sc by pft, rest are zero
      summarise_all(.funs =list(nzmean)) 
    # this is the average across SCPF. Could have been weighted by biomass of stem and leaves
    write.table(
      moisture,
      file = file.path(outdir, paste0("livefuel.moisture.txt")),
      row.names = FALSE
    )
    return(TRUE)
  }
