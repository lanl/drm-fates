### Restart file outputs -------
extract_litter <-
  function(sam.start,
           sam.end,
           outdir,
           runroot,
           filebase,
           var.vec.re,
           var.vec.re_2,
           filterFile,
           finalyear,
           cell.side,
           fates_c2b) {
    library(ncdf4)
    library(tidyverse)
    filter.arr <-
      read.table(file.path(outdir, filterFile), header = F) # on server
    sam.vec <-
      c(sam.start:sam.end)[filter.arr$V1]
    nsam <- length(sam.vec)
    # base lists
    all.sam.list <- rep(list(), nsam)
    all.sam.list.moisture <- rep(list(), nsam)
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
      #For litter
      res.arr <- vector("list", length = length(var.vec.re))  
      for (v in 1:length(var.vec.re)) {
        var.name <- var.vec.re[v]
        names(res.arr)[v] <- var.name
        res.arr[[v]] <- ncvar_get(nc, var.name)
      }
      res.all.df <- do.call(cbind.data.frame, res.arr)
      res.all.df$nsam <- i
      all.sam.list[[i]] <- res.all.df

      # For litter moisture
      res.arr.moisture <- vector("list", length = length(var.vec.re_2))
      for (v in 1:length(var.vec.re_2)) {
        var.name <- var.vec.re_2[v]
        names(res.arr.moisture)[v] <- var.name
        res.arr.moisture[[v]] <- ncvar_get(nc, var.name)
      }
      res.all.df.moisture <- do.call(cbind.data.frame, res.arr.moisture)
      res.all.df.moisture$nsam <- i
      all.sam.list.moisture[[i]] <- res.all.df.moisture

      nc_close(nc)
    }
    all.sam.var <- do.call(rbind, all.sam.list)
    all.sam.var.moisture <- do.call(rbind, all.sam.list.moisture)

    litter <- all.sam.var %>% 
      group_by(nsam) %>%
      #Values are only present per patch, rest are zero
      summarise_all(list(sum), na.rm = TRUE) %>%
      mutate(total.litter = rowSums(.)*fates_c2b) %>% # From kgC/m2 to kg/m2
      select(nsam, total.litter)
    
      ## Assign x-y location. Assume simulation index, nsam, increases from East to West, then increases northwards.
    nsam_side = sam.end^0.5
    cell.xy <- data.frame(nsam = c(1:sam.end),
      xmin = rep(cell.side*c(1:nsam_side-1) + 1, times = nsam_side),
      ymin = rep(cell.side*c(1:nsam_side-1) + 1, each = nsam_side)) %>%
      mutate(xmax = xmin + cell.side - 0.1, ymax = ymin + cell.side - 0.01)
    litter <- litter %>%
    left_join(cell.xy, by = "nsam") %>%
    rowwise() %>%
    mutate(x = runif(1, min = xmin, max = xmax),
           y = runif(1, min = ymin, max = ymax)) %>%
    select(-xmin, -xmax, -ymin, -ymax)

    litter.moisture <- all.sam.var.moisture %>% 
    left_join(cell.xy, by = "nsam") %>%
    rowwise() %>%
    mutate(x = runif(1, min = xmin, max = xmax),
           y = runif(1, min = ymin, max = ymax)) %>%
    select(-xmin, -xmax, -ymin, -ymax)

    write.table(
      litter,
      file = file.path(outdir, paste0("litter.txt")),
      row.names = FALSE
    )
    write.table(
      litter.moisture,
      file = file.path(outdir, paste0("litter.moisture.txt")),
      row.names = FALSE
    )
    return(TRUE)
  }
