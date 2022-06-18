### Restart file outputs -------
extract_treelist <-
  function(sam.start,
           sam.end,
           outdir,
           runroot,
           filebase,
           var.vec.re,
           filterFile,
           finalyear,
           cell.side,
           fates_CWD_frac_twig,
           fates_c2b,
           leafdensity,
           wooddensity) {
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
      # Write all available vars
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
      res.all.df <- do.call(cbind.data.frame, res.arr)
      res.all.df$nsam <- i
      all.sam.list[[i]] <- res.all.df
    }
    all.sam.var <- do.call(rbind, all.sam.list)

    all.sam.var <- all.sam.var %>%
      mutate(finalyear = finalyear,
             # Crown dia per plant
      fates_crown_dia = 2*(fates_cohort_area/fates_nplant/pi)^0.5,
    # Find nplants per cell, from nplants per ha
      fates_nplant.cell = fates_nplant*cell.side^2/10000, # All patch areas together amount to 1 ha or 10000 m2
      fates_height_cbb = fates_height*(1-fates_crown_depth))

    all.sam.var <- all.sam.var %>%
    # From cohort to individual scale biomass
    mutate(fates_bagw_twig = sum(c(sapw_c_val_001, store_c_val_001, repro_c_val_001, struct_c_val_001))*fates_CWD_frac_twig*fates_c2b,
      fates_bleaf = leaf_c_val_001*fates_c2b) %>%
      select(-c(fates_crown_depth, fates_cohort_area, leaf_c_val_001, sapw_c_val_001, store_c_val_001,repro_c_val_001, struct_c_val_001))
   # Biomass weighted bulk density
   all.sam.var <- all.sam.var %>%
      mutate(leaf_twig_bulkd = c(leafdensity*fates_bleaf + wooddensity*fates_bagw_twig)/c(fates_bleaf + fates_bagw_twig))

    ## From nplants per size-pft to individual plants
    trees.whole <- all.sam.var %>%  
      # Selecting only whole trees # This could be made complex by grouping trees into size class bins and checkign if the size class makes it to a whole tree (nplant.400 = 1)
      mutate(fates_nplant.cell = as.integer(round(fates_nplant.cell, 0))) %>% 
      subset(fates_nplant.cell >= 1) %>% 
    # repeat each plant number of plants 
      uncount(fates_nplant.cell) %>%
      select(-fates_nplant)
    trees.whole$treeid <- 1:nrow(trees.whole)

    ## Assign x-y location. Assume simulation index, nsam, increases from East to West, then increases northwards.
    nsam_side = sam.end^0.5
    cell.xy <- data.frame(nsam = c(1:sam.end), 
      xmin = rep(cell.side*c(1:nsam_side-1) + 1, times = nsam_side), 
      ymin = rep(cell.side*c(1:nsam_side-1) + 1, each = nsam_side)) %>% 
      mutate(xmax = xmin + cell.side - 0.1, ymax = ymin + cell.side - 0.01)
    trees.whole <- trees.whole %>%
    left_join(cell.xy, by = "nsam") %>%
    rowwise() %>%  
    mutate(x = runif(1, min = xmin, max = xmax),
           y = runif(1, min = ymin, max = ymax)) %>%
    select(c(fates_pft, x, y, fates_height, fates_height_cbb, fates_crown_dia, leaf_twig_bulkd, fates_bagw_twig, fates_bleaf, treeid))

    write.table(
      trees.whole,
      file = file.path(outdir, paste0("treelist.txt")),
      row.names = FALSE
    )
    return(TRUE)
  }
