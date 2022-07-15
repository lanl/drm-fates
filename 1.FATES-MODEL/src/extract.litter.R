### Restart file outputs -------
extract_litter <-
  function(sam.start,
           sam.end,
           outdir,
           VDM2FM,
	   runroot,
           filebase,
           var.vec.re,
           filterFile,
           finalyear,
           fire_res,
           fates_res,
           fates_c2b,
           cycle_index) {
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
      filetag <- paste0("elm.r.", finalyear + 1, "-", "01-01-00000.nc") # a restart file at the end of final year has timestamp for the beginning of next year
      filename <-
        paste0(runroot,
               "/",
               casename,
               "/run/",
               casename,
               ".",
               filetag)
      nc <- nc_open(filename, write = TRUE)
      res.arr <- vector("list", length = length(var.vec.re))  
      for (v in 1:length(var.vec.re)) {
        var.name <- var.vec.re[v]
        names(res.arr)[v] <- var.name
        res.arr[[v]] <- ncvar_get(nc, var.name)

        #--------------
        # Reset litter to 0: Satisfies current full fire events, but when fire spread is partial through the area
        # this may be modified in post-fire src/restart.update.treelist.py script
        #--------------
        ncdf4::ncvar_put(nc, var.name, rep(0, length(res.arr[[v]])))
      }
      res.all.df <- do.call(cbind.data.frame, res.arr)
      res.all.df$nsam <- i
      all.sam.list[[i]] <- res.all.df
      nc_close(nc)
    }
    all.sam.var <- do.call(rbind, all.sam.list)

    ## Assign x-y location. Assume simulation index, nsam, increases from East to West, then increases northwards.
    nsam_side = sam.end^0.5
    cell.xy <- data.frame(nsam = c(1:sam.end),
      xmin = rep(fates_res*c(1:nsam_side-1) + 1, times = nsam_side),
      ymin = rep(fates_res*c(1:nsam_side-1) + 1, each = nsam_side)) %>%
      mutate(xmax = xmin + fates_res - 0.1, ymax = ymin + fates_res - 0.01)

    litter_by_nsam <- all.sam.var %>% 
      select(c(contains("leaf"), contains("cwd"), nsam)) %>%
      group_by(nsam) %>%
      #Values are only present per patch, rest are zero
      summarise_all(list(sum), na.rm = TRUE) %>%
      mutate(litter = rowSums(.)*fates_c2b) %>% # From kgC/m2 to kg/m2
      select(nsam, litter) %>%
      left_join(cell.xy, by = "nsam")

    grass_by_nsam <- all.sam.var %>%
      select(contains("grass"), nsam) %>%
      group_by(nsam) %>%
      #Values are only present per patch, rest are zero
      summarise_all(list(grass = sum), na.rm = TRUE) %>%
      mutate(grass = grass*fates_c2b) %>% # From kgC/m2 to kg/m2
      left_join(cell.xy, by = "nsam")

    # fates litter is given for fates_res scale, but needs to be averaged for 2m resolution
    # Creating trees_res x trees_res coordinate structure and finding correspondence with low-res coordinates in cell.xy 
    trees_res <- 2
    trees_seq <- seq(from = trees_res, to = fire_res, by = trees_res)
    coord_for_TREES <- data.frame(x = rep(trees_seq, times = fire_res/trees_res), 
                       y = rep(trees_seq, each = fire_res/trees_res)) %>%
      mutate(xmax = as.numeric(as.character(cut(x, breaks = c(0, unique(cell.xy$xmax)), labels = unique(cell.xy$xmax)))),
             ymax = as.numeric(as.character(cut(y, breaks = c(0, unique(cell.xy$ymax)), labels = unique(cell.xy$ymax))))) 

    litter <- coord_for_TREES %>%
      left_join(litter_by_nsam %>% select(litter, xmax, ymax), by = c("xmax", "ymax")) %>%
      select(-xmax, -ymax) %>%
      arrange(desc(y), x) %>%
      # Make sure litter doen't have zeroes
      mutate(litter = ifelse(litter == 0, 0.001, litter))

    grass <- coord_for_TREES %>%
      left_join(grass_by_nsam %>% select(grass, xmax, ymax), by = c("xmax", "ymax"))  %>%
      select(-xmax, -ymax) %>%
      arrange(desc(y), x)  %>%
      # Make sure litter doen't have zeroes
      mutate(grass = ifelse(grass == 0, 0.001, grass))
          
    litter.m <- matrix(litter$litter, nrow = length(trees_seq), byrow = TRUE)
    grass.m <- matrix(grass$grass, nrow = length(trees_seq), byrow = TRUE)

    if(!dir.exists(VDM2FM)) {dir.create(VDM2FM)}
    write.table(
      litter.m,
      file = file.path(VDM2FM, "VDM_litter_trees.dat"),
      row.names = FALSE, col.names = FALSE
    )
    write.table(
      grass.m,
      file = file.path(VDM2FM, "VDM_litter_WG.dat"),
      row.names = FALSE, col.names = FALSE
    )
    write.table(
      litter.m,
      file = file.path(VDM2FM, paste0("VDM_litter_trees.", cycle_index,".dat")),
      row.names = FALSE, col.names = FALSE
    )
    write.table(
      grass.m,
      file = file.path(VDM2FM, paste0("VDM_litter_WG.", cycle_index,".dat")),
      row.names = FALSE, col.names = FALSE
    )
    return(TRUE)
 } 
