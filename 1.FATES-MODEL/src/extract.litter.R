### Restart file outputs -------
extract_litter <-
  function(sam.start,
           sam.end,
           outdir,
           VDM2FM,
	   runroot,
           filebase,
           var.vec.re.R,
           var.vec.re.H,
           filterFile,
           finalyear,
           fire_res,
           fates_res,
           fates_c2b,
           cycle_index,
	   finaltag_month_ci0) {
    library(ncdf4)
    library(tidyverse)
    filter.arr <-
      read.table(file.path(outdir, filterFile), header = F) # on server
    sam.vec <-
      c(sam.start:sam.end)[filter.arr$V1]
    nsam <- length(sam.vec)
    # base lists
    all.sam.list.R <- rep(list(), nsam)
    all.sam.list.H <- rep(list(), nsam)
    for (i in 1:nsam) {
      sample <- sam.vec[i]
      casename <- paste0(filebase, ".", sample)
      #--------------
      # Extracting from Restart file (should have only .R suffix):-------
      #--------------
      if (cycle_index==0 && finaltag_month_ci0 != 12) {
	finalhour <- "00000"
        filetag.R <- paste0("elm.r.", finalyear, "-01-01-", finalhour,".nc") # at the end of the nmonth
	filename.R <- paste0(runroot, "/", casename, "/run/", casename, ".", filetag.R)
	if (file.exists(filename.R) == FALSE) {
	  finalhour <- "01800" 
	  filetag.R <- paste0("elm.r.", finalyear, "-01-01-", finalhour,".nc")
	}	
      } else {
	finalhour <- "00000"
        filetag.R <- paste0("elm.r.", finalyear + 1, "-01-01-", finalhour, ".nc") # a restart file at the end of final year has timestamp for the beginning of next year
      }  
      filename.R <-
        paste0(runroot,
               "/",
               casename,
               "/run/",
               casename,
               ".",
               filetag.R)
      nc.R <- nc_open(filename.R, write = TRUE)
      res.arr.R <- vector("list", length = length(var.vec.re.R))  
      for (v in 1:length(var.vec.re.R)) {
        var.name.R <- var.vec.re.R[v]
        names(res.arr.R)[v] <- var.name.R
        res.arr.R[[v]] <- ncvar_get(nc.R, var.name.R)

        #--------------
        # Reset litter to 0: Satisfies current full fire events, but when fire spread is partial through the area
        # this may be modified in post-fire src/restart.update.treelist.py script
        #--------------
        ncdf4::ncvar_put(nc.R, var.name.R, rep(0, length(res.arr.R[[v]])))
      }
      nc_close(nc.R)
      res.all.df.R <- do.call(cbind.data.frame, res.arr.R)
      res.all.df.R$nsam <- i
      all.sam.list.R[[i]] <- res.all.df.R
      #--------------
      # Extracting from History file (should have only .H suffix):-------
      #--------------
      if (cycle_index==0) {
        filetag.H <- paste0("elm.h0.", finalyear, "-", "01.nc") # at the end of the first month
      } else {
	filetag.H <- paste0("elm.h0.", finalyear, "-", "12.nc") # for full years
      } # History file h1 doesn't have these variables by default
      filename.H <-
        paste0(runroot,
               "/",
               casename,
               "/run/",
               casename,
               ".",
               filetag.H)
      nc.H <- nc_open(filename.H, write = TRUE)
      res.arr.H <- vector("list", length = length(var.vec.re.H))   
      for (v in 1:length(var.vec.re.H)) {
        var.name.H <- var.vec.re.H[v]
        names(res.arr.H)[v] <- var.name.H
        # "FATES_LITTER_CWD_ELDC" has ncwd = 4 size classes.
        # Summing only first two size classes (twigs and small branches), excluding large branches and trunk
        res.arr.H[[v]] <- sum(ncvar_get(nc.H, var.name.H)[1:2], na.rm = TRUE) 
      }
      nc_close(nc.H)
      res.all.df.H <- do.call(cbind.data.frame, res.arr.H)
      res.all.df.H$nsam <- i
      all.sam.list.H[[i]] <- res.all.df.H
    }
    all.sam.var.R <- do.call(rbind, all.sam.list.R)
    all.sam.var.H <- do.call(rbind, all.sam.list.H)

    ## Assign x-y location. Assume simulation index, nsam, increases from East to West, then increases northwards.
    nsam_side = sam.end^0.5
    cell.xy <- data.frame(nsam = c(1:sam.end),
      xmin = rep(fates_res*c(1:nsam_side-1) + 1, times = nsam_side),
      ymin = rep(fates_res*c(1:nsam_side-1) + 1, each = nsam_side)) %>%
      mutate(xmax = xmin + fates_res - 0.1, ymax = ymin + fates_res - 0.01)

    litter_by_nsam <- all.sam.var.R %>%
      select(c(contains("leaf"), nsam)) %>%
      group_by(nsam) %>%
      #Values are only present per patch, rest are zero
      summarise_all(list(sum), na.rm = TRUE) %>%
      # Joining all.sam.var.H here, since values are per site
      left_join(all.sam.var.H, by = "nsam") %>%
      mutate(litter = rowSums(select(., -nsam))) %>% # Already in kg/m2
      select(nsam, litter) %>%
      left_join(cell.xy, by = "nsam")

    grass_by_nsam <- all.sam.var.R %>%
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
