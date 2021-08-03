## works both for for vars like H2OSOI that have a vector of values (by depth) and also vars with one with  a single value

extractres_h0 <-
  function(sam.start,
           sam.end,
           outdir,
  	   runroot,
	   filebase,
           var.vec.h0,
           scale.vec.h0,
           filterFile,
           start.year,
           end.year) {
    library(ncdf4)
    nmonth <- 12
    nyears <- end.year - start.year + 1
    ncol <- nmonth * nyears
    filter.arr <-
       read.table(file.path(outdir, filterFile), header = F) # on server
    sam.vec <-
      c(sam.start:sam.end)[filter.arr$V1]

    nsam <- length(sam.vec)
    
    ## to get the length of value vector (for H2OSOI the number of depths)#-----
    casename <- paste(filebase, sam.vec[1], sep = "")
    filetag <- paste0("clm2.h0.", start.year, "-")
    monstr <- sprintf("%02d", 1)
    filename <-
      paste0(runroot,
             "/",
             casename,
             "/run/",
             casename,
             ".",
             filetag,
             monstr,
             ".nc")
    nc <- nc_open(filename, write = F)
    res.arr <- vector("list", length = length(var.vec.h0))
    var.length <- vector()
    for (v in 1:length(var.vec.h0)) {
      var.name <- var.vec.h0[v]
      val <- ncvar_get(nc, var.name)
      var.length[v] <- length(val)
      for (k in 1:var.length[v]) {
        res.arr[[v]][[k]] <- matrix(NA, nsam, ncol)
      }
    }
    
    ##----
    pb <- txtProgressBar(min = 0, max = nsam, style = 3)
    cnames <- NULL
    for (yr in start.year:end.year) {
      for (j in 1:12) {
        cnames <- c(cnames, paste0(yr, "-", sprintf("%02d", j)))
      }
    }
    for (i in 1:nsam) {
      sample <- sam.vec[i]
      setTxtProgressBar(pb, i)
      casename <- paste(filebase, sample, sep = "")
      for (yr in start.year:end.year) {
        filetag <- paste("clm2.h0.", yr, "-", sep = "")
        for (j in 1:nmonth) {
          #months
          monstr <- sprintf("%02d", j)
          filename <-
            paste0(runroot,
                   "/",
                   casename,
                   "/run/",
                   casename,
                   ".",
                   filetag,
                   monstr,
                   ".nc")
          nc <- nc_open(filename, write = F)
          for (v in 1:length(var.vec.h0)) {
            var.name <- var.vec.h0[v]
            scale <- scale.vec.h0[v]
            val <- ncvar_get(nc, var.name)
            yri <- yr - start.year
            indx <- j + nmonth * yri
            for (k in 1:var.length[v]) {
              res.arr[[v]][[k]][i, indx] <- val[k] * scale
            }
          }
          nc_close(nc)
        }
      }
      if(i %% 50 == 0) print(i)
    }
    for (v in 1:length(var.vec.h0)) {
      var.name <- var.vec.h0[v]
      for (k in 1:var.length[v]) {
        res.arr[[v]][[k]] <- as.data.frame(res.arr[[v]][[k]])
        colnames(res.arr[[v]][[k]]) <- cnames
        rownames(res.arr[[v]][[k]]) <- sam.vec
      }
      # data is stored as a list
      var.res.arr <- res.arr[[v]]
      out.file.name <- paste0(var.name, ".h0.extract.Rdata")
      if(!dir.exists(file.path(outdir, "extract"))) {dir.create(file.path(outdir, "extract"))}
      save(var.res.arr, file = file.path(outdir, "extract", out.file.name)) # on server
    }
    return("TRUE")
  }

### FOR DAILY output -------
extractres_h1 <-
  function(sam.start,
           sam.end,
           outdir,
   	   runroot,
           filebase,
           var.vec.h1,
           scale.vec.h1,
           filterFile,
           start.year,
           end.year) {
    library(ncdf4)
    cnames <- seq(from = as.Date(paste0(start.year, "-01-01")),
                       to = as.Date(paste0(end.year, "-12-31")),
                       by = 1)
    ## removing Feb 29 because CLM only produces output for 365 days not 366
    cnames <- cnames[format(cnames, "%m-%d") != "02-29"]
    ncol <- length(cnames)
    
    filter.arr <- read.table(file.path(outdir, filterFile), header = F) # on server
    sam.vec <-
      c(sam.start:sam.end)[filter.arr$V1] 
    nsam <- length(sam.vec)
    
    ## to get the length of value vector (for H2OSOI the number of depths)#-----
    casename <- paste(filebase, sam.vec[1], sep = "")
    filetag <- paste0("clm2.h1.", start.year, "-")
    filename <- paste0(runroot,
                       "/",
                       casename,
                       "/run/",
                       casename,
                       ".",
                       filetag,
                       "01-01-00000.nc")
    nc <- nc_open(filename, write = F)
    res.arr <- vector("list", length = length(var.vec.h1))
    var.dim <- vector()
    for (v in 1:length(var.vec.h1)) {
      var.name <- var.vec.h1[v]
      val <- ncvar_get(nc, var.name)
      if (length(dim(val)) == 1) {
        var.dim[v] <- 1
      } else {
        var.dim[v] <- dim(val)[1]
      }
      for (k in 1:var.dim[v]) {
        res.arr[[v]][[k]] <- matrix(NA, nsam, ncol)
      }
    }
    ##----
    #pb <- txtProgressBar(min = 0, max = nsam, style = 3)
    for (i in 1:nsam) {
      sample <- sam.vec[i]
      #setTxtProgressBar(pb, i)
      casename <- paste(filebase, sample, sep = "")
      for (yr in start.year:end.year) {
        filetag <- paste("clm2.h1.", yr, "-", sep = "")
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
        for (v in 1:length(var.vec.h1)) {
          var.name <- var.vec.h1[v]
          scale <- scale.vec.h1[v]
          val <- ncvar_get(nc, var.name)
          index.start <- (yr - start.year) * 365 + 1
          index.end <- index.start + 365 - 1
          if (var.dim[v] == 1) {
            res.arr[[v]][[1]][i, index.start:index.end] <- val * scale
          } else {
            for (k in 1:var.dim[v]) {
              res.arr[[v]][[k]][i, index.start:index.end] <- val[k,] * scale
            }
          } 	
        }
        nc_close(nc)
      }
      #print(sample)
      # if (i %% 200 == 0)
      #  print(sample)
    }
    for (v in 1:length(var.vec.h1)) {
      var.name <- var.vec.h1[v]
      for (k in 1:var.dim[v]) {
        res.arr[[v]][[k]] <- as.data.frame(res.arr[[v]][[k]])
        colnames(res.arr[[v]][[k]]) <- cnames
        rownames(res.arr[[v]][[k]]) <- sam.vec
      }
      # data is stored as a list
      var.res.arr <- res.arr[[v]]
      out.file.name <- paste0(var.name, ".h1.extract.Rdata")
      if(!dir.exists(file.path(outdir, "extract"))) {dir.create(file.path(outdir, "extract"))}
      save(var.res.arr, file = file.path(outdir, "extract", out.file.name)) # on server
    }
      return(TRUE)
  }
# load(file.path(outdir, "data-raw", "extract", out.file.name))
### FOR HOURLY output -------
extractres_h2 <-
  function(sam.start,
           sam.end,
           outdir,
  	   runroot,
           filebase,
           var.vec.h2,
           scale.vec.h2,
           start.year,
           end.year) {
    library(ncdf4)
    cnames <-
      seq(from = as.POSIXct(paste0(start.year, "-01-01 00:00:00", tz = "UTC")),
               to = as.POSIXct(paste0(end.year, "-12-31 23:00:00", tz = "UTC")),
               by = "hour")
    cnames <- cnames[format(cnames, "%m-%d") != "02-29"]
    ncol <- length(cnames)
    filter.arr <- read.table(file.path(outdir, filterFile), header = F) # on server
    sam.vec <-
      c(sam.start:sam.end)[filter.arr$V1]
    nsam <- length(sam.vec)
    
    ## to get the length of value vector (for H2OSOI the number of depths)#-----
    casename <- paste(filebase, sam.vec[1], sep = "")
    filetag <- paste0("clm2.h2.", start.year, "-")
    filename <- paste0(runroot,
                       "/",
                       casename,
                       "/run/",
                       casename,
                       ".",
                       filetag,
                       "01-01-00000.nc")
    nc <- nc_open(filename, write = F)
    val <- ncvar_get(nc, var.name)
    res.arr <- list()
    if (length(dim(val)) == 1) {
      var.dim <- 1
    } else {
      var.dim <- dim(val)[1]
    }
    for (k in 1:var.dim) {
      res.arr[[k]] <- as.data.frame(matrix(NA, nsam, ncol))
    }
    
    ##----
    pb <- txtProgressBar(min = 0, max = nsam, style = 3)
    # days_in_year <- function(year) {
    #   365 + (year %% 4 == 0) - (year %% 100 == 0) + (year %% 400 == 0)
    # }
    for (i in 1:nsam) {
      sample <- sam.vec[sample]
      setTxtProgressBar(pb, i)
      casename <- paste(filebase, sample, sep = "")
        for (yr in start.year:end.year) {
          filetag <- paste("clm2.h2.", yr, "-", sep = "")
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
          val <- ncvar_get(nc, var.name)
          index.start <- 24 * (yr - start.year) * 365 + 1
          index.end <- index.start + 365 * 24 - 1
          if (var.dim == 1) {
            res.arr[[1]][i, index.start:index.end] <- val * scale
          } else {
            for (k in 1:var.dim) {
              res.arr[[k]][i, index.start:index.end] <- val[k,] * scale
            }
          }
          nc_close(nc)
        }
      if (i %% 50 == 0) print(i)
    }
    for (k in 1:var.dim) {
      colnames(res.arr[[k]]) <- cnames
      rownames(res.arr[[k]]) <- sam.vec
      # data is stored as a list
      var.res.arr <- res.arr[[v]]
      out.file.name <- paste0(var.name, ".h2.extract.Rdata")
      if(!dir.exists(file.path(outdir, "extract"))) {dir.create(file.path(outdir, "extract"))}
      save(var.res.arr, file = file.path(outdir, "extract", out.file.name)) # on server
    }
    return("TRUE")
}

