#-------------------------------------------------
# Generate parameter data files with varying parameters 
# Modified by Rutuja Chitra-Tarak
# July 15, 2021
#-------------------------------------------------

#-------------------------------------------------
# Create the files
generate.parameter.files <- function(Param.tab, PARAM_PATH, CLONE_TYPE, clone_base, file_to_clone) {
  #*******************
  ## Load packages----
  #*******************
  # package ncdf4 and tidyverse are needed and are installed wit hthe conda env
 
  #*******************
  ## Process inputs----
  #******************* 
  n.sam <- nrow(Param.tab)
  
  #*******************
  ## Create clones of parameter data file to hold changed parameter----
  #******************* 
  param.dir <- paste0(PARAM_PATH, "/", CLONE_TYPE ,".params")
  unlink(param.dir, recursive=TRUE)
  if(!dir.exists(param.dir)) {dir.create(param.dir)}
  
  for (i in 1:n.sam) {
    basefile <- paste0(PARAM_PATH,"/", file_to_clone)
    clonefile <- paste0(param.dir, "/", clone_base, i, ".nc")
    file.copy(basefile, clonefile, overwrite = T)
  }
  
  #-------------------------------------------------
  #Change the parameter value for texture and organic matter
  #-------------------------------------------------

  par.names <- colnames(Param.tab)
  n.par <- length(par.names)
  
  pb <- txtProgressBar(min = 0, max = n.sam, style = 3)
  for (i in 1:n.sam) {
    setTxtProgressBar(pb, i)
    clonefile2 <- paste0(param.dir, "/", clone_base, i, ".nc")
    param_para <- ncdf4::nc_open(clonefile2, write = T)
    
    for (j in 1:n.par) {
      var.name <- par.names[j]
      if (var.name != "NA") {
        current.param.para.values <- ncdf4::ncvar_get(param_para, var.name)
        num_values <- length(current.param.para.values)
        if (var.name == "aveDTB") {
        param.para.values <- round(Param.tab[i, par.names[j]], 0)
	} else {
	param.para.values <- Param.tab[i, par.names[j]]
	}
        # A temperory solution for multiple PFTs for which both values are same. Ideally separate values should be provided for each PFT.
        ncdf4::ncvar_put(param_para, var.name, rep(param.para.values, num_values))
      } 
    } #j
    ncdf4::nc_close(param_para)
  } #i
  ## Check whether a file has been converted as desired:
  clonefile.test <- paste0(param.dir, "/", clone_base, 1, ".nc")
  para.test <- ncdf4::nc_open(clonefile.test, write = T)
  ## value in the file
  file.value <- ncdf4::ncvar_get(para.test, par.names[1])
  ## value that was to be substituted
  desired.value <- Param.tab[1, par.names[1]]
  ncdf4::nc_close(para.test)
  ## make sure the two do match
  
  # for parametetr like aveDTB, are rounded off to a meter, even though higher res is the desired value
  if (par.names[1] == "aveDTB") {
    match <- all(round(file.value, 0) == round(desired.value, 0))
    } else {
    match <- all(file.value == desired.value)
  }
  return(match)
}

