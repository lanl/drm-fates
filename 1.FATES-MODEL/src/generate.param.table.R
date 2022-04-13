#-------------------------------------------------
# Generate parameter table for sensitivity analysis 
# Modified by Rutuja Chitra-Tarak
# December 21, 2021
#-------------------------------------------------

generate.param.table <- function(df_r, SLICES, outdir, PARAM_Table) {
  # Package slhd is needed and is installed with the conda env
  param.info <- list(par.names = df_r$par.names,
                    min.param = df_r$min.param,
                    max.param = df_r$max.param)
  all.info <- list(par.names = c(param.info$par.names),
                   min.param = c(param.info$min.param),
                   max.param = c(param.info$max.param))

  # Generated parameters should be saved and not changed; unless a new parameter set is required.
  n.sam <- as.numeric(SLICES)
  n.par <- length(all.info$par.names)
  out.slhd <- SLHD::maximinSLHD(t = 1, m = n.sam, k = n.par)
  grid <- out.slhd$StandDesign
  all.params <- matrix(0, n.sam, n.par)

  ## generating ensembles
  for (i in 1: n.sam) {
    for (j in 1: n.par) {
      all.params[i, j] <- qunif(grid[i, j], min = all.info$min.param[j], max = all.info$max.param[j])
    }
  }
  # 
  all.params.df <- data.frame(all.params); colnames(all.params.df) <- all.info$par.names
  if(!dir.exists(outdir)) {dir.create(outdir)}
  write.csv(all.params.df, file = PARAM_Table, row.names = FALSE)
  return(file.exists(PARAM_Table))
}
