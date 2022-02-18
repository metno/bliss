#+
write_yenv_rdata <- function( argv, y_env) {
#------------------------------------------------------------------------------

  save( file=argv$off_yenv_rdata, y_env)
  print(paste("output saved on file",argv$off_yenv_rdata))

}
