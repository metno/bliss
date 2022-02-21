#+
write_rdata <- function( file, argv, env) {
#------------------------------------------------------------------------------

  save( file=file, env)
  print(paste("output saved on file",file))

}
