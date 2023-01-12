#+
write_rdata <- function( file, arguments, environment) {
#------------------------------------------------------------------------------

  save( file=file, arguments, environment)
  print(paste("output saved on file",file))

}
