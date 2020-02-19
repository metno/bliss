  if ( argv$mode=="ensigap" ) {
    save( file= argv$off_rdata,
          # cv
          yo_cv, ya_cv, ya_cv_var, yb_cv, yidi_cv, ya_cv_alpha, ya_cv_pdf_par, Yb_cv
          )
  }
  print(paste("output saved on file",argv$off_rdata))
