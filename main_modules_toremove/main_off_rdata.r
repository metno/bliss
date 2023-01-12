  if ( argv$mode=="ensigap" ) {
    cat( paste( "rmse_b=", round(sqrt(mean((yo_cv-yb_cv)**2)),3),
                "rmse_a=", round(sqrt(mean((yo_cv-ya_cv)**2)),3),"\n"))
    save( file= argv$off_rdata,
          # cv
          yo_cv, ya_cv, ya_cv_var, yb_cv, yidi_cv, ya_cv_henoi_varu, 
          ya_cv_pdf_par, Yb_cv, ya_cv_henoi_Dh, 
          VecS_cv, VecZ_cv,VecX_cv, VecY_cv, VecXorig_cv, VecYorig_cv
          )
  }
  print(paste("output saved on file",argv$off_rdata))
