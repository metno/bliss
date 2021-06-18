  if (!exists("yidi"))  yidi  <- rep( -999, length(yo))
  if (!exists("yidiv")) yidiv <- rep( -999, length(yo))
  if (argv$mode=="ensigap") {
    cat("date;sourceId;prid;x;y;z;yo;yb;ya;ya_var;yav;yav_var;yidi;yidiv;dqc;ya_gamma_shape;ya_gamma_rate;ya_henoi_varu;yav_gamma_shape;yav_gamma_rate;yav_henoi_varu;\n",
        file=argv$off_y_table,append=F)
    cat(paste(argv$date_out,
              formatC(VecS,format="f",digits=0),
              formatC(prId,format="f",digits=0),
              formatC(VecX,format="f",digits=0),
              formatC(VecY,format="f",digits=0),
              formatC(VecZ,format="f",digits=0),
              formatC(yo,format="f",digits=2),
              formatC(yb,format="f",digits=2),
              formatC(ya,format="f",digits=2),
              formatC(ya_var,format="f",digits=4),
              formatC(ya_cv,format="f",digits=2),
              formatC(ya_cv_var,format="f",digits=4),
              formatC(yidi,format="f",digits=4),
              formatC(yidi_cv,format="f",digits=4),
              rep(0,length(VecS)),
              formatC(ya_pdf_par[1],format="f",digits=4),
              formatC(ya_pdf_par[2],format="f",digits=4),
              formatC(ya_henoi_varu,format="f",digits=4),
              formatC(ya_cv_pdf_par[1],format="f",digits=4),
              formatC(ya_cv_pdf_par[2],format="f",digits=4),
              formatC(ya_cv_henoi_varu,format="f",digits=4),
              "\n",sep=";"),
        file=argv$off_y_table,append=T)
  } else {
    cat("date;sourceId;x;y;z;yo;yb;ya;yav;yidi;yidiv;dqc;\n",
        file=argv$off_y_table,append=F)
    cat(paste(argv$date_out,
              formatC(VecS,format="f",digits=0),
              formatC(VecX,format="f",digits=0),
              formatC(VecY,format="f",digits=0),
              formatC(VecZ,format="f",digits=0),
              formatC(yo,format="f",digits=2),
              formatC(yb,format="f",digits=2),
              formatC(ya,format="f",digits=2),
              formatC(yav,format="f",digits=2),
              formatC(yidi,format="f",digits=4),
              formatC(yidiv,format="f",digits=4),
              rep(0,length(VecS)),
              "\n",sep=";"),
        file=argv$off_y_table,append=T)
  }
  print(paste("output saved on file",argv$off_y_table))
