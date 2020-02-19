  cat("date;sourceId;x;y;z;yo_cv;yb_cv;ya_cv;yidi_cv;dqc;\n",
      file=argv$off_cv_table,append=F)
  cat(paste(argv$date_out,
            formatC(VecS_cv,format="f",digits=0),
            formatC(VecX_cv,format="f",digits=0),
            formatC(VecY_cv,format="f",digits=0),
            formatC(VecZ_cv,format="f",digits=0),
            formatC(yo_cv,format="f",digits=2),
            formatC(yb_cv,format="f",digits=2),
            formatC(ya_cv,format="f",digits=2),
            formatC(yidi_cv,format="f",digits=4),
            rep(0,length(VecS_cv)),
            "\n",sep=";"),
      file=argv$off_cv_table,append=T)
  print(paste("output saved on file",argv$off_cv_table))
