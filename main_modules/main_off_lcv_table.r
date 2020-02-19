  cat("date;sourceId;x;y;z;yo;yb;ya_lcv;yidi_lcv;yidiv_lcv;dqc;\n",
      file=argv$off_lcv_table,append=F)
  cat(paste(argv$date_out,
            formatC(VecS,format="f",digits=0),
            formatC(VecX,format="f",digits=0),
            formatC(VecY,format="f",digits=0),
            formatC(VecZ,format="f",digits=0),
            formatC(yo,format="f",digits=2),
            formatC(yb,format="f",digits=2),
            formatC(ya_lcv,format="f",digits=2),
            formatC(yidi_lcv,format="f",digits=4),
            formatC(yidiv_lcv,format="f",digits=4),
            rep(0,length(VecS)),
            "\n",sep=";"),
      file=argv$off_lcv_table,append=T)
  print(paste("output saved on file",argv$off_lcv_table))
