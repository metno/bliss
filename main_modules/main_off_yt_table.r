  cat("date;sourceId;x;y;z;yto;ytb;yta;ytav;yidi;yidiv;dqc;\n",
      file=argv$off_yt_table,append=F)
  cat(paste(argv$date_out,
            formatC(VecS,format="f",digits=0),
            formatC(VecX,format="f",digits=0),
            formatC(VecY,format="f",digits=0),
            formatC(VecZ,format="f",digits=0),
            formatC(yto,format="f",digits=2),
            formatC(ytb,format="f",digits=2),
            formatC(yta,format="f",digits=2),
            formatC(ytav,format="f",digits=2),
            formatC(yidi,format="f",digits=4),
            formatC(yidiv,format="f",digits=4),
            rep(0,length(VecS)),
            "\n",sep=";"),
      file=argv$off_yt_table,append=T)
  print(paste("output saved on file",argv$off_yt_table))
