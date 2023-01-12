  cat("date;sourceId;x;y;z;yto_lcv;ytb_lcv;yta_lcv;yidi_lcv;yidiv_lcv;dqc;\n",
      file=argv$off_lcvt_table,append=F)
  cat(paste(argv$date_out,
            formatC(VecS,format="f",digits=0),
            formatC(VecX,format="f",digits=0),
            formatC(VecY,format="f",digits=0),
            formatC(VecZ,format="f",digits=0),
            formatC(yto_lcv,format="f",digits=2),
            formatC(ytb_lcv,format="f",digits=2),
            formatC(yta_lcv,format="f",digits=2),
            formatC(yidi_lcv,format="f",digits=4),
            formatC(yidiv_lcv,format="f",digits=4),
            rep(0,length(VecS)),
            "\n",sep=";"),
      file=argv$off_lcvt_table,append=T)
  print(paste("output saved on file",argv$off_lcvt_table))
