
analyzeInclIndex <- function(plate.setup="splicing",plate.format="96-well"){
  
  ### Checking arguments ###
  # setting plate format
  if (plate.format == "96-well"){
    nwells <- 97
  }
  # Setting sample naming patter
  if(plate.setup == "splicing"){
    sample.pattern <- c("Cells","Gene","Primer","Treatment")
  }
  
  ## reseting graphics device
  graphics.off()
  
  ## choose file  and create dir ###
  exp.link <- file.choose()
  
  if (length(exp.link) < 1){
    stop("Something is wrong with the file choosen")
  }
  ## enter experiment name (maybe different from plate name)
  cat("\n","Enter ExpNr","\n") # prompt
  exp.nr=scan(n=1,what="character") # read 1 line from console
  # change working dir
  exp.dir <- dirname(exp.link)
  exp.dir <- file.path(exp.dir,"analysis")
  dir.create(exp.dir)
  setwd(exp.dir)
  
  ### reading data ###
  
  raw.df <- read.delim(file=exp.link,skip=10,nrows=nwells,stringsAsFactors=F)
  # removing trailing information
  raw.df <- raw.df[c(1:(grep(pattern="Slope",raw.df$Well)-1)),]
  # remove empty wells
  raw.df <- raw.df[grep(pattern=".*_.*",raw.df$Sample.Name),]
  # split sample name column
  df <- data.frame(raw.df$Well
                   ,colsplit(raw.df$Sample.Name,"_",sample.pattern)
                   ,raw.df$Detector.Name
                   ,as.numeric(raw.df$Ct)
                   ,raw.df$Threshold
  )
  
  ## remove primer effiency wells
  df <- df[df$Treatment != "PE",]
  header <- read.delim(file=exp.link,nrows=9,header=T)
  plate.name <- as.character(header[1,2])
  # workaround colnames
  aa <- colnames(df)
  aa[7] <- "Ct"
  colnames(df) <- aa
  
  
  # workaround CT as numeric
  df$Ct <- as.numeric(as.character(df$Ct))
  # H2O Kontrollen entfernen
  df <- df[df$Cells != "H2O",]  
  # alle roh daten in langes format bringen
  df.m <- melt(data=df
               ,id.vars=c("raw.df.Well","Cells","Gene","Primer","Treatment")
               ,measure.vars="Ct"
  )
  # first draft of results
  gg <- ggplot(df.m,aes(x=Treatment,y=value,fill=Primer))
  gg <- gg + geom_boxplot(stat="boxplot")
  gg <- gg + ylab("Ct") + xlab(NULL) + ggtitle(paste("EP No",exp.nr,"\n","Plate No",plate.name))# + coord_trans(y="log2")
  gg <- gg + scale_y_continuous(trans = "reverse") #+ coord_trans(y="log2")#scale_y_continuous(trans=log2_trans())
  gg <- gg + facet_grid(~Gene)
  #devSVGTips(file=paste(plate.name,"_plot1_rawCt.svg",sep=""))
  png(file=paste(exp.nr,"--",plate.name,"_plot1_rawCt.png",sep=""),width=600,height=600)
  plot(gg)  
  dev.off()
  
  # calculating means
  res.df <- ddply(df.m
                  ,c("Cells","Gene","Primer","Treatment")
                  ,summarize
                  ,mean = mean(value,na.rm=T)
                  ,sd = sd(value,na.rm=T)
  )
  
  
  
  # remove genes from plate, which are only checked for expression
  df.m <- df.m[df.m$Primer != "exp",]
  
  
  #treatment <- as.vector(rle(sort(as.character(res.df$Treatment)))[[2]])
  treatment <- c("NI","P12")
  
  # sortieren nach Treatment  
  res.ca <- dcast(data = df.m
                  ,formula = raw.df.Well + Cells + Gene + Primer ~ Treatment
                  ,value.var = "value"
  )
  # defining house keeping control and removing it
  # entfernen von HK und exp aus der einzeltabelle 
  for(i in treatment){
    HK <- (res.df[res.df$Primer == "HK" & res.df$Treatment == i,"mean"])[1]
    # dCT durch abziehen der GAPDH
    inx <- colnames(res.ca) == i
    res.ca[,inx] <- res.ca[,inx] - HK
  }
  res.df <- res.df[res.df$Primer != "HK",]
  df.m <- df.m[df.m$Primer != "HK" & df.m$Primer != "exp",]
  res.m <- melt(res.ca,id.vars=c("raw.df.Well","Cells","Gene","Primer"),measure.vars=treatment)
  res.m <- res.m[res.m$Primer != "HK",]
  #erstmal nur begrenzt auf NI und P12. P12dPAI wird entfernt
  res.m <- res.m[res.m$variable == "NI" | res.m$variable == "P12",]
  res.dCt <- na.omit(res.m)
  
  
  ####
  #dCt calculated in res.dCt
  ####
  
  # ddCt for each condition and primer pair
  res.mean <- ddply(res.dCt
                    ,c("Cells","Gene","Primer","variable")
                    ,summarize
                    ,mean = mean(value)
                    ,sd = sd(value,na.rm=T)
  )    
  
  
  # create index of Gene and treatment
  genes <- unique(res.mean$Gene)
  primer <- unique(res.mean$Primer)
  cells <- unique(res.mean$Cells) # später mal auflösen mit faktoren allgemein
  # eventuelle warnung bei mehr als zwei zuständen
  #   if (length(primer) > 2){
  #     print("There are more primers than expected")
  #   }
  
  out.df <- NULL
  res <- NULL
  for(c in cells){
    for(g in genes){
      res1 <- res.dCt[res.dCt$Gene == g,]
      for(tr in treatment) {
        res2 <- res1[res1$variable == tr,]
        res3 <- 2^-((res2[res2$Primer == "alt","value"])-(res.mean[res.mean$Gene == g 
                                                                   & res.mean$variable == tr 
                                                                   & res.mean$Primer == "con","mean"]))
        res4 <- data.frame(c,g,tr,res3) 
        res <- rbind(res,res4)
      }
      if (all(res$res3 > 1)){ # invertieren des ratios, wenn alt form die dominante ist
        res$res3 <- 1/res$res3
        res$g <- paste(res$g,"(inv)",sep="_")
      }
      
      out.df <- rbind(out.df,res)
      res <- NULL
      #### Achtung, Alt und Con könnten wechseln !!!
      
    }
  }
  colnames(out.df) <- c("Cells","Gene","Treatment","alt.Iso.frac")
  genes <- as.vector(rle(sort(as.character(out.df$Gene)))[[2]]) # recreate because of inversion
  ddCt.df <- out.df
  
  # Plotting
  gg <- ggplot(ddCt.df,aes(x=Gene,y=alt.Iso.frac,fill=Treatment))
  gg <- gg + geom_boxplot(position="dodge",stat="boxplot",show_guide = T)
  gg <- gg + ylab("Alternative isoform fraction") + xlab(NULL) + ggtitle(paste("EP No: ",exp.nr,"// // Cells: ",cells,"\n","Plate No: ",plate.name)) #+ scale_y_continuous(limits=c(0,1))
  gg <- gg + facet_grid(~Gene,scales="free_x")
  #gg <- gg + geom_errorbar(limits, position="dodge", width=0.25)
  png(file=paste(exp.nr,"--",plate.name,"_plot2_altIsoform_fraction.png",sep=""),width=600,height=600)
  plot(gg)  
  dev.off()
  
  
  # statistical testing
  
  results <- ddply(ddCt.df
                   ,c("Cells","Gene","Treatment")
                   ,summarize
                   ,mean.iso.fraction = mean(alt.Iso.frac)
                   ,sd.iso.fraction = sd(alt.Iso.frac,na.rm=T)
                   )
  ## berechnen der delta alternative isoform fraction dAIF
  ## Achtung NI und P12 fest eingebaut
  results$dAIF <- NA
  for (i in levels(results$Gene)){
    results$dAIF[results$Gene == i] <- results[results$Gene == i & results$Treatment == "P12" ,"mean.iso.fraction"] - results[results$Gene == i & results$Treatment == "NI" ,"mean.iso.fraction"]
  }
  
  for(g in genes){
    tt <-t.test(alt.Iso.frac ~ Treatment,data = ddCt.df,subset = (Gene == g))
    results[results$Gene == g,"p.value"] <- round(x=tt$p.value,digits=5)
    results[results$Gene == g,"conf.interval"] <- tt$conf.int
    results[results$Gene == g,"means.tt"] <- tt$estimate
  }
  
  # combine single observations and t.test results
  ddCt.df <- merge(x = ddCt.df,y = results,by = c("Cells","Gene","Treatment"),all.x=T)
  
  
  gg <- ggplot(ddCt.df,aes(x=Treatment,y=alt.Iso.frac,fill=Treatment))
  gg <- gg + geom_boxplot(position="dodge",stat="boxplot",show_guide = T)
  gg <- gg + ylab("Alternative isoform fraction") + xlab(NULL) + ggtitle(paste("EP No: ",exp.nr,"// // Cells: ",cells,"\n","Plate No: ",plate.name))#+ scale_y_continuous(limits=c(0,1))
  gg <- gg + facet_grid(Cells~Gene,scales="free_x",labeller= ) + geom_text(aes(x=1.5,y = -0.2,label=p.value))
  #gg <- gg + geom_hline(yintercept=1)
  #gg <- gg + geom_errorbar(limits, position="dodge", width=0.25)
  png(file=paste(exp.nr,"--",plate.name,"_plot3_altIsoform_fraction_pval.png",sep=""),width=600,height=600)
  plot(gg)  
  dev.off()
  
  ## create ouput table
  results$EP.No <- exp.nr
  results$plate.name <- plate.name
  results$Analysis.date <- Sys.time()
  results$Package.version <- packageVersion("kangoo")
  # write output
  write.table(results
              , paste(exp.nr,"--", plate.name,"_report_analysis.csv",sep="")
              , row.names = F
              , sep = ";")
  ## fürs ebook
  write.table(results
              , paste(exp.nr,"--", plate.name,"_report_analysis.xls",sep="")
              , row.names = F
              , sep = "\t")
  # 
  #write.csv(results,file.path("","","FILESERVER","ebook-mey-02","F-J","Glowinski-FG","summary_RTresults.csv"),append=T,row.names=F,col.names=F,quote=T,sep="\t")
  # hard link
  write.table(results
              ,"//FILESERVER/ebook-mey-02/F-J/Glowinski-FG/summary_RTresults.csv"
              , append = TRUE
              , row.names = FALSE
              , col.names = TRUE
              , quote = TRUE
              , sep = ";"
  )
  sink(paste(exp.nr, "--", plate.name, "_report_sessionInfo.txt", sep=""))
  print(Sys.time())
  print(sessionInfo())
  sink()
}

