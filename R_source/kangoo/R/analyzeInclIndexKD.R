
analyzeInclIndexKD <- function(plate.setup="splicing_knockdown",plate.format="96-well"){
  
  ### Checking arguments ###
  # setting plate format
  if (plate.format == "96-well"){
    nwells <- 97
  }
  if (plate.format == "unknown"){
    nwells <- 999
  }
  # Setting sample naming patter
  if(plate.setup == "splicing"){
    sample.pattern <- c("Cells","Gene","Primer","Treatment")
  }
  if (plate.setup == "splicing_knockdown") {
    sample.pattern <- c("Cells", "Gene","Knockdown", "Primer", "Treatment")
  }
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
  if(file.exists(exp.dir)){
    cat("Folder already exists. Using the existing folder for output.")
  } else {
  dir.create(exp.dir)
  }
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
  header <- read.delim(file = exp.link, nrows = 9, header = T)
  plate.name <- as.character(header[1, 2])
  # workaround colnames
  aa <- colnames(df)
  aa[8] <- "Ct"
  colnames(df) <- aa
  
  
  # workaround CT as numeric
  df$Ct <- as.numeric(as.character(df$Ct))
  # H2O Kontrollen entfernen
  df <- df[df$Cells != "H2O", ]
  df.m <- melt(data = df, id.vars = c("raw.df.Well", "Cells", 
                                      "Gene", "Primer", "Treatment","Knockdown"), measure.vars = "Ct")
  # Plot mit CT werten
  gg <- ggplot(df.m, aes(x = Treatment, y = value, fill = Primer))
  gg <- gg + geom_boxplot(stat = "boxplot")
  gg <- gg + ylab("Ct") + xlab(NULL) + ggtitle(paste("EP No",exp.nr,"\n","Plate No",plate.name))# + coord_trans(y="log2")
  gg <- gg + scale_y_continuous(trans = "reverse") #+ coord_trans(y="log2")#scale_y_continuous(trans=log2_trans())
  gg <- gg + facet_grid( Gene ~ Knockdown)
  png(file = paste(exp.nr,"--",plate.name, "_plot1_rawCt.png", sep = ""), 
      width = 600, height = 600)
  plot(gg)
  dev.off()
  
  # calculating means
  res.df <- ddply(df.m
                  , c("Cells", "Gene", "Primer", "Treatment","Knockdown")
                  , summarize
                  , mean = mean(value, na.rm = T)
                  , sd = sd(value, na.rm = T))
  
  df.m <- df.m[df.m$Primer != "exp", ] ## OK bis hier
  treatment <- c("NI", "P12")
  knockdown.genes <- unique(res.df$Knockdown)
  res.df2 <- NULL
  for (k in knockdown.genes) {
    df.kd <- df.m[df.m$Knockdown == k,]
    res.ca <- dcast(data = df.kd, fun.aggregate = mean, formula = raw.df.Well + Cells + 
                      Gene + Primer  ~ Treatment , value.var = "value")
    
    # abziehen der Referenz (GAPDH, HK)
    for (i in treatment) {
      HK <- (res.df[res.df$Primer == "HK" & 
                      res.df$Treatment == i & 
                      res.df$Knockdown == k
                    , "mean"])[1]
      inx <- colnames(res.ca) == i
      res.ca[, inx] <- res.ca[, inx] - HK
    }
    res.ca <- cbind(res.ca,k)
    res.df2 <- rbind(res.df2,res.ca) ## OK bis hier
  }
  
  res.df <- res.df[res.df$Primer != "HK", ]
  df.m <- df.m[df.m$Primer != "HK" & df.m$Primer != "exp", ] ## OK, bis hier
  res.m <- melt(res.df2, id.vars = c("raw.df.Well", "Cells", 
                                     "Gene", "Primer","k"), measure.vars = treatment)
  res.m <- res.m[res.m$Primer != "HK", ]
  res.m <- res.m[res.m$variable == "NI" | res.m$variable == 
                   "P12", ]
  res.dCt <- na.omit(res.m) ## ok, bis hier
  res.mean <- ddply(res.dCt, c("Cells", "Gene", "Primer", "k","variable"), 
                    summarize, mean = mean(value), sd = sd(value, na.rm = T))
  genes <- unique(res.mean$Gene)
  primer <- unique(res.mean$Primer)
  cells <- unique(res.mean$Cells)
  out.df <- NULL
  res.g <- NULL
  for (c in cells) {
    for (k in knockdown.genes){
      for (g in genes) {
        res1 <- res.dCt[res.dCt$Cells == c &
                          res.dCt$k == k &
                          res.dCt$Gene == g
                        , ]
        
        for (tr in treatment) {
          if (nrow(res1[res1$variable == tr, ]) > 0){
          res2 <- res1[res1$variable == tr, ]
          res3 <- 2^-((res2[res2$Primer == "alt", "value"]) - 
                        (res.mean[res.mean$Gene == g & res.mean$variable == 
                                    tr & res.mean$Primer == "con" & res.mean$k == k , "mean"]))
          res <- data.frame(c, g, k, tr, res3)
          res.g <- rbind(res.g,res)
          res <- NULL
        }
        }
        if(plate.format != "unknown"){
          if (all(res.g$res3 > 1)){ # invertieren des ratios, wenn alt form die dominante ist
            res.g$res3 <- 1/res.g$res3
            res.g$g <- paste(res.g$g,"(inv)",sep="_")
            
          }
        }
        out.df <- rbind(out.df,res.g)
        res.g <- NULL
        #### Achtung, Alt und Con könnten wechseln !!!
      }
      
    }
  }
  
  colnames(out.df) <- c("Cells", "Gene", "Knockdown", "Treatment", "alt.Iso.frac")
  ddCt.df <- out.df
  
  # Plotting
  gg <- ggplot(ddCt.df,aes(x = Treatment, y=alt.Iso.frac,fill=Treatment))
  gg <- gg + geom_boxplot(position="dodge",stat="boxplot",show_guide = T)
  gg <- gg + ylab("Alternative isoform fraction") + xlab(NULL) + ggtitle(paste("EP No: ",exp.nr,"// // Cells: ",cells,"\n","Plate No: ",plate.name)) #+ scale_y_continuous(limits=c(0,1))
  gg <- gg + facet_grid(Gene~Knockdown,scales="free_y")
  #gg <- gg + geom_errorbar(limits, position="dodge", width=0.25)
  png(file=paste(exp.nr,"--",plate.name,"_plot2_altIsoform_fraction.png",sep=""),width=600,height=600)
  plot(gg)  
  dev.off()
  
  
  # statistical testing
  
  results <- ddply(ddCt.df
                   ,c("Cells","Gene","Treatment","Knockdown")
                   ,summarize
                   ,mean.iso.fraction = mean(alt.Iso.frac)
                   ,sd.iso.fraction = sd(alt.Iso.frac,na.rm=T)
  )
  ## berechnen der delta alternative isoform fraction dAIF
  ## Achtung NI und P12 fest eingebaut
  results$dAIF <- NA
  for (i in unique(results$Gene)){
    for (k in knockdown.genes){
      results$dAIF[results$Gene == i & results$Knockdown == k] <- results[results$Gene == i & results$Knockdown == k & results$Treatment == "P12" ,"mean.iso.fraction"] - results[results$Gene == i & results$Knockdown == k & results$Treatment == "NI" ,"mean.iso.fraction"]
    }
  }
  for(g in unique(results$Gene)){
    for (k in knockdown.genes){
      tt <- t.test(alt.Iso.frac ~ Treatment, data = ddCt.df, 
                   subset = (Gene == g & Knockdown == k))
      results[results$Gene == g & results$Knockdown == k, "p.value"] <- round(x = tt$p.value
                                                                              , digits = 5)
      results[results$Gene == g & results$Knockdown == k, "conf.interval"] <- tt$conf.int
      results[results$Gene == g & results$Knockdown == k, "means.tt"] <- tt$estimate
    }
  }
  ## Falls der p Wert kleiner als 1 * 10^-5 ist, wird "p < 0,00001" ausgegeben
  p.small <- results$p.value < 0.00001
  results$p.value[p.small] <- "p < 0,00001"
  
  
  ddCt.df <- merge(x = ddCt.df, y = results, by = c("Cells", 
                                                    "Gene", "Treatment", "Knockdown"), all.x = T)
  gg <- ggplot(ddCt.df, aes(x = Treatment, y = alt.Iso.frac, 
                            fill = Treatment))
  gg <- gg + geom_boxplot(position = "dodge", stat = "boxplot", 
                          show_guide = T)
  gg <- gg + ylab("Alternative isoform fraction") + xlab(NULL)
  gg <- gg + facet_grid(Gene ~ Knockdown, scales = "free_y") 
    geom_text(aes(x = 1.5, y = 0.05 + max(ddCt.df$alt.Iso.frac), 
                  label = as.factor(p.value)))
  png(file=paste(exp.nr,"--",plate.name,"_plot3_altIsoform_fraction_pval.png",sep=""),width=600,height=600)  
  plot(gg)
  dev.off()
  
  ## create output table
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

