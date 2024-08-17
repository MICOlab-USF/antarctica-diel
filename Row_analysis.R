create_plot = function(x_vec,y_vec,lm_vec){
  plot(x_vec,y_vec)
  abline(lm_vec)
}

Eukaryota_short <- Eukaryota.merged[1:10000,]

for(i in 1:nrow(Eukaryota_short)){
  if(i == 1){
    time1 <- Sys.time()
    names(Eukaryota_short) <- unlist(lapply(names(Eukaryota_short),
                                            function(x){unlist(strsplit(x,
                                                                        split = "_",
                                                                        fixed = TRUE))[1]}))
    df_meta <- meta.data[,c("SampleID","par","AssemblyGroup")]
    SampleID <- names(Eukaryota_short[i,2:25])
    
    seq_vec <- NA
    n <- 1
  }
  vec <- as.numeric(Eukaryota_short[i,2:25])
  
  
  df.temp <- merge(cbind.data.frame(vec,SampleID),df_meta,by = "SampleID")
  
  df.summ <- summarySE(data = df.temp, "vec",
                       groupvars = c("AssemblyGroup"))
  
  dfOut <- merge(df.summ,df_meta,by = "AssemblyGroup") %>% 
    distinct(AssemblyGroup, .keep_all = TRUE)
  
  lm_vec <- lm(vec~par,dfOut)
  
  if(summary(lm_vec)$r.squared > .90 & sum(dfOut$vec > 0) > 2){
    seq_vec[n] <- Eukaryota_short[i,1]
    n <- n+1
    
    create_plot(dfOut$par,dfOut$vec,lm_vec)
  }
  
  if(i == nrow(Eukaryota_short)){
    print(Sys.time() - time1)
  }
}
