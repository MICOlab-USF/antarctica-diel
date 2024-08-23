create_plot = function(x_vec,y_vec,lm_vec){
  plot(x_vec,y_vec)
  abline(lm_vec)
}


Input_data <- Eukaryota.merged[1:1000,]
# Input_data <- read.csv("CSV_files/Bacteria.subsection.TPM.csv")

for(i in 1:nrow(Input_data)){
  if(i == 1){
    time1 <- Sys.time()
    names(Input_data) <- unlist(lapply(names(Input_data),
                                       function(x){unlist(strsplit(x,
                                                                   split = "_",
                                                                   fixed = TRUE))[1]}))
    df_meta <- meta.data[,c("SampleID","par","AssemblyGroup")]
    SampleID <- names(Input_data[i,2:25])

    seq_vec <- NA
    seq_df <- cbind.data.frame(unique(df_meta$AssemblyGroup),unique(df_meta$par))
    names(seq_df) <- c("AssemblyGroup","par")
    n <- 1
  }
  vec <- as.numeric(Input_data[i,2:25])
  
  
  df.temp <- merge(cbind.data.frame(vec,SampleID),df_meta,by = "SampleID")
  
  df.summ <- summarySE(data = df.temp, "vec",
                       groupvars = c("AssemblyGroup"))
  
  dfOut <- merge(df.summ,df_meta,by = "AssemblyGroup") %>% 
    distinct(AssemblyGroup, .keep_all = TRUE)
  
  lm_vec <- lm(vec~par,dfOut)
  
  if(summary(lm_vec)$r.squared > .95 & sum(dfOut$vec > 0) > 2){
    seq_vec[n] <- Input_data[i,1]
    
    vec1 <- dfOut$vec
    
    seq_df <- data.frame(cbind.data.frame(seq_df,vec1))
    names(seq_df)[names(seq_df) == "vec1"] <- seq_vec[n]
    
    n <- n+1
    
    create_plot(dfOut$par,dfOut$vec,lm_vec)
  }
  
  if(i%%1000 == 0){
    print(paste("Row", i))
    print(Sys.time() - time1)
  }
  
  if(i == nrow(Input_data)){
    print("Finished")
    print(Sys.time() - time1)
  }
}



seq_df <- read.csv("CSV_files/Eukaryota_row_reg.csv")
filter <- Eukaryota.merged$SequenceID %in% names(seq_df)
Eukaryota_light <- Eukaryota.merged[filter,]
dfOut <- Eukaryota_light



# seq_df <- read.csv("CSV_files/Bacteria_row_reg.csv")
# Bacteria.TFA <- read.csv("CSV_files/Bacteria.TFA_NoRep.csv")
# 
# filter <- Bacteria.TFA$SequenceID %in% names(seq_df)
# Bacteria_light <- Bacteria.TFA[filter,]
# 
# dfOut <- Bacteria_light

# Create a bar plot
ggplot(dfOut, aes(x = Class))+
  geom_bar(aes(y = after_stat(count)/sum(after_stat(count))))+
  labs(title = "Class Distribution",
       x = "Class",
       y = "Count")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# ggsave("Light_correlated_Class.png",dpi = 300,height = 5, width = 6)

# Create a bar plot
ggplot(Eukaryota.merged, aes(x = Class))+
  geom_bar(aes(y = after_stat(count)/sum(after_stat(count))))+
  labs(title = "Class Distribution",
       x = "Class",
       y = "Count")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# ggsave("All_Class.png",dpi = 300,height = 5, width = 6)
