create_plot = function(x_vec,y_vec,lm_vec){
  plot(x_vec,y_vec)
  abline(lm_vec)
}


Input_data <- Eukaryota.merged[1:1000,]
# Input_data <- read.csv("CSV_files/Bacteria.subsection.TPM.csv")

vec_assembly <- unique(meta.data$AssemblyGroup)

for(i in 1:length(vec_assembly)){
  if(i == 1){
    names(Input_data) <- unlist(lapply(names(Input_data),
                                       function(x){unlist(strsplit(x,
                                                                   split = "_",
                                                                   fixed = TRUE))[1]}))
    
    df_meta <- meta.data[,c("SampleID","par","AssemblyGroup")]
    
    df_means <- Input_data$SequenceID
    par_vec <- NA
    time1 <- Sys.time()
  }
  
  df_temp <- Input_data[,meta.data$SampleID[meta.data$AssemblyGroup == vec_assembly[i]]]
  
  meta_temp <- df_meta %>% 
    distinct(AssemblyGroup,.keep_all = TRUE) %>% 
    filter(AssemblyGroup == vec_assembly[i])
  
  par_vec[i] <- meta_temp$par
  
  vec <- rowMeans(df_temp)
  df_means <- cbind.data.frame(df_means,vec)
  
  names(df_means)[i+1] <- vec_assembly[i]
  
  if(i == length(vec_assembly)){
    # df_means <- as.data.frame(t(df_means))
    names(df_means)[1] <- "SequenceID"
    
    df_means$r.squared <- apply(df_means[,2:9],1,function(row){summary(lm(row~par_vec))$r.squared})
    df_means$n_nonzero <- apply(df_means[,2:9],1,function(row){sum(row > 0)})
    
    seq_df <- df_means[df_means$r.squared >.95 & df_means$n_nonzero > 2,]
    
    print(Sys.time() - time1)
  }
}










# rowMeans(Input_data[,meta.data$SampleID[meta.data$AssemblyGroup == vec_assembly[1]]])

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



seq_df2 <- read.csv("CSV_files/Eukaryota_row_reg.csv")
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

