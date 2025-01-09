library(DESeq2)
library(ggplot2)

load("CSV_files/Subsections_numreads.RData",verbose = TRUE)

#####
Input_data <- Bact.Sub
vec_assembly <- unique(meta.data$AssemblyGroup)

for(i in 1:length(vec_assembly)){
  if(i == 1){
    names(Input_data) <- unlist(lapply(names(Input_data),
                                      function(x){unlist(strsplit(x,
                                                                  split = "_",
                                                                  fixed = TRUE))[1]}))
    
    df_meta <- meta.data[,c("SampleID","AssemblyGroup","Day.num","ToD")]
    
    df_means <- Input_data$SequenceID
    par_vec <- NA
    time1 <- Sys.time()
  }
  
  df_temp <- Input_data[,meta.data$SampleID[meta.data$AssemblyGroup == vec_assembly[i]]]
  
  meta_temp <- df_meta %>% 
    distinct(AssemblyGroup,.keep_all = TRUE) %>% 
    filter(AssemblyGroup == vec_assembly[i])
  
  vec <- rowMeans(df_temp)
  df_means <- cbind.data.frame(df_means,vec)
  
  names(df_means)[i+1] <- vec_assembly[i]
}

df_meta <- df_meta %>% 
  distinct(AssemblyGroup,.keep_all = TRUE)

#####

countData <- round(df_means[,-1])
countData <- cbind.data.frame("SequenceID" = Bact.Sub$SequenceID,countData)


head(countData)

dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=df_meta, 
                              design=~ToD, tidy = TRUE)

dds

# Now we’re ready to run DESEQ function ---------------------------------------

dds <- DESeq(dds)

# Take a look at the results table --------------------------------------------

res <- results(dds)

head(results(dds, tidy = TRUE))

# Summary of differential gene expression -------------------------------------

summary(res) #summary of results

# Sort summary list by p-value ------------------------------------------------

res <- res[order(res$padj),]
head(res)

# Volcano Plot ----------------------------------------------------------------

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))