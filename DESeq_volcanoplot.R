library(DESeq2)
library(ggplot2)

load("CSV_files/Subsections_numreads.RData",verbose = TRUE)

#####
# Input_data <- Bact.Sub

Input_data <- df.count # I used an artificial numread vector; count was achieved from TPM that I multiplied with 10^6.
names(Input_data)[25] <- "SequenceID"

# # For loop that averages the each assembly group into 8 time points rather that the 8 time points with replicates
# for(i in 1:24){
#   names(Input_data)[i] <- paste(meta.data$SampleID[meta.data$SampleName == names(Input_data)[i]],"quant",sep = "_")
# }
# 
# vec_assembly <- unique(meta.data$AssemblyGroup)
# 
# for(i in 1:length(vec_assembly)){
#   if(i == 1){
#     names(Input_data) <- unlist(lapply(names(Input_data),
#                                       function(x){unlist(strsplit(x,
#                                                                   split = "_",
#                                                                   fixed = TRUE))[1]}))
#     
#     df_meta <- meta.data[,c("SampleID","SampleName","AssemblyGroup","Day.num","ToD")]
#     
#     df_means <- Input_data$SequenceID
#     par_vec <- NA
#     time1 <- Sys.time()
#   }
#   
#   df_temp <- Input_data[,meta.data$SampleID[meta.data$AssemblyGroup == vec_assembly[i]]]
#   
#   meta_temp <- df_meta %>% 
#     distinct(AssemblyGroup,.keep_all = TRUE) %>% 
#     filter(AssemblyGroup == vec_assembly[i])
#   
#   vec <- rowMeans(df_temp)
#   # vec <- rowSums(df_temp)
#   df_means <- cbind.data.frame(df_means,vec)
#   
#   names(df_means)[i+1] <- vec_assembly[i]
# }
# 
# # creating factors for the ToD in order to analyse the time points chronologically
# df_meta <- df_meta %>% 
#   #distinct(AssemblyGroup,.keep_all = TRUE) %>% 
#   mutate(ToD = factor(ToD,levels = c("morning","afternoon","evening","night")))

#####

# countData <- round(df_means[,-1])
countData <- cbind.data.frame("SequenceID" = df_means[,1],countData)

countData <- round(df.count[,-25])
countData <- cbind.data.frame("SequenceID" = df.count[,25],countData)
# countData <- cbind.data.frame("SequenceID" = Bact.Sub$SequenceID,countData)

head(countData)

# can filter the data in order to perform a differential expression analysis on two time points of choice (or all).
GroupOI <- c("morning","evening","afternoon","night")

dds <- DESeqDataSetFromMatrix(countData=countData[,c("SequenceID",df_meta$SampleName[df_meta$ToD %in% GroupOI])],
                              colData=df_meta[df_meta$ToD %in% GroupOI,], 
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

# par(mfrow=c(2,3))
# 
# plotCounts(dds, gene="Sequence_15334987", intgroup="ToD")
# plotCounts(dds, gene="Sequence_14032611", intgroup="ToD")
# plotCounts(dds, gene="Sequence_1192986", intgroup="ToD")
# plotCounts(dds, gene="Sequence_12328187", intgroup="ToD")
# plotCounts(dds, gene="Sequence_10000657", intgroup="ToD")
# plotCounts(dds, gene="Sequence_10000667", intgroup="ToD")

# par(mfrow=c(2,3))
# 
# plotCounts(dds, gene="Sequence_10746420", intgroup="Day.num")
# plotCounts(dds, gene="Sequence_11903479", intgroup="Day.num")
# plotCounts(dds, gene="Sequence_14032611", intgroup="Day.num")
# plotCounts(dds, gene="Sequence_14169545", intgroup="Day.num")
# plotCounts(dds, gene="Sequence_16277791", intgroup="Day.num")
# plotCounts(dds, gene="Sequence_16116212", intgroup="Day.num")

par(mfrow=c(2,3))

plotCounts(dds, gene="ko:K03883", intgroup="ToD")
plotCounts(dds, gene="ko:K07810", intgroup="ToD")
plotCounts(dds, gene="ko:K12976", intgroup="ToD")
plotCounts(dds, gene="ko:K02037", intgroup="ToD")
plotCounts(dds, gene="ko:K07127", intgroup="ToD")
plotCounts(dds, gene="ko:K03046", intgroup="ToD")

# par(mfrow=c(2,3))
# 
# plotCounts(dds, gene="ko:K00404", intgroup="Day.num")
# plotCounts(dds, gene="ko:K09459", intgroup="Day.num")
# plotCounts(dds, gene="ko:K01113", intgroup="Day.num")
# plotCounts(dds, gene="ko:K09930", intgroup="Day.num")
# plotCounts(dds, gene="ko:K21594", intgroup="Day.num")
# plotCounts(dds, gene="ko:K02902", intgroup="Day.num")

# Volcano Plot ----------------------------------------------------------------
pAdj_thresh <- 0.01

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res,
            padj<pAdj_thresh ),
     points(log2FoldChange,-log10(pvalue),
            pch=19,
            col="blue"))

with(subset(res,
            padj<pAdj_thresh & abs(log2FoldChange)>2 & !is.na(padj<pAdj_thresh)),
     points(log2FoldChange, -log10(pvalue),
            pch=19,
            col="red"))

res[res$padj<pAdj_thresh & !is.na(res$padj<pAdj_thresh),]
