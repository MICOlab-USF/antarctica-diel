library(SpiecEasi)
library(qiime2R)
library(tidyverse)
library(phyloseq)

#Meta-data----------------------------------------------------------------------
samplemap<- read.csv("UrbanOKI_SpiecEasi/meta.csv")
samplemap$date <- dmy(samplemap$date)
samplemap$Month <- as.factor(as.numeric(samplemap$Month))

samplemap <- samplemap %>% 
  mutate(Day=day(date)) %>% 
  unite(Name, c(LandUse, Number, Position, Month, Day), sep = "_", remove = FALSE)

samplemap$Position <- factor(samplemap$Position, levels =c("S", "C", "N"))
samplemap$LandUse <- factor(samplemap$LandUse, levels = c("U", "R"))

row.names(samplemap)<- samplemap$SampleID
META <- sample_data(samplemap)

#Import data -------------------------------------------------------------------

ps <- qza_to_phyloseq(features="UrbanOKI_SpiecEasi/merged_table.qza")

taxonomy <- read.csv("UrbanOKI_SpiecEasi/16taxonomy.csv", stringsAsFactors = FALSE, header = FALSE)
names(taxonomy) <- c("row", "tax", "Confidence") #change the headers (column names)
row.names(taxonomy) <-taxonomy[[1]] #move the feature ID column to become the row names
taxonomy <- taxonomy[,(-1)] #delete the feature ID  column 
taxonomy <-  separate(taxonomy, tax, c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species", "D7", "D8", "D9", "D10", "D11", "D12", "D13", "D14"), sep = ";", fill = "right")
taxonomy <- taxonomy[,c(1:7)]

taxonomy$D0 <- with(taxonomy, ifelse(Order == "D_3__Chloroplast", "Chloroplast", "Bacteria"))

col_order<- c("D0", "Domain","Phylum", "Class", "Order", "Family", "Genus", "Species" )
taxonomy<- taxonomy[, col_order]

taxmat <- as.matrix(taxonomy)
TAX = tax_table(taxmat)

ps = merge_phyloseq(ps, TAX, META) 

#Test SpiecEasi

UrbanOki <- spiec.easi(ps, method='mb', lambda.min.ratio=1e-2,
                           nlambda=20, pulsar.params=list(rep.num=50))
ig2.mb <- adj2igraph(getRefit(UrbanOki),  vertex.attr=list(name=taxa_names(ps)))
plot_network(ig2.mb, UrbanOki, type='taxa')
