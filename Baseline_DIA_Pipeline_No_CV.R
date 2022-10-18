''' 
For this pipeline you need to make a meta file with the samples and conditions
as seen in the proteus vignette. Load the libraries before anything. gg repel
overlaps is for the datapoints in the ggplots.
Proteus and limma refuse to flip the log FC so that has to be done manually.

IMPORTANT:
Any write functions MUST be changed so that it paths to your own computer.
!!!

Proteus only accepts the file if it is reading it as a txt file.
The reason that the file is uploaded into R before this is because of the CV step,
but that step is removed because we did not need it.
''' 
library(proteus)
library(limma)
library(tidyverse)
library(ggrepel)
library(dplyr)

options(ggrepel.max.overlaps = Inf)

#upload your proteins file here, this should be the output from the previous 
#python step. 

DIA_proteins <- read.csv("C:/Users/ASUS/Desktop/export_DIA_table.csv")
#get rid of the labels here, need to relabel.
DIA_proteins <- DIA_proteins[-1,]

proteus_pipeline_no_CV <- function(proteins.file, meta.file) {
  
  conditions <- c("Neurotypical", "Autism")
  
  colnames(proteins.file) <- gsub(x = colnames(proteins.file), pattern = "\\.", replacement = " ")
  
  #IMPORTANT: rename these to write to your own computer!!
  write.table(proteins.file, file = "C:/Users/ASUS/Desktop/proteins_file.txt", sep = "\t", row.names = F, col.names = T)
  main.proteins <- readProteinGroups("C:/Users/ASUS/Desktop/proteins_file.txt", meta.file)
  
  count <- plotCount(main.proteins)
  print(count)
  PCA <- plotPCA(main.proteins)
  print(PCA)
  
  main.proteins.med <- normalizeData(main.proteins)
  
  main.proteins.med.res <- limmaDE(main.proteins.med, sig.level = 0.05, conditions = conditions)
  volcano <- plotVolcano(main.proteins.med.res)
  
  print(volcano)
  
  main.proteins.med.res$impFC <- "NO"
  main.proteins.med.res$impFC[main.proteins.med.res$logFC > 0.6 | main.proteins.med.res$logFC < -0.6] <- "YES"
  
  main.proteins.med.res$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "DOWN" because we flip these.
  main.proteins.med.res$diffexpressed[main.proteins.med.res$logFC > 0.6 & main.proteins.med.res$P.Value < 0.05] <- "DOWN"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "UP" because we flip these.
  main.proteins.med.res$diffexpressed[main.proteins.med.res$logFC < -0.6 & main.proteins.med.res$P.Value < 0.05] <- "UP"
  
  colnames(main.proteins.med.res)[1] <- "Majority protein IDs"
  joined.main.proteins.med.res <- left_join(main.proteins.med.res, proteins.file, by = "Majority protein IDs", copy = T)
  
  joined.main.proteins.med.res$delabel <- NA
  joined.main.proteins.med.res$delabel[joined.main.proteins.med.res$diffexpressed != "NO"] <- joined.main.proteins.med.res$`Genes`[joined.main.proteins.med.res$diffexpressed != "NO"]
  joined.main.proteins.med.res$delabel[joined.main.proteins.med.res$impFC == "YES"] <- joined.main.proteins.med.res$`Genes`[joined.main.proteins.med.res$impFC == "YES"]
  joined.main.proteins.med.res$delabel[joined.main.proteins.med.res$P.Value <= 0.05] <- joined.main.proteins.med.res$`Genes`[joined.main.proteins.med.res$P.Value <= 0.05]
  
  "Log2 FC" <- -1*(joined.main.proteins.med.res$logFC)
  
  final.plot <- ggplot(data=joined.main.proteins.med.res, aes(x=`Log2 FC`, y=-log10(P.Value), label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_hline(yintercept=-log10(0.05), col="red")
  
  print(final.plot)
  
  #IMPORTANT: rename this file to path to your own computer!!
  write.csv(joined.main.proteins.med.res, file = "c:/Users/ASUS/Desktop/placeholder_results.csv")
}

'''
Plots: can generate violin plots with this if you take the results file, delete
everything except the labels, LFQ intensities and genes, then transpose, so
the genes are the column names, the LFQ intensities are the row names. Can do
this in excel.
'''

THE_VIOLIN <- function(boxplot_file) {
  for(i in colnames(boxplot_file)) {
    print(i)
    title <- paste(i, "Expression", sep = " ")
    X = "Label"
    if(i != "Genes" & i != "Label") {
      g <- ggplot(data = boxplot_file, aes_string(x = X, y = i, fill = X)) +
        ggtitle(title) +
        ylab( "log2(LFQ Intensity)") + 
        geom_violin() +
        geom_boxplot(width = 0.1, fill = "white") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
      print(g)
    }
  }
}