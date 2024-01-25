Create Voicano plots from differential expressed genes generated from RNA seq
Cutoff values of fold change need to be defined first. Either FDR or p values may be used.

############################ Volcano plot (Fold Change >= 2 or <= -2) #################################
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualization. dplyr, for data manipulation.
library(RColorBrewer) # for a colorful plot
library(ggrepel) # for nice annotations

# Set input path
path <- "The path having 'igenomes.genes.combined.txt'"
setwd(path)

# Import DGE results
df <- read.delim(paste0(path, 'igenomes.genes.combined.txt'))

# Check files in the path
list.files(path)

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
df$diffexpressed <- "NO"

# if log2Foldchange > 1 and FDR < 0.05, set as "UP"
df$diffexpressed[df$IgG_Control.VS.Lep_ETP.logFC > 1 & df$IgG_Control.VS.Lep_ETP.FDR < 0.05] <- "UP"
# if log2Foldchange < -1 and FDR < 0.05, set as "DOWN"
df$diffexpressed[df$IgG_Control.VS.Lep_ETP.logFC < -1 & df$IgG_Control.VS.Lep_ETP.FDR < 0.05] <- "DOWN"

head(df[order(df$IgG_Control.VS.Lep_ETP.FDR) & df$diffexpressed == 'DOWN', ])

# Sort the data frame based on -log10(FDR)
df <- df[order(-log10(df$IgG_Control.VS.Lep_ETP.FDR)), ]

# Create a new column "delabel" that will contain the name of the top 45 differential expressed genes (NA in case they are not)
# These top 45 genes were further filtered by IgG_Control.VS.Lep_ETP.logFC (>=1 or <=-1)
df$delabel <- ifelse(df$FEATURE_NAME %in% head(df[order(df$IgG_Control.VS.Lep_ETP.FDR), "FEATURE_NAME"], 45) & abs(df$IgG_Control.VS.Lep_ETP.logFC) >= 1, df$FEATURE_NAME, NA)

# Biostatsquid theme
theme_set(theme_classic(base_size = 10) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,10,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(10,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

# Create a basic volcano plot
ggplot(data = df, aes(x = IgG_Control.VS.Lep_ETP.logFC, y = -log10(IgG_Control.VS.Lep_ETP.FDR), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "red"), # to set the colors of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 80), xlim = c(-8, 6)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customize the breaks in the x axis
  geom_text_repel(max.overlaps = Inf) # To show all labels

##################################################################################################
# export all genes with FDR < 0.05 and log2FC >=1 or <=-1
library(openxlsx)
df_filtered <- filter(df, IgG_Control.VS.Lep_ETP.FDR < 0.05 & abs(IgG_Control.VS.Lep_ETP.logFC) >= 1)
write.xlsx(df_filtered, 'filtered_gene_list.xlsx')