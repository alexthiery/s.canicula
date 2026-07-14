###########################################################################
# Differential expression, SOM clustering, GO enrichment and regulatory
# network analysis of the Scyliorhinus canicula dental lamina RNA-seq.
#
# Run from the repository root (e.g. Rscript DE_analysis.R).
# Outputs are written to suppl_files/dea_output/.
###########################################################################

# Load packages
library(tidyverse)
library(kohonen)
library(pheatmap)
library(DESeq2)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(gridExtra)
library(GOstats)
library(GSEABase)


###########################################################################
# Setup
###########################################################################

# Run this script with the repository root as the working directory.
output_path = "./suppl_files/dea_output/"
dir.create(output_path, recursive = T)

# Tissue colours used throughout (BHTB, TTJ, SL, ET, LT).
tissue_colours <- c('BHTB' = '#288fc9', 'TTJ' = '#cf29ef', 'SL' = '#bc1212', 'ET' = '#f07c19', 'LT' = '#faea7b')


###########################################################################
# Load and prepare data
###########################################################################

# Gene annotations
annotations <- read.delim("./data/annotations/protein_annotations_de_novo.tab")
annotations$ID <- make.unique(as.character(annotations$ID))

# add column for shortened gene name
annotations$gene.name <- make.unique(sub("_.*", "", annotations$ID))

# replace gene names according to replace_names dataframe
replace_names <- data.frame(from = 'CTNB1', to = 'CTNNB1')
annotations$gene.name <- sapply(annotations$gene.name, function(x) ifelse(x %in% replace_names[['from']],
                                                                          as.character(dplyr::filter(replace_names, from == x) %>% dplyr::pull(to)),
                                                                          x))

# readcounts
read.counts <- read.delim("./data/expression/de_novo.gene.counts.matrix")
colnames(read.counts) <- c("gene.name", "BHTB_1", "BHTB_2", "BHTB_3", "TTJ_1", "TTJ_2", "TTJ_3", "SL_1", "SL_2", "SL_3", "ET_1", "ET_2", "ET_3", "LT_1", "LT_2", "LT_3")

# Filter read counts by protein annotations
read.counts <- read.counts[read.counts[,1] %in% annotations$trinity_id,]

# reorder annotations to match the TRINITY ID order in read.counts, then swap the TRINITY ID for the gene name
annotations <- annotations[match(read.counts[,1], annotations$trinity_id),]
read.counts[,1] <- annotations$gene.name

# set rownames and delete column
rownames(read.counts) <- read.counts[,1]
read.counts[,1] <- NULL

#Add group and pool to metadata
colData <- data.frame(Group = factor(x = sub("_.*", "", colnames(read.counts)), levels = c('BHTB', 'TTJ', 'SL', 'ET', 'LT')),
                      row.names = colnames(read.counts),
                      Pool = sub(".*_", "", colnames(read.counts)))

# TRINITY outputs fractional counts for multi-mapped reads so I round the reads for downstream DE with DEseq.
read.counts <- round(read.counts)

# Make deseq object and filter low count genes
deseq <- DESeqDataSetFromMatrix(read.counts, design = ~ 0 + Group + Pool, colData = colData)

deseq$Group <- droplevels(deseq$Group)
deseq$Pool <- droplevels(deseq$Pool)

# normalise using median of ratios method
deseq <- estimateSizeFactors(deseq)

### Remove genes which do not have 5 readcounts in at least 2 samples
deseq <- deseq[rowSums(counts(deseq) > 5) > 2]

# rlog transform the data. The regularised log stops the most highly expressed
# genes from dominating downstream clustering; fitType 'local' estimates dispersion.
rld <- rlog(deseq, blind=FALSE, fitType='local')


###########################################################################
# QC plots
###########################################################################

# Check distribution of data after transformation
png(paste0(output_path, "rlog_Hist.png"), height = 15, width = 15, units = "cm", res = 200)
hist(x = rowMeans(assay(rld)),
     main = "Histogram of average rlog values for each gene",
     xlab = "rlog")
graphics.off()

##### Plot sample-sample distance and PCA

# sample-sample distance
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(rld))
colnames(sampleDistMatrix) <- paste(colnames(rld))
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(paste0(output_path, "SampleDist.png"), height = 12, width = 15, units = "cm", res = 400)
pheatmap(sampleDistMatrix, color = colours)
graphics.off()

# PCA
png(paste0(output_path, "SamplePCA.png"), height = 12, width = 12, units = "cm", res = 200)
plotPCA(rld, intgroup = "Group") +
  scale_colour_manual(values = tissue_colours) +
  theme(aspect.ratio=1,
        panel.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_rect(fill = NA, color = NA))
graphics.off()


###########################################################################
# LRT differential expression and SOM clustering
###########################################################################
# SOMs are calculated across all replicates and visualised as the per-tissue
# average; each cluster is assigned to its highest-expressing tissue.


# Now we carry out DEA using likelihood ratio test (LRT) in order to identify genes which vary significantly between the grouping variable (tissue)
# LRT compares two models, one full model and one reduced one where the grouping variable is removed
# This is to identify genes which vary across the tissues, which will then be clustered by SOM
deseq <- DESeq(deseq, test = 'LRT', reduced = ~ 0 + Pool, fitType = 'local')

res_LRT <- results(deseq)

# Although fold changes present they are not directly associated with the actual hypothesis test therefore genes are filtered by padj (FDR) < 0.0001 not FC
sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  dplyr::filter(padj < 0.0001)

# 3275 genes vary significantly as a function of tissue
length(sig_res_LRT$gene)

# get counts for DE genes
DE.counts <- assay(rld)[rownames(rld) %in% sig_res_LRT$gene,]

# Plot heatmap of differentially expressed genes
png(paste0(output_path, "DE_LRT_HM.png"), height = 30, width = 20, units = "cm", res = 200)
pheatmap(DE.counts, cluster_rows=T, show_rownames=FALSE, scale = "row",
         cluster_cols=T, annotation_col=as.data.frame(colData(deseq)["Group"]) %>% rename(Group = 'Tissue'),
         annotation_colors = list(Tissue = tissue_colours),
         treeheight_row = 20, treeheight_col = 40)
graphics.off()


# cluster using SOM
set.seed(2)

# calculate SOM based on tissue average for each gene
scaled.DE.counts <- t(scale(t(DE.counts)))

# calculate average SOM in order to get cluster identity
som.out <- som(scaled.DE.counts, grid = somgrid(5, 5, "hexagonal"))
som.out$grid$pts[,2] <- rev(som.out$grid$pts[,2])

# Get gene cluster IDs from SOM
gene.cluster <- som.out$unit.classif
names(gene.cluster) <- row.names(as.data.frame(som.out$data))

# get gene counts for som clusters for labelling violin plots
labels = data.frame(cluster = c(1:25), label = paste0('n = ', as.vector(table(gene.cluster))))
labels <- labels %>% mutate(label = paste0('C', cluster, '; ', label))

######## plot SOM for all replicates

# Add cluster IDs to expression data
group_melt <- reshape2::melt(as.matrix(scaled.DE.counts))
group_melt[['cluster']] <- gene.cluster[match(group_melt$Var1, names(gene.cluster))]

# label clusters by highest median expression for each cluster -> first average expression for each tissue
group_melt[["tissue_expressed"]] <- sub("_.*", "", group_melt$Var2)

max_med <- group_melt %>% dplyr::group_by(cluster, tissue_expressed) %>% dplyr::summarise(median = median(value, na.rm = TRUE)) %>% dplyr::group_by(cluster) %>% dplyr::filter(median == max(median)) %>% dplyr::ungroup()
max_med <- as.data.frame(max_med)
group_melt[["colour"]] <- unlist(lapply(group_melt$cluster, function(x) max_med$tissue_expressed[max_med$cluster %in% x]))
group_melt$colour <- factor(group_melt$colour, levels = c("BHTB", "TTJ", "SL", "ET", "LT"))

png(paste0(output_path, "SOM_DEgenes_violin.png"), height = 30, width = 35, units = "cm", res = 200)
ggplot(group_melt, aes(x = Var2, y = value)) +
  geom_violin(aes(group = Var2, fill = colour)) +
  scale_fill_manual(values = tissue_colours) +
  ylim(-2, 2.5) +
  facet_wrap(~cluster, ncol = 5, scales = "free_y", labeller = labeller(cluster = c(1:25))) +
  geom_text(data = labels, aes(x = 3.5, y = 2.3, label = label), size = 4) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1), axis.line = element_line(colour = "black"), strip.text.x = element_blank())
graphics.off()


# save list of genes in each cluster
cluster_gene_list <- group_melt %>% dplyr::select(c(Var1, cluster, colour)) %>% unique() %>% arrange(cluster) %>% dplyr::rename(Gene = Var1, Cluster = cluster, `Assigned Tissue` = colour)
write.csv(cluster_gene_list, file = paste0(output_path, 'cluster_gene_list.csv'), quote = FALSE, row.names = F)

######## Plot SOM for average tissue expression

# calculate average tissue expression for SOM visualisation
group.mean <- sapply(c("BHTB", "TTJ", "SL", "ET", "LT"), function(x) rowMeans(DE.counts[,grep(x, colnames(DE.counts))]))
group.mean <- t(scale(t(group.mean)))

# Add cluster IDs to expression data
group_melt <- reshape2::melt(as.matrix(group.mean))
group_melt$cluster <- gene.cluster[match(group_melt$Var1, names(gene.cluster))]

# label clusters by highest median expression for each cluster
max_med <- group_melt %>% group_by(cluster, Var2) %>% summarise(median = median(value, na.rm = TRUE)) %>% group_by(cluster) %>% filter(median == max(median)) %>% ungroup()
max_med <- as.data.frame(max_med)
group_melt[["Tissue"]] <- unlist(lapply(group_melt$cluster, function(x) max_med$Var2[max_med$cluster %in% x]))

png(paste0(output_path, "SOM_DEgenes_violin_mean.png"), height = 25, width = 25, units = "cm", res = 200)
ggplot(group_melt, aes(x = Var2, y = value)) +
  geom_violin(aes(group = Var2, fill = Tissue)) +
  ylim(-2, 2.5) +
  scale_fill_manual(values = tissue_colours) +
  facet_wrap(~cluster, ncol = 5, scales = "free_y", labeller = labeller(cluster = c(1:25))) +
  geom_text(data = labels, aes(x = 2, y = 2.3, label = label), size = 4) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1), axis.line = element_line(colour = "black"), strip.text.x = element_blank())
graphics.off()


###########################################################################
# GO enrichment (per SOM tissue cluster)
###########################################################################

# Extract gene lists from SOM clusters
som_clusters_genes <- group_melt %>%
  arrange(cluster) %>%
  dplyr::mutate(gene_id = Var1) %>%
  dplyr::select(c(gene_id, Tissue)) %>%
  split(.$Tissue)

som_clusters_genes <- lapply(som_clusters_genes, function(x) {
  x <- as.character(x$gene_id)
  return(x[!duplicated(x)])
})

# remove genes without GO blast hits
filt_annotations <- annotations %>%
  filter(grepl("GO", GO_blast))

# split GO annotations
go_terms <- apply(filt_annotations, 1, function(x){
  sub(unlist(strsplit(as.character(x['GO_blast'][1]), split = "`")), pattern = "\\^.*", replacement = "")
})

# name genes with gene name or trinity ID depending on preference
names(go_terms) <- filt_annotations$gene.name
go_terms <- reshape2::melt(go_terms)
goFrame = data.frame(frame.go_id = go_terms$value, frame.Evidence = 'TAS', frame.gene_id = go_terms$L1)
goFrame=GOFrame(goFrame, organism="s.canicula")
goAllFrame=GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

# Prepare GSEA data for hyperGTests - set background as all genes in the RNAseq - run one hyperGTest per som tissue group
go_enrichment <- lapply(names(som_clusters_genes), function(x){
  GSEAGOHyperGParams(name=paste0("GO_SOM_", x, "_clusters"),
                     geneSetCollection=gsc,
                     geneIds = som_clusters_genes[[x]],
                     universeGeneIds = unique(go_terms$L1),
                     ontology = "BP",
                     pvalueCutoff = 0.01,
                     conditional = FALSE,
                     testDirection = "over")
})

names(go_enrichment) <- names(som_clusters_genes)

# Run hyperGTest
go_enrichment <- lapply(go_enrichment, hyperGTest)

# Extract summary **and add FDR column**
go_enrichment_summary <- lapply(go_enrichment, function(hg) {
  out <- summary(hg)                      # original GOstats table
  out$FDR <- p.adjust(out$Pvalue,         # BH correction
                      method = "BH")
  out
})

# Write csv for GO enrichment for each tissue
curr_output_path <- paste0(output_path, 'functional_enrichment/')
dir.create(curr_output_path)

for(tissue in names(go_enrichment_summary)){
  write.csv(go_enrichment_summary[[tissue]], file = paste0(curr_output_path, 'top_go_', tissue, '.csv'), quote = TRUE, row.names = FALSE)
}


########### Manually select GO terms to plot

GO_oi <- list(
  BHTB = c(
    'GO:0050909' = 'sensory perception of taste',
    'GO:0008283' = 'cell population proliferation',
    'GO:0070887' = 'cellular response to chemical stimulus',
    'GO:0032103' = 'positive regulation of response to external stimulus',
    'GO:0014070' = 'response to organic cyclic compound'
  ),
  TTJ = c(
    'GO:0016477' = 'cell migration',
    'GO:1903034' = 'regulation of response to wounding',
    'GO:0043066' = 'negative regulation of apoptotic process',
    'GO:0010942' = 'positive regulation of cell death',
    'GO:0030855' = 'epithelial cell differentiation'
  ),
  SL = c(
    'GO:0042127' = 'regulation of cell population proliferation',
    'GO:0008283' = 'cell population proliferation',
    'GO:0001667' = 'ameboidal-type cell migration',
    'GO:0060070' = 'canonical Wnt signaling pathway',
    'GO:0061138' = 'morphogenesis of a branching epithelium',
    'GO:0050678' = 'regulation of epithelial cell proliferation',
    'GO:0007492' = 'endoderm development',
    'GO:0048870' = 'cell motility',
    'GO:0050680' = 'negative regulation of epithelial cell proliferation',
    'GO:0014033' = 'neural crest cell differentiation',
    'GO:0006935' = 'chemotaxis',
    'GO:0042060' =  'wound healing',
    'GO:0035107' =  'appendage morphogenesis',
    'GO:0031099' =  'regeneration',
    'GO:0048864' =  'stem cell development'
  ),  
  ET = c(
    'GO:0001503' = 'ossification',
    'GO:0030154' = 'cell differentiation',
    'GO:0016477' = 'cell migration',
    'GO:0060485' = 'mesenchyme development',
    'GO:0048870' = 'cell motility',
    'GO:0042476' = 'odontogenesis',
    'GO:0042127' = 'regulation of cell population proliferation',
    'GO:0045596' = 'negative regulation of cell differentiation',
    'GO:0030509' = 'BMP signaling pathway',
    'GO:0060070' = 'canonical Wnt signaling pathway',
    'GO:0030855' = 'epithelial cell differentiation',
    'GO:0048863' = 'stem cell differentiation'
  ),
  LT = c(
    'GO:0030198' = 'extracellular matrix organization',
    'GO:0048870' = 'cell motility',
    'GO:0048771' = 'tissue remodeling',
    'GO:0061448' = 'connective tissue development',
    'GO:0010941' = 'regulation of cell death'
  )
)


names(GO_oi) <- NULL
GO_oi <- unlist(GO_oi)
GO_oi <- GO_oi[!duplicated(GO_oi)]


go_enrichment_summary_subset <- lapply(go_enrichment_summary, function(x) x[,c('GOBPID', 'FDR', 'OddsRatio', 'Term')])


filtered_go <- lapply(go_enrichment_summary_subset, function(x){
  go_hits <- x[x$GOBPID %in% names(GO_oi),]
  if(nrow(go_hits) < length(GO_oi)){
    missing_go <- GO_oi[!names(GO_oi) %in% x$GOBPID]
    missing_go <- data.frame(GOBPID = names(missing_go), FDR = 1, OddsRatio = 0, Term = missing_go, row.names = NULL)
    return(rbind(go_hits, missing_go))
  } else {
    return(go_hits)
  }
})


# combine tissues
filtered_go <- bind_rows(filtered_go, .id = "Tissue")

# transform FDR → –log10
filtered_go$`-log10(FDR)` <- -log10(filtered_go$FDR)

# clamp values for plotting
filtered_go$`-log10(FDR)`[filtered_go$`-log10(FDR)` > 5] <- 5
filtered_go$OddsRatio[filtered_go$OddsRatio > 100]         <- 100

# factor ordering as before
filtered_go$Tissue <- factor(filtered_go$Tissue,
                             levels = c("BHTB","TTJ","SL","ET","LT"))
filtered_go$Term   <- factor(filtered_go$Term, levels = rev(GO_oi))


# Plot
curr_plot_path <- paste0(output_path, 'dotplots/')
dir.create(curr_plot_path)

png(paste0(curr_plot_path, "selected_go.png"),
    width = 25, height = 25, units = "cm", res = 400)
ggplot(filtered_go, aes(x = Term, y = Tissue)) +
  geom_point(aes(size = OddsRatio,
                 fill = `-log10(FDR)`),          # updated legend label
             colour = "black", shape = 21) +
  scale_size("Gene Count/Expected", range = c(0, 15)) +
  scale_fill_gradientn(
    colours = viridisLite::magma(100),
    guide   = guide_colorbar(ticks.colour = "black",
                             frame.colour = "black"),
    name    = "-log10(FDR)"                   # updated legend title
  ) +
  ylab("Tissue") + xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title  = element_text(size = 14)) +
  coord_flip()
graphics.off()


###########################################################################
# RNEA regulatory network
###########################################################################

curr_output_path = paste0(output_path, 'RNEA/')
dir.create(curr_output_path)


# ---------------------------------------------------------------------------
# RNEA regulatory-network inference is run OUTSIDE this repository.
#
# RNEA (Chouvardas, Kollias & Nikolaou, "Inferring active regulatory networks
# from gene expression data using a combination of prior knowledge and
# enrichment analysis", BMC Bioinformatics 2016, doi:10.1186/s12859-016-1040-7;
# https://sites.google.com/a/fleming.gr/rnea/) is third-party software. We do
# NOT redistribute it in this repository because (a) it is not released under a
# licence that permits redistribution, and (b) its reference files are derived
# from external databases with their own terms - notably KEGG, which restricts
# redistribution. The original hosted tarball is also no longer downloadable
# (the Google Sites URL now returns a login page rather than the archive).
#
# To reproduce this step: obtain RNEA.R and its ReferenceFiles/ from the authors,
# source('RNEA.R'), and run RNEA() on the RNEA_input.txt written below (the exact
# call is retained, commented out, further down). This repository instead ships
# RNEA's network output (suppl_files/RNEA_output/Network.csv) so that
# the downstream network analysis runs end-to-end from that checkpoint.
#
# Original (now non-functional) download + load, kept for reference only:
# if (!file.exists('./RNEA.R')) {
#   download.file(url = 'https://sites.google.com/a/fleming.gr/rnea/downloads/RNEA.tar.gz',
#                 destfile = 'RNEA.tar.gz')
#   untar('RNEA.tar.gz')
#   file.copy(c("./RNEA/ReferenceFiles", "./RNEA/RNEA.R"), getwd(), recursive = TRUE)
#   unlink(c('./RNEA.tar.gz', './RNEA'), recursive = TRUE)
# }
# source('RNEA.R')
# ---------------------------------------------------------------------------



# In order to establish a predictive network for dental initiation within the SL, we need to identify which genes we wish to input into RNEA

# First we run DEA on SL vs all other tissues in order to identify genes Up/Down regulated in the SL

# DEA for SL samples vs everything else
deseq_res <- results(deseq, contrast=list("GroupSL", c("GroupBHTB", "GroupTTJ", "GroupET", "GroupLT")), listValues=c(1, -1/4))

# Generate summary table of DEA
dea_summary <- as.data.frame(deseq_res[abs(deseq_res$log2FoldChange) > 0.75 & deseq_res$padj < 0.001,])
dea_summary <- arrange(dea_summary, desc(log2FoldChange))
dea_summary <- rownames_to_column(dea_summary, var = 'gene')
sum(dea_summary$log2FoldChange > 0.75) # 748 genes upregulated in SL relative to other tissues
sum(dea_summary$log2FoldChange < 0.75) # 642 genes downregulated in SL relative to other tissues

write.table(dea_summary, paste0(output_path, 'SLvsALL_dea_res.tsv'), quote = FALSE, sep = '\t', row.names = FALSE)

####

deseq_res <- as.data.frame(deseq_res) %>%
  rownames_to_column(var = "geneName") %>%
  dplyr::select(geneName, log2FoldChange, padj) %>%
  mutate(RNEA_geneName = sub("\\..*", "", geneName)) %>% # remove variant information from RNEA input for database searching
  mutate(sig = ifelse(abs(log2FoldChange) > 0.75 & padj < 0.001, TRUE, FALSE))

## Run RNEA
# prepare data frame for RNEA
RNEA_input <- deseq_res[,c("RNEA_geneName", "log2FoldChange", "padj")]
colnames(RNEA_input) <- c("geneName", "LogFC", "pvalue")

# We filter only DE genes (abs(log2FC) > 0.75 & FDR < 0.001) - however we also include SOX2 and BMP4 even though they are not significantly differentially expressed as we know they are important in the regulation of stem cell fate in the SL.
additional_genes <- filter(RNEA_input, geneName %in% c('SOX2', 'BMP4'))
RNEA_input <- filter(RNEA_input, abs(LogFC) > 0.75 & pvalue < 0.001) %>% bind_rows(additional_genes)

write_delim(RNEA_input, paste0(curr_output_path, "/RNEA_input.txt"), col_names = T, delim = "\t")

# RNEA was run externally on the RNEA_input.txt written above (see note at the
# top of this section). The exact call used was:
# RNEA(filename = paste0(curr_output_path, "RNEA_input.txt"), identifier = "GeneName", species = "Human", FC_threshold = 0, PV_threshold = 1,
#      output = curr_output_path, network = "regulatory", type_of_output = "csv")

# Read RNEA's network output shipped with this repository (reproducibility checkpoint).
network <- read.csv("./suppl_files/RNEA_output/Network.csv", stringsAsFactors = F)

# remove rows with empty strings (there is an empty row added by RNEA)
network <- network[!apply(network, 1, function(x) any(x=="")),] 

# # First: recover corresponding variants when reading back in RNEA network data for DE genes
de_deseq_res <- deseq_res %>% filter(sig == TRUE)
de_network_genes <- de_deseq_res$geneName[de_deseq_res$RNEA_geneName %in% unique(unlist(network))]

# Next: recover corresponding variants when reading back in RNEA network data for RNEA additions
rnea_additions <- unique(unlist(network))[unlist(lapply(unique(unlist(network)), function(x) {
  !any(grepl(x, de_network_genes, ignore.case = TRUE))}
))]

# Rbind the two together to get all network isoforms
full_network_genes <- c(de_network_genes, rnea_additions)

# remove genes from network which have been added by RNEA and are not in our dataset
missing_genes <- unique(unlist(network))[unlist(lapply(unique(unlist(network)), function(x) {
  !any(grepl(x, rownames(read.counts), ignore.case = TRUE))}
))]


if(length(missing_genes) > 0){
  network <- network[!apply(network, 1, function(x) any(x %in% missing_genes)),] 
}

# We subsequently filtered our putative network by removing genes which are not within two degrees of separation of our central nodes.

# We also iteratively remove any markers which interact with only a single other gene, as these genes exhibit low network connectivity and are therefore less likely to play a key role in our TF regulatory network.

# subset no degrees of separation
central_nodes <-  c("CTNNB1", "PITX1", "PITX2", "SOX2", "LEF1", "BMP4", "SHH")
no_sep <- filter(network, Source %in% central_nodes | Target %in% central_nodes)

write_csv(no_sep, paste0(curr_output_path, "no_sep.csv"), col_names = T)

# subset one degree of separation
secondary_gene_list <- unique(unlist(no_sep))
one_sep <- filter(network, Source %in% secondary_gene_list | Target %in% secondary_gene_list)
write_csv(one_sep, paste0(curr_output_path, "one_sep.csv"), col_names = T)

# function to identify genes with a single interacting partner
single_genes <- function(network_table) {
  # split interactions into list for each gene
  filtered_network <- lapply(unique(unlist(network_table)), function(x) network_table %>% filter_all(any_vars(. %in% x)))
  names(filtered_network) <- unique(unlist(network_table))
  return(names(filtered_network[sapply(filtered_network, function(x) length(unique(unlist(x))) < 3)]))
}

# iteratively remove genes which only interact with one gene until all genes have more than one partner
filtered_network <- one_sep

while(length(single_genes(filtered_network)) > 0){
  filt <- single_genes(filtered_network)
  filtered_network <- filtered_network[!apply(filtered_network, 1, function(x) any(x %in% filt)),]
}

# Next we identify transcription factors based on gene ontology terms GO:0003700, GO:0043565, and GO:0000981
tf_annotations <- filt_annotations %>%
  dplyr::filter(grepl('GO:0003700|GO:0043565|GO:0000981', GO_pfam)) %>%
  dplyr::filter(gene.name %in% full_network_genes)


# Finally we also generate a metadata table which defines the differential expression status of a given gene for visualisation in cytoscape
grn_metadata <- deseq_res %>%
  filter(geneName %in% full_network_genes) %>%
  mutate(centralNode = ifelse(geneName %in% central_nodes, TRUE, FALSE)) %>%
  mutate(
    DE = case_when(
      sig & log2FoldChange > 0.75 ~ 'upregulated',
      sig & log2FoldChange <= 0.75 ~ 'downregulated',
      TRUE ~ 'RNEA_addition'
    )) %>%
  mutate(TF = ifelse(geneName %in% tf_annotations$gene.name, TRUE, FALSE)) %>%
  mutate(filteredNetwork = ifelse(RNEA_geneName %in% unique(unlist(filtered_network)), TRUE, FALSE))


# Collapse the metadata by isoform and keep the highest abs(log2FC)
grn_metadata_iso_collapsed <- grn_metadata %>%
  dplyr::select(-geneName) %>%
  group_by(RNEA_geneName) %>%
  # Keep the row with the highest absolute value of 'log2FoldChange'
  filter(abs(log2FoldChange) == max(abs(log2FoldChange))) %>%
  ungroup()



# Gene names in grn_metadata don't match the network directly, because each RNEA
# gene can correspond to multiple transcript variants. Match RNEA gene names back
# to all potential variants in our DEA to build the final network.

# Remove var information from gene names
grn_metadata$geneNameShort <- sapply(strsplit(grn_metadata$geneName, "\\."), function(x) x[1])

# Join the dataframes on the keys
network_join1 <- merge(network, grn_metadata, by.x = 'Source', by.y = 'geneNameShort') %>% dplyr::select(c(Source,Target,geneName)) %>% dplyr::rename(Source_long = geneName)
network_join2 <- merge(network, grn_metadata, by.x = 'Target', by.y = 'geneNameShort') %>% dplyr::select(c(Source,Target,geneName)) %>% dplyr::rename(Target_long = geneName)


# Joining the dataframes based on Source and Target columns
network_long_isoform <- left_join(network_join1, network_join2, by = c("Source", "Target"), relationship = "many-to-many") %>%
  dplyr::select(c(Source_long, Target_long)) %>%
  dplyr::rename(Source = Source_long, Target = Target_long) %>%
  na.omit() %>%
  distinct()

network_long_collapsed_iso <- left_join(network_join1, network_join2, by = c("Source", "Target"), relationship = "many-to-many") %>%
  dplyr::select(c(Source, Target)) %>%
  na.omit() %>%
  distinct()

write_csv(network_long_isoform, paste0(curr_output_path, "network_renamed_iso.csv"), col_names = T)
write_csv(network_long_collapsed_iso, paste0(curr_output_path, "network_renamed_iso_collapsed.csv"), col_names = T)

# Subset filtered genes and save
filtered_network_genes_iso <- grn_metadata$geneName[grn_metadata$filteredNetwork]
filtered_network_long_iso <- network_long_isoform[apply(network_long_isoform, 1, function(x) all(x %in% filtered_network_genes_iso)),]
rownames(filtered_network_long_iso) <- NULL
write_csv(filtered_network_long_iso, paste0(curr_output_path, "filtered_network_iso.csv"), col_names = T)

# Filter the final network for just direct interactions with central nodes
filtered_network_long_iso_nosep <- filter(filtered_network_long_iso, Source %in% central_nodes | Target %in% central_nodes)
write_csv(filtered_network_long_iso_nosep, paste0(curr_output_path, "filtered_network_iso_nosep.csv"), col_names = T)


# Add no sep filtered network in grn metadata
grn_metadata <- grn_metadata %>%
  mutate(filteredNetwork_nosep = ifelse(geneName %in% unique(unlist(filtered_network_long_iso_nosep)), TRUE, FALSE))


# Subset filtered genes and save after collapsing isoforms
filtered_network_genes_collapsed_iso <- grn_metadata$geneNameShort[grn_metadata$filteredNetwork]
filtered_network_long_collapsed_iso <- network_long_collapsed_iso[apply(network_long_collapsed_iso, 1, function(x) all(x %in% filtered_network_genes_collapsed_iso)),]
rownames(filtered_network_long_collapsed_iso) <- NULL
write_csv(filtered_network_long_collapsed_iso, paste0(curr_output_path, "filtered_network_iso_collapsed.csv"), col_names = T)

# Filter the final network for just direct interactions with central nodes
filtered_network_long_collapsed_iso_nosep <- filter(filtered_network_long_collapsed_iso, Source %in% central_nodes | Target %in% central_nodes)
write_csv(filtered_network_long_collapsed_iso_nosep, paste0(curr_output_path, "filtered_network_iso_collapsed_nosep.csv"), col_names = T)


write_csv(arrange(grn_metadata, geneName), paste0(curr_output_path, "grn_metadata.csv"), col_names = T)
write_csv(arrange(grn_metadata_iso_collapsed, RNEA_geneName), paste0(curr_output_path, "grn_metadata_iso_collapsed.csv"), col_names = T)



###########################################################################
# Network heatmaps
###########################################################################

# Heatmap for network genes
plot_metadata <- grn_metadata %>%
  filter(DE != 'RNEA_addition') %>%
  rename(DE = 'Differentially expressed')

# get counts for DEA genes
DE.counts <- assay(rld)
HM.counts <- DE.counts[rownames(DE.counts) %in% plot_metadata$geneName,]

png(paste0(output_path, "network_hm.png"), height = 60, width = 30, units = "cm", res = 200)
pheatmap(HM.counts, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, scale = "row",
         cluster_cols=FALSE, annotation_col=as.data.frame(colData(deseq)["Group"]) %>% rename(Group = 'Tissue'),
         annotation_row = plot_metadata %>% dplyr::select(geneName, `Differentially expressed`) %>% column_to_rownames('geneName'),
         annotation_names_row = F, annotation_colors = list(Tissue = tissue_colours, `Differentially expressed` = c('upregulated' = '#66BD63', 'downregulated' = '#9970AB', 'non-significant' = '#BABABA')),
         treeheight_row = 20, treeheight_col = 40)
graphics.off()

# Heatmap for filtered network genes
plot_metadata <- plot_metadata %>%
  filter(filteredNetwork == TRUE)


HM.counts <- DE.counts[rownames(DE.counts) %in% plot_metadata$geneName,]

png(paste0(output_path, "filtered_network_hm.png"), height = 25, width = 28, units = "cm", res = 200)
pheatmap(HM.counts, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, scale = "row",
         cluster_cols=FALSE, annotation_col=as.data.frame(colData(deseq)["Group"]) %>% rename(Group = 'Tissue'),
         annotation_row = plot_metadata %>% dplyr::select(geneName, `Differentially expressed`) %>% column_to_rownames('geneName'),
         annotation_names_row = F, annotation_colors = list(Tissue = tissue_colours, `Differentially expressed` = c('upregulated' = '#66BD63', 'downregulated' = '#9970AB', 'non-significant' = '#BABABA')),
         treeheight_row = 20, treeheight_col = 40)
graphics.off()



# Heatmap for filtered network no sep (same as mini GRN)
plot_metadata <- grn_metadata %>%
  mutate(DE = ifelse(DE == "RNEA_addition", "non-significant", DE)) %>% # Modify DE
  rename(DE = 'Differentially expressed') %>% # Rename DE
  filter(filteredNetwork_nosep == TRUE) # Filter rows

HM.counts <- DE.counts[rownames(DE.counts) %in% plot_metadata$geneName,]

png(paste0(output_path, "filtered_nosep_network_hm.png"), height = 25, width = 25, units = "cm", res = 200)
pheatmap(HM.counts, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, scale = "row",
         cluster_cols=FALSE, annotation_col=as.data.frame(colData(deseq)["Group"]) %>% rename(Group = 'Tissue'),
         annotation_row = plot_metadata %>% dplyr::select(geneName, `Differentially expressed`) %>% column_to_rownames('geneName'),
         annotation_names_row = F, annotation_colors = list(Tissue = tissue_colours, `Differentially expressed` = c('upregulated' = '#66BD63', 'downregulated' = '#9970AB', 'non-significant' = '#BABABA')),
         treeheight_row = 20, treeheight_col = 40)
graphics.off()

# Heatmap for filtered network no sep TFs
plot_metadata <- grn_metadata %>%
  mutate(DE = ifelse(DE == "RNEA_addition", "non-significant", DE)) %>% # Modify DE
  rename(DE = 'Differentially expressed') %>% # Rename DE
  filter(filteredNetwork_nosep == TRUE) %>% # Filter rows
  filter(TF == TRUE)

HM.counts <- DE.counts[rownames(DE.counts) %in% plot_metadata$geneName,]

png(paste0(output_path, "filtered_nosep_network_hm_TF.png"), height = 10, width = 25, units = "cm", res = 200)
pheatmap(HM.counts, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, scale = "row",
         cluster_cols=FALSE, annotation_col=as.data.frame(colData(deseq)["Group"]) %>% rename(Group = 'Tissue'),
         annotation_row = plot_metadata %>% dplyr::select(geneName, `Differentially expressed`) %>% column_to_rownames('geneName'),
         annotation_names_row = F, annotation_colors = list(Tissue = tissue_colours, `Differentially expressed` = c('upregulated' = '#66BD63', 'downregulated' = '#9970AB', 'non-significant' = '#BABABA')),
         treeheight_row = 20, treeheight_col = 40)
graphics.off()

#### Plot SOM with labels for genes from filtered network

# extract key genes in filtered GRN for plotting on SOM
selected_network_genes <- c("BMP4", "PITX3", "PITX2", "ISL1", "SOX2", "CTNNB1", "KLF4.1", "GLI2", "BMP2.1", "MYCN.2", "RUNX2", "DKK1")

selected_network_data <- filter(group_melt, Var1 %in% selected_network_genes) %>%
  group_by(cluster, Var2) %>%
  mutate(position = 1:n())

# calculate mean profile for each cluster
tissue.mean <- group_melt %>% 
  group_by(group_melt$cluster,Var2) %>% 
  summarise(value = mean(value))
colnames(tissue.mean)[1] <- "cluster"

######## plot SOMs as violins with GOI ##########
colours <- c("#003f5c", "#bc5090", "#ffa600")
selected_network_data$plot_colour <-  colours[selected_network_data$position]
colours <- selected_network_data$plot_colour
names(colours) <- selected_network_data$Var1

plots <- list()
plots[['violins']] <- ggplot(group_melt, aes(x = Var2, y = value)) +
  geom_violin(aes(group = Var2, fill = Tissue)) +
  scale_fill_manual(values = tissue_colours) +
  geom_line(data = selected_network_data, aes(group=Var1, color = Var1), size=0.7) +
  geom_text(data = filter(selected_network_data, Var2 == 'BHTB') , aes(x = 3, y = -1.8, label = Var1, color = Var1), position = position_dodge(width = 5), size = 3.3, check_overlap = FALSE) +
  scale_colour_manual(values = colours) +
  ylim(-2, 2.5) +
  facet_wrap(~cluster, ncol = 5, scales = "free_y", labeller = labeller(cluster = c(1:25))) +
  geom_text(data = labels, aes(x = 2.5, y = 2.3, label = label), size = 4) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  guides(colour = FALSE) +
  theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1), axis.line = element_line(colour = "black"), strip.text.x = element_blank()) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))

legend_colours <- c('BHTB' = 'gray70', 'TTJ' = 'gray70', 'SL' = 'gray70', 'ET' = '#f07c19', 'LT' = 'gray70')
plots[['legend_1']] <- ggplot(group_melt %>% filter(Tissue == 'ET'), aes(x = Var2, y = value)) +
  geom_violin(aes(group = Var2, fill = Var2)) +
  scale_fill_manual(values = tissue_colours) +
  geom_text(aes(x = Var2, y = 3, label=Var2, color=Var2), size = 5) +
  scale_colour_manual(values = tissue_colours) +
  ylim(-2, 7) +
  theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(0,1,0,0, "cm"),
        strip.text.x = element_blank(),
        legend.position = "none")

plots[['legend_2']] <- ggplot(group_melt %>% filter(Tissue == 'ET'), aes(x = Var2, y = value)) +
  geom_violin(aes(group = Var2, fill = Var2)) +
  scale_fill_manual(values = legend_colours) +
  ylim(-4, 5) +
  theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(0,1,0,0, "cm"),
        strip.text.x = element_blank(),
        legend.position = "none")

plots[['legend_3']] <- ggplot(group_melt %>% filter(Tissue == 'ET'), aes(x = Var2, y = value)) +
  geom_violin(aes(group = Var2, fill = Tissue)) +
  scale_fill_manual(values = tissue_colours) +
  ylim(-6, 3) +
  theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(0,1,0,0, "cm"),
        strip.text.x = element_blank(),
        legend.position = "none")


lay <- rbind(c(2,1,1,1,1),
             c(3,1,1,1,1),
             c(4,1,1,1,1))

png(paste0(output_path, "SOM_candidate_overlay.png"), height = 23, width = 35, units = "cm", res = 200)
grid.arrange(grobs = plots, layout_matrix = lay)
graphics.off()


###########################################################################
# Candidate gene heatmap (MYC family)
###########################################################################
# MYCN is identified as the novel regulator from the GRN analysis above; here we
# plot rlog expression of the MYC-family transcripts and SHH across tissues.
candidate_hm_genes <- c('MYC1', 'SHH', 'MYCN.1', 'MYCN.2')

# per-gene differential expression status (SL vs all), as a row annotation
candidate_de <- deseq_res %>%
  dplyr::filter(geneName %in% candidate_hm_genes) %>%
  dplyr::mutate(`Differentially expressed` = dplyr::case_when(
    sig & log2FoldChange > 0.75  ~ 'upregulated',
    sig & log2FoldChange <= 0.75 ~ 'downregulated',
    TRUE                         ~ 'non-significant')) %>%
  dplyr::select(geneName, `Differentially expressed`) %>%
  column_to_rownames('geneName')

png(paste0(output_path, "goi_heatmap.png"), height = 10, width = 25, units = "cm", res = 200)
pheatmap(assay(rld)[candidate_hm_genes,], cluster_rows = TRUE, show_rownames = TRUE, show_colnames = FALSE, scale = "row",
         cluster_cols = FALSE, annotation_col = as.data.frame(colData(deseq)["Group"]) %>% rename(Group = 'Tissue'),
         annotation_row = candidate_de, annotation_names_row = FALSE,
         annotation_colors = list(Tissue = tissue_colours,
                                  `Differentially expressed` = c('upregulated' = '#66BD63', 'downregulated' = '#9970AB', 'non-significant' = '#BABABA')),
         treeheight_row = 20, treeheight_col = 40)
graphics.off()


###########################################################################
# Candidate regulator bar plots (Figure 7)
###########################################################################

# Bar plots of DESeq2-normalised counts for candidate regulators (Figure 7)
# Bars = mean normalised count per tissue, error bars = +/- SEM, points = replicates.
barplot_genes <- c('LEF1', 'BMP4', 'DKK1', 'PITX1', 'PITX3', 'MYCN.1')

# DESeq2 median-of-ratios normalised counts ("Normalised counts (DESeq2)")
norm_counts <- counts(deseq, normalized = TRUE)[barplot_genes,]

barplot_data <- reshape2::melt(as.matrix(norm_counts), varnames = c('gene', 'sample'), value.name = 'norm_count')
barplot_data$Tissue <- factor(sub("_.*", "", barplot_data$sample), levels = c('BHTB', 'TTJ', 'SL', 'ET', 'LT'))
barplot_data$gene   <- factor(barplot_data$gene, levels = barplot_genes)

# per-gene, per-tissue mean and standard error of the mean
barplot_summary <- barplot_data %>%
  dplyr::group_by(gene, Tissue) %>%
  dplyr::summarise(mean = mean(norm_count),
                   se   = sd(norm_count) / sqrt(dplyr::n()),
                   .groups = 'drop')

# one separate bar plot per gene
dir.create(paste0(output_path, "candidate_gene_barplot/"))
for (g in barplot_genes) {
  gene_summary <- dplyr::filter(barplot_summary, gene == g)
  gene_points  <- dplyr::filter(barplot_data, gene == g)
  
  set.seed(42) # deterministic point jitter
  png(paste0(output_path, "candidate_gene_barplot/", g, ".png"), height = 10, width = 10, units = "cm", res = 200)
  print(
    ggplot(gene_summary, aes(x = Tissue, y = mean, fill = Tissue)) +
      geom_col(width = 0.9) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3) +
      geom_point(data = gene_points, aes(x = Tissue, y = norm_count),
                 inherit.aes = FALSE, size = 1,
                 position = position_jitter(width = 0.12, height = 0)) +
      scale_fill_manual(values = tissue_colours) +
      ggtitle(g) +
      ylab("Normalised counts (DESeq2)") + xlab("") +
      theme_classic() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            axis.text  = element_text(colour = "black"))
  )
  graphics.off()
}
