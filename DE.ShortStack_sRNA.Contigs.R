# Author: This script was written by Jose Roberto Bermudez Barrientos with suggestions by Cei Abreu
# Journal entry: Week 9/24, Semester 6/8, Oct 29 - Nov 2, 2018
# Description: this script performs a differential expression analysis of ShortStack clusters between parariste-host interacting libraries vs host libraries
# Running instructions: hit the Source button in Rstudio
rm(list = ls())
library(edgeR)
library(RColorBrewer) # to use Dark2 color palette
library(Rsamtools) # to readFasta, also to scanBam
library(VennDiagram) # to use venn.diagram() which takes a list directly and produces a venn diagram
library(Sushi) # to use opaque(colors, transparency = 0.3) to add transparency to our set of choice of colors

# Parameters {BEGIN}
DEA_units <- "contigs" # "indiv_seqs", "contigs", "clusters"
read_assignment_method_label <- "informed0.2" # 6 Feb 2019 Obed's method to distribute reads, "SppC" stands for supporting counts
data_size <- "real" # either "real" or "trial"
first_run <- FALSE # the first run is slower, but it is necessary, if you are only adjusting images use first_run <- FALSE if you already have run it once
home_dir <- "/Users/beto"
md <- paste(home_dir, "Hp/II.Hp-sRNAs_in_mouse_cells.contigs_extra.informed0.2", sep = "/")
org_map_dir <- paste(md, "01b.org_map_table", sep = "/") # mapping information
counts_dir <- paste(md, "count_table_Obed", sep = "/")
#SS_dir <- paste(home_dir, "Hp/II.Hp-sRNAs_in_mouse_cells.clusters/02.ShortStack_sRNA_clusters.real", sep = "/")
genome_dir <- paste(home_dir, "Hp/II.Hp-sRNAs_in_mouse_cells.clusters/00.merged_index", sep = "/")
map_dir <- paste(md, "01d.contigs_to_merged_reference_tables", sep = "/")
#bam_file <- paste(SS_dir, "merged_alignments.bam", sep = "/")
de_dir <- paste(md, paste("03.DE.ShortStack_sRNA.Contigs", data_size, sep = "."), sep = "/")
# DEA key parameters {BEGIN}
logFC_thres <- 0
FDR_thres <- 0.1
min_cpm <- 1
min_libs <- 2
# DEA key parameters {END}
parasite_label <- "p" # 29 Oct 2018
host_label <- "h" # 29 Oct 2018
ambiguous_label <- "a" # 21 Nov 2018
unmapped_label <- "n" # 21 Nov 2018
parasite_color <-  "#27AE60" # Nephritis, green color will contrast with other darkest colors
host_color <- "#A9A9A9" # darker color for a high number of sequences
ambiguous_color <- "#8E44AD"
unmapped_color <- "#2980B9"
bin_number <- 30 # stacked barplots parameter # 29 Oct 2018
non_counts_columns <- c("origin", "mmhp", "len") # 5 Feb 2019 I use this variable to be able to subset only those columns that have libray read counts used for DEA
# Parameters {END}

# functions setting zone {BEGIN}
# this function is used to parse large numbers in easier to read letter representations
f2si2<-function (number,rounding=F) 
{
  lut <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06, 
           0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21, 
           1e+24)
  pre <- c("y", "z", "a", "f", "p", "n", "u", "m", "", "k", 
           "M", "G", "T", "P", "E", "Z", "Y")
  ix <- findInterval(number, lut)
  if (lut[ix]!=1) {
    if (rounding==T) {
      sistring <- paste(signif(number/lut[ix], digits = 1), pre[ix])
    }
    else {
      sistring <- paste(number/lut[ix], pre[ix])
    } 
  }
  else {
    sistring <- as.character(number)
  }
  return(sistring)
}
# functions setting zone {END}

dir.create(de_dir, recursive = T, showWarnings = F)

## 15 Nov 2018 {BEGIN}
#chr <- paste(sapply(strsplit(df_category$Locus, ":"), "[[", 1), sapply(strsplit(df_category$Locus, ":"), "[[", 2), sep = ":")
#coords_from_to <- sapply(strsplit(df_category$Locus, ":"), "[[", 3)
#coords_from <- sapply(strsplit(coords_from_to, "-"), "[[", 1)
#coords_to <- sapply(strsplit(coords_from_to, "-"), "[[", 2)
#df_bed <- data.frame(chr = chr, coords_from = coords_from, coords_to = coords_to, id = rownames(df_category), strand = df_category$Strand, org = df_category$org, stringsAsFactors = F)
## 15 Nov 2018 {BEGIN}
  
## 29 Oct 2018 {BEGIN}
# set the colors for each org
org_categories <- c(parasite_label, host_label, ambiguous_label, unmapped_label)
df_color <- data.frame(row.names = org_categories, color = c(parasite_color, host_color, ambiguous_color, unmapped_color), stringsAsFactors = FALSE)
## 29 Oct 2018 {END}

# loading annotations 22 Nov 2018
annot_file <- paste(home_dir, "Hp/Annotation/cei", "nHp.2.0.segmentation_annotation_GR.Rds", sep = "/")
GR_annot <- readRDS(annot_file)
seqlevels(GR_annot) <- paste("Hp", seqlevels(GR_annot), sep = ":")

# mapping files 22 Nov 2018
#GR_unique_file <- paste(map_dir, "GR_parasite_unique.RDS", sep = "/") # here unique refers contigs that map to a unique site in the parasite genome
#GR_parasite_unique <- readRDS(file = GR_unique_file)
#GR_multi_file <- paste(map_dir, "GR_parasite_multi.RDS", sep = "/") # here multi refers contigs that map to multiple sites in the parasite genome
#GR_parasite_multi <- readRDS(file = GR_multi_file)

# SECTION I: reading count table and performing filterings {BEGIN}
# IMPORTANT: this table has additional columns other than counts, I shall call it table_full and use the name m_count for only those columns involving read counts
counts_file <- paste(counts_dir, "c_alg_unitigs_extra.tbl", sep = "/")

RDS_counts_file <- paste(de_dir, "m_count.RDS", sep = "/") # a fast loading RDS file to load the filtered matrix faster for Rmarkdown
#RDS_bam_file <- paste(de_dir, "GR_bam.RDS", sep = "/") 
print_first_run <- paste("first_run flag:", first_run, sep = " ")
print(print_first_run)

if (first_run == TRUE){
  print(system.time({
  
    # declare an empty matrix, we will fill it later
    table_full <- read.table(counts_file, header = TRUE)
  }))
  print("time taken to read counts file file")
  
  # exclude libraries
  libraries_to_exclude <- c("M_neg_24_1", "vesicle_poly_1", "vesicle_poly_2")
  table_full <- table_full[,!colnames(table_full) %in% libraries_to_exclude]
  
  # reading and parsing ShortStack counts file {END}
  
  saveRDS(table_full, file = RDS_counts_file)
  
  # reading BAM file, contigs section {BEGIN} 22 Nov 2018
  print("time taken to read BAM file and convert it to GRanges object")
  print(system.time({
    what_to_extract <- c("rname", "qwidth", "pos", "strand")
    #l_bam <- scanBam(bam_file, param = ScanBamParam(what = what_to_extract, flag = scanBamFlag(isUnmappedQuery = FALSE) ))
    #df_bam <- do.call("data.frame", l_bam)
    #rm(l_bam)
    
    #df_bam$org <- sapply(strsplit(as.character(df_bam$rname), ":"), "[[", 1)
    
    #df_bam <- df_bam[df_bam$org == parasite_label,]
    
    #GR_bam <- GRanges(seqnames = as.character(df_bam$rname), ranges = IRanges(start = df_bam$pos, width = df_bam$qwidth), strand = df_bam$strand)
    #saveRDS(object = GR_bam, file =  RDS_bam_file) 
  }))
  # reading BAM file, contigs section {END} 22 Nov 2018
  
}else{
  
  print(system.time({
  table_full <- readRDS(file = RDS_counts_file)
  #GR_bam <- readRDS(file = RDS_bam_file)
  }))
  print("time taken to read RDS table_full and bam objects")
}

m_count_full <- table_full[,colnames(table_full)[!colnames(table_full) %in% non_counts_columns]]

## WARNING: contigs exclusive? {BEGIN}
# are ambiguous sequences shorter?
l_split <- split(table_full$len, table_full$origin)
c <- 1
ecdf_plot_file <- paste(de_dir, paste(DEA_units, "cumulative_length_distrib", "pdf", sep = "."), sep = "/")
pdf(file = ecdf_plot_file)
for(origin in names(l_split)){
  if(c == 1 ){
    plot(ecdf(l_split[[origin]]), xlim=c(18, 30), col=c, main="cumulative length distrib de novo assembly", ylab="Fraction", xlab="length")
  }else{
    plot(ecdf(l_split[[origin]]), add=T, col=c)
  }
  c <- c +1
}
legend("bottomright", legend = names(l_split), col = 1:c, pch = 19)
dev.off()
## WARNING: contigs exclusive? {END}

# then set the color for each clusters by appending an extra column to df_category
table_full$color <- df_color[table_full$origin,]

print(paste("number of reads in full matrix for", DEA_units, sep = " ")) # 30 Dic 2018
print(sum(m_count_full))

# reading genome to get all levels 22 Nov 2018
genome_file <- paste(genome_dir, "Hp_Mm_index.fa", sep = "/")
fa_info <- scanFaIndex(genome_file)
fa_info <- fa_info[grepl(pattern = parasite_label, seqnames(fa_info))]
seq_info <- Seqinfo(seqnames = as.character(seqnames(fa_info)), seqlengths = width(fa_info))
#seqlevels(GR_parasite_unique) <- seqlevels(fa_info)
#seqlevels(GR_parasite_multi) <- seqlevels(fa_info)

# reading contigs fasta file to get all contigs lengths
contigs_fa_file <- paste(counts_dir, "contigs_unitigs.fa", sep = "/")
#contigs_fa_index <- scanFaIndex(contigs_fa_file)

print(head(cpm(m_count_full)))
print(head(m_count_full))


### low expression filtering
print("Filtering low expresseds elements")
m_count <- m_count_full[rowSums(cpm(m_count_full) >= min_cpm) >= min_libs,]

# which libraries where more affected by the filtering?
print("before filtering")
print(round(colSums((m_count_full))/1e06,1))
print("after filtering")
print(round(colSums((m_count))/1e06,1))

print("m_count rows before filtering")
print(nrow(m_count_full))
print("m_count rows after filtering")
print(nrow(m_count))

print(paste("number of reads in matrix filtered for", DEA_units, sep = " ")) # 30 Dic 2018
print(sum(m_count))

# SECTION I: reading count table and performing filterings {BEGIN}

#### SECTION II: defining relevant comparisons
# preparing to to edgeR, setting experimental design and grouping
factors_setting <- sub("_.$","",colnames(m_count)) ## This line may vary according to the naming convention of samples
grp <- factor(factors_setting)
design <- model.matrix(~0+grp)
colnames(design) <- levels(grp)

# define contrasts of interest # IMPORTANT: contrast will vary between different symbiont systems and data sets 
contrasts <- makeContrasts(
  "MODE-K_vesicle_vs_no_treatment" = "(M_ves_24 + M_ves_4)/2 - (M_neg_24 + M_neg_4)/2",
  "MODE-K_vesicle_vs_no_treatment_4hrs" = "M_ves_4 - M_neg_4",
  levels=design)

# now do edgeR
print("building DGEList object")
de <- DGEList(counts=m_count, group=grp)

# normalization
par(mfrow = c(1,2))
plotMDS(de, main="prior to normalization")
de = calcNormFactors(de)
plotMDS(de, main="after normalization")

## Principal Component Analysis {BEGIN}
# We define colors according to factors
pca_colors <- brewer.pal(n = length(levels(grp)), name = "Paired")
names(pca_colors) <- levels(grp)
interquantrs <- apply(de$counts, 1, IQR)
OrderedExp = de$counts[order(interquantrs, decreasing=TRUE), ]
exprsPCA = log1p(OrderedExp[1:100,])
pca = prcomp(t(exprsPCA))#, retx = T, center = T)
percVar = round(100*(pca$sdev ** 2)/sum(pca$sdev ** 2))

pdf(file = paste(de_dir, "PCA_100_highest_var_C1_C2.pdf", sep = "/"), width = 5, height = 5)
par(mar=c(5,5,2,1))
plot(pca$x[,c(1,2)], main="PCA of mouse cells small RNA-Seq data", col=pca_colors[factors_setting], pch= 18, cex=2, xlab = paste("PC1 (", percVar[1],"%)", sep = ""), ylab = paste("PC2 (", percVar[2],"%)", sep = "")) #, xlim=c(-30,30), ylim=c(-6,6))
text(pca$x[,c(1,2)], rownames(pca$x), pos=3, cex=0.8)
dev.off()
## Principal Component Analysis {END}

# estimating data dispersion
de <- estimateGLMCommonDisp(de, design, verbose=T)
de <- estimateGLMTrendedDisp(de, design) # as we are comparing very different conditions it is better to use the common dispersion only
de <- estimateGLMTagwiseDisp(de, design) # it is interesting to calculate both trended and tagwise dispersion first just to visualize it, eventhough you may not use it in the end

plotBCV_file <- paste(de_dir, "Biological_Coefficient_of_Variation.png", sep = "/")
png(filename = plotBCV_file, width = 960, height = 720)
plotBCV(de)
dev.off()

venn_pos_l <- list()
venn_neg_l <- list()
# estimating the dispersion in all contrasts
print(paste("common dispersion", de$common.dispersion))
de$common.dispersion <- 1.62627  # basal common dispersion (reads) for reads, contigs, clusters # 26 March 2019
print(paste("common dispersion", de$common.dispersion, "after using basal common dispersion")) # 26 March 2019
fit <- glmFit(de, design, dispersion=de$common.dispersion) # IMPORTANT: fix a common dispersion, instead of using more complex trended or tagwise dispersion

# the counts per million matrix is useful for manual comparisons
# NOT to be used by edgeR
m_cpm <- cpm(de$counts)

l_pos_de <- list()
l_neg_de <- list()

#### SECTION III: performing DE analyses
par(mfrow = c(1,1))
for (contrast in colnames(contrasts)) {
  print("######################################")
  print(contrast)
  # create output directory
  contrast_dir <- paste(de_dir, contrast, sep = "/")
  dir.create(contrast_dir, showWarnings = F)
  
  lrt <- glmLRT(fit, contrast=contrasts[,contrast])
  venn_col <- as.data.frame(decideTestsDGE(lrt, adjust.method = "BH", p.value = FDR_thres))
  colnames(venn_col) <- contrast
  deTab <- topTags(lrt, n=Inf)$table
  deTab <- signif(deTab, digits = 3)
  deGenes <- rownames(deTab)[deTab$FDR < FDR_thres & abs(deTab$logFC) > logFC_thres]
  
  # include cluster category and overlapping reference cluster {BEGIN} # 14 August 2018
  deTab$org <- table_full[rownames(deTab),]$origin
  # include cluster category and overlapping reference cluster {END} # 14 August 2018
  
  # next line was used for example purposes
  #deGenes = rownames(deTab)[deTab$FDR < FDRthres]
  posdeTab <- deTab[(deTab$FDR < FDR_thres & abs(deTab$logFC) >= logFC_thres & deTab$logFC>0),]
  negdeTab <- deTab[(deTab$FDR < FDR_thres & abs(deTab$logFC) >= logFC_thres & deTab$logFC<0),]
  
  nondeTab <- deTab[(deTab$FDR >= FDR_thres),] ## 14 Nov 2018
  comdeTab <- deTab[!(deTab$FDR < FDR_thres & abs(deTab$logFC) >= logFC_thres & deTab$logFC>0),] # com stands for complement, complement set relative to enriched 26 Nov 2018

  # host noise mapped as parasite {BEGIN} # 19 Nov 2018 
  # make a subset of non-DE clusters with only parasite ones
  parasite_nondeTab <- nondeTab[nondeTab$org == parasite_label,]
  # sort from most expressed to less expressed
  parasite_nondeTab <- parasite_nondeTab[order(parasite_nondeTab$logCPM, decreasing = T),]
  top_expressed_parasite_nondeTab <- table_full[rownames(head(parasite_nondeTab, n=10)),]$origin
  # host noise mapped as parasite {END} # 19 Nov 2018
  
  parasite_posdeTab <- posdeTab[posdeTab$org == parasite_label,]
  
  posDE <- nrow(posdeTab)
  print("high expression in treatment:")
  print(posDE)
  negDE <- nrow(negdeTab)
  print("high expression in negative control:")
  print(negDE)
  
  # keep track of the NUMBER of DE genes in each contrast
  l_pos_de[[contrast]] <- posDE
  l_neg_de[[contrast]] <- negDE
  
  # we keep track of the sequences found positively DE in all contrasts
  venn_pos_l[[contrast]] <- rownames(posdeTab)
  venn_neg_l[[contrast]] <- rownames(negdeTab)
  
  # we merge counts and expansion/losses statistical information
  posdeTab <- cbind(posdeTab, signif(m_cpm[rownames(posdeTab),], digits = 3))
  negdeTab <- cbind(negdeTab, signif(m_cpm[rownames(negdeTab),], digits = 3))
  
  # define and create threshold directories
  thres_dir <- paste(contrast_dir, paste("FDR", FDR_thres, "logFC", logFC_thres, sep = "_"), sep = "/")
  dir.create(thres_dir, showWarnings = F)
  
  # writing ouput tables
  write.table(deTab, file=paste(thres_dir, "DE_sRNA_edgeR.txt", sep = "/"), sep="\t", quote=F)
  write.table(posdeTab, file=paste(thres_dir, paste(contrast, "treatment_sRNAs_edgeR.tab", sep = "_"), sep = "/"), sep="\t", quote=F)
  write.table(negdeTab, file=paste(thres_dir, paste(contrast, "control_sRNAs_edgeR.tab", sep = "_"), sep = "/"), sep="\t", quote=F)
  
  # # assembled contigs section {BEGIN} 22 Nov 2018
  # GR_parasite_unique$contig <- as.character(GR_parasite_unique$contig)
  # posdeTab_parasite <- posdeTab[posdeTab$org == parasite_label,]
  # GR_parasite_unique_DE <- GR_parasite_unique[GR_parasite_unique$contig %in% rownames(posdeTab_parasite),]
  # 
  # nondeTab_parasite <- nondeTab[nondeTab$org == parasite_label,]
  # GR_parasite_unique_NDE <- GR_parasite_unique[GR_parasite_unique$contig %in% rownames(nondeTab_parasite),]
  # 
  # GR_parasite_multi$contig <- as.character(GR_parasite_multi$contig)
  # #GR_parasite_multi$counts <- countOverlaps(query = GR_parasite_multi, subject=GR_bam)
  # 
  # # order possible mapping sites from most to less expressed 
  # GR_parasite_multi <- GR_parasite_multi[order(GR_parasite_multi$contig, GR_parasite_multi$counts, decreasing = T)]
  # # now remove duplicated entries
  # GR_parasite_multi <- GR_parasite_multi[!duplicated(GR_parasite_multi$contig)]
  # 
  # GR_parasite_multi_DE <- GR_parasite_multi[GR_parasite_multi$contig %in% rownames(posdeTab_parasite),]
  # GR_parasite_multi_NDE <- GR_parasite_multi[GR_parasite_multi$contig %in% rownames(nondeTab_parasite),]
  # mcols(GR_parasite_multi_DE)$counts <- NULL
  # mcols(GR_parasite_multi_NDE)$counts <- NULL
  # 
  # # merge multiple mapping and unique mapping GRanges objects for parasite into a single one
  # GR_parasite_DE <- c(GR_parasite_unique_DE) #, GR_parasite_multi_DE)
  # GR_parasite_NDE <- c(GR_parasite_unique_NDE) #, GR_parasite_multi_NDE)
  # 
  # df_bed_DE <- data.frame(row.names = NULL, chrom = as.character(seqnames(GR_parasite_DE)), chromStart=start(GR_parasite_DE), chromEnd = end(GR_parasite_DE), contig = GR_parasite_DE$contig, strand= strand(GR_parasite_DE), org= parasite_label)
  # df_bed_NDE <- data.frame(row.names = NULL, chrom = as.character(seqnames(GR_parasite_NDE)), chromStart=start(GR_parasite_NDE), chromEnd = end(GR_parasite_NDE), contig = GR_parasite_NDE$contig, strand= strand(GR_parasite_NDE), org= parasite_label)
  # # here ends the usage of mappings with their original size for later extraction of reads
  # # assembled contigs section {END} 22 Nov 2018
  # 
  # # resize {BEGIN} 25 Nov 2018
  # GR_parasite_multi <- readRDS(file = GR_multi_file) # read again GR_parasite multi, as I modified it in a way that can't be recovered
  # GR_parasite_unique <- resize(GR_parasite_unique, width=1, fix="center")
  # GR_parasite_multi <- resize(GR_parasite_multi, width=1, fix="center")
  # 
  # GR_parasite_unique_DE <- GR_parasite_unique[GR_parasite_unique$contig %in% rownames(posdeTab_parasite),]
  # GR_parasite_unique_NDE <- GR_parasite_unique[GR_parasite_unique$contig %in% rownames(nondeTab_parasite),]
  # 
  # GR_parasite_unique_DE$logCPM <- posdeTab_parasite[match(GR_parasite_unique_DE$contig, rownames(posdeTab_parasite)),]$logCPM # 25 Nov 2018
  # GR_parasite_unique_DE$type <- GR_annot[subjectHits(findOverlaps(GR_parasite_unique_DE, GR_annot))]$type # 25 Nov 2018
  # GR_parasite_unique_NDE$logCPM <- nondeTab_parasite[match(GR_parasite_unique_NDE$contig, rownames(nondeTab_parasite)),]$logCPM # 25 Nov 2018
  # GR_parasite_unique_NDE$type <- GR_annot[subjectHits(findOverlaps(GR_parasite_unique_NDE, GR_annot))]$type # 25 Nov 2018
  # 
  # GR_parasite_unique_NDE$cpm <- 2^GR_parasite_unique_NDE$logCPM # 25 Nov 2018  
  # GR_parasite_unique_DE$cpm <- 2^GR_parasite_unique_DE$logCPM # 25 Nov 2018
  # type_nondeTab <- aggregate(x = GR_parasite_unique_NDE$cpm, by=list(type = GR_parasite_unique_NDE$type), FUN= sum)
  # type_posdeTab <- aggregate(x = GR_parasite_unique_DE$cpm, by=list(type = GR_parasite_unique_DE$type), FUN= sum)
  # 
  # # dealing with multimapping contigs {BEGIN}
  # GR_parasite_multi_DE <- GR_parasite_multi[GR_parasite_multi$contig %in% rownames(posdeTab_parasite),]
  # GR_parasite_multi_DE <- GR_parasite_multi_DE[order(GR_parasite_multi_DE$contig)]
  # GR_parasite_multi_NDE <- GR_parasite_multi[GR_parasite_multi$contig %in% rownames(nondeTab_parasite),]
  # GR_parasite_multi_NDE <- GR_parasite_multi_NDE[order(GR_parasite_multi_NDE$contig)]
  # # getting annotations for all mapping sites
  # GR_parasite_multi_DE$type <- GR_annot[subjectHits(findOverlaps(GR_parasite_multi_DE, GR_annot))]$type # 25 Nov 2018
  # GR_parasite_multi_NDE$type <- GR_annot[subjectHits(findOverlaps(GR_parasite_multi_NDE, GR_annot))]$type # 25 Nov 2018
  # print("time taken to assign counts to DE multiple mapping contigs")
  # print(system.time({
  #   GR_parasite_multi_DE$cpm <- 0
  #   for(contig in unique(GR_parasite_multi_DE$contig)){
  #     current_cpm <- 2^posdeTab_parasite[contig,]$logCPM
  #     GR_parasite_multi_DE_sub <- GR_parasite_multi_DE[GR_parasite_multi_DE$contig == contig]
  #     num_mapping_sites <- length(GR_parasite_multi_DE_sub)
  #     cpm_per_site <- current_cpm/num_mapping_sites
  #     GR_parasite_multi_DE[GR_parasite_multi_DE$contig == contig]$cpm <- cpm_per_site
  #   }
  # }))
  # print("time taken to assign counts to non-DE multiple mapping contigs")
  # print(system.time({
  #   GR_parasite_multi_NDE$cpm <- 0
  #   for(contig in unique(GR_parasite_multi_NDE$contig)){
  #     current_cpm <- 2^nondeTab_parasite[contig,]$logCPM
  #     GR_parasite_multi_NDE_sub <- GR_parasite_multi_NDE[GR_parasite_multi_NDE$contig == contig]
  #     num_mapping_sites <- length(GR_parasite_multi_NDE_sub)
  #     cpm_per_site <- current_cpm/num_mapping_sites
  #     GR_parasite_multi_NDE[GR_parasite_multi_NDE$contig == contig]$cpm <- cpm_per_site
  #   }
  # }))
  # 
  # type_nondeTab_multi <- aggregate(x = GR_parasite_multi_NDE$cpm, by=list(type = GR_parasite_multi_NDE$type), FUN= sum)
  # type_posdeTab_multi <- aggregate(x = GR_parasite_multi_DE$cpm, by=list(type = GR_parasite_multi_DE$type), FUN= sum)
  # # dealing with multimapping contigs {END}
  # 
  # rownames(type_posdeTab_multi) <- type_posdeTab_multi$type # 25 Nov 2018
  # tab_type_posdeTab_multi <- type_posdeTab_multi[,"x"]
  # names(tab_type_posdeTab_multi) <- rownames(type_posdeTab_multi)
  # rownames(type_nondeTab_multi) <- type_nondeTab_multi$type # 25 Nov 2018
  # tab_type_nondeTab_multi <- type_nondeTab_multi[,"x"]
  # names(tab_type_nondeTab_multi) <- rownames(type_nondeTab_multi)
  # 
  # rownames(type_posdeTab) <- type_posdeTab$type # 25 Nov 2018
  # tab_type_posdeTab <- type_posdeTab[,"x"]
  # names(tab_type_posdeTab) <- rownames(type_posdeTab)
  # rownames(type_nondeTab) <- type_nondeTab$type # 25 Nov 2018
  # tab_type_nondeTab <- type_nondeTab[,"x"]
  # names(tab_type_nondeTab) <- rownames(type_nondeTab)
  # 
  # # Gathering annotation counts for both unique and multi-mapping annotations for non-DE contigs
  # de_annot_ids <- unique(names(c(tab_type_nondeTab, tab_type_nondeTab_multi)))
  # tab_type_nondeTab_all <- numeric(0)
  # c <- 1
  # for(de_annot_id in de_annot_ids){
  #   current_val <- 0
  #   if(de_annot_id %in% names(tab_type_nondeTab)){
  #     current_val <- current_val + tab_type_nondeTab[de_annot_id]
  #   }
  #   if(de_annot_id %in% names(tab_type_nondeTab_multi)){
  #     current_val <- current_val + tab_type_nondeTab_multi[de_annot_id]
  #   }
  #   tab_type_nondeTab_all[de_annot_id] <- current_val
  #   c <- c +1
  # }
  # 
  # # Gathering annotation counts for both unique and multi-mapping annotations for DE contigs
  # de_annot_ids <- unique(names(c(tab_type_posdeTab, tab_type_posdeTab_multi)))
  # tab_type_posdeTab_all <- numeric(0)
  # c <- 1
  # for(de_annot_id in de_annot_ids){
  #   current_val <- 0
  #   if(de_annot_id %in% names(tab_type_posdeTab)){
  #     current_val <- current_val + tab_type_posdeTab[de_annot_id]
  #   }
  #   if(de_annot_id %in% names(tab_type_posdeTab_multi)){
  #     current_val <- current_val + tab_type_posdeTab_multi[de_annot_id]
  #   }
  #   tab_type_posdeTab_all[de_annot_id] <- current_val
  #   c <- c +1
  # }
  # 
  # l_annot <- list()
  # l_annot[["DE parasite"]] <- tab_type_posdeTab_all
  # l_annot[["non-DE parasite"]] <- tab_type_nondeTab_all
  # 
  # # format as matrix
  # all_categories <- unique(unlist(lapply(l_annot, names)))
  # m_annot <- matrix(data = NA, nrow = length(l_annot), ncol = length(all_categories))
  # rownames(m_annot) <- names(l_annot)
  # colnames(m_annot) <- all_categories
  # 
  # for(row in rownames(m_annot)){
  #   m_annot[row,] <- l_annot[[row]][all_categories]
  # }
  # 
  # out_file <- paste(thres_dir, paste("genome_categories", parasite_label, "txt", sep = "."), sep = "/") # 25 Nov 2018
  # write.table(m_annot, out_file, quote = FALSE, sep = "\t")
  # # resize {END} 25 Nov 2018
  # 
  # # write bed files of DE clusters
  # write.table(df_bed_DE, file=paste(thres_dir, paste("DE", paste(parasite_label, "bed", sep = "."), sep = "_"), sep = "/"), sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
  # write.table(df_bed_NDE, file=paste(thres_dir, paste("non-DE", paste(parasite_label, "bed", sep = "."),  sep = "_"), sep = "/"), sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
  # 
  # # write bed files of DE clusters
  # idx_tmp <- match(rownames(posdeTab[posdeTab$org == parasite_label,]), seqnames(contigs_fa_index))
  # df_ids_DE_parasite <- data.frame(rownames(posdeTab[posdeTab$org == parasite_label,]), from = start(contigs_fa_index[idx_tmp]), to = width(contigs_fa_index[idx_tmp]))
  # 
  # idx_tmp <- match(rownames(comdeTab[comdeTab$org == parasite_label,]), seqnames(contigs_fa_index))
  # df_ids_NDE_parasite <- data.frame(rownames(comdeTab[comdeTab$org == parasite_label,]), from = start(contigs_fa_index[idx_tmp]), to = width(contigs_fa_index[idx_tmp]))
  # 
  # idx_tmp <- match(rownames(posdeTab[posdeTab$org == host_label,]), seqnames(contigs_fa_index))
  # df_ids_DE_host <- data.frame(rownames(posdeTab[posdeTab$org == host_label,]), from = start(contigs_fa_index[idx_tmp]), to = width(contigs_fa_index[idx_tmp]))
  # 
  # idx_tmp <- match(rownames(comdeTab[comdeTab$org == host_label,]), seqnames(contigs_fa_index))
  # df_ids_NDE_host <- data.frame(rownames(comdeTab[comdeTab$org == host_label,]), from = start(contigs_fa_index[idx_tmp]), to = width(contigs_fa_index[idx_tmp]))
  # 
  # write.table(df_ids_DE_parasite, file=paste(thres_dir, paste("ids", paste("DE", parasite_label, sep = "_"), sep = "."), sep = "/"), sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
  # write.table(df_ids_NDE_parasite, file=paste(thres_dir, paste("ids", paste("non-DE", parasite_label,  sep = "_"), sep = "."), sep = "/"), sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
  # write.table(df_ids_DE_host, file=paste(thres_dir, paste("ids", paste("DE", host_label, sep = "_"), sep = "."), sep = "/"), sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
  # write.table(df_ids_NDE_host, file=paste(thres_dir, paste("ids", paste("non-DE", host_label,  sep = "_"), sep = "."), sep = "/"), sep="\t", quote=F, row.names = FALSE, col.names = FALSE)
  # 
  
  # getting IDs and raw counts average {BEGIN} # 30 Dic 2018
  ids_DE_parasite <- rownames(posdeTab[posdeTab$org == parasite_label,])
  ids_NDE_parasite <- rownames(comdeTab[comdeTab$org == parasite_label,])
  ids_DE_host <- rownames(posdeTab[posdeTab$org == host_label,])
  ids_NDE_host <- rownames(comdeTab[comdeTab$org == host_label,])
  ids_DE_ambiguous <- rownames(posdeTab[posdeTab$org == ambiguous_label,])
  ids_NDE_ambiguous <- rownames(comdeTab[comdeTab$org == ambiguous_label,])
  ids_DE_unmap <- rownames(posdeTab[posdeTab$org == unmapped_label,])
  ids_NDE_unmap <- rownames(comdeTab[comdeTab$org == unmapped_label,])
  
  parasite_cols_regex <- "ves"
  idx_parasite_cols <- grepl(pattern = parasite_cols_regex, colnames(m_count))
  
  avg_raw_counts_DE_parasite <- sum(apply(m_count[ids_DE_parasite,idx_parasite_cols], 1, mean))
  avg_raw_counts_NDE_parasite <- sum(apply(m_count[ids_NDE_parasite,idx_parasite_cols], 1, mean))
  avg_raw_counts_DE_host <- sum(apply(m_count[ids_DE_host,idx_parasite_cols], 1, mean))
  avg_raw_counts_NDE_host <- sum(apply(m_count[ids_NDE_host,idx_parasite_cols], 1, mean))
  # particular to contigs and reads
  avg_raw_counts_DE_ambiguous <- sum(apply(m_count[ids_DE_ambiguous,idx_parasite_cols], 1, mean))
  avg_raw_counts_NDE_ambiguous <- sum(apply(m_count[ids_NDE_ambiguous,idx_parasite_cols], 1, mean))
  avg_raw_counts_DE_unmap <- sum(apply(m_count[ids_DE_unmap,idx_parasite_cols], 1, mean))
  avg_raw_counts_NDE_unmap <- sum(apply(m_count[ids_NDE_unmap,idx_parasite_cols], 1, mean))
  
  avg_cpm_DE_parasite <- sum(apply(m_cpm[ids_DE_parasite,idx_parasite_cols], 1, mean))
  avg_cpm_NDE_parasite <- sum(apply(m_cpm[ids_NDE_parasite,idx_parasite_cols], 1, mean))
  avg_cpm_DE_host <- sum(apply(m_cpm[ids_DE_host,idx_parasite_cols], 1, mean))
  avg_cpm_NDE_host <- sum(apply(m_cpm[ids_NDE_host,idx_parasite_cols], 1, mean))
  # particular to contigs and reads
  avg_cpm_DE_ambiguous <- sum(apply(m_cpm[ids_DE_ambiguous,idx_parasite_cols], 1, mean))
  avg_cpm_NDE_ambiguous <- sum(apply(m_cpm[ids_NDE_ambiguous,idx_parasite_cols], 1, mean))
  avg_cpm_DE_unmap <- sum(apply(m_cpm[ids_DE_unmap,idx_parasite_cols], 1, mean))
  avg_cpm_NDE_unmap <- sum(apply(m_cpm[ids_NDE_unmap,idx_parasite_cols], 1, mean))
  # getting IDs and raw counts average {END} # 30 Dic 2018
  
  # reporting counts comparable with Figure 5 {BEGIN} # 29 Mar 2019
  raw_counts_DE_parasite <- sum(m_count[ids_DE_parasite,idx_parasite_cols])
  raw_counts_NDE_parasite <- sum(m_count[ids_NDE_parasite,idx_parasite_cols])
  raw_counts_DE_host <- sum(m_count[ids_DE_host,idx_parasite_cols])
  raw_counts_NDE_host <- sum(m_count[ids_NDE_host,idx_parasite_cols])
  raw_counts_DE_ambiguous <- sum(m_count[ids_DE_ambiguous,idx_parasite_cols])
  raw_counts_NDE_ambiguous <- sum(m_count[ids_NDE_ambiguous,idx_parasite_cols])
  mmhp_counts <- matrix(data = 
                          c(raw_counts_DE_parasite, raw_counts_NDE_parasite, 
                            raw_counts_DE_host, raw_counts_NDE_host, 
                            raw_counts_DE_ambiguous, raw_counts_NDE_ambiguous), nrow = 2, ncol = 3)
  rownames(mmhp_counts) <- c("DE", "NDE")
  colnames(mmhp_counts) <- c(parasite_label, host_label, ambiguous_label)
  raw_counts_file <- paste(thres_dir, "mmhp_raw_counts_per_org.txt", sep = "/")
  write.table(mmhp_counts, file = raw_counts_file, quote = FALSE, sep = "\t")
  stop("miau")
  # reporting counts comparable with Figure 5 {END} # 29 Mar 2019
  
  
  paired_colors <- brewer.pal(n=10, "Paired") ## 14 Nov 2018 
  colors <- ifelse((deTab$logFC > logFC_thres & deTab$FDR < FDR_thres & deTab$org == parasite_label), paired_colors[10],  #### 27 Oct 2018
                   ifelse((deTab$logFC > logFC_thres & deTab$FDR < FDR_thres & deTab$org == host_label), paired_colors[4], 
                          ifelse(deTab$FDR > FDR_thres & deTab$org == parasite_label, paired_colors[9],
                                 paired_colors[3])))
  colors_transparent <- paste(colors, "4C", sep = "")
  
  plot_main <- paste(paste(DEA_units, read_assignment_method_label), contrast, paste(f2si2(nrow(deTab)), DEA_units, "tested", sep = " "), paste("FDR", FDR_thres, "min_cpm",min_cpm, "min_libs", min_libs, sep = " "), sep = "\n")
  plot_main_simplified <- paste(DEA_units, contrast, paste(f2si2(nrow(deTab)), DEA_units, "tested", sep = " "), paste("FDR", FDR_thres, "min_cpm",min_cpm, "min_libs", min_libs, sep = " "), sep = "\n")
  
  # MA plot showing enriched (orange) and depleted (deepskyblue) clusters
  print("MA plot showing enriched (orange) and depleted (deepskyblue) clusters")
  png(filename = paste(thres_dir, paste(contrast, "DE_logFC_and_logCPM.png", sep = "_"), sep = "/"), width = 720, height = 720)
  # bottom, left, top and right
  par(mar= c(5, 5, 5, 2))
  plot(deTab$logCPM, deTab$logFC, col=colors_transparent, main=plot_main, ylab = "log2FC", xlab="log2CPM", pch=1, lwd=3, cex=1.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
  legend("topright", legend = c(   #### 27 Oct 2018
    paste("enriched", parasite_label, "(",  table(posdeTab$org)[parasite_label], ")", sep= " "),
    paste("enriched", host_label, "(",  table(posdeTab$org)[host_label], ")", sep= " "),
    paste("non-DE", parasite_label, "(", table(nondeTab$org == parasite_label)[["TRUE"]], ")", sep = " "), ## 14 Nov 2018 
    paste("non-DE", host_label, "(", table(nondeTab$org == host_label)[["TRUE"]], ")", sep = " ")),## 14 Nov 2018 
    col = c(paired_colors[10], paired_colors[4], paired_colors[9], paired_colors[3]), pch=1, cex=1.5, pt.lwd =3, pt.cex = 2)
  #abline(a = NULL, b=NULL, h = 0, v = NULL, col = "yellow", lwd=4, lty = "dashed")
  dev.off()
  
  orgs_order <- c("h","p", "a") # 4 Dic 2018 {BEGIN}
  # MA plot showing enriched (orange) and depleted (deepskyblue) clusters
  print("MA plot showing parasite last")
  png(filename = paste(thres_dir, paste(contrast, "DE_logFC_and_logCPM_parasite_on_top.png", sep = "_"), sep = "/"), width = 720, height = 720)
  # bottom, left, top and right
  par(mar= c(5, 5, 6, 2))
  plot(NULL, xlim=c(min(deTab$logCPM),max(deTab$logCPM)), ylim=c(min(deTab$logFC),max(deTab$logFC)), ylab="log2FC", xlab="log2CPM", main=plot_main_simplified, cex=1.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
  for(org in orgs_order){
    
    deTab_sub <- deTab[deTab$org == org,]
    
    colors_sub <- ifelse(test = (deTab_sub$logFC > logFC_thres & deTab_sub$FDR < FDR_thres & deTab_sub$org == parasite_label), yes = "#386CB0",  #### 27 Oct 2018
                    no = ifelse(test = (deTab_sub$logFC > logFC_thres & deTab_sub$FDR < FDR_thres & deTab_sub$org == host_label), yes = "#33A02C",
                      ifelse(test = (deTab_sub$logFC > logFC_thres & deTab_sub$FDR < FDR_thres & deTab_sub$org == ambiguous_label), yes = "#8C6BB1",
                        ifelse(test = deTab_sub$FDR > FDR_thres & deTab_sub$org == ambiguous_label, yes = paired_colors[9], 
                          ifelse(test = deTab_sub$FDR > FDR_thres & deTab_sub$org == parasite_label, yes = paired_colors[1],
                            no = paired_colors[3])))))
    colors_transparent_sub <- paste(colors, "4C", sep = "")
    
    points(deTab_sub$logCPM, deTab_sub$logFC, col=colors_sub, pch=1, lwd=3, cex=2)
    
  } 
  legend("bottomright", legend = c(   #### 27 Oct 2018
    paste("Up", parasite_label, "(",  table(posdeTab$org)[parasite_label], "seqs,", sum(m_count[ids_DE_parasite,]), "counts", ")", sep= " "),
    paste("Up", host_label, "(",  table(posdeTab$org)[host_label], "seqs,", sum(m_count[ids_DE_host,]), "counts", ")", sep= " "),
    paste("Up", ambiguous_label, "(",  table(posdeTab$org)[ambiguous_label], "seqs,", sum(m_count[ids_DE_ambiguous,]), "counts", ")", sep= " "),
    paste("non-Up", parasite_label, "(", table(nondeTab$org == parasite_label)[["TRUE"]], "seqs,", sum(m_count[ids_NDE_parasite,]), "counts", ")", sep = " "), ## 14 Nov 2018 
    paste("non-Up", host_label, "(", table(nondeTab$org == host_label)[["TRUE"]], "seqs,", sum(m_count[ids_NDE_host,]), "counts", ")", sep = " "),
    paste("non-Up", ambiguous_label, "(", table(nondeTab$org == ambiguous_label)[["TRUE"]], "seqs,", sum(m_count[ids_NDE_ambiguous,]), "counts", ")", sep = " ")),## 14 Nov 2018 
    col = c("#386CB0", "#33A02C", "#8C6BB1", paired_colors[1], paired_colors[3], paired_colors[9]), pch=1, cex=1.5, pt.lwd =3, pt.cex = 2)
  #abline(a = NULL, b=NULL, h = 0, v = NULL, col = "yellow", lwd=4, lty = "dashed")
  dev.off()
  # 4 Dic 2018 {BEGIN}
  
  # MA plot showing enriched (orange) and depleted (deepskyblue) clusters zoom # 3 dic 2018
  pdf(file = paste(thres_dir, paste(contrast, "DE_logFC_and_logCPM_zoom.pdf", sep = "_"), sep = "/"), width = 7, height = 4)
  # bottom, left, top and right
  par(mar= c(5, 5, 5, 2))
  plot(deTab$logCPM, deTab$logFC, col=colors_transparent, main=plot_main, ylab = "log2FC", xlab="log2CPM", pch=1, lwd=3, cex=1, cex.axis=1, cex.lab=1, cex.main=1, ylim = c(0, max(deTab$logFC)))
  legend("topright", legend = c(   #### 27 Oct 2018
    paste("enriched", parasite_label, "(",  table(posdeTab$org)[parasite_label], ")", sep= " "),
    paste("enriched", host_label, "(",  table(posdeTab$org)[host_label], ")", sep= " "),
    paste("non-DE", parasite_label, "(", table(nondeTab$org == parasite_label)[["TRUE"]], ")", sep = " "), ## 14 Nov 2018 
    paste("non-DE", host_label, "(", table(nondeTab$org == host_label)[["TRUE"]], ")", sep = " ")),## 14 Nov 2018 
    col = c(paired_colors[10], paired_colors[4], paired_colors[9], paired_colors[3]), pch=1, cex=1, pt.lwd =3, pt.cex = 2)
  #abline(a = NULL, b=NULL, h = 0, v = NULL, col = "yellow", lwd=4, lty = "dashed")
  dev.off()
  
  # getting tables {BEGIN} # 26 Nov 2018
  # number of DEA features (contigs, clusters, reads)
  posdeTab$c <- 1
  comdeTab$c <- 1
  df_posdeTab_org_units <- aggregate(x = posdeTab$c, by=list(org = posdeTab$org), FUN= sum)
  df_comdeTab_org_units <- aggregate(x = comdeTab$c, by=list(org = comdeTab$org), FUN= sum)
  df_org_DE_units <- data.frame(row.names = df_posdeTab_org_units$org, enriched = df_posdeTab_org_units$x, complement = df_comdeTab_org_units$x )
  
  # actual counts in each set set: DE/NDE, parsite/host
  df_posdeTab_org_reads <- aggregate(x = posdeTab$logCPM, by=list(org = posdeTab$org), FUN= sum)
  df_comdeTab_org_reads <- aggregate(x = comdeTab$logCPM, by=list(org = comdeTab$org), FUN= sum)
  df_org_DE_reads <- data.frame(row.names = df_posdeTab_org_reads$org, enriched = df_posdeTab_org_reads$x, complement = df_comdeTab_org_reads$x )
  
  out_table <- paste(thres_dir, "df_org_DE_units.txt", sep = "/")
  write.table(df_org_DE_units, file = out_table, quote = FALSE, sep = "\t")
  out_table <- paste(thres_dir, "df_org_DE_reads.txt", sep = "/")
  write.table(df_org_DE_reads, file = out_table, quote = FALSE, sep = "\t")
  # getting tables {END} # 26 Nov 2018
  
  legend_names <- c(parasite_label, host_label)
  # MA plot with all clusters categories: poly-P, mono-P, other, outside_hit
  col_cluster_type_transparency <- paste(table_full[rownames(deTab),]$color, "4C", sep = "") # adding transparency # 1 Aug 2018
  png(filename = paste(thres_dir, paste(contrast, "cluster_type_DE_logFC_and_logCPM.png", sep = "_"), sep = "/"), width = 720, height = 720)
  # bottom, left, top and right
  par(mar= c(5, 5, 5, 2))
  plot(deTab$logCPM, deTab$logFC, col=col_cluster_type_transparency, main=plot_main, ylab = "log2FC", xlab="log2CPM", pch=1, lwd=6, cex=1.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
  #points(miR_deTab$logCPM, miR_deTab$logFC, col=miR_deTab$color, pch=1, cex=1.5, lwd=6)
  legend("topright", legend = legend_names, col = df_color[legend_names,], pch=16, cex=2)
  abline(a = NULL, b=NULL, h = 0, v = NULL, col = "yellow", lwd=4, lty = "dashed")
  dev.off()

  xlims <- c(min(deTab$logCPM), max(deTab$logCPM))  
  ylims <- c(min(deTab$logFC), max(deTab$logFC))
  # plots for individual organisms {BEGIN}
  # MA plot for each individual org clusters: poly-P, mono-P, other, outside_hit
  for (cat in unique(table_full$origin)){
    org_ids <- rownames(table_full[table_full$origin == cat,])
    
    # build an index for those clusters under current organism
    idx_org <- match(org_ids, rownames(deTab))
    # remove NAs
    idx_org <- idx_org[!is.na(idx_org)]
    
    if(length(idx_org) != 0){
      org_deTab <- deTab[idx_org,]
      
      plot_cat_main <- paste(DEA_units, contrast, paste(f2si2(nrow(org_deTab)), "clusters for organism", sep = " "), sep = "\n")
      
      png(filename = paste(thres_dir, paste(contrast, "cluster_type", cat, "DE_logFC_and_logCPM.png", sep = "_"), sep = "/"), width = 720, height = 720)
      # bottom, left, top and right
      par(mar= c(5, 5, 5, 2))
      plot(org_deTab$logCPM, org_deTab$logFC, col=table_full[rownames(org_deTab),]$color, main=plot_cat_main, ylab = "log2FC", xlab="log2CPM", pch=1, cex=1.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, ylim = ylims, xlim= xlims)
      #points(miR_deTab$logCPM, miR_deTab$logFC, col=miR_deTab$color, pch=1, cex=1.5, lwd=3)
      legend("topright", legend = rownames(df_color), col = df_color$color, pch=16, cex=2)
      abline(a = NULL, b=NULL, h = 0, v = NULL, col = "yellow", lwd=4, lty = "dashed")
      dev.off()
      
      # zoom for individual org categories # 3 dic 2018
      pdf(file = paste(thres_dir, paste(contrast, "cluster_type", cat, "DE_logFC_and_logCPM_zoom.pdf", sep = "_"), sep = "/"), width = 7, height = 4)
      # bottom, left, top and right
      par(mar= c(5, 5, 5, 2))
      plot(org_deTab$logCPM, org_deTab$logFC, col=table_full[rownames(org_deTab),]$color, main=plot_cat_main, ylab = "log2FC", xlab="log2CPM", pch=1, cex=1, cex.axis=1, cex.lab=1, cex.main=1, ylim = c(0, max(ylims)), xlim= xlims)
      #points(miR_deTab$logCPM, miR_deTab$logFC, col=miR_deTab$color, pch=1, cex=1.5, lwd=3)
      legend("topright", legend = rownames(df_color), col = df_color$color, pch=16, cex=1)
      abline(a = NULL, b=NULL, h = 0, v = NULL, col = "yellow", lwd=4, lty = "dashed")
      dev.off()
      
      out_table <- paste(thres_dir, paste(cat, DEA_units, "ids", sep = "."), sep = "/")
      write.table(org_ids, file = out_table, quote = FALSE, sep = "\t", row.names = F, col.names = F)
      
      up_reg_org_ids <- org_ids[org_ids %in% rownames(posdeTab)]  # 22 Feb 2019
      out_up_reg_file <- paste(thres_dir, paste(cat, DEA_units, "up_reg", sep = "."), sep = "/")  # 22 Feb 2019
      write.table(up_reg_org_ids, file = out_up_reg_file, quote = FALSE, sep = "\t", row.names = F, col.names = F)  # 22 Feb 2019
    }
  }
  # plots for individual organisms {END}
  
  # bins is a factor of the different logFC values to be grouped together given a numer of bins, this is where discretization step occurs
  bins <- cut(deTab$logFC, bin_number)
  
  # add the corresponding bin to each row
  deTab <- cbind(deTab, bins)
  deTab$CPM <- 2^deTab$logCPM # fixing CPM MA bar plots sum # 11 Dic 2018
  
  # correction made on 06 April 2018 {BEGIN}
  # See journal entry: 2018 Week 14/52 Apr 2 - Apr 6. second aggregate, we lost information regarding logCPMs 
  # merge all the rows with the same mapping information and the exact same logFC (a thousands of different possible values)
  deTab_binned <- aggregate(x = deTab$CPM, by= list(org = deTab$org, bins= deTab$bins), FUN= sum) # fixing CPM MA bar plots sum # 11 Dic 2018
  
  # fraction of Ambiguous reads {BEGIN} # 2 August 2018 
  df_org_frac <- aggregate(x = deTab$CPM, by= list(org = deTab$org), FUN= sum) # fixing CPM MA bar plots sum # 11 Dic 2018
  rownames(df_org_frac) <- df_org_frac$org
  colnames(df_org_frac) <- c("org", "CPM")
  
  barplot(df_org_frac[c(parasite_label, host_label),]$CPM/ sum(df_org_frac$CPM), names.arg = c(parasite_label, host_label), col = df_color[c(parasite_label, host_label),], ylab= "% CPM")
  # fraction of Ambiguous reads {END} # 2 August 2018 
  
  # you may split a data.frame by a character or factor column
  l_binned <- split(deTab_binned, deTab_binned$bins)
  # now I need to re-format l_binned to a matrix-like object
  
  # we need to then build a matrix m with rows equal to mapping organisms and columns equal to the number of bins
  rows <- unique(as.character(deTab_binned$org))
  cols <- names(l_binned)
  m_binned <- matrix(data = NA, nrow = length(rows), ncol = length(cols))
  rownames(m_binned) <- rows
  colnames(m_binned) <- cols
  
  # for loop to fill m_binned
  for(col in cols){
    
    for(row in rows){
      if(row %in% l_binned[[col]][,"org"]){
        
        # where_is is a logic vector that help us subset the appropiate row in the data.frame
        where_is <- row == l_binned[[col]][, "org"]
        
        m_binned[row, col] <- l_binned[[col]][where_is, "x"]
      }
      
    }
    
  }
  m_binned[is.na(m_binned)] <- 0
  
  # change the order of rows in m_binned to change the plotting order
  m_binned <- m_binned[c(parasite_label, host_label),]
  # correction made on 6 April 2018 {END}  
  
  # raw stacked barplot
  stacked_barplot_file <- paste(thres_dir, "stacked_barplot.pdf", sep = "/")
  pdf(file = stacked_barplot_file, width = 7, height = 7)
  # make barplot of org & logFC
  # bottom, left, top and right
  par(mar= c(5, 7, 5, 1))
  barplot_main <- paste(DEA_units, paste(parasite_label, host_label, sep = "-"), sep = "\n") # 10 Dic 2018
  barplot(m_binned, las=2, col = df_color[rownames(m_binned),], horiz = TRUE, main = barplot_main)
  legend("top", legend = rownames(m_binned), col = df_color[rownames(m_binned),], pch = 15, pt.cex = 2)
  dev.off()
  
  # Parasite only reads barplot
  stacked_barplot_file <- paste(thres_dir, paste(paste("stacked_barplot", parasite_label, sep="_"), "pdf", sep="."), sep = "/")
  pdf(file = stacked_barplot_file, width = 3.5, height = 7)
  # make barplot of org & logFC
  # bottom, left, top and right
  par(mar= c(5, 7, 5, 1))
  barplot(m_binned[parasite_label,], las=2, col = df_color[parasite_label,], horiz = TRUE)
  legend("topright", legend = parasite_label, col = df_color[parasite_label,], pch = 15, pt.cex = 2)
  dev.off()
  
  # WARNING: i sum absolute values to avoid negative numbers in fraction bar plot, is this wrong?
  m_binned_frac <- t(t(abs(m_binned))/colSums(abs(m_binned)))
  
  # fraction barplot
  stacked_barplot_file <- paste(thres_dir, "stacked_barplot_frac.pdf", sep = "/")
  pdf(file = stacked_barplot_file, width = 2.5, height = 7)
  # make barplot of org & logFC
  # bottom, left, top and right
  par(mar= c(5, 7, 5, 1))
  barplot(m_binned_frac, las=2, col = df_color[rownames(m_binned_frac),], horiz = TRUE)
  legend("bottom", legend = rownames(m_binned), col = df_color[rownames(m_binned),], pch = 15, pt.cex = 2, bg = "white")
  dev.off()
  ## stacked barplot with organism mapping information {END}
  
  
  if(contrast == "MODE-K_vesicle_vs_no_treatment"){
    stop("miau")
  }
  
}

df_de_num <- data.frame(row.names = names(l_pos_de), treatment = unlist(l_pos_de), control = unlist(l_neg_de), stringsAsFactors = FALSE)

df_de_num_file <- paste(de_dir, paste(paste("contrasts_DE_nums_FDR", FDR_thres, "logFC", logFC_thres, sep = "_"), "txt", sep=".") , sep = "/")
write.table(df_de_num, file = df_de_num_file, quote = FALSE, sep = "\t", col.names = TRUE)

