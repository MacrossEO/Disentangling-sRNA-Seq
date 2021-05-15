#03/ago/19
#trim contigs after assembly

rm(list=ls())
suppressMessages(library(protr)) 
suppressMessages(library(dplyr))
suppressMessages(library(Rsamtools))

computer="production"

if(computer=="test"){
  wd="/home/macross/compartida/data/miRNA/mm_hp"
  setwd(wd)
  contig_file = "dea_contigs/inchworm.fasta" #fasta contigs
  all_reads = "dea_contigs/hostPar.fa"
  output = "dea_contigs/contigs.fa" #fasta contigs
  threads = 7

  bowtie_build_bin = "~/Documentos/bin/bowtie-1.2.2/bowtie-build"
  bowtie_bin = "~/Documentos/bin/bowtie-1.2.2/bowtie"
  tally_to_fasta_bin = "~/Dropbox/proyectos/miRNA/bin/tally_to_fasta.pl"
  samtools_bin = "~/Documentos/bin/samtools-1.5/samtools"
}

if(computer=="production"){
  args = commandArgs(trailingOnly=TRUE)
  contig_file = args[1] 
  all_reads = args[2]
  output = args[3]
  threads = args[4]

  bowtie_build_bin = args[5]
  bowtie_bin = args[6]
  tally_to_fasta_bin = args[7]
  samtools_bin = args[8]
}

dir = paste(dirname(contig_file),"tmp_dir",sep="/")#switch when necessary
dir.create(dir)

#file names
sam_file = paste(dir,"all.sam",sep="/")
bam_file = paste(dir,"all.bam",sep="/")
contig_db = paste(dir,"contig_db",sep="/")
  
#bowtie index
bowtie_build_call = paste(bowtie_build_bin,"--threads",threads,"-f",contig_file,contig_db);bowtie_build_call
system(bowtie_build_call,intern = T,ignore.stdout=T,ignore.stderr=T)

#bowtie mapping
bowtie_call = paste(bowtie_bin,"-S -f -v 0 -a --norc -p",threads,contig_db,all_reads,sam_file);bowtie_call
system(bowtie_call,intern = T,ignore.stdout=T,ignore.stderr=T)
  
#sam to bam 
view_call = paste(samtools_bin, "view", sam_file, "-Sb >", bam_file, sep = " ");view_call
system(view_call)
sort_call = paste(samtools_bin, "sort", bam_file, "-o", bam_file, sep = " ");sort_call
system(sort_call)
index_call = paste(samtools_bin, "index", bam_file,sep = " ");index_call
system(index_call)
  
#library(Rsamtools)
# set the index file name
index_file <- paste(bam_file, "bai", sep = ".")
  
# contig counts
what_to_extract <- c("rname","pos","seq")
l_bam <- scanBam(bam_file, param = ScanBamParam(what = what_to_extract, flag = scanBamFlag(isUnmappedQuery = FALSE)))
df_bam <- do.call("data.frame", l_bam)
tmp_table = table(df_bam$rname)
tmp_column <- data.frame(row.names = names(tmp_table), lib = as.numeric(tmp_table), stringsAsFactors = FALSE)

#contigs seqs
contigs = readFASTA(contig_file)
  
#fasta to dataframe
seqs = data.frame(unlist(contigs),stringsAsFactors = F)
  
#check rownames from bam and fasta 
if(nrow(tmp_column) == sum(rownames(seqs)==rownames(tmp_column))){print("rownames are OK, keep working")}
  
#################
##trim contigs###
#################
len_seqs = nchar(seqs[,1])
contig_tbl = cbind(seqs,tmp_column,len_seqs)
  
#get start-end position alignment of each read to contigs 
df_bam$len = nchar(df_bam$seq)
df_bam$end = (df_bam$pos-1) + df_bam$len

#get start-end position aligment of reads to each contig
df_bam %>%
  group_by(rname) %>%
  summarise(start=min(pos),end=max(end)) ->
  df_bam_contigs
  
#trimm those contigs with extra nucleotides
contig_tbl$rname = rownames(contig_tbl)
contig_tbl_start_end = inner_join(contig_tbl,df_bam_contigs,by="rname")
contig_tbl_start_end$trimmed_seqs =  substr(contig_tbl_start_end$unlist.contigs.,contig_tbl_start_end$start,contig_tbl_start_end$end)
contig_tbl_start_end$len_seqs = nchar(contig_tbl_start_end$trimmed_seqs)

#dta frame to tally
tally_out_f = sub(pattern = ".fasta",replacement = "_f.tally",contig_file);tally_out_f
contig_tbl = subset(contig_tbl_start_end,select=c(trimmed_seqs,lib,len_seqs))
write.table(contig_tbl, file = tally_out_f, quote = FALSE, sep = "\t", col.names = F,row.names = F)
  
#tally to fasta
tally_to_fasta_call = paste(tally_to_fasta_bin,"-f",tally_out_f,"-i c >",output);tally_to_fasta_call
system(tally_to_fasta_call)
  
#rm temp files
system(paste("rm -r ",dir))
system(paste("rm",tally_out_f))