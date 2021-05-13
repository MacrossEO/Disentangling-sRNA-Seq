#13/jan/21
#Dr. Obed Ramirez
#map reads to uniqueReads
#print counts
rm(list=ls())
suppressMessages(library(Rsamtools))

computer="production"
if(computer=="test"){
  setwd("/home/macross/compartida/data/miRNA/mm_hp")
  sample = "mapped_data/hostPar/M_ves_24_1.reduced.fa.gz"
  uniqueReads_db = "mapped_data/bowtieIndex/host_hostPar_tag"
  output = "dea_reads/hostPar/M_ves_24_1.counts"
  threads = 1
  
  bowtie_bin = "~/Documentos/bin/bowtie-1.2.2/bowtie"
  samtools_bin = "~/Documentos/bin/samtools-1.5/samtools"
}

if(computer=="production"){
  args = commandArgs(trailingOnly=TRUE)
  sample = args[1]
  uniqueReads_db = args[2]
  output = args[3]
  threads = args[4]

  bowtie_bin = args[5]
  samtools_bin = args[6]
}

tmp_dir = paste(dirname(output),"tmp_dir",sep="/")
dir.create(tmp_dir)

map_reads_to_reads = function(query_fasta,ref_fasta_db,params){
  df_bam = map_bam_df(query_fasta,ref_fasta_db,params)
  
  #===GET PERFECT alignments un_fa (actual_lib) vs (global) unitigs 
  #length from un_fa and length from unitigs(rlen)
  df_bam$rname = as.character(df_bam$rname)
  ref_ids = split_ids(df_bam$rname,",")
  colnames(ref_ids) = c("rid","rcounts","rlength","rtag")
  
  df_bam$qname = as.character(df_bam$qname)
  query_ids = split_ids(df_bam$qname,",")
  colnames(query_ids) = c("qid","qcounts","qlength")
  
  df_bam = cbind(df_bam,query_ids,ref_ids)
  em = df_bam[ref_ids$rlength==query_ids$qlength,]
  print(head(em))
  
  #check point
  u = length(unique(df_bam$qname))
  print(u);print(nrow(em))
  if(u==nrow(em)){print("counts were correctly assingend, keep working :D")}
  
  return(em)
}

map_bam_df = function(query,reference,params){
  #get params
  bowtie_params = params[["bowtie_params"]];bowtie_params
  query_type = params[["query_type"]] ;query_type
  if(query_type == "fasta"){bowtie_params=paste(bowtie_params,"-f")}
  tmp_dir = params[["tmp_dir"]]
  what_to_extract = params[["what_to_extract"]]
  unmapped = params[["unmapped"]]
  
  #file names
  sam_file = paste(tmp_dir,paste0(basename(query),".sam"),sep="/")
  bam_file = paste(tmp_dir,paste0(basename(query),".bam"),sep="/")
  
  #bowtie mapping
  if(unmapped==""){
    bowtie_call = paste(bowtie_bin,bowtie_params,reference,query,sam_file);print(bowtie_call)
  }else{
    unmapped = paste(tmp_dir,unmapped,sep="/")
    bowtie_call = paste(bowtie_bin,bowtie_params,reference,query,sam_file,"--un",unmapped);print(bowtie_call)
  }
  
  system(bowtie_call)
  
  #sam to bam 
  view_call = paste(samtools_bin, "view", sam_file, "-Sb >", bam_file, sep = " ");print(view_call)
  system(view_call)
  sort_call = paste(samtools_bin, "sort", bam_file, "-o", bam_file, sep = " ");sort_call
  system(sort_call)
  index_call = paste(samtools_bin, "index", bam_file,sep = " ");print(index_call)
  system(index_call)
  index_file <- paste(bam_file, "bai", sep = ".")
  
  l_bam <- scanBam(bam_file, param = ScanBamParam(what = what_to_extract, flag = scanBamFlag(isUnmappedQuery = FALSE)))
  df_bam <- do.call("data.frame", l_bam)
  
  system(paste("rm ",bam_file, sam_file))
  return(df_bam)
}

split_ids=function(vector_ids,by){
  ids = strsplit(vector_ids,by)
  ids = data.frame(matrix(unlist(ids), nrow = length(ids), byrow = T),stringsAsFactors = F)
  return(ids)
}


params = list(bowtie_params = paste("-S -v 0 -p",threads,"-a --norc"), 
              query_type="fasta", 
              tmp_dir=tmp_dir, 
              what_to_extract=c("qname","rname"),
              unmapped="")
reads_map = map_reads_to_reads(sample,uniqueReads_db,params)
write.table(reads_map[,c("rid","qcounts")],output,col.names=F,row.names=F,quote=F,sep="\t")
