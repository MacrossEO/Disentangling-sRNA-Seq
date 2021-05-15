#21/dec/2020
#Dr. Obed Ramirez
#get support counts (uniquely mapped reads)

rm(list=ls())
suppressMessages(library(Rsamtools))
suppressMessages(library(dplyr))

server="production"

if(server=="test"){
  wd = "/home/macross/compartida/data/miRNA/mm_hp"
  setwd(wd)
  contigs_db = "dea_contigs/bowtieIndex/contigs"
  samples = "mapped_data/hostPar.fa.gz"
  output = "dea_contigs/hostPar.sppc"
  threads= 7
  
  bowtie_bin = "/home/macross/Documentos/bin/bowtie-1.2.2/bowtie"
  tally_to_fasta_bin = "/home/macross/Dropbox/proyectos/miRNA/bin/tally_to_fasta.pl"
  samtools_bin = "/home/macross/Documentos/bin/samtools-1.10/samtools"
}
if(server=="production"){
  args = commandArgs(trailingOnly=TRUE)
 
  contigs_db= args[1]
  samples= args[2]
  output= args[3]
  threads = args[4] 
  
  bowtie_bin = args[5]
  tally_to_fasta_bin = args[6]
  samtools_bin = args[7]
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
  sam_file = paste0(basename(query),".sam")
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
  
  system(paste("rm",sam_file,bam_file,unmapped))
  return(df_bam)
}

split_ids=function(vector_ids,by){
  ids = strsplit(vector_ids,by)
  ids = data.frame(matrix(unlist(ids), nrow = length(ids), byrow = T),stringsAsFactors = F)
  return(ids)
}

tmp_dir = paste(dirname(output),"tmp_dir",sep="/")
dir.create(tmp_dir)

#===map reads to ref contigs/clusters===#
mbd_params = list(bowtie_params = paste("-S -v 0 -p",threads,"-a --norc"), 
                  query_type="fasta", 
                  tmp_dir=tmp_dir, 
                  what_to_extract=c("qname","rname"),
                  unmapped="")
df_bam = map_bam_df(samples,contigs_db,mbd_params)

#=== contigs-reads table===#
df_bam$qname = as.character(df_bam$qname)
df_bam$rname = as.character(df_bam$rname)
query_ids = split_ids(df_bam$qname,",")
ref_ids = split_ids(df_bam$rname,",")
contigs_reads = data.frame(contig_id=ref_ids$X1,read_id=query_ids$X1,read_counts=as.numeric(query_ids$X2),read_length=as.numeric(query_ids$X3),stringsAsFactors = F)

#=== Unique-mapping reads (umr)===#
df_tab = table(contigs_reads$read_id)
umr = names(df_tab)[df_tab == 1]

#=== base counts===#
umr_contigs_reads = contigs_reads[contigs_reads$read_id %in% umr,]
base_counts = tapply(umr_contigs_reads$read_counts,factor(umr_contigs_reads$contig_id),sum)

contigs_counts = left_join(data.frame(contig_id=unique(contigs_reads$contig_id),stringsAsFactors=F), 
                           data.frame(contig_id=names(base_counts),base_counts=base_counts,stringsAsFactors=F),by="contig_id")
contigs_counts$base_counts[is.na(contigs_counts$base_counts)] = 0

write.table(contigs_counts,output,row.names = F,col.names = T,quote = F,sep="\t")
print("global counts were correctly assigned, :D")
