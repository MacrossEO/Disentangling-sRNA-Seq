#08/jan/2020
#Dr. Obed Ramirez
#read support counts file(uniquely mapped reads)
#distritube multi-mapping reads

rm(list=ls())
suppressMessages(library(Rsamtools))
suppressMessages(library(dplyr))

server="production"

if(server=="test"){
  wd = "/home/macross/compartida/data/miRNA/mm_hp"
  setwd(wd)
  contigs_db = "dea_contigs/bowtieIndex/contigs"
  samples = "mapped_data/hostPar/M_ves_24_1.reduced.fa.gz"
  uniqueReads_db = "mapped_data/bowtieIndex/host_hostPar_tag"
  sppc_file = "dea_contigs/hostPar.sppc"
  output = "dea_contigs/hostPar/M_ves_24_1.counts"
  threads= 7
  
  bowtie_bin = "/home/macross/Documentos/bin/bowtie-1.2.2/bowtie"
  tally_to_fasta_bin = "/home/macross/Dropbox/proyectos/miRNA/bin/tally_to_fasta.pl"
  samtools_bin = "/home/macross/Documentos/bin/samtools-1.10/samtools"
}
if(server=="production"){
  args = commandArgs(trailingOnly=TRUE)
  
  contigs_db= args[1]
  samples= args[2]
  uniqueReads_db = args[3]
  sppc_file = args[4]
  output= args[5]
  threads = args[6] 
  
  bowtie_bin = args[7]
  tally_to_fasta_bin = args[8]
  samtools_bin = args[9]
}

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

tmp_dir = paste(dirname(output),"tmp_dir",sep="/")
dir.create(tmp_dir)

#===map reads to ref contigs/clusters===#
unmapped_reads_fa = paste0(basename(samples),"_un_to_contigs.fa")
params = list(bowtie_params = paste("-S -v 0 -p",threads,"-a --norc"), 
                  query_type="fasta", 
                  tmp_dir=tmp_dir, 
                  what_to_extract=c("qname","rname","seq"),
                  unmapped=unmapped_reads_fa)
df_bam = map_bam_df(samples,contigs_db,params)

#=== contigs-reads table===#
df_bam$qname = as.character(df_bam$qname)
df_bam$rname = as.character(df_bam$rname)
query_ids = split_ids(df_bam$qname,",")
ref_ids = split_ids(df_bam$rname,",")
contigs_reads = data.frame(contig_id=ref_ids$X1,read_id=query_ids$X1,read_counts=as.numeric(query_ids$X2),read_length=as.numeric(query_ids$X3),stringsAsFactors = F)

#=== Unique-mapping reads(umr) / multi-mapping reads(mmr)===#
df_tab = table(contigs_reads$read_id)
umr = names(df_tab)[df_tab == 1]
mmr = names(df_tab)[df_tab > 1]
umr_contigs_reads = contigs_reads[contigs_reads$read_id %in% umr,]
mmr_contigs_reads = contigs_reads[contigs_reads$read_id %in% mmr,]

#=== read spp, assign sppc on mmr_contigs_readsc ===#
informed_contigs_counts = read.table(sppc_file,header=T,stringsAsFactors = F)
informed_contigs_counts = informed_contigs_counts[informed_contigs_counts$contig_id %in% unique(mmr_contigs_reads$contig_id),]
mmr_contigs_reads = merge(mmr_contigs_reads,informed_contigs_counts,by="contig_id")

#=== contigs with very small sppc are discarted (DEPRECATED) ===#
#min_counts = get_min_counts(threshold,base_counts_file)
#mmr_contigs_reads$base_counts[mmr_contigs_reads$base_counts<=min_counts] = 0

#=== Get support counts sum by each mmr_contig ===#
bcs = tapply(mmr_contigs_reads$base_counts,factor(mmr_contigs_reads$read_id),sum)#base_counts sum for each mmr
mmr_contigs_reads = left_join(mmr_contigs_reads,data.frame(bcs=bcs,read_id=names(bcs)),by="read_id")

#===reads whose contigs have non zero supported_reads_sum, otherwise will remain being ambiguous reads===#
nz = mmr_contigs_reads[mmr_contigs_reads$bcs>0,] # 0 or min_counts
nz = nz[order(nz$read_id,nz$base_counts,decreasing = c(T,F)),]

#===calculate the supported counts fraction and SPLIT COUNTS according to===#
nz$bcf = nz$base_counts/nz$bcs
nz$shared = round(nz$bcf * nz$read_counts,0)

#=== Checkpoint ===#
#head(nz,100)
#head(table(nz$base_counts),50)
#===#

#=== base counts===#
base_counts = tapply(umr_contigs_reads$read_counts,factor(umr_contigs_reads$contig_id),sum)
contigs_counts = left_join(data.frame(contig_id=unique(contigs_reads$contig_id),stringsAsFactors=F), 
                           data.frame(contig_id=names(base_counts),base_counts=base_counts,stringsAsFactors=F),by="contig_id")
contigs_counts$base_counts[is.na(contigs_counts$base_counts)] = 0

#=== add the shared reads to base_counts ===#
sr = tapply(nz$shared,factor(nz$contig_id),sum)# get the sum of shared reads
if(length(sr)>0){
  contigs_counts = left_join(contigs_counts, data.frame(contig_id=names(sr),sr_counts=sr,stringsAsFactors=F), by="contig_id")
  contigs_counts$sr_counts[is.na(contigs_counts$sr_counts)] = 0
  contigs_counts$final_counts = contigs_counts$base_counts + contigs_counts$sr_counts 
}else{
  contigs_counts$final_counts = contigs_counts$base_counts
}
final_cc = contigs_counts[,c("contig_id","final_counts")]

#=== get ambiguous reads ===#
amb_reads = mmr_contigs_reads[mmr_contigs_reads$bcs<=0,c("read_id","read_counts")] #0 bcs or min_counts
amb_reads = amb_reads[!duplicated(amb_reads),]

#=== ambiguous reads to fasta ===#
if(nrow(amb_reads)>0){
  df_bam = cbind(df_bam,query_ids)
  amb_reads = df_bam[df_bam$X1 %in% amb_reads$read_id,]
  amb_reads = amb_reads[!duplicated(amb_reads$X1),c("seq","X2","X3")]
  
  amb_tally = paste(tmp_dir,paste0(basename(samples),".tally"),sep="/")
  amb_fa = paste(tmp_dir,paste0(basename(samples),"_amb.fa"),sep="/")
  write.table(amb_reads,amb_tally,row.names=F,col.names=F,quote=F,sep="\t")
  tally2fasta_call = paste(tally_to_fasta_bin,"-f",amb_tally,">",amb_fa);print(tally2fasta_call)
  system(tally2fasta_call)
  system(paste("rm",amb_tally))
}

#=== map amb-reads/unalg-reads to host_hostPar uniqueReads  ===#
unmapped_reads_fa = paste(tmp_dir,unmapped_reads_fa,sep="/")
if(!file.exists(unmapped_reads_fa)){unmapped_reads_fa=""}
if(nrow(amb_reads)==0){amb_fa=""}
if(unmapped_reads_fa!="" | amb_fa!=""){
  unamb_fa = paste(tmp_dir,paste0(basename(samples),"_unamb.fa"),sep="/")
  cat_call = paste("cat",unmapped_reads_fa,amb_fa,">",unamb_fa);print(cat_call)
  system(cat_call)
  
  params = list(bowtie_params = paste("-S -v 0 -p",threads,"-a --norc"), 
                    query_type="fasta", 
                    tmp_dir=tmp_dir, 
                    what_to_extract=c("qname","rname"),
                    unmapped="")
  unamb_map = map_reads_to_reads(unamb_fa,uniqueReads_db,params)
  
  final_cc = rbind(final_cc,data.frame(contig_id=unamb_map$rid,final_counts=unamb_map$qcounts))
  system(paste("rm",unmapped_reads_fa,amb_fa))
}

write.table(final_cc,output,sep="\t",quote=F,row.names=F,col.names=F)

