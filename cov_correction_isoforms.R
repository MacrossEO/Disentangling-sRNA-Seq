#03/oct/18
#second contig coverage correction
#ambigous reads are asigned proportionally to contigs
rm(list=ls())
library(protr)
suppressMessages(library(Rsamtools))
suppressMessages(library(dplyr))

###input vars###
server="mac"

if(server=="mac"){
  args = commandArgs(trailingOnly=TRUE)
  md = args[1]#output_dir
  lib_file = args[2]#fastq_fasta_file
  alg_file = args[3]#alg.fa
  unitigs_file = args[4]
  contig_file <- args[5]
  threshold = as.double(args[6])
  tmp_dir = args[7]
  option = args[8]#check_options
  base_counts_file = args[9]#"base_counts_hostpar.counts"
  #base_counts_file = paste(md,base_counts_file,sep="/")
}
if(server=="linux"){
  md = "/Users/macross/Documents/data/sRNAs_assembly/assembly_eval/mm_hp/mapped/diSrna/DEA_contigs_informed_0.1"#output_dir
  setwd(md)
  lib_file = "mmhp_mm.fa"#fastq_fasta_file
  alg_file = "alg.fa"#alg.fa
  unitigs_file = "unitigs.fa"
  contig_file <- "contigs.fa"
  threshold = "0.1"
  tmp_dir = "tmp_dir"
  option = "write_counts_hostpar"#check_options
  base_counts_file = "empty"#"base_counts_hostpar.counts"
  #base_counts_file = paste(md,base_counts_file,sep="/")
}



bowtie_build_bin = "/Users/macross/Documents/bin/bowtie-1.2.2/bowtie-build"
bowtie_bin = "/Users/macross/Documents/bin/bowtie-1.2.2/bowtie"
tally_to_fasta_bin = "/Users/macross/Dropbox/proyectos/miRNA/bin/tally_to_fasta.pl"
samtools_bin = "/Users/macross/Documents/bin/samtools-1.5/samtools"


check_options = function(option){
  if(option=="informed_SppC"  || option=="informed_SpliT" 
     || option=="write_counts_host" || option == "write_counts_hostpar" || option == "write_counts_par"
     || option == "library_SppC" || option == "library_SpliT"){}
  else{
    print("possible options: 
          informed_SppC
          informed_SpliT
          write_counts_host
          write_counts_hostpar
          write_counts_par
          library_SppC
          library_SpliT")
       stop("incorrect option")
     }
}

tbl_to_fasta = function(tbl_data,fasta_name){
  #tally to fasta
  out_tally = sub(".fa",".tally",fasta_name)
  write.table(tbl_data,out_tally,row.names = F,col.names = F,quote = F,sep="\t")
  command = paste(tally_to_fasta_bin,"-f",out_tally,">",fasta_name);command
  system(command)
}

map_reads_to_reads = function(query_fasta,ref_fasta,dir){
  mbd_params = list(bowtie_params = "-S -v 0 -p 20 -a --norc", 
                    query_type="fasta", 
                    tmp_dir=dir, 
                    what_to_extract=c("qname","rname","seq"))
  df_bam = map_bam_df(query_fasta,ref_fasta,mbd_params)
  
  #===GET PERFECT alignments un_fa (actual_lib) vs (global) unitigs 
  #length from un_fa and length from unitigs(rlen)
  rid = strsplit(as.character(df_bam$rname),",")
  rid = data.frame(matrix(unlist(rid), nrow = length(rid), byrow = T),stringsAsFactors = F)
  df_bam$rlen = as.numeric(rid[,3])
  df_bam$rid = rid[,1]
  df_bam$rcounts = as.numeric(rid[,2])
  

  qid = strsplit(as.character(df_bam$qname),",")
  qid = data.frame(matrix(unlist(qid), nrow = length(qid), byrow = T),stringsAsFactors = F)
  df_bam$qid = qid[,1]
  df_bam$qcounts = as.numeric(qid[,2])
  df_bam$qlen = nchar(df_bam$seq)
  
  em = df_bam[df_bam$qlen==df_bam$rlen,]#exact match
  print(head(em))
  
  #check point
  u = length(unique(df_bam$qname))
  print(u);print(nrow(em))
  if(u==nrow(em)){print("counts were correctly assingend, keep working :D")}
  
  return(em)
}

map_bam_df = function(query,reference,mbd_params){
  #get params
  bowtie_params = mbd_params[["bowtie_params"]];bowtie_params
  query_type = mbd_params[["query_type"]] ;query_type
  if(query_type == "fasta"){bowtie_params=paste(bowtie_params,"-f")}
  tmp_dir = mbd_params[["tmp_dir"]]
  what_to_extract = mbd_params[["what_to_extract"]]
  unmmaped = mbd_params[["unmmaped"]]
  
  #file names
  sam_file = paste(tmp_dir,"tmp.sam",sep="/")
  bam_file = paste(tmp_dir,"tmp.bam",sep="/")
  tmp_index = paste(tmp_dir,"tmp_index",sep="/")
  unmmaped = paste(tmp_dir,unmmaped,sep="/")
  
  #construct bowtie index
  bowtie_build_call = paste(bowtie_build_bin,reference,tmp_index,"--threads 20");bowtie_build_call
  system(bowtie_build_call,intern = T,ignore.stdout=T,ignore.stderr=T)
  
  #bowtie mapping
  if(unmmaped==""){
    bowtie_call = paste(bowtie_bin,bowtie_params,tmp_index,query,sam_file);bowtie_call
  }else{
    bowtie_call = paste(bowtie_bin,bowtie_params,tmp_index,query,sam_file,"--un",unmmaped);bowtie_call
  }
  
  system(bowtie_call)
  
  #sam to bam 
  view_call = paste(samtools_bin, "view", sam_file, "-Sb >", bam_file, sep = " ");view_call
  system(view_call)
  sort_call = paste(samtools_bin, "sort", bam_file, "-o", bam_file, sep = " ");sort_call
  system(sort_call)
  index_call = paste(samtools_bin, "index", bam_file,sep = " ");index_call
  system(index_call)
  index_file <- paste(bam_file, "bai", sep = ".")
  
  l_bam <- scanBam(bam_file, param = ScanBamParam(what = what_to_extract, flag = scanBamFlag(isUnmappedQuery = FALSE)))
  df_bam <- do.call("data.frame", l_bam)
  
  rm_call = paste(sam_file,bam_file);rm_call
  system(paste("rm",rm_call))
  rm_call = paste(tmp_dir,"*ebwt",sep="/");rm_call
  system(paste("rm",rm_call))
  return(df_bam)
}
split_ids=function(vector_ids,by){
  ids = strsplit(vector_ids,by)
  ids = data.frame(matrix(unlist(ids), nrow = length(ids), byrow = T),stringsAsFactors = F)
  return(ids)
}

fasta_file_to_dataframe=function(contig_file,by=","){
  contigs = readFASTA(contig_file)
  contigs = data.frame(unlist(contigs),stringsAsFactors = F)
  ids = split_ids(rownames(contigs),by)
  contigs = cbind(ids,contigs) 
}

get_min_counts=function(threshold,base_counts_file){
  if(threshold==0){
    min_counts = 0
  }else{
    data = read.table(base_counts_file, header=T)
    data = data[order(data$base_counts,decreasing=T),]
    nrow(data)
    top_percent = threshold
    idx = ceiling(nrow(data) * (top_percent/100) )
    min_counts = data[idx,2]; 
    print(paste("min_counts at ",top_percent,"are:",min_counts))
    print(paste("contigs/clusters ",nrow(data[data$base_counts>min_counts,]))) 
  }
  return(min_counts)
}

print_control=function(option){
  
  #unique reads
  x1 = unique(contigs_reads$read_id)
  x2 = unique(final$read_id)#22g
  print(paste("unique reads in contigs/clusters:",length(x1)))
  print(paste("unique reads in .contigs/clusters file:",length(x2)))
  print(paste("missing unique reads = ", length(x1)  - length(x2)))
  lost_reads = setdiff(x1,x2)
  lost_reads = contigs_reads[contigs_reads$read_id %in% lost_reads,]
  print(head(lost_reads,20))
  write.table(lost_reads,paste(md,"missing_reads.tab",sep="/"),col.names = T,row.names = F,quote = F,sep="\t")
  
  
  #reads counts
  l1 = sum(reads$counts[reads$read_id %in% x1])
  l2 = sum(final$final_counts)      
  l3 = sum(ccounts$final_counts)
  print(paste("reads counts mapped to contigs/clusters:",l1))
  print(paste("reads counts in .counts file:",l3))
  print(paste("mising counts (mapped-.counts):",l1-l3))
  if(length(grep("SpliT",option)) == 1){
    lost_counts = sum(nz[!duplicated(nz$read_id),"read_counts"]) - sum(sr)
    print(paste("missing_counts due rounding :",lost_counts))
    print(paste("mapped - (.counts + missing_counts) :",l1-(l3+lost_counts)))
    print(paste("reads counts in .cr file:",l2))  
  }
  print(".counts and .cr must be equal")

  ##
  print("mmr to contigs");print(length(mmr))##
  print("reads to contigs");print(summary(as.data.frame(df_tab)$Freq))##
  print(paste("total contigs",length(unique(contigs_reads$contig_id))))##
  
  if(length(grep("SpliT",option)) == 1){
    print(paste("contigs having non-cero counts:",sum(contigs_counts$final_counts>0)))
    print(paste("amb reads:",nrow(amb_reads)))
    sum_amb_counts = sum(amb_reads$final_counts)
    print(paste("amb reads counts:",sum_amb_counts))
  }
  
  #global checkpoint
  nreads = nrow(reads)
  cr_reads = length(unique(final_f$read_id))
  print(paste("TOTAL unique reads in lib:",nreads))
  print(paste("TOTAL unique reads in .cr file:",  cr_reads))
  print(paste("diff = ",nreads-cr_reads))
  
  l4 = sum(as.numeric(final_f$final_counts))      
  l5 = sum(as.numeric(ccounts_f$final_counts))
  print(paste("TOTAL COUNTS in lib:",sum(reads$counts)))
  print(paste("TOTAL COUNTS in .counts:",l5))
  print(paste("TOTAL COUNTS in .cr:",l4))
}


#################################MAIN PROGRAM#####################################
dir = paste(md,tmp_dir,sep="/")#switch when necessary
dir.create(dir)
check_options(option)

#===calculate min_counts in threshold===#
#lib = fasta_file_to_dataframe(lib_file,",")
#colnames(lib) = c("read_id","counts","sequence")
#N = sum(as.numeric(lib$counts))
#min_counts = floor ((threshold * N) / 1000000 )
#print(paste("min_counts:",min_counts))
#min_counts = threshold

#===map library to alg+unitigs reference===#
#alg are reads that have alignment to contigs/clusters
#unitigs are reads that do not map to contigs/clusters
#unitigs do not exist in hostpar libraries but in host libraries
alg_unitigs_fasta = paste(dir,"alg_unitigs.fasta",sep="/")
cat_call = paste("cat",alg_file,unitigs_file,">",alg_unitigs_fasta);cat_call
system(cat_call)
mapped_reads = map_reads_to_reads(lib_file,alg_unitigs_fasta,dir)

tmp_fa = paste(dir,"tmp.fa",sep="/")
tbl_to_fasta(mapped_reads[,c("seq","rid","qcounts","qlen")],tmp_fa)


#===map reads to ref contigs/clusters===#
mbd_params = list(bowtie_params = "-S -v 0 -p 20 -a", 
                  query_type="fasta", 
                  tmp_dir=dir, 
                  what_to_extract=c("qname","rname"),
                  unmmaped="un_to_contigs.fa")
df_bam = map_bam_df(tmp_fa,contig_file,mbd_params)

#=== contigs-reads table===#
df_bam$qname = as.character(df_bam$qname)
df_bam$rname = as.character(df_bam$rname)
query_ids = split_ids(df_bam$qname,",")
ref_ids = split_ids(df_bam$rname,",")
contigs_reads = data.frame(contig_id=ref_ids$X1,read_id=query_ids$X2,read_counts=as.numeric(query_ids$X3),read_length=as.numeric(query_ids$X4),stringsAsFactors = F)
contigs_reads$contig_id = gsub("s","c",contigs_reads$contig_id)

#=== Unique-mapping reads (umr)/Multi-mapping reads (mmr) to contigs/clusters===#
df_tab = table(contigs_reads$read_id)
umr = names(df_tab)[df_tab == 1]#umr
mmr = names(df_tab)[df_tab > 1]#mmr

#=== base counts===#
if(option=="informed_SppC" || option == "informed_SpliT"){
  informed_contigs_counts = read.table(base_counts_file,header=T,stringsAsFactors = F)
  rownames(informed_contigs_counts) = informed_contigs_counts$contig_id
}
#=== base counts===#
umr_contigs_reads = contigs_reads[contigs_reads$read_id %in% umr,]
base_counts = tapply(umr_contigs_reads$read_counts,factor(umr_contigs_reads$contig_id),sum)
contigs_counts = left_join(data.frame(contig_id=unique(contigs_reads$contig_id),stringsAsFactors=F), 
                           data.frame(contig_id=names(base_counts),base_counts=base_counts,stringsAsFactors=F),by="contig_id")
contigs_counts$base_counts[is.na(contigs_counts$base_counts)] = 0

#=== write experiment base counts and exit===#
if(option=="write_counts_hostpar" || option=="write_counts_host" || option=="write_counts_par"){
  base_name = sub("write_counts","base_counts",option)
  contigs_counts_file = paste(base_name,"txt",sep=".")
  contigs_counts_file = paste(md,contigs_counts_file,sep="/")
  write.table(contigs_counts,contigs_counts_file,row.names = F,col.names = T,quote = F,sep="\t")
  print("global counts were correctly assigned, please ignore this error messge")
  stop(paste("counts of experiment:",base_name))
}
#======#

#===get contigs/clusters ids===#
contigs = fasta_file_to_dataframe(contig_file)[1]
names(contigs) = c("contig_id")
contigs$contig_id = gsub("s","c",contigs$contig_id)

#===get aligned seqs===#
alg = fasta_file_to_dataframe(alg_file,",")
names(alg) = c("read_id","counts","length","sequence")#alg$contig_id = ids[,1]

#===get reads (library) sequences===#
reads = fasta_file_to_dataframe(tmp_fa,",")
reads = reads[,-1]
colnames(reads) = c("read_id","counts","length","sequence")
reads$counts = as.numeric(reads$counts)
reads$length = as.numeric(reads$length)

#===subset mmr on contigs_reads ===#
mmr_contigs_reads = contigs_reads[contigs_reads$read_id %in% mmr,]

#=========================#
if(length(grep("SppC",option)) == 1){
  if(option == "informed_SppC"){
    contigs_counts = informed_contigs_counts[unique(mmr_contigs_reads$contig_id),]
  }

  #join contigs with contigs_counts
  ccounts = left_join(contigs,contigs_counts,by="contig_id")
  ccounts$base_counts[is.na(ccounts$base_counts)] = 0
  ccounts$contig_id = gsub("s","c",ccounts$contig_id)
  
  #===1) ##write counts to a dataframe===#
  #get mmr reads seqs
  mmr_reads = tapply(mmr_contigs_reads$read_counts,factor(mmr_contigs_reads$read_id),function(x)sample(x,1))
  
  #map mmr to alg
  alg = left_join(alg[,c("read_id","sequence")], data.frame(read_id=names(mmr_reads),counts=mmr_reads,stringsAsFactors=F), by="read_id")
  alg$counts[is.na(alg$counts)] = 0

  ccounts = rbind(ccounts[,c("contig_id","base_counts")],data.frame(contig_id=alg$read_id, base_counts=alg$counts))
  mmr_alg = merge(data.frame(read_id=names(mmr_reads),counts2=mmr_reads,stringsAsFactors=F), reads, by="read_id")####used in 22g
  
  #2)
  #===umr reads can be easily tracked unsing the umr table===#
  cr = contigs_reads[contigs_reads$read_id %in% umr,c("contig_id","read_id","read_counts")]
  colnames(cr) = c("contig_id","read_id","final_counts")
  
}else if(length(grep("SpliT",option)) == 1){
  #===Get support counts from contigs by each mmmr ===#
  #add contigs_counts to mmr_contigs_reads #(the trick) to get the support counts sum
  if(option == "informed_SpliT"){
    mmr_contigs_counts = informed_contigs_counts[unique(mmr_contigs_reads$contig_id),]
    min_counts = get_min_counts(threshold,base_counts_file)
  }else{
    mmr_contigs_counts = contigs_counts[unique(mmr_contigs_reads$contig_id),]
    #min_counts = eo #this parameter remains undefinide in non_split algorithm
  }
  mmr_contigs_reads = merge(mmr_contigs_reads,mmr_contigs_counts,by="contig_id")
  mmr_contigs_reads$base_counts[mmr_contigs_reads$base_counts<=min_counts] = 0##explain this
  head(mmr_contigs_reads[order(mmr_contigs_reads$read_id),],200)
  
  #===Get support counts sum by each mmr===#
  bcs = tapply(mmr_contigs_reads$base_counts,factor(mmr_contigs_reads$read_id),sum)#base_counts sum for each mmr
  mmr_contigs_reads = merge(mmr_contigs_reads,data.frame(bcs=bcs,read_id=names(bcs)),by="read_id")
  head(mmr_contigs_reads[order(mmr_contigs_reads$read_id),],200)
  
  #===reads whose contigs have non zero supported_reads_sum, otherwise will remain being ambiguous reads===#
  #zero means: more than "min_counts" support counts, miminal for sifnificance in poisson test
  nz = mmr_contigs_reads[mmr_contigs_reads$bcs>min_counts,] #mmr_contigs_reads, non_zero
  nz = nz[order(nz$read_id,nz$base_counts,decreasing = c(T,F)),]
  
  #===calculate the supported counts fraction and SPLIT COUNTS according to===#
  nz$bcf = nz$base_counts/nz$bcs
  nz$shared = floor(nz$bcf * nz$read_counts)
  sr = tapply(nz$shared,factor(nz$contig_id),sum)# get the sum of shared reads
  head(nz,100)
  #===#
  
  #===get first round of ambiguous reads===#
  z = mmr_contigs_reads[mmr_contigs_reads$bcs<=min_counts,] #cero bcs
  z = z[order(z$read_id),]
  amb_reads = tapply(z$read_counts,factor(z$read_id),function(x)sample(x,1))
  amb_reads = data.frame(read_id=names(amb_reads),final_counts=amb_reads,stringsAsFactors=F)
  head(amb_reads)

  #add the shared reads to base_counts
  if(length(sr)>0){
    contigs_counts = left_join(contigs_counts, data.frame(contig_id=names(sr),sr_counts=sr,stringsAsFactors=F), by="contig_id")
    contigs_counts$sr_counts[is.na(contigs_counts$sr_counts)] = 0
    contigs_counts$final_counts = contigs_counts$base_counts + contigs_counts$sr_counts 
  }else{
    contigs_counts$final_counts = contigs_counts$base_counts
  }
  
  #===1) ##write counts to a dataframe===#
  #join contigs with contigs_counts
  ccounts = left_join(contigs,contigs_counts,by="contig_id")
  ccounts$final_counts[is.na(ccounts$final_counts)] = 0
  
  #map amb to alg
  if(nrow(amb_reads)>0){
    alg = left_join(alg,amb_reads,by="read_id")
    alg$final_counts[is.na(alg$final_counts)] = 0
    ccounts = rbind(ccounts[,c("contig_id","final_counts")],data.frame(contig_id=alg$read_id,final_counts=alg$final_counts))
    mmr_alg = reads[reads$read_id %in% amb_reads$read_id,]####used in 22g 
  }else{
    alg$final_counts = 0
    ccounts = rbind(ccounts[,c("contig_id","final_counts")],data.frame(contig_id=alg$read_id,final_counts=alg$final_counts))
    mmr_alg = data.frame()
  }
  
  
  #2)
  #===umr reads can be easily tacked unsing the umr table===#
  cr = contigs_reads[contigs_reads$read_id %in% umr,c("contig_id","read_id","read_counts")]
  colnames(cr) = c("contig_id","read_id","final_counts")
  
  #===mmr reads were asigned in nz steps===#
  cr = rbind(cr, data.frame(contig_id=nz$contig_id,read_id=nz$read_id,final_counts=nz$shared))
  cr = cr[cr$final_counts>0,]
  cr = cr[order(cr$contig_id),]
}

##
colnames(ccounts) = c("contig_id","final_counts")

#add unitigs
unitigs_ref = fasta_file_to_dataframe(unitigs_file,",")
colnames(unitigs_ref) = c("ref_id","counts","length","sequence")
if(length(list.files(dir,"un_to_contigs.fa"))==1){
  #get unitigs_query
  unitigs_query_fasta = paste(dir,"un_to_contigs.fa",sep="/")
  unitigs_query = fasta_file_to_dataframe(unitigs_query_fasta,",")
  colnames(unitigs_query) = c("old_id","ref_id","counts","length","sequence")
 
  #map unitigs_query to unitigs_ref 
  unitigs_ref = left_join(unitigs_ref,data.frame(ref_id=unitigs_query$ref_id,final_counts=unitigs_query$counts,stringsAsFactors=F),by="ref_id")
  unitigs_ref$final_counts[is.na(unitigs_ref$final_counts)] = 0
  #unitigs_ref$final_counts[unitigs_ref$ref_id %in% unitigs_query$ref_id] = unitigs_query$counts
}else{
  unitigs_ref$final_counts = 0
}
ccounts_f = rbind(ccounts,data.frame(contig_id=unitigs_ref$ref_id, final_counts=unitigs_ref$final_counts))
tally_out_f = paste(md,paste(gsub(".fa","",basename(lib_file)),".counts",sep=""),sep="/")
write.table(ccounts_f, file = tally_out_f, quote = FALSE, sep = "\t", col.names = F,row.names = F)

#===Prepare data to 22G enrichment analysis===#
cr = left_join(cr,reads[,c("read_id","sequence")],by="read_id")
mmr_un_22g = data.frame(contig_id=mmr_alg$read_id, read_id=mmr_alg$read_id, final_counts=mmr_alg$counts, sequence=mmr_alg$sequence)
if(length(list.files(dir,"un_to_contigs.fa"))==1){
  ur = unitigs_ref[unitigs_ref$final_counts>0,]
  ur = data.frame(contig_id=ur$ref_id, read_id=ur$ref_id, final_counts=ur$final_counts, sequence=ur$sequence)
  final = rbind(cr,mmr_un_22g)
  final_f = rbind(final,ur)
}else{
  final = rbind(cr,mmr_un_22g)
  final_f = final
}
contigs_reads_file = paste(md,paste(gsub(".fa","",basename(lib_file)),".cr",sep=""),sep="/")
write.table(final_f,contigs_reads_file,quote = F,row.names = F,col.names = T,sep="\t")

print_control(option)
#===rm temp files===#
rm_files = paste(dir,"*",sep="/")
system(paste("rm",rm_files))





