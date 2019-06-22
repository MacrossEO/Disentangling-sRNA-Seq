#tag ids DEA table
rm(list=ls())
library(protr)

computer="mazorka"
if(computer=="mac"){
  md = "/Users/macross/Documents/data/sRNAs_assembly/assembly_eval/mu_ls/mapped/muls/DEA_contigs_informed_0.1"
  amb = "/Users/macross/Documents/data/sRNAs_assembly/assembly_eval/mu_ls/mapped/muls/amb_contigs"
  counts_table = "unitigs.tbl"
  counts_table = paste(md,counts_table,sep="/")
  par_ids = "contigs_alg_unitigs_mPar.fa"
  par_ids = paste(amb,par_ids,sep="/")
  host_ids = "contigs_alg_unitigs_mHost.fa"
  host_ids = paste(amb,host_ids,sep="/")
  amb_ids = "contigs_alg_unitigs_amb.txt"
  amb_ids = paste(amb,amb_ids,sep="/")
  contigs_reads_unitigs_file = "contigs_alg_unitigs.fa"
  contigs_reads_unitigs_file = paste(amb,contigs_reads_unitigs_file,sep="/")
  bt2_file = "contigs_alg_unitigs_un_f.fa"
  bt2_file = paste(amb,bt2_file,sep="/")
  only_reads = "T"#T or F

}
if(computer=="mazorka"){
  args = commandArgs(trailingOnly=TRUE)
  counts_table = args[1]
  par_ids = args[2]
  host_ids = args[3]
  amb_ids = args[4]
  contigs_reads_unitigs_file = args[5]
  bt2_file = args[6]
  only_reads = args[7]
  amb = args[8]
}

counts = read.table(counts_table,header=T,sep="\t")
counts$origin = NA
counts$len = 0

#mmhp counts and seq length
cru = readFASTA(contigs_reads_unitigs_file)
ids = strsplit(names(cru),",")
ids = data.frame(matrix(unlist(ids), nrow = length(ids), byrow = T),stringsAsFactors = F)

#checkpoint
all(rownames(counts)==ids$X1)

#
counts[ids$X1,"len"] = ids$X3

#host
host_fasta = readFASTA(host_ids)
ids = strsplit(names(host_fasta),",")
ids = data.frame(matrix(unlist(ids), nrow = length(ids), byrow = T),stringsAsFactors = F)
counts[ids$X1,"origin"] = "h"

#par
par_fasta = readFASTA(par_ids)
ids = strsplit(names(par_fasta),",")
ids = data.frame(matrix(unlist(ids), nrow = length(ids), byrow = T),stringsAsFactors = F)
counts[ids$X1,"origin"] = "p"

#amb
amb_txt = read.table(amb_ids,sep=",",header=F,stringsAsFactors = F)
counts[amb_txt$V1,"origin"] = "a"

#bt2
if(length(grep(basename(bt2_file),list.files(amb))) != 0){
  bt2_fasta = readFASTA(bt2_file)
  if(length(bt2_fasta)>1){
    ids = strsplit(names(bt2_fasta),",")
    ids = data.frame(matrix(unlist(ids), nrow = length(ids), byrow = T),stringsAsFactors = F)
    ids$contig_id = gsub("[hpa]","c",ids$X2)
    ids$origin = gsub("[0-9]","",ids$X2)
    counts[ids$contig_id,"origin"] = ids$origin
  }
}

print(unique(counts$origin))
print(paste("amb rows txt-dtframe",nrow(amb_txt),sum(counts$origin=="a")))
#counts[head(amb_txt)[,1],"origin"]

if(only_reads=="T"){
  idx = grep("c",rownames(counts))
  counts = counts[-idx,]
}
out_file = gsub(".tbl","_extra.tbl",counts_table);out_file
write.table(counts,out_file,row.names = T,col.names = T,sep="\t",quote = F)


