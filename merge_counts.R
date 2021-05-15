#10/jan/21
#Dr. Obed Ramirez
#Tabulate individual contigs_counts

rm(list=ls())
computer="production"
if(computer=="test"){
  setwd("/home/macross/compartida/data/miRNA/mm_hp")
  contigs_tag = "NA" #NA or dea_contigs/contigs_tag.fa
  uniqueReads_tag = "mapped_data/host_hostPar_tag.fa"
  hostPar_dir = "dea_reads/hostPar"
  host_dir = "dea_reads/host"
  output = "dea_reads/contigs_tag.rds"
}
if(computer=="production"){
  args = commandArgs(trailingOnly=TRUE)
  contigs_tag = as.character(args[1])
  uniqueReads_tag = args[2]
  hostPar_dir = args[3]
  host_dir = args[4]
  output = args[5]
}

#---prepare big Matrix---#
if(contigs_tag!="NA"){
  contigs_ids = system(paste("grep \">\"",contigs_tag),intern=T)
  contigs_ids = sub(">","",contigs_ids)
  contigs_ids = data.frame(matrix(unlist(strsplit(contigs_ids,",")),ncol=4,byrow=T))
}else{
  contigs_ids = data.frame()
}

reads_ids = system(paste("grep \">\"",uniqueReads_tag),intern=T)
reads_ids = sub(">","",reads_ids)
reads_ids = data.frame(matrix(unlist(strsplit(reads_ids,",")),ncol=4,byrow=T))

system(paste("rm -r",paste(hostPar_dir,"tmp_dir",sep="/")))
system(paste("rm -r",paste(host_dir,"tmp_dir",sep="/")))
samples = c(list.files(hostPar_dir,full.names=T),
            list.files(host_dir,full.names=T))

nrows = nrow(contigs_ids) + nrow(reads_ids)
ncols = length(samples)
tcounts = data.frame(matrix(0,nrow=nrows,ncol=ncols))

n = length(unlist(strsplit(basename(samples[1]),"[.]")))
colnames(tcounts) = matrix(unlist(strsplit(basename(samples),"[.]")),ncol=n,byrow=T)[,1]
rownames(tcounts) = c(contigs_ids$X1,reads_ids$X1)

for(i in 1:ncols){
  print(paste("merging count files:",samples[i]))
  scounts = read.table(samples[i],header=F)
  tcounts[scounts[,1],i] = scounts[,2]
}

tcounts$tag = c(contigs_ids$X4,reads_ids$X4)
save(tcounts,file=output)
